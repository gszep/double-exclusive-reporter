from numpy import array,unique,log,exp,ones,zeros,mean,inf,gradient,linspace,sqrt,amax,append,full,prod,matmul
from numpy.random import uniform

from sympy.functions.combinatorial.factorials import binomial
from scipy.ndimage import laplace

from collections import defaultdict
from re import findall,sub,search

from .utils import isnumber
from .colors import colors

from .roots import roots_parallel
from scipy import optimize

class Model(object) :
    '''Contains all we need to calculate steady states and simulate
    a chemical reaction system parsed from a crn file.'''

    def __init__(self,**kwargs) :

        # default timegrid
        self.time = 0.0
        self._dt = 0.01

        self._parse_args(kwargs)
        self._set_parameters()

        self._set_states(spatial=hasattr(self,'_spatial'))
        self._set_rates()

        self._set_stoichiometry()
        self._set_propensity()

        if hasattr(self,'_spatial') :
            print(colors.yellow,'[model]',colors.reset,'setting spatial directive')

            self._set_boundary()
            self.space = linspace(0,self._xmax,self._nx)
            self.dx = self.space[1]-self.space[0]

        else :
            print(colors.yellow,'[model]',colors.reset,'non-spatial directive')


    def _parse_args(self,kwargs):
        '''parse keyword arguments to private variables'''

        for key,value in kwargs.items() :
            attr_name = '_'+key

            value = self._parse_type(value)
            setattr(self,attr_name,value)


    def _parse_type(self,value):
        '''parse all data types'''
        dtype = type(value)

        if dtype is dict :
            value = self._parse_dict(value)

        elif dtype is list :
            value = self._parse_list(value)

        elif dtype is str :
            value = self._parse_str(value)

        elif dtype is type(None) :
            value = None

        elif not isnumber(value) :
            raise Exception('type {} not parsable'.format(dtype))

        return value


    def _parse_dict(self,value):
        '''parse dictionary types recursivesly'''

        self._parse_args(value)
        return value


    def _parse_list(self,values):
        '''iteratively parse through list types'''

        for i,value in enumerate(values) :

            value = self._parse_type(value)
            values[i] = value

        return values


    def _parse_str(self,value):
        '''parse string types as expressions involving class attributes'''

        # extract unique attribute names from value string
        for attr_name in unique(findall('\w+',value)) :
            if not isnumber(attr_name) :

                # prepend 'self' to reference attributes
                value = sub('\\b'+attr_name+'\\b','self.'+attr_name,value)

        return value


    def _set_property(self,attr_name,value=None,setter=False):
        '''create getter/setter for given attribute name,value pair'''

        # ensure attribute names are private
        assert type(attr_name) is str, 'attr_name must be string'
        attr_name = attr_name if attr_name[0] == '_' else '_'+attr_name

        if value is None : # get value from existing attribute
            def fget(self):
                return getattr(self,attr_name)

        else : # or evaulate given mathematical expression

            # precompile statements for speed boost
            value = compile(str(value),'<string>','eval')

            def fget(self):
                return eval(value)

        if setter : # optionally define setter
            def fset(self,value):
                setattr(self,attr_name,value)
        else :
            fset = None

        # create property for private variable
        setattr(self.__class__,attr_name[1:],property(fget,fset))


    def _set_boundary(self):
        '''parseing in boundary conditions'''

        if hasattr(self, '_boundary'):

            if 'ZeroFlux' in self._boundary :
                self.boundary = 'nearest'

            elif 'Periodic' in self._boundary :
                self.boundary = 'wrap'

            else :
                raise Exception('boundary condition "{}" not recognised'.format(self._boundary))
        else :
            raise Exception('no boundary condition provided')

    def _set_parameters(self):
        '''set attribute for each parameter'''

        for parameter in self._parameters:
            setattr(self,parameter,getattr(self,'_'+parameter))


    def _set_states(self,spatial=True):
        '''create propery that returns current state variables'''

        # set initial conditions
        for state,value in self._init.items():

            # substitute parameter values
            if type(value) is str :
                value = getattr(self,value)

            init = full(self._dimensions*(self._nx,),value) if spatial else value
            setattr(self,state,init)

        # identify state variables
        self.names = list(self._init)
        fget = lambda self : array([ getattr(self,name) for name in self.names ])
        fset = lambda self,values : [ setattr(self,name,values[k]) for k,name in enumerate(self.names) ]
        self.__class__.states = property(fget,fset)

        # identify diffusables
        self.diffusibles = defaultdict(lambda : 0.0)
        if spatial :
            self.diffusibles.update({ key:value for param in self._spatial['diffusibles']
                                      for key,value in param.items() })


    def _set_rates(self):
        '''create propery that calcuates rates from current states'''

        for rate in self._rates:
            expression = getattr(self,'_'+rate)
            self._set_property(rate,expression)


    def _set_stoichiometry(self):
        '''create property that returns stoichiometric coefficient matrix'''

        self.stoichiometry = zeros((len(self.names),len(self._reactions)))

        pattern = r'(?<=->)[ ]*[\[\{].*[\]\}]'
        for i,reaction in enumerate(self._reactions) :

            match = search(pattern,reaction)
            rate = match.group() if match is not None else '{1.0}'

            reaction = reaction.replace(rate,'')
            reactants,products = reaction.split('->')

            self.stoichiometry.T[i] = array([ products.count('self.'+name)- \
                                               reactants.count('self.'+name) \
                                               for name in self.names ])


    def _set_propensity(self):
        '''create property that calculates propensities from current states'''

        self.n_reactions = len(self._reactions)
        self._propensity = self.n_reactions*[None]
        self._nullcines = self.n_reactions*[None]

        pattern = r'(?<=->)[ ]*[\[\{].*[\]\}]'
        for i,reaction in enumerate(self._reactions) :
            match = search(pattern,reaction)

            rate = match.group().strip() if match is not None else '{1.0}'
            action_type = rate[0]+rate[-1]

            # linear rates
            rate = rate.replace('{','').replace('}','').replace('[','').replace(']','')
            self._propensity[i] = str(rate)

            # mass action rule
            if action_type == '{}':
                reactants,products = reaction.split('->')

                state_product = '*'.join([
                    ('('+str(binomial('state',reactants.count('self.'+name))
                        .expand(func=True))+')').replace('state','self.'+name)
                             for name in self.names if reactants.count('self.'+name) != 0 ])

                self._propensity[i] += '*'+state_product

            # precompile statements for speed boost
            self._nullcines[i] = self._propensity[i]


        self._nullcines = [ '+'.join(['0.0']+[ '*'.join([str(c),rate])
                            for c,rate in zip(cs,self._nullcines) if c != 0 ])
                            for cs in self.stoichiometry ]

        # filtering out trivial nullcines
        self.nontrivials = [ self.names[k] for k,nullcine in enumerate(self._nullcines) if nullcine != '0.0' ]
        self._nullcines = [ nullcine for k,nullcine in enumerate(self._nullcines) if self.names[k] in self.nontrivials ]
        self.n_nontrivials = len(self.nontrivials)

        # evaluate remaining expressions TODO(gszep) this is crap.. need refactoring
        for attr_name,value in vars(self).items():
            if type(value) is str :

                for k,nullcine in enumerate(self._nullcines) :
                    expression = '('+getattr(self,attr_name)+')'
                    self._nullcines[k] = sub('self.'+attr_name[1:],expression,nullcine)

        for attr_name,value in vars(self).items():
            if type(value) is str :

                for k,nullcine in enumerate(self._nullcines) :
                    expression = '('+getattr(self,attr_name)+')'
                    self._nullcines[k] = sub('self.'+attr_name[1:],expression,nullcine)

        for attr_name,value in vars(self).items():
            if type(value) is str :

                for k,nullcine in enumerate(self._nullcines) :
                    expression = '('+getattr(self,attr_name)+')'
                    self._nullcines[k] = sub('self.'+attr_name[1:],expression,nullcine)


        for i,_ in enumerate(self._nullcines) :
            self._nullcines[i] = compile(self._nullcines[i],'<string>','eval')

        for i,_ in enumerate(self._propensity) :
            self._propensity[i] = compile(self._propensity[i],'<string>','eval')

        fget = lambda self : array([ eval(self._propensity[i]) for i in range(self.n_reactions) ])
        self.__class__.propensity = property(fget)

    def reaction(self):
        '''return reactive term'''

        # TODO(gszep) matmul each time inefficient?
        return matmul(self.stoichiometry,self.propensity)


    def diffusion(self):
        '''return diffusive term'''

        if hasattr(self,'_spatial') : # TODO(gszep) checking each time inefficient?
            return array([ self.diffusibles[name]*laplace(getattr(self,name), \
                           mode=self.boundary)/self.dx**2 for name in self.names ])

        else : # no diffusion for non-spatial setting
            return array([ 0.0 for name in self.names ])


    def time_step(self) :
        '''propagate reaction-diffusion system using forward euler method'''

        self.states += ( self.reaction() + self.diffusion() )*self._dt
        self.time += self._dt


    def nullcines(self,states,parameters) :
        '''evaluate nontrivial nullcines at given state and parameter values'''

        # updating parameters and state variable
        for k,name in enumerate(self.nontrivials):
            setattr(self,name,states[k])

        for name,value in parameters.items():
            setattr(self,name,value)

        # evaluating nullcines
        return [ eval(self._nullcines[i]) for i in range(self.n_nontrivials) ]

        # NOTE: unexpected behaviour in joblib.Parallel leads to
        # NameError: name 'self' is not defined
        # when the above return statement is rewritten as
        # return [ eval(nullcine) for nullcine in self._nullcines ]

    def get_steady_state(self,c,clip=None) :
        '''solves f(x,c)=0 for steady state concentration x,
        given input signal c, which can be a grid of values'''

        # flattening input
        input_shape = c.shape
        c.shape = (-1,2)

        # parallel solve for steady states
        args = [ { 'c6':c6, 'c12':c12 } for c6,c12 in c ]

        self._set_states(spatial=False)
        steady_states = roots_parallel(self.nullcines, interval=[-5,5], args=args, nvar=self.n_nontrivials)
        self._set_states(spatial=True)

        # return original shaped array
        output_shape = input_shape[:-1] + (self.n_nontrivials,)
        steady_states.shape = output_shape
        return steady_states
