from numpy import array,unique,log,exp,ones,zeros,mean,inf,gradient,linspace,sqrt,amax,append,full,prod,matmul
from numpy.random import uniform

from sympy.functions.combinatorial.factorials import binomial
from scipy.ndimage import laplace

from collections import defaultdict
from re import findall,sub,search
from .utils import isnumber

from .roots import roots_parallel
from scipy import optimize

class Model(object) :
    '''Contains all we need to calculate steady states and simulate
    a chemical reaction system parsed from a crn file.'''

    def __init__(self,**kwargs) :
        self._parse_args(kwargs)

        self._set_boundary()
        self._set_parameters()

        self._set_states()
        self._set_rates()

        self._set_stoichiometry()
        self._set_propensity()

        self.space = linspace(0,self._xmax,self._nx)
        self.dx = self.space[1]-self.space[0]
        self.time = 0.0


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

        if self._boundary == 'self.ZeroFlux' :
            self.boundary = 'nearest'

        elif self._boundary == 'self.Periodic' :
            self.boundary = 'wrap'

        else :
            _,condition = self._boundary.split('.')
            raise Exception('boundary condition "{}" not recognised'.format(condition))


    def _set_parameters(self):
        '''set attribute for each parameter'''

        for parameter in self._parameters:
            setattr(self,parameter,getattr(self,'_'+parameter))


    def _set_states(self,spatial=True):
        '''create propery that returns current state variables'''

        # set initial conditions
        for state,value in self._init.items():

            init = full(self._dimensions*(self._nx,),value) if spatial else value
            setattr(self,state,init)

        # identify state variables
        self.names = list(self._init)
        fget = lambda self : array([ getattr(self,name) for name in self.names ])
        fset = lambda self,values : [ setattr(self,name,values[k]) for k,name in enumerate(self.names) ]
        self.__class__.states = property(fget,fset)

        # identify diffusables
        self.diffusibles = defaultdict(lambda : 0.0)
        self.diffusibles.update(self._spatial['diffusibles'])


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
            rate = match.group()

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

            rate = match.group().strip()
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

        # evaluate remaining expressions
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
        '''return reactive part of system'''

        return matmul(self.stoichiometry,self.propensity)


    def diffusion(self):
        '''return diffusive part of system'''

        return array([ self.diffusibles[name]*laplace(getattr(self,name), \
                       mode=self.boundary)/self.dx**2 for name in self.names ])


    def time_step(self) :
        '''propagate reaction-diffusion system using forward euler method'''

        self.states += ( self.reaction() + self.diffusion() )*self._dt
        self.time += self._dt

    def get_diffusives(self,c) :
        '''return the concentrations of c6 and c12 in nM for given
        dimensionless signal inputs c and cdash '''
        return ( self.k / (self.alpha**((1-c)/self.n)-1) ).T

    def get_inputs(self,diffusives) :
        '''return the dimensionless signal inputs c and cdash for
        a given concentration of c6 and c12 in nM'''
        return (1+self.n*log(diffusives/(self.k+diffusives))/log(self.alpha) ).T

    def nullcines(self,states,parameters) :
        '''return parametrised nontrivial nullcine values'''

        self._set_states(spatial=False)
        for k,name in enumerate(self.nontrivials):
            setattr(self,name,states[k])

        for name,value in parameters.items():
            setattr(self,name,value)

        return [ eval(self._nullcines[i]) for i in range(self.n_nontrivials) ]

    def get_steady_state(self,c,clip=None) :
        '''solves f(x,c)=0 for steady state concentration x,
        given input signal c, which can be a grid of values

        --- parameters ---
        c : <ndarray>
            ndarray whose last axis has two numbers
            between -1 and 1 representing the input
            signal range for c6 and c12
        clip : <float>
            Clipping threshold x**clip, below which
            the 'off' state is defined

        --- returns ---
        x : <ndarray>
            steady state values for LacI and TetR in
            an array of the same shape as c
        '''

        # flattening input
        input_shape = c.shape
        c.shape = (-1,2)

        # parallel solve for steady states
        args = [ { 'c6':c6, 'c12':c12 } for c6,c12 in c ]
        steady_states = roots_parallel(self.nullcines, interval=[-5,5], args=args, nvar=self.n_nontrivials)

        # return original shaped array
        output_shape = input_shape[:-1] + (self.n_nontrivials,)
        steady_states.shape = output_shape
        return steady_states
