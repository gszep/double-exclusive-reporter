from numpy import array,unique,log,exp,ones,zeros,mean,inf,gradient,linspace,sqrt,amax,append,full,prod,matmul
from numpy.random import uniform

from sympy.functions.combinatorial.factorials import binomial
from scipy.ndimage import laplace

from collections import defaultdict
from re import sub,search,finditer,findall
from yaml import load

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

        ###### converting to dimensionless parameters ######

        # # inihibitions
        # self.alpha = array([(self.capacity*self.aR33)/(self.dR + self.growth),(self.capacity*self.aS175)/(self.dS+ self.growth)]) **2
        # self.nu = array([self.dR,self.dS]) + self.growth
        #
        # # activations
        # self.beta = self.capacity * array([(self.aL*self.a1R*self.KGR_76)/(self.dL + self.growth),(self.aT*self.a1S*self.KGS_81)/(self.dT + self.growth)])
        # self.mu = array([self.dL,self.dT]) + self.growth
        #
        # # baseline inhibitor production and  saturation
        # self.omega = array([self.a0_76/(self.a1R*self.KGR_76),self.a0_81/(self.a1S*self.KGS_81)])
        # self.Omega = array([self.KGR_76,self.KGS_81])
        #
        # # signalling hill functions
        # self.n = array([self.nR,self.nS])
        # self.exponents = array([self.nT,self.nL])
        #
        # # signalling dissociation constants
        # self.k = array([1.0/self.KR6,1.0/self.KS12])
        #
        # # diffusion coefficients
        # self.D = array([self.c6,self.c12])
        #
        # # one dimensional lattice of states
        # self.epsilon = 1e-5
        # self.state = { 'diffusables':full((self._nx,2),self.epsilon),
        #                'activators':full((self._nx,2),self.epsilon),
        #                'inhibitors':full((self._nx,2),self.epsilon) }

        self.space = linspace(0,self._xmax,self._nx)
        self.time = 0.0
        self.dx = self.space[1]-self.space[0]


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
                value = sub(attr_name,'self.'+attr_name,value)

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
            #self._propensity[i] = compile(self._propensity[i],'<string>','eval')

        fget = lambda self : array([ eval(self._propensity[i]) for i in range(self.n_reactions) ])
        self.__class__.propensity = property(fget)

        self._nullcines = [ '+'.join(['0.0']+[ '*'.join([str(c),rate])
                            for c,rate in zip(cs,self._nullcines) if c != 0 ])
                            for cs in self.stoichiometry ]

        # filtering out trivial nullcines
        self.nontrivials = [ self.names[k] for k,nullcine in enumerate(self._nullcines) if nullcine != '0.0' ]
        self._nullcines = [ nullcine for k,nullcine in enumerate(self._nullcines) if self.names[k] in self.nontrivials ]
        self.n_nontrivials = len(self.nontrivials)


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
        steady_states = roots_parallel(self.nullcines, interval=[-3,3], args=args, nvar=self.n_nontrivials)

        # # applying clipping threshold
        # if clip is not None :
        #
        #     x,y = c.T
        #     for i in xrange(n_roots) :
        #
        #         mask = (L<clip)[:,i]
        #         mask = (T<clip)[:,i]
        #
        #         L[:,i][(x<mean(x))*(y<mean(y))*mask] = -inf
        #         T[:,i][(x<mean(x))*(y<mean(y))*mask] = -inf

        # return original shaped array
        output_shape = input_shape[:-1] + (self.n_nontrivials,)
        steady_states.shape = output_shape
        return steady_states


def fromcrn(file_path) :
    '''Creates model object from parameters in a given crn file

    --- parameters ---
    file_path : <str>
        path to crn file to read parameters from

    --- returns ---
    model : <Model>
        model object initialised with parameters
    '''

    # open and read file as string
    kwargs = { }
    with open(file_path, 'r') as file:

        # import as string without comments
        file_string = sub(r'(//).*',' ',file.read())

        # parse reactions
        kwargs['reactions'] = []
        pattern = r'[\w\[\]\{\}\(\)*/\-+\^ ]*->[\w\[\]\{\}\(\)*/\-+\^ ]*'
        for match in finditer(pattern,file_string) :

                reaction = match.group()
                reaction = reaction.replace('^','**')
                kwargs['reactions'] += [ reaction ]

        # convert to yaml compatible syntax
        file_string = file_string.replace(']','}')
        file_string = file_string.replace('[','{')
        file_string = file_string.replace(';',',')
        file_string = file_string.replace('=',':')

        # find assigned parameters
        pattern = r'([0-9|\w]+)[ ]*:[ ]*[\w0-9.eE+-]+'
        for match in finditer(pattern,file_string) :

                item = match.group()
                name,value = item.split(':')

                # format as yaml
                value = value.strip() if isnumber(value) else '"'+value.strip()+'"'
                yaml_format = '"'+name.strip()+'":'+value
                file_string = file_string.replace(item,yaml_format)

        # parse each directive to dictionary
        pattern = r'(?<=directive ).*?(?=(directive|init))'
        for match in finditer(pattern,file_string.replace('\n',' ')) :
            item = match.group()

            # attempt to parse automagically
            if 'simulator' not in item :
                key,value = item.split(" ",1)
                try : kwargs[key] = load(value)

                # otherwise parse manually
                except :

                    kwargs[key] = {}
                    for s in value[1:-1].split(',') :
                        name,value = s.strip().split(':')

                        # convert to python compatible syntax
                        value = value.replace('^','**')

                        value = value.replace('"','').replace('{','').replace('}','')
                        name = name.replace('"','').replace('{','').replace('}','')

                        value = float(value) if isnumber(value) else value.strip()
                        kwargs[key][name.strip()] = value

        # parse initial conditions
        kwargs['init'] = {}
        pattern = r'init[ ]*\w+[ ]*[0-9.eE+-]+'
        for match in finditer(pattern,file_string) :

                item = match.group()
                _,name,value = item.split(' ')
                kwargs['init'][name.strip()] = float(value)

    return Model(**kwargs)


def isnumber(value):
    try:
        float(value)
        return True

    except (ValueError,TypeError):
        return False
