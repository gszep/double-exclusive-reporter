from numpy import array,log,exp,ones,zeros,mean,inf,gradient,linspace,sqrt,amax,append,full
from numpy.linalg import norm
from scipy.ndimage import laplace
<<<<<<< HEAD

from .roots import roots_parallel
from re import sub,search,finditer
from yaml import load


class Model(object) :
    '''Contains all we need to calculate steady states and simulate
    a chemical reaction system parsed from a crn file.'''
=======
from roots import roots_parallel
>>>>>>> parent of 7452ce7... making get_scaffold.py compatible with python3

    def __init__(self,**kwargs) :

        # parse all set parameters to attributes
        for key in kwargs :
            setattr(self, key, kwargs[key])

        # parse out parameters and rates
        for param in self.parameters :
            setattr(self, param, self.parameters[param])

        # substitute in parameter values
        for rate in self.rates :
            if self.rates[rate] in self.parameters :
                setattr(self, rate, self.parameters[self.rates[rate]])
            else :
                setattr(self, rate, self.rates[rate])

        # convert boundary condition keyword
        if self.spatial['boundary'] == 'ZeroFlux' : self.boundary = 'nearest'
        if self.spatial['boundary'] == 'Periodic' : self.boundary = 'wrap'

        ###### converting to dimensionless parameters ######

        # inihibitions
        self.alpha = array([(self.capacity*self.aR33)/(self.dR + self.growth),(self.capacity*self.aS175)/(self.dS+ self.growth)]) **2
        self.nu = array([self.dR,self.dS]) + self.growth

        # activations
        self.beta = self.capacity * array([(self.aL*self.a1R*self.KGR_76)/(self.dL + self.growth),(self.aT*self.a1S*self.KGS_81)/(self.dT + self.growth)])
        self.mu = array([self.dL,self.dT]) + self.growth

        # baseline inhibitor production and  saturation
        self.omega = array([self.a0_76/(self.a1R*self.KGR_76),self.a0_81/(self.a1S*self.KGS_81)])
        self.Omega = array([self.KGR_76,self.KGS_81])

        # signalling hill functions
        self.n = array([self.nR,self.nS])
        self.exponents = array([self.nT,self.nL])

        # signalling dissociation constants
        self.k = array([1.0/self.KR6,1.0/self.KS12])

        # diffusion coefficients
        self.D = array([self.c6,self.c12])

        # one dimensional lattice of states
        self.epsilon = 1e-5
        self.state = { 'diffusables':full((self.nx,2),self.epsilon),
                       'activators':full((self.nx,2),self.epsilon),
                       'inhibitors':full((self.nx,2),self.epsilon) }

        self.space = linspace(0,self.xmax,self.nx)
        self.time = array([0.0])
        self.dx = self.space[1]-self.space[0]

    def laplacian(self,states):
        return array([ laplace(state,mode=self.boundary)/self.dx**2 for state in states.T ]).T

    def time_step(self) :

        self.time = append(self.time,self.time[-1]+self.dt)
        self.c = ( self.state['diffusables']/(self.k+self.state['diffusables']) )**self.n

        self.state['inhibitors'] += self.mu*(
            self.beta*(self.c*self.state['activators']**2+self.omega)\
                     /(self.c*self.state['activators']**2*self.Omega+1)\
                            - self.state['inhibitors']) * self.dt

        self.state['activators'] += self.nu*(
            sqrt(self.alpha)/(1+self.state['inhibitors'][:,::-1]**self.exponents)\
                              - self.state['activators'] ) * self.dt

        self.state['diffusables'] += self.D*self.laplacian(self.state['diffusables']) * self.dt
        return True

    def get_diffusives(self,c) :
        '''return the concentrations of c6 and c12 in nM for given
        dimensionless signal inputs c and cdash '''
        return ( self.k / (self.alpha**((1-c)/self.n)-1) ).T

    def get_inputs(self,diffusives) :
        '''return the dimensionless signal inputs c and cdash for
        a given concentration of c6 and c12 in nM'''
        return (1+self.n*log(diffusives/(self.k+diffusives))/log(self.alpha) ).T

    def get_nullcines(self,x,c,state=None) :
        '''return implicit the value of nullcines f(x,c), such
        that f(x,c)=0 gives the steady state concentration x,
        given input signal c

        --- parameters ---
        x : <1darray>
            array that has two numbers representing
            concentrations of LacI and TetR
        c : <1darray>
            array that has two numbers between -1 and 1
            representing the input signal range for
            c6 and c12

        --- returns ---
        nullcines : <1darray>
            contains the values of f(x,c) for LacI and TetR
            in array of the same shape as x and c
        '''


        # transform to logspace
        C = self.alpha**c
        x = 10**(x*ones(C.shape))

        # common calculations
        offset = (self.beta-x*self.Omega)*C
        scaling = self.beta*self.omega-x
        rates = zeros(C.shape)

        # calculate nullcines
        rates[0] = (self.beta[1]*(C[1]+self.omega[1]*(1+x[0])**2)/ \
                                 (C[1]*self.Omega[1]+(1+x[0])**2))**4
        rates[1] =  self.beta[0]*(C[0]+self.omega[0]*(1+x[1]**4)**2)/ \
                                 (C[0]*self.Omega[0]+(1+x[1]**4)**2)

        nullcines = offset+scaling*(1+rates)**2
        if state is not None :
            return nullcines[state]
        else :
            return nullcines

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
        n_roots = 3

        # parallel solve for steady states
        args_list = zip(c,c.shape[0]*[0])
        L = roots_parallel(self.get_nullcines,interval=[-3,3],args=args_list,n_roots=n_roots,epsilon=0.01)
        args_list = zip(c,c.shape[0]*[1])
        T = roots_parallel(self.get_nullcines,interval=[-3,3],args=args_list,n_roots=n_roots,epsilon=0.01)

        # applying clipping threshold
        if clip is not None :

            x,y = c.T
            for i in xrange(n_roots) :

                mask = (L<clip)[:,i]
                mask = (T<clip)[:,i]

                L[:,i][(x<mean(x))*(y<mean(y))*mask] = -inf
                T[:,i][(x<mean(x))*(y<mean(y))*mask] = -inf

        # return original shaped array
        output_shape = input_shape[:-1] + (3,)
        L.shape,T.shape = output_shape,output_shape
        return L,T


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
                value = value.strip() if isfloat(value) else '"'+value.strip()+'"'
                json_format = '"'+name.strip()+'":'+value
                file_string = file_string.replace(item,json_format)

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

                        value = float(value) if isfloat(value) else value.strip()
                        kwargs[key][name.strip()] = value

    return Model(**kwargs)


def isfloat(value):
    try:
        float(value)
        return True

    except ValueError:
        return False
