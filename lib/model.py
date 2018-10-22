from numpy import array,log,exp,ones,zeros,mean,inf
from roots import roots_parallel

from re import search
from inspect import getargspec

class DoubleExclusive(object) :
    '''Contains all we need to calculate steady states of the
    double exclusive reporter system; assumptions are currently
    that there is no cross-talk between reporters.'''

    def __init__(self,a0_76=1.92,a0_81=0.822,
                 aL=0.001170,a1R=1.86e3,KGR_76=0.00145301678863517,dL=0.000117,
                 aT=0.000951,a1S=704,KGS_81=1.00086836995962000e-5,dT=0.666,
                 aR33=9.1,dR=0.267,aS175=3.58,dS=0.319,
                 KR6=3.50695213045775e-8,nR=0.699,
                 KS12=0.0095893786931253,nS=1.25,
                 growth=1,capacity=67.5) :

        ###### converting to dimensionless parameters ######

        # inihibitions
        self.alpha = array([(capacity*aR33)/(dR + growth),(capacity*aS175)/(dS+ growth)]) **2

        # activations
        self.beta = capacity * array([(aL*a1R*KGR_76)/(dL + growth),(aT*a1S*KGS_81)/(dT + growth)])

        # baseline inhibitor production
        self.omega = array([a0_76/(a1R*KGR_76),a0_81/(a1S*KGS_81)])

        # inhibitor saturation
        self.Omega = array([KGR_76,KGS_81])

        # signalling hill functions
        self.n = array([nR,nS])

        # signalling dissociation constants
        self.k = array([1.0/KR6,1.0/KS12])

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
    '''Creates DoubleExclusive model object from parameters in a given crn file

    --- parameters ---
    file_path : <str>
        path to crn file to read parameters from

    --- returns ---
    model : <DoubleExclusive>
        model object initialised with parameters
    '''

    # open and read file as string
    parameters = {}
    with open(file_path, 'r') as file:
        file_string=file.read().replace('\n', ' ')

        # iterate through model object constructor arguments
        for arg in getargspec(DoubleExclusive.__init__)[0] :
            if arg is not 'self' :

                # use regex to match arguments in file string
                pattern = r'[ ]*=[ ]*[0-9.eE+-]+'
                match = search(arg+pattern,file_string)
                if match is not None :

                    # parse argument names and values
                    name,value = match.group().split('=')
                    name = name.strip()
                    parameters[name] = float(value)

                else :
                    raise Exception('parameter "{}" not found in crn file'.format(arg))

    return DoubleExclusive(**parameters)
