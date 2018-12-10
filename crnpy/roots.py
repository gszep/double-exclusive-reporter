from numpy import array,log10,ones
from numpy.random import uniform

from scipy import optimize
from joblib import Parallel,delayed


def root(function,interval=[-5,5],args=(),nvar=1, logspace=False) :

    variable = uniform(*tuple(interval),size=nvar)
    initial_guess = 10.0**variable if logspace else variable
    optimization = optimize.root(function,initial_guess,args=args)

    if optimization.success :
        return optimization.x
    else :
        return root(function,interval,args,nvar)


def roots_parallel(function, interval, args=[()], nvar=1, logspace=False ):
    '''roots() function parallelised across a list of arguments

    ---parameters---
    function : <function>
        scalar function of one variable x, with possible arguments
    interval : [<float>,<float>]
        a list of two numbers; interval to search for roots
    args : [...<tuple>...]
        list of tuples for additional arguments to parallelise across
    epsilon : <float>
        the precision to use when searching the interval
    n_roots : <int>
        maximum number of roots to expect; defaults to no assumptions

    ---returns---
    roots : ...<1darray>...
        tuples of array of roots of f(x)
    '''
    root_list = Parallel(n_jobs=-1)(
                delayed(root)(function,interval=interval,args=arg,nvar=nvar, logspace=logspace)
                for arg in args)

    return array(root_list)
