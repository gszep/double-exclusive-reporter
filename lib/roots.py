from numpy import array,log10,ones
from scipy.optimize import bisect
from joblib import Parallel, delayed


def rootsearch(f,interval,dx=1e-3,args=()):
    '''Searches an interval from below in increments of dx for
    locations where f(x) and f(x+dx) have opposite signs, suggesting
    the existance of a root

    ---parameters---
    f : <function>
        scalar function of one variable x, with possible arguments
    interval : [<float>,<float>]
        a list of two numbers; interval to search for roots
    dx : <float>
        the increment to use when searching the interval
    args : <tuple>
        additional arguments to use in f(x,args)

    ---returns---
    (x,x+dx) : <tuple>
        an infinitessimal interval that contains a root of f(x)
    '''
    x,xmax = interval

    # find interval where function crosses zero
    while f(x,*args)*f(x+dx,*args) > 0.0:

        # stop once whole interval has been covered
        if xmax <= x:
            return None

        # move along interval
        x += dx

    return x,x+dx


def roots(function, interval, args=(), epsilon=1e-3, n_roots=None ):
    '''Returns an array of roots of function within the given interval

    ---parameters---
    function : <function>
        scalar function of one variable x, with possible arguments
    interval : [<float>,<float>]
        a list of two numbers; interval to search for roots
    args : <tuple>
        additional arguments to use in f(x,args)
    epsilon : <float>
        the precision to use when searching the interval
    n_roots : <int>
        maximum number of roots to expect; defaults to no assumptions

    ---returns---
    roots : <1darray>
        array of roots of f(x)
    '''

    precision = -int(log10(epsilon))
    lower,upper = interval
    roots = []

    while True: # find interval around a root
        interval = rootsearch(function,[lower,upper],epsilon,args)
        if interval is not None:

            # return root within interval
            root = bisect(function,*interval,args=args)
            roots += [round(root,precision)]

            if len(roots) is n_roots :
                return ones(n_roots)*array(roots)


            # exclude root from new search interval
            _,lower = interval

        else: # break loop once search interval closes
            return ones(n_roots)*array(roots)


def roots_parallel(function, interval, args=[()], epsilon=1e-3, n_roots=None ):
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
                delayed(roots)(function,interval=interval,args=arg,
                n_roots=n_roots,epsilon=epsilon) for arg in args)

    return array(root_list)
