from lib.model import fromcrn
from lib.colors import cyan,yellow

from sys import argv
from argparse import ArgumentParser

from numpy import linspace,logspace,meshgrid,vstack,dstack,inf,mean,ones,array,log10
from matplotlib.pyplot import *

def get_args() :
    '''parse arguments from command line'''
    parser = ArgumentParser(description='creates hysteresis plot from crn file parameters')

    parser.add_argument('--crn_path', type=str,
                        help='path to crn file')
    parser.add_argument('--c12_const', type=float, default=120,
                        help='constant value of c12 in nM')
    parser.add_argument('--c6_range', nargs=2, type=float, default=[1e-1,1e6],
                        help='input range for c6 in nM')
    parser.add_argument('--N', type=int, default=250,
                        help='number of grid points per dimension to use')

    return vars(parser.parse_args())


def get_hysteresis(crn_path,N,c6_range,c12_const):
    '''Calculate hysteresis diagram for double exclusive reporter
    for a given range of diffusives c6 and a constant value of c12.

    ---parameters---
    crn_path : <str>
        path to crn file
    N : <int>
        number of grid points per dimension to use

    c6_range : [<float>,<float>]
        input range for c6 in nM
    c12_const : <float>
        input constant for c12 in nM
    '''

    c6 = logspace(*tuple(log10(c6_range)),num=N)
    c12 = c12_const*ones(N)
    model = fromcrn(crn_path)

    diffusives = vstack([c6,c12]).T
    c_grid = model.get_inputs(diffusives).T
    L,T = model.get_steady_state(c_grid)

    return c6,c12,L,T


def generate_figure(c6,c12,L,T,c12_const):
    '''main program figure display'''

    figure(figsize=(10,10))
    plot(c6,10**L,'cyan',marker='.',linestyle='')
    plot(c6,10**T,'gold',marker='.',linestyle='')

    xscale('log')
    yscale('log')

    xlabel(r'Diffusive Signal $C_{6}$ / nM',fontsize=16)
    ylabel(r'Inhibitors $L,T$ / nM',fontsize=16)
    text(1e1,1e2,r'Constant $C_{12}$ = '+str(c12_const)+'nM',fontsize=16)
    show()


def main(crn_path,N=250,c12_const=120,c6_range=[1e-1,1e6]) :
    '''parametrisation of main program'''

    c6,c12,L,T = get_hysteresis(crn_path,N,c6_range,c12_const)
    generate_figure(c6,c12,L,T,c12_const)


# execute main program
if __name__ == '__main__' :

    args = get_args()
    main(**args)
