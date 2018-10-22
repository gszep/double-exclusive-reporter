from lib.model import fromcrn
from lib.colors import cyan,yellow

from sys import argv
from argparse import ArgumentParser

from numpy import linspace,logspace,meshgrid,vstack,dstack,inf,mean,ones,array
from matplotlib.pyplot import *

def get_args() :
    '''parse arguments from command line'''

    parser = ArgumentParser(description='creates hysteresis plot from crn file parameters')
    parser.add_argument('--crn_path', type=str,
                        help='path to crn file')
    parser.add_argument('--c12_const', type=float, default=120,
                        help='constant value of c12')
    parser.add_argument('--N', type=int, default=250,
                        help='number of grid points per dimension to use')

    return vars(parser.parse_args())


def main(crn_path,N=250,c12_const=120) :
    '''parametrisation of main program'''

    c6 = logspace(-1,6,N)
    c12 = c12_const*ones(N)
    model = fromcrn(crn_path)

    diffusives = vstack([c6,c12]).T
    c_grid = model.get_inputs(diffusives).T
    L,T = model.get_steady_state(c_grid)

    figure(figsize=(10,10))
    plot(c6,10**L,'cyan',marker='.',linestyle='')
    plot(c6,10**T,'gold',marker='.',linestyle='')

    xscale('log')
    yscale('log')

    xlabel(r'Diffusive Signal $C_{6}$ / nM',fontsize=16)
    ylabel(r'Inhibitors $L,T$ / nM',fontsize=16)
    text(1e1,1e2,r'Constant $C_{12}$ = '+str(c12_const)+'nM',fontsize=16)
    show()


# execute main program
if __name__ == '__main__' :

    args = get_args()
    main(**args)
