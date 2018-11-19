from lib.model import fromcrn
from lib.colors import cyan,yellow

from sys import argv
from argparse import ArgumentParser

from numpy import linspace,logspace,meshgrid,vstack,dstack,inf,mean,ones,array,log10
from matplotlib.pyplot import *

def get_args() :
    '''parse arguments from command line'''
    parser = ArgumentParser(description='creates bifrucation plot from crn file parameters')

    parser.add_argument('crn_path', type=str,
                        help='path to crn file')
    parser.add_argument('--N', type=int, default=50, metavar='gridpoints',
                        help='number of grid points per dimension to use')
    parser.add_argument('--c6_range', nargs=2, type=float, default=[-6,8],
                        help='log10 input range for c6 such that 10**input is in nM',metavar=('min','max'))
    parser.add_argument('--c12_range', nargs=2, type=float, default=[-0.5,5],
                        help='log10 input range for c12 such that 10**input is in nM',metavar=('min','max'))
    parser.add_argument('--clip', type=float,default=-0.5,metavar='value',
                        help='threshold concentration 10**clip below which system is off')
    parser.add_argument('--eps', type=float,default=1e-3,metavar='value',
                        help='precision to use when computing roots')
    return vars(parser.parse_args())


def get_bifurcations(crn_path,N,c6_range,c12_range,clip,eps):
    '''Calculate bifrucation diagram for double exclusive reporter
    for a given range of diffusives c6 and c12.

    ---parameters---
    crn_path : <str>
        path to crn file
    N : <int>
        number of grid points per dimension to use

    c6_range : [<float>,<float>]
        input range for c6 in nM
    c12_range : [<float>,<float>]
        input range for c12 in nM

    clip : <float>
        threshold concentration 10**clip below which system is off
    eps : <float>
        precision to use when computing roots
    '''

    # initialisation of model
    model = fromcrn(crn_path)

    c = logspace(*c6_range,num=N)
    cdash = logspace(*c12_range,num=N)

    c6,c12 = meshgrid(c,cdash,copy=False)
    c_grid = dstack([c6,c12])

    # calculation of steady states
    steady_state = model.get_steady_state(c_grid,clip=clip)

    L = steady_state[:,:,[s.lower() for s in model.nontrivials ].index('laci')]
    T = steady_state[:,:,[s.lower() for s in model.nontrivials ].index('tetr')]

    return c6,c12,L,T


def generate_figure(c6,c12,L,T):
    '''main program figure display'''

    figure(figsize=(10,10))

    contourf(c6,c12,log10(L),cmap='cyan',alpha=0.25)
    contourf(c6,c12,log10(T),cmap='yellow',alpha=0.25)

    xlabel(r'Diffusive Signal $C_{6}$ / nM',fontsize=16)
    ylabel(r'Diffusive Signal $C_{12}$ / nM',fontsize=16)
    xscale('log')
    yscale('log')

    # labelling regions
    text(1e-4, 1.5,r'$Off$ State',fontsize=16)
    text(120, 3000,r'$Bistable$ Region',fontsize=16)
    text(1e-5, 100,r'$Monostable$ Region',fontsize=16)
    text(1e3, 10,r'$Monostable$ Region',fontsize=16)
    show()


def main(crn_path,N=50,c6_range=[1e-6,1e8],c12_range=[10**-0.5,1e5],clip=-0.5,eps=1e-3) :
    '''parametrisation of main program'''

    print('Calculating steady states...')
    c6,c12,L,T = get_bifurcations(crn_path,N,c6_range,c12_range,clip,eps)
    print('Done')

    generate_figure(c6,c12,L,T)


# execute main program
if __name__ == '__main__' :
    args = get_args()
    main(**args)
