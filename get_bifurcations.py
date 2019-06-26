from crnpy.parser import fromcrn
from crnpy.colors import cyan,yellow

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
    parser.add_argument('--c6_range', nargs=2, type=float, default=[-0.5,5],
                        help='log10 input range for c6 such that 10**input is in nM',metavar=('min','max'))
    parser.add_argument('--c12_range', nargs=2, type=float, default=[-0.5,5],
                        help='log10 input range for c12 such that 10**input is in nM',metavar=('min','max'))
    parser.add_argument('--atc', type=float, default=0.0,
                        help='level of atc in nM',metavar='value')
    parser.add_argument('--iptg', type=float, default=0.0,
                        help='level of iptg in nM',metavar='value')
    parser.add_argument('--clip', type=float,default=-0.5,metavar='value',
                        help='threshold concentration 10**clip below which system is off')
    parser.add_argument('--eps', type=float,default=1e-3,metavar='value',
                        help='precision to use when computing roots')
    return vars(parser.parse_args())


def get_bifurcations(crn_path,N=50,c6_range=[-0.5,5],c12_range=[-0.5,5],atc=0.0,iptg=0.0,clip=-0.5,eps=1e-3):
    '''Calculate bifrucation diagram for double exclusive reporter
    for a given range of diffusives c6 and c12.

    ---parameters---
    crn_path : <str>
        path to crn file
    N : <int>
        number of grid points per dimension to use

    c6_range : [<float>,<float>]
        log10 input range for c6 such that 10**input is in nM
    c12_range : [<float>,<float>]
        log10 input range for c12 such that 10**input is in nM

    atc : <float>
        level of atc in nM
    iptg : <float>
        level of iptg in nM

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

    model.ATC = atc
    model.IPTG = iptg

    # calculation of steady states
    steady_state = model.get_steady_state(c_grid,clip=clip,logspace=True)

    cfp = steady_state[:,:,model.nontrivials.index('lacI')]
    yfp = steady_state[:,:,model.nontrivials.index('tetR')]

    return c6,c12,cfp,yfp


def generate_figure(c6,c12,cfp,yfp,atc,iptg):
    '''main program figure display'''

    figure(figsize=(10,10))
    title('ATC = {} nM     IPTG = {} nM'.format(atc,iptg),fontsize=16,y=1.02)

    contourf(c6,c12,cfp,cmap='cyan',alpha=0.5)
    contourf(c6,c12,yfp,cmap='yellow',alpha=0.5)
    contour(c6,c12,cfp-yfp,levels=[0.0],colors=['k'],alpha=0.5)

    xlabel(r'Diffusive Signal $C_{6}$ / nM',fontsize=16)
    ylabel(r'Diffusive Signal $C_{12}$ / nM',fontsize=16)

    xscale('log'); yscale('log')
    show()


def main(crn_path,N=50,c6_range=[-0.5,5],c12_range=[-0.5,5],atc=0.0,iptg=0.0,clip=-0.5,eps=1e-3) :
    '''parametrisation of main program'''

    print('Calculating steady states...')
    c6,c12,cfp,yfp = get_bifurcations(crn_path,N,c6_range,c12_range,atc,iptg,clip,eps)
    print('Done')

    generate_figure(c6,c12,cfp,yfp,atc,iptg)


# execute main program
if __name__ == '__main__' :
    args = get_args()
    main(**args)
