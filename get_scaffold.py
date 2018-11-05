from lib.model import fromcrn
from get_bifurcations import get_bifurcations
from lib.colors import cyan,yellow

from sys import argv
from argparse import ArgumentParser
from os import system

from numpy import linspace,logspace,meshgrid,vstack,dstack,inf,mean,ones,array,argmin
from matplotlib.pyplot import *

def get_args() :
    '''parse arguments from command line'''
    parser = ArgumentParser(description='simulates scaffold from crn file parameters for a given amount of time')

    parser.add_argument('--crn_path', type=str,
                        help='path to crn file')
    parser.add_argument('--N', type=int, default=150,
                        help='number of grid points per dimension to use')
    parser.add_argument('--c6_range', nargs=2, type=float, default=[1e-6,1e8],
                        help='input range for c6 in nM')
    parser.add_argument('--c12_range', nargs=2, type=float, default=[10**-0.5,1e5],
                        help='input range for c12 in nM')
    parser.add_argument('--clip', type=float,default=-0.5,
                        help='threshold concentration 10**clip below which system is off')
    parser.add_argument('--eps', type=float,default=1e-3,
                        help='precision to use when computing roots')
    return vars(parser.parse_args())


def get_scaffold(model,c6,c12,L,T):
    '''Simulate concentrations and scaffold for double exclusive reporter
    in a one dimensional space for a given initialisation of c6 and c12.

    ---parameters---
    model : <str>
        path to crn file

    c6 : [<float>,<float>]
        input range for c6 in nM
    c12 : [<float>,<float>]
        input range for c12 in nM

    L : <float>
        threshold concentration 10**clip below which system is off
    T : <float>
        precision to use when computing roots
    '''


    Lattractor = array([
        L[argmin(abs(c12[:,0]-c12x)),argmin(abs(c6[0]-c6x))]
        for c6x,c12x in model.state['diffusables'] ]).reshape(-1)

    Tattractor = array([
        T[argmin(abs(c12[:,0]-c12x)),argmin(abs(c6[0]-c6x))]
        for c6x,c12x in model.state['diffusables'] ]).reshape(-1)

    space = linspace(0,model.xmax,len(Lattractor))
    return space,Lattractor,Tattractor


def generate_frame(j,model,space,L,T):
    '''main program figure display'''

    figure(figsize=(10,10))

    # plot scaffolds
    plot(space,10**L,'.',color='darkcyan')
    plot(space,10**T,'.',color='gold')
    plot(-1,-1,'k.',label='Local Steady State')

    # plot system state
    plot(model.space,model.state['inhibitors'].T[0],color='darkcyan')
    plot(model.space,model.state['inhibitors'].T[1],color='gold')
    plot(-1,-1,'k',label='Concentration at $t = {}h$'.format(model.time[-1]))

    legend(fontsize=16)
    xlim(0,model.xmax)

    ylabel(r'Inhibitors $L(x,t),T(x,t)$ / nM',fontsize=16);
    xlabel(r'Space $x$ / cm',fontsize=16);
    yscale('log')

    savefig(str(j).zfill(4)+'.png')
    close()

    return j+1


def create_animation(C6,C12):
    '''execute system commands to create animation from frames'''
    system('convert -delay 10 -loop 0 *.png c6-{}_c12-{}.gif'.format(C6,C12))
    system('rm *.png')


def main(crn_path,C6=120,C12=120,width=0.01,t_final=750.0,N=150,c6_range=[1e-6,1e8],c12_range=[10**-0.5,1e5],clip=-0.5,eps=1e-3) :
    '''parametrisation of main program'''

    # import model from file
    model = fromcrn(crn_path)

    # get scaffold from bifrucation diagram
    c6,c12,L,T = get_bifurcations(crn_path,N,c6_range,c12_range,clip,eps)

    # initial condition
    model.state['diffusables'][:,0][model.space<width] = model.xmax * C6 / width
    model.state['diffusables'][:,1][model.space>(model.xmax-width)] = model.xmax * C12 / width

    # evolve system
    i,j = 0,0

    n = int(t_final/model.dt)
    while model.time[-1] < t_final :

        model.time_step()
        space,Lattractor,Tattractor = get_scaffold(model,c6,c12,L,T)

        if i % int(n/300) == 0 :
            j = generate_frame(j,model,space,Lattractor,Tattractor)
        i += 1

    create_animation(C6,C12)


# execute main program
if __name__ == '__main__' :
    args = get_args()
    main(**args)
