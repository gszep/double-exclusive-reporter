from crnpy.parser import fromcrn
from get_bifurcations import get_bifurcations
from crnpy.colors import cyan,yellow

from sys import stdout
from argparse import ArgumentParser

from numpy import linspace,logspace,meshgrid,vstack,dstack,inf,mean,ones,array,argmin,log10
from matplotlib.pyplot import *

from os import system
from pickle import dump,load
from os.path import isfile

def get_args() :
    '''parse arguments from command line'''
    parser = ArgumentParser(description='simulates scaffold from crn file parameters for a given amount of time')

    parser.add_argument('crn_path', type=str,
                        help='path to crn file')
    parser.add_argument('--C6', type=float,default=120,metavar='value',
                        help='equilibrium value of c6 in nM')
    parser.add_argument('--C12', type=float,default=120,metavar='value',
                        help='equilibrium value of c12 in nM')
    parser.add_argument('--width', type=float,default=0.01,metavar='value',
                        help='initial width of c6/c12 regions')
    parser.add_argument('--n_timepoints', type=int, default=30, metavar='timepoints',
                        help='number of time points for animation')

    parser.add_argument('--N', type=int, default=50, metavar='gridpoints',
                        help='number of grid points per dimension to use')
    parser.add_argument('--c6_range', nargs=2, type=float, default=[-6,8],
                        help='input range for c6 in nM',metavar=('min','max'))
    parser.add_argument('--c12_range', nargs=2, type=float, default=[-0.5,5],
                        help='input range for c12 in nM',metavar=('min','max'))
    parser.add_argument('--clip', type=float,default=None,metavar='value',
                        help='threshold concentration 10**clip below which system is off')
    parser.add_argument('--eps', type=float,default=1e-3,metavar='value',
                        help='precision to use when computing roots')
    return vars(parser.parse_args())


def get_scaffold(model,c6,c12,L,T):
    '''Return current scaffold for given state of the model
    and the precalculated bifrucation diagram for L and T

    ---parameters---
    model : <model>
        model object to extract state from

    c6 : <2darray>
        c6/c12 meshgrid in nM
    c12 : <2darray>
        c6/c12 meshgrid in nM

    L : <3darray>
        bifrucation grid
    T : <3darray>
        bifrucation grid
    '''

    Lattractor = array([
        L[argmin(abs(c12[:,0]-c12x)),argmin(abs(c6[0]-c6x))]
        for c6x,c12x in zip(model.c6,model.c12) ]).reshape(-1)

    Tattractor = array([
        T[argmin(abs(c12[:,0]-c12x)),argmin(abs(c6[0]-c6x))]
        for c6x,c12x in zip(model.c6,model.c12) ]).reshape(-1)

    space = linspace(0,model._xmax,len(Lattractor))
    return space,Lattractor,Tattractor


def import_bifurcations(crn_path,N,c6_range,c12_range,clip,eps):
    '''imports stored bifrucation diagram; otherwise calculates it'''

    # get scaffold from bifrucation diagram
    if isfile('bifurcations') : # from file if precalculated
        with open('bifurcations','r') as file :
            c6,c12,L,T = load(file)

    else :
        with open('bifurcations','w') as file : # otherwise calculate

            print('Calculating steady states...')
            c6,c12,L,T,atc,iptg = get_bifurcations(crn_path,N,c6_range,c12_range,clip,eps)
            print('Done')

            dump((c6,c12,L,T),file)

    return c6,c12,L,T


def generate_frame(j,model,space,L,T):
    '''main program figure display'''

    figure(figsize=(10,10))

    # plot scaffolds
    plot(space,L,'.',color='darkcyan')
    plot(space,T,'.',color='gold')
    plot(-1,-1,'k.',label='Local Steady State')

    # plot system state
    n,m = list(model._plots)
    plot(model.space,getattr(model,n),color='darkcyan')
    plot(model.space,getattr(model,m),color='gold')
    plot(-1,-1,'k',label='Concentration at $t = {}h$'.format(model.time))

    legend(fontsize=16)
    xlim(0,model._xmax)

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


def main(crn_path,C6,C12,width,n_timepoints,N,c6_range,c12_range,clip,eps) :
    '''parametrisation of main program'''

    # import model from file
    model = fromcrn(crn_path)
    c6,c12,L,T = import_bifurcations(crn_path,N,c6_range,c12_range,clip,eps)

    # initial condition
    model.c6[model.space<width] = model._xmax * C6 / width
    model.c12[model.space>(model._xmax-width)] = model._xmax * C12 / width

    # evolve system
    i,j = 0,0
    n = int(model._final/model._dt)

    try :

        # run simulation
        while model.time < model._final :
            progress(model.time, model._final, status='Simulating')
            model.time_step()

            # output results at given temporal resolution
            if i % int(n/n_timepoints) == 0 :
                space,Lattractor,Tattractor = get_scaffold(model,c6,c12,L,T)
                j = generate_frame(j,model,space,Lattractor,Tattractor)

            i += 1
        print '\n'

    finally : # create animation on exit
        create_animation(C6,C12)


def progress(count, total, status=''):
    '''commandline output to visualise progress'''

    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
    stdout.flush()


# execute main program
if __name__ == '__main__' :
    args = get_args()
    main(**args)
