from lib.model import fromcrn
from lib.colors import cyan,yellow

from sys import argv
from argparse import ArgumentParser

from numpy import linspace,logspace,meshgrid,vstack,dstack,inf,mean,ones,array
from matplotlib.pyplot import *

def get_args() :
    '''parse arguments from command line'''

    parser = ArgumentParser(description='creates bifrucation plot from crn file parameters')
    parser.add_argument('--crn_path', type=str,
                        help='path to crn file')
    parser.add_argument('--N', type=int, default=50,
                        help='number of grid points per dimension to use')
    parser.add_argument('--c_range', nargs=2, type=float, default=[-1,1],
                        help='input range of non-dimensionalised c6')
    parser.add_argument('--cdash_range', nargs=2, type=float, default=[0.5,1],
                        help='input range of non-dimensionalised c12')
    parser.add_argument('--clip', type=float,default=-0.5,
                        help='threshold concentration x**clip below which system is off')
    parser.add_argument('--eps', type=float,default=1e-3,
                        help='precision to use when computing roots')
    return vars(parser.parse_args())


def main(crn_path,N=50,c_range=[-1,1],cdash_range=[0.5,1],clip=-0.5,eps=1e-3) :
    '''parametrisation of main program'''

    eps = array([eps,-eps])
    c_range = tuple(array(c_range)+eps)
    cdash_range = tuple(array(cdash_range)+eps)

    # initialisation of model
    model = fromcrn(crn_path)
    c = linspace(*c_range,num=N)
    cdash = linspace(*cdash_range,num=N)
    x,y = meshgrid(c,cdash,copy=False)
    c_grid = dstack([x,y])

    # calculation of steady states
    L,T = model.get_steady_state(c_grid,clip=clip)
    c6,c12 = model.get_diffusives(c_grid).reshape(2,N,N)

    # create bifrucation figure
    figure(figsize=(10,10))

    contour(c6,c12,T[:,:,-1]-L[:,:,0],levels=[0],colors=['k'])
    contour(c6,c12,L[:,:,-1]-T[:,:,0],levels=[0],colors=['k'])
    contourf(c6,c12,L[:,:,0],cmap='cyan',alpha=0.25)
    contourf(c6,c12,L[:,:,-1],cmap='cyan',alpha=0.25)
    contourf(c6,c12,T[:,:,0],cmap='yellow',alpha=0.25)
    contourf(c6,c12,T[:,:,-1],cmap='yellow',alpha=0.25)

    xlabel(r'Diffusive Signal $C_{6}$ / nM',fontsize=16)
    ylabel(r'Diffusive Signal $C_{12}$ / nM',fontsize=16)
    xscale('log')
    yscale('log')

    # labelling regions
    text(1e-5, 4,r'$Off$ State',fontsize=16)
    text(120, 3000,r'$Bistable$ Region',fontsize=16)
    text(1e-5, 100,r'$Monostable$ Region',fontsize=16)
    text(1e3, 10,r'$Monostable$ Region',fontsize=16)
    show()


# execute main program
if __name__ == '__main__' :

    args = get_args()
    main(**args)
