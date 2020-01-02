from crnpy.parser import fromcrn
from crnpy.colors import cyan,yellow

from os import system
from os.path import isfile

from matplotlib.pyplot import *
from matplotlib.legend_handler import HandlerTuple
from numpy import *

from sys import argv
from argparse import ArgumentParser


def get_args() :
	'''parse arguments from command line'''
	parser = ArgumentParser(description='creates bifrucation plot from crn file parameters')

	parser.add_argument('crn_path', type=str, help='path to crn file')
	parser.add_argument('--C6', type=float, default=100.0, help='equilibrium concentration nM')
	parser.add_argument('--C12', type=float, default=300.0, help='equilibrium concentration nM')
	parser.add_argument('--t_final', type=float, default=24.0,help='final time in hours')
	parser.add_argument('--n_timepoints', type=int, default=50,help='number of time points for movie')
	parser.add_argument('--initial', type=str, default='both',help='initial condition = "c6", "c12" or "both"')
	return vars(parser.parse_args())


def save_frame(model,region,j=0) :

    fig, (ax1, ax3) = subplots(1,2, figsize=(14,9))
    C6,C12 = mean(model.c6),mean(model.c12)
    fig.suptitle(r'$t$ = {:.3} hours, $C6$ = {:.0f} nM, $C12$ = {:.0f} nM'.format(model.time,C6,C12), fontsize=16, y=0.15, x=0.69)
    ax2 = twinx(ax1)

    p1 = ax1.fill_between(100*model.space,model.cfp,color='#00b0f0',linewidth=0,alpha=0.5)
    p2 = ax1.fill_between(100*model.space,model.yfp,color='#ffc000',linewidth=0,alpha=0.5)

    p3, = ax2.plot(100*model.space,model.c6,color='#000099',linewidth=5)
    p4, = ax2.plot(100*model.space,model.c12,color='#ff6600',linewidth=5)

    # cgrid = vstack([model.c6,model.c12]).T
    # local_equilibria = model.get_local_equilibria(cgrid)
    # local_cfp = local_equilibria[:,model.nontrivials.index('cfp')]
    # local_yfp = local_equilibria[:,model.nontrivials.index('yfp')]

    # p5, = ax1.plot(100*model.space,local_cfp,'o',color='#00b0f0')
    # p6, = ax1.plot(100*model.space,local_yfp,'o',color='#ffc000')

    ax1.set_xlabel('space, $x$ / cm',fontsize=16)
    ax1.set_ylabel('concentrations / nM',fontsize=16)
    
    ax1.set_ylim(0,); ax2.set_ylim(0,); ax1.set_ylim(0,)
    ax1.set_yticks([]); ax2.set_yticks([])
    xlim(0,)
    
    mask = model.cfp > model.yfp
    ax3.plot(model.c12[mask],model.c6[mask],'.',color='gray',ms=12)
    ax3.plot(model.c12[~mask],model.c6[~mask],'.',color='gray',ms=12)

    region_c6,region_c12 = region
    p7, = ax3.plot(C12,C6,'x',color='black')
    ax3.plot(10**region_c12,10**region_c6,color='black',linewidth=3)

    p5 = ax3.quiver( model.c12, model.c6,
        model.capacity*model.P76*model.kC12/model.dlasI,
        model.capacity*model.P81*model.kC6/model.dluxI, color='red' )

    ax3.set_xlabel(r'morphogen, $C_{12}$ / nM',fontsize=16)
    ax3.set_ylabel(r'morphogen, $C_{6}$ / nM',fontsize=16)

    ax3.text(10**2.6,10**2.2,'bistable region',fontsize=16)
    ax3.set_xscale('log'); ax3.set_yscale('log')
    ax3.set_xlim(0.04,75000); ax3.set_ylim(0.04,75000) 
    
    fig.legend([(p1, p2), (p5), (p3,p4), p7], [r'cell responses $CFP$,$YFP$', r'relay reaction', r'morphogens $C_{6}$,$C_{12}$', r'homogenous equilibrium'],
               handler_map={tuple: HandlerTuple(ndivide=None)}, fontsize=16, loc=(0.18,0.06), borderaxespad=0.1)
    subplots_adjust(bottom=0.35)
    savefig(str(j).zfill(4)+'.png')
    close()

    return j+1


def create_animation(model):
    '''execute system commands to create animation from frames'''
    C6,C12 = mean(model.c6),mean(model.c12)
    system('convert -delay 10 -loop 0 *.png c6-{:.0f}_c12-{:.0f}.gif'.format(C6,C12))
    system('rm *.png')


def main(crn_path,C6,C12,n_timepoints,t_final,initial) :
    region = load('./bifurcations.npy')
    model = fromcrn(crn_path)
    width = model._xmax / 4.0

    # initial condition
    model.c6[:] =  1e-1
    model.c12[:] = 1e-1

    if initial == 'both' :
        model.c6[model.space<width] = model._xmax * C6 / width
        model.c12[model.space>=(model._xmax-width)] = model._xmax * C12 / width

    elif initial == 'c6' :
        model.c6[model.space<width] = model._xmax * C6 / width

    elif initial == 'c12' :
        model.c12[model.space>=(model._xmax-width)] = model._xmax * C12 / width
    
    else :
        raise Exception('invalid initial condition. initial = "c6", "c12" or "both"')

    j = 0
    while model.time < t_final :

        record = int(model.time/model._dt) % int(t_final/model._dt/n_timepoints) == 0
        if record & (model.time > 0.16) : j = save_frame(model,region,j)
        model.time_step()

    create_animation(model)

# execute main program
if __name__ == '__main__' :
	args = get_args()
	main(**args)