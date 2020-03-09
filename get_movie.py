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
	parser.add_argument('--C12', type=float, default=400.0, help='equilibrium concentration nM')
	parser.add_argument('--t_final', type=float, default=24.0,help='final time in hours')
	parser.add_argument('--n_timepoints', type=int, default=50,help='number of time points for movie')
	parser.add_argument('--initial', type=str, default='both',help='initial condition = "c6", "c12" or "both"')
	return vars(parser.parse_args())


def save_frame(model,region,j=0) :

    fig, (ax1, ax3) = subplots(1,2, figsize=(14,7))
    C6,C12 = mean(model.c6),mean(model.c12)
    fig.text(0.15,0.15,r'$t$ = {:.3} hours'.format(model.time), fontsize=16)
    fig.text(0.57,0.15,r'$\overline{C6}$'+' = {:.0f} nM, '.format(C6)+r'$\overline{C12}$'+' = {:.0f} nM'.format(C12), fontsize=16)
    ax2 = twinx(ax1)

    p1 = ax1.fill_between(100*model.space,log10(model.cfp),color='#00b0f0',linewidth=0,alpha=0.5)
    p2 = ax1.fill_between(100*model.space,log10(model.yfp),color='#ffc000',linewidth=0,alpha=0.5)

    p3, = ax2.plot(100*model.space,log10(model.c6),color='#000099',linewidth=5)
    p4, = ax2.plot(100*model.space,log10(model.c12),color='#ff6600',linewidth=5)

    ax1.set_xlabel('space, $x$ / cm',fontsize=16)
    ax1.set_ylabel('concentrations',fontsize=16)

    ax1.set_ylim(3,5); ax2.set_ylim(-1,5); ax1.set_ylim(0,)
    ax1.set_yticks([]); ax2.set_yticks([]); ax1.set_xlim(0,100*amax(model.space))
    xlim(0,)

    mask = model.cfp > model.yfp
    ax3.plot(model.c12[mask],model.c6[mask],'.',color='#00b0f0',ms=12)
    ax3.plot(model.c12[~mask],model.c6[~mask],'.',color='#ffc000',ms=12)

    region_c12,region_c6 = region
    p7, = ax3.plot(C12,C6,'x',color='black')
    ax3.plot(10**region_c12,10**region_c6,color='black',linewidth=3)

    try :
        p5 = ax3.quiver( model.c12, model.c6,
            model.capacity*model.P76*model.kC12/model.dlasI,
            model.capacity*model.P81*model.kC6/model.dluxI, color='red' )

        if model.kC12 == 0.0 and model.kC6 == 0 :
            fig.legend([p7], [r'morphogen spatial average'],
                handler_map={tuple: HandlerTuple(ndivide=None)}, fontsize=16, loc=(0.55,0.81), borderaxespad=0.1)
        else :
            fig.legend([p5, p7], [r'relay reaction', r'morphogen spatial average'],
                handler_map={tuple: HandlerTuple(ndivide=None)}, fontsize=16, loc=(0.55,0.76), borderaxespad=0.1)
    except :
        p5 = ax3.quiver( model.c12, model.c6,
            model.capacity*model.P76*0.0,
            model.capacity*model.P81*0.0, color='red' )

        fig.legend([p7], [r'morphogen spatial average'],
            handler_map={tuple: HandlerTuple(ndivide=None)}, fontsize=16, loc=(0.55,0.76), borderaxespad=0.1)

    ax3.set_xlabel(r'morphogen, C12 / nM',fontsize=16)
    ax3.set_ylabel(r'morphogen, C6 / nM',fontsize=16)

    ax3.set_xscale('log'); ax3.set_yscale('log')
    ax3.set_xlim(1.6,2.5e4); ax3.set_ylim(1.6,2.5e4)

    fig.legend([(p1, p2), (p3,p4) ], [r'cell responses $CFP$,$YFP$', r'morphogens C6, C12' ],
               handler_map={tuple: HandlerTuple(ndivide=None)}, fontsize=16, loc=(0.13,0.76), borderaxespad=0.1)
    savefig(str(j).zfill(4)+'.png')
    close()

    return j+1


def create_animation(model):
    '''execute system commands to create animation from frames'''
    C6,C12 = mean(model.c6),mean(model.c12)
    system('convert -delay 10 -loop 0 *.png c6-{:.0f}_c12-{:.0f}.gif'.format(C6,C12))
    system('rm *.png')


def main(crn_path,C6,C12,n_timepoints,t_final,initial) :
    region = load('./bifurcations.npz').T
    model = fromcrn(crn_path)
    width = 0.004#model._xmax / 4.0

    # initial condition
    model.c6[:] = 0.0
    model.c12[:] = 0.0

    if initial == 'both' :
        model.c6[model.space<width] = model._xmax * C6 / width
        model.c12[model.space>=(model._xmax-width)] = model._xmax * C12 / width

    elif initial == 'c6' :
        model.c6[abs(model.space-model._xmax/2)<width/2] = 4000

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
