from get_bifurcations import get_bifurcations
from scipy.interpolate import UnivariateSpline,interp1d
from scipy.signal import convolve2d

from crnpy.parser import fromcrn
from crnpy.colors import cyan,yellow
from matplotlib.pyplot import *

from numpy import *
from scipy.interpolate import interp2d
from json import loads,dumps

from os import system
from os.path import isfile
from re import finditer,sub,search

from matplotlib.pyplot import *
from matplotlib.legend_handler import HandlerTuple
from numpy import *
from numpy.random import uniform
from scipy.interpolate import interp2d

from numpy.random import normal
from scipy.interpolate import UnivariateSpline
from yaml import load


def save_frame(model,j=0) :

    fig, (ax1, ax3) = subplots(1,2, figsize=(14,8))
    fig.suptitle(r'$t$ = {:.3} hours'.format(model.time), fontsize=16, y=0.15, x=0.63)
    ax2 = twinx(ax1)

    p1 = ax1.fill_between(100*model.space,model.cfp,color='#00b0f0',linewidth=3,alpha=0.5)
    p2 = ax1.fill_between(100*model.space,model.yfp,color='#ffc000',linewidth=3,alpha=0.5)

    p3, = ax2.plot(100*model.space,model.c6,color='#000099',linewidth=5)
    p4, = ax2.plot(100*model.space,model.c12,color='#ff6600',linewidth=5)

    ax1.set_xlabel('space, $x$ / cm',fontsize=16)
    ax1.set_ylabel('concentrations / nM',fontsize=16);
    
    ax1.set_ylim(0,); ax2.set_ylim(0,)
    ax1.set_yticks([]); ax2.set_yticks([])
    xlim(0,1.6)
    
    mask = model.cfp > model.yfp
    ax3.plot(model.c6[mask]+model.c12[mask],model.c6[mask]/model.c12[mask],'.',color='#00b0f0',ms=12)
    ax3.plot(model.c6[~mask]+model.c12[~mask],model.c6[~mask]/model.c12[~mask],'.',color='#ffc000',ms=12)

    ax3.set_xlabel(r'density, $C_{6}$+$C_{12}$ / nM',fontsize=16)
    ax3.set_ylabel(r'ratio, $C_{6}$/$C_{12}$',fontsize=16);

    ax3.set_xscale('log'); ax3.set_yscale('log')
    ax3.set_xlim(0.04,75000); ax3.set_ylim(0.04,75) 
    
    fig.legend([(p1, p2), (p3,p4)], [r'cell responses $CFP$,$YFP$', r'morphogens $C_{6}$,$C_{12}$'],
               handler_map={tuple: HandlerTuple(ndivide=None)}, fontsize=16, loc=(0.31,0.08), borderaxespad=0.1)
    subplots_adjust(bottom=0.3)
    savefig(str(j).zfill(4)+'.png')
    close()

    return j+1


def create_animation(C6,C12):
    '''execute system commands to create animation from frames'''
    system('convert -delay 10 -loop 0 *.png c6-{}_c12-{}.gif'.format(C6,C12))
    system('rm *.png')


crn_path = 'models/direct-constant.crn'

model = fromcrn(crn_path)
t_final = 30

C6=100; C12=300
width = model._xmax / 4.0
n_timepoints = 50

# initial condition
model.c6[model.space<width] = model._xmax * C6 / width
model.c12[model.space>=(model._xmax-width)] = model._xmax * C12 / width

j = 0
while model.time < t_final :
    model.time_step() 

    record = int(model.time/model._dt) % int(t_final/model._dt/n_timepoints) == 0
    if record : j = save_frame(model,j)

create_animation(C6,C12)