from sys import argv
from argparse import ArgumentParser
from lib.utils import str2bool

from lib.model import Model
from lib.parsers import crn_parameters
from models.doubleExclusive import system_specifications,parameters

from pandas import read_csv
from lib.colors import cyan,yellow

from numpy import meshgrid,log10
from matplotlib.pyplot import plot,scatter,figure,xlim,ylim,xlabel,ylabel,xscale,yscale,colorbar,pcolor,fill_between,show,legend

def get_args() :
	'''parse arguments from command line'''
	parser = ArgumentParser(description='creates bifrucation plot from crn file parameters')

	parser.add_argument('crn_path', type=str, help='path to crn file')
	parser.add_argument('--model', type=str2bool, default=True, help='model predictions')
	parser.add_argument('--data_path', type=str, default='./data/liquid/char_ExRep_1_R33S175ExRepTet33AAVLac300ND.csv',help='liquid culture dataset')
	return vars(parser.parse_args())


def get_bifurcations(crn_path,model,data_path):
	'''Calculate bifrucation diagram for double exclusive reporter
	for a given range of diffusives c6 and c12'''

	if model :

		parameters.update(crn_parameters(crn_path))
		parameters['c6'],parameters['c12'] = 4.88,-1.39
		model = Model(pars = parameters , **system_specifications)

		model.integrate()
		model.get_cusp('c12','c6')

		# parameters.update(crn_parameters(crn_path))
		# parameters['c6'],parameters['c12'] = -1.39,4.88
		# model = Model(pars = parameters , **system_specifications)

		# model.integrate()
		# model.get_cusp('c6','c12')
	
	liquid_data = read_csv(data_path)
	cfp,yfp = liquid_data.pivot('C6','C12','P(ECFP/mRFP1)'), liquid_data.pivot('C6','C12','P(EYFP/mRFP1)')
	c12,c6 = meshgrid(cfp.columns.values,cfp.index.values)

	return model,c6,c12,cfp,yfp


def generate_figure(model,c6,c12,cfp,yfp):
	'''main program figure display'''

	figure(figsize=(9,7))

	c12shift,c6shift = 2*cfp.columns.values, 2*cfp.index.values
	pcolor(c12+c12shift,c6+c6shift[:,None],log10(cfp.values),cmap='cyan',vmin=1,vmax=2.1)
	colorbar(pad=-0.1,ticks=[1,2,3]).ax.set_yticklabels(['$10^{1}$','$10^{2}$','$10^{3}$'])
	
	pcolor(c12+c12shift,c6+c6shift[:,None],log10(yfp.values),cmap='yellow',vmin=1,vmax=2.1)
	colorbar(ticks=[]).ax.set_yticklabels([])
	
	if model :
		region = model.bifurcations['LC1']
		region_c6,region_c12 = region.curve[:-1,region.params].T
		fill_between(10**region_c12,10**region_c6,facecolor='white',hatch='///', edgecolor='k',linewidth=0)
	
	# flow cytometry points
	scatter(1e-1,100,s=200,facecolor='#00b0f0',marker='s',edgecolors='k')
	scatter(1e-1,1e-1,s=200,facecolor='white',marker='s',edgecolors='k')
	scatter(200,1e-1,s=200,facecolor='#ffc000',marker='s',edgecolors='k')

	scatter(25000,8333,s=200,facecolor='#00b0f0',marker='s',edgecolors='k')
	plot(25000,2777,linewidth=0,fillstyle='left',marker='s',markersize=15,markeredgewidth=1,markeredgecolor='k',markerfacecolor='#ffc000',markerfacecoloralt='#00b0f0')
	plot(25000,926,linewidth=0,fillstyle='left',marker='s',markersize=15,markeredgewidth=1,markeredgecolor='k',markerfacecolor='#ffc000',markerfacecoloralt='#00b0f0')
	scatter(25000,309,s=200,facecolor='#ffc000',marker='s',edgecolors='k')

	scatter(5000,2777,s=200,facecolor='#00b0f0',marker='s',edgecolors='k')
	plot(5000,926,linewidth=0,fillstyle='left',marker='s',markersize=15,markeredgewidth=1,markeredgecolor='k',markerfacecolor='#ffc000',markerfacecoloralt='#00b0f0')
	plot(5000,308,linewidth=0,fillstyle='left',marker='s',markersize=15,markeredgewidth=1,markeredgecolor='k',markerfacecolor='#ffc000',markerfacecoloralt='#00b0f0')
	plot(5000,102,linewidth=0,fillstyle='left',marker='s',markersize=15,markeredgewidth=1,markeredgecolor='k',markerfacecolor='#ffc000',markerfacecoloralt='#00b0f0')
	scatter(5000,34,s=200,facecolor='#ffc000',marker='s',edgecolors='k')

	scatter(1000,926,s=200,facecolor='#00b0f0',marker='s',edgecolors='k')
	plot(1000,308,linewidth=0,fillstyle='left',marker='s',markersize=15,markeredgewidth=1,markeredgecolor='k',markerfacecolor='#ffc000',markerfacecoloralt='#00b0f0')
	plot(1000,102,linewidth=0,fillstyle='left',marker='s',markersize=15,markeredgewidth=1,markeredgecolor='k',markerfacecolor='#ffc000',markerfacecoloralt='#00b0f0')
	plot(1000,34,linewidth=0,fillstyle='left',marker='s',markersize=15,markeredgewidth=1,markeredgecolor='k',markerfacecolor='#ffc000',markerfacecoloralt='#00b0f0')
	scatter(1000,11,s=200,facecolor='#ffc000',marker='s',edgecolors='k')

	scatter(200,308,s=200,facecolor='#00b0f0',marker='s',edgecolors='k')
	plot(200,100,linewidth=0,fillstyle='left',marker='s',markersize=15,markeredgewidth=1,markeredgecolor='k',markerfacecolor='#ffc000',markerfacecoloralt='#00b0f0')
	plot(200,34,linewidth=0,fillstyle='left',marker='s',markersize=15,markeredgewidth=1,markeredgecolor='k',markerfacecolor='#ffc000',markerfacecoloralt='#00b0f0')
	scatter(200,11,s=200,facecolor='#ffc000',marker='s',edgecolors='k')

	scatter(40,102,s=200,facecolor='#00b0f0',marker='s',edgecolors='k')
	plot(40,34,linewidth=0,fillstyle='left',marker='s',markersize=15,markeredgewidth=1,markeredgecolor='k',markerfacecolor='#ffc000',markerfacecoloralt='#00b0f0')
	plot(40,11,linewidth=0,fillstyle='left',marker='s',markersize=15,markeredgewidth=1,markeredgecolor='k',markerfacecolor='#ffc000',markerfacecoloralt='#00b0f0')
	scatter(40,3.8,s=200,facecolor='#ffc000',marker='s',edgecolors='k')

	scatter(8,34,s=200,facecolor='#00b0f0',marker='s',edgecolors='k')
	scatter(8,11,s=200,facecolor='#ffc000',marker='s',edgecolors='k')

	yscale('log'); xscale('log')
	xlim(0.04,75000); ylim(0.04,75000)

	xlabel('Morphogen $C_{12}$ / nM',fontsize=16)
	ylabel('Morphogen $C_{6}$ / nM',fontsize=16)

	plot([None],[None],linewidth=0,fillstyle='left',marker='s',markersize=15,markeredgewidth=1,markeredgecolor='k',markerfacecolor='#ffc000',markerfacecoloralt='#00b0f0',label='flow cytometry')
	plot([None],[None],linewidth=0,fillstyle='left',marker='s',markersize=15,markeredgewidth=0,markerfacecolor='#ffc000',markerfacecoloralt='#00b0f0',label='plate reader')
	scatter(None,None,s=200,facecolor='white',marker='s',edgecolors='none',hatch='///',label='bistable prediction')
	legend(fontsize=16,loc=2)
	show()


def main(crn_path,model,data_path) :
	'''parametrisation of main program'''

	model,c6,c12,cfp,yfp = get_bifurcations(crn_path,model,data_path)
	generate_figure(model,c6,c12,cfp,yfp)


# execute main program
if __name__ == '__main__' :
	args = get_args()
	main(**args)