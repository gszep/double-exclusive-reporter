from sys import argv
from argparse import ArgumentParser
from lib.utils import str2bool

from lib.model import Model
from lib.parsers import crn_parameters
import models.doubleExclusive as model1
import models.doubleExclusive_v2 as model2

from pandas import read_csv
from lib.colors import cyan,yellow

from numpy import meshgrid,log10,save
from matplotlib.pyplot import plot,scatter,figure,xlim,ylim,xlabel,ylabel,xscale,yscale,colorbar,pcolor,fill_between,show,legend

def get_args() :
	'''parse arguments from command line'''
	parser = ArgumentParser(description='creates bifrucation plot from crn file parameters')

	parser.add_argument('crn_path', type=str, help='path to crn file')
	parser.add_argument('--predictions', type=str2bool, default=True, help='model predictions')
	parser.add_argument('--version', type=int, default=1, help='model version')
	parser.add_argument('--data_path', type=str, default='./data/liquid/char_ExRep_1_R33S175ExRepTet33AAVLac300ND.csv',help='liquid culture dataset')
	parser.add_argument('--save_path', type=str, default='',help='save bifrucation only')
	return vars(parser.parse_args())

def define_model(version, crn_path=None):
	print('Generating bifurcation diagram for version %d'%version)
	if version == 1:
		parameters = model1.parameters
		system_specifications = model1.system_specifications
	elif version == 2:
		parameters = model2.parameters
		system_specifications = model2.system_specifications
	else:
		raise Exception('Unknown model version')

	if crn_path is not None:
		parameters.update(crn_parameters(crn_path))
	
	return parameters, system_specifications

def get_bifurcations(parameters, system_specifications, **kwargs):
	'''Calculate bifurcation diagram for double exclusive reporter
	for a given range of concentrations of c6 and c12'''

	try : 
		print('######### search along c6 axis #########')
		model = Model(pars = parameters , **system_specifications)
		model.get_cusp(['c6','c12'], **kwargs)
	except : 
		print('######### search along c12 axis #########')
		model = Model(pars = parameters , **system_specifications)
		model.get_cusp(['c12','c6'], **kwargs)

	return model

def load_data(data_path):
	liquid_data = read_csv(data_path)
	cfp,yfp = liquid_data.pivot('C6','C12','P(ECFP/mRFP1)'), liquid_data.pivot('C6','C12','P(EYFP/mRFP1)')
	c12,c6 = meshgrid(cfp.columns.values,cfp.index.values)
	return c6,c12,cfp,yfp

def generate_figure(model, data_path):
	'''main program figure display'''

	figure(figsize=(9,7))
	if data_path != '' :
		c6,c12,cfp,yfp = load_data(data_path)
		c12shift,c6shift = 2*cfp.columns.values, 2*cfp.index.values
		pcolor(c12+c12shift,c6+c6shift[:,None],log10(cfp.values),cmap='cyan',vmin=1,vmax=2.1)
		colorbar(pad=-0.1,ticks=[1,2,3]).ax.set_yticklabels(['$10^{1}$','$10^{2}$','$10^{3}$'])
		
		pcolor(c12+c12shift,c6+c6shift[:,None],log10(yfp.values),cmap='yellow',vmin=1,vmax=2.1)
		colorbar(ticks=[]).ax.set_yticklabels([])
	
	if model :
		region_forward = model.bifurcations['LC1for']
		region_c6, region_c12 = region_forward.curve[:-1,region_forward.params].T
		plot(10**region_c12,10**region_c6,color='black',linewidth=3)
		region_backward = model.bifurcations['LC1back']
		region_c6, region_c12 = region_backward.curve[:-1,region_backward.params].T
		plot(10**region_c12,10**region_c6,color='black',linewidth=3)
	
	if data_path != '' :
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

	if data_path != '' :
		plot([None],[None],linewidth=0,fillstyle='left',marker='s',markersize=15,markeredgewidth=1,markeredgecolor='k',markerfacecolor='#ffc000',markerfacecoloralt='#00b0f0',label='flow cytometry')
		plot([None],[None],linewidth=0,fillstyle='left',marker='s',markersize=15,markeredgewidth=0,markerfacecolor='#ffc000',markerfacecoloralt='#00b0f0',label='plate reader')
	
	plot([None],[None],color='black',linewidth=3, label='saddle-node bifurcation')
	legend(fontsize=16, loc=2)
	show()


def main(crn_path,predictions,version,data_path,save_path) :
	'''parametrisation of main program'''

	parameters, system_specifications = define_model(version, crn_path)
	model = get_bifurcations(parameters, system_specifications, predictions)
	if save_path != '':
		region = model.bifurcations['LC1']
		save(save_path,region.curve[:-1,region.params].T)
	else:
		generate_figure(model,data_path,c6,c12,cfp,yfp)


# execute main program
if __name__ == '__main__' :
	args = get_args()
	main(**args)