from sys import argv
from argparse import ArgumentParser

from lib.model import Model
from lib.parsers import crn_parameters
from models.doubleExclusive import system_specifications,parameters

from pandas import read_csv
from lib.colors import cyan,yellow

from numpy import meshgrid,log10
from matplotlib.pyplot import plot,figure,xlim,ylim,xlabel,ylabel,xscale,yscale,colorbar,pcolor,fill_between,show

def get_args() :
	'''parse arguments from command line'''
	parser = ArgumentParser(description='creates bifrucation plot from crn file parameters')

	parser.add_argument('crn_path', type=str, help='path to crn file')
	parser.add_argument('--model', type=bool, default=True, help='model predictions')
	parser.add_argument('--data_path', type=str, default='./data/liquid/char_ExRep_1_R33S175ExRepTet33AAVLac300ND.csv',help='liquid culture dataset')
	return vars(parser.parse_args())


def get_bifurcations(crn_path,model,data_path):
	'''Calculate bifrucation diagram for double exclusive reporter
	for a given range of diffusives c6 and c12'''

	if model :

		parameters.update(crn_parameters(crn_path))
		model = Model(pars = parameters , **system_specifications)

		model.integrate()
		model.get_cusp('c6','c12')
	
	liquid_data = read_csv(data_path)
	cfp,yfp = liquid_data.pivot('C6','C12','P(ECFP/mRFP1)'), liquid_data.pivot('C6','C12','P(EYFP/mRFP1)')
	c12,c6 = meshgrid(cfp.columns.values,cfp.index.values)

	return model,c6,c12,cfp,yfp


def generate_figure(model,c6,c12,cfp,yfp):
	'''main program figure display'''

	figure(figsize=(9,7))
	
	pcolor(c12,c6,log10(cfp.values),cmap='cyan',vmin=1,vmax=2.1)
	colorbar(pad=-0.1,ticks=[1,2,3]).ax.set_yticklabels(['$10^{1}$','$10^{2}$','$10^{3}$'])
	
	pcolor(c12,c6,log10(yfp.values),cmap='yellow',vmin=1,vmax=2.1)
	colorbar(ticks=[]).ax.set_yticklabels([])
	
	if model :
		region = model.bifurcations['LC1']
		region_c6,region_c12 = region.curve[:-1,region.params].T
		fill_between(10**region_c12,10**region_c6,facecolor="none",hatch="///", edgecolor="k",linewidth=0)
	
	yscale('log'); xscale('log')
	xlim(0.04,2500); ylim(0.04,2500)
	
	plot(200,100,'kx',ms=10,mew=5)
	plot(1e-1,100,'kx',ms=10,mew=5)
	plot(200,1e-1,'kx',ms=10,mew=5)
	plot(1e-1,1e-1,'kx',ms=10,mew=5)
	
	xlabel('Morphogen $C_{12}$ / nM',fontsize=16)
	ylabel('Morphogen $C_{6}$ / nM',fontsize=16)
	show()


def main(crn_path,model,data_path) :
	'''parametrisation of main program'''

	model,c6,c12,cfp,yfp = get_bifurcations(crn_path,model,data_path)
	generate_figure(model,c6,c12,cfp,yfp)


# execute main program
if __name__ == '__main__' :
	args = get_args()
	main(**args)