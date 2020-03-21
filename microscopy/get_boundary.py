from numpy import *
from matplotlib.pyplot import *
from skimage import io

from colors import cyan,yellow
from sklearn.preprocessing import FunctionTransformer
from sklearn.linear_model import BayesianRidge
from sklearn.pipeline import Pipeline

from glob import glob
from os import system

from sys import argv
from argparse import ArgumentParser

# basis function
def sigmoid(x) :
	return 1.0/(1.0+exp(-x))

# feature vector
def features(x,degree=3) :
	x = squeeze(x)
	fx = [ sigmoid(sigma*(x-mu)) for mu in linspace(0.25,.75,num=degree) for sigma in [-20,20] ]
	return stack(fx,axis=-1)


def fit(data,threshold = 25000,t_final=24,x_width=1.0) :

	# create coordinate meshes
	n_frames,height,width,n_channels = data.shape
	x,t = linspace(0,x_width,num=width),linspace(0,t_final,num=n_frames)

	xx,tt = meshgrid(x,t)
	X,_ = meshgrid(x,linspace(0,1,num=height))
	cfp,yfp,rfp = 0,1,2

	# intensity normalisation
	data /= array([5519.0,4994.0,2296.0])

	# setting values outside grid squares to nan
	mask = data[-1,...,rfp] < 0.5*amax(data[-1,...,rfp])
	mask = stack([[ mask for _ in range(n_frames) ] for _ in range(n_channels) ],axis=-1)
	data[mask] = NaN

	# create regressors
	yfp_predictor = Pipeline([("feature-map", FunctionTransformer(features)), ("regressor", BayesianRidge() )])
	cfp_predictor = Pipeline([("feature-map", FunctionTransformer(features)), ("regressor", BayesianRidge() )])

	# fit both channels for each timepoint
	yfp_predictions,cfp_predictions = [],[]
	for t in range(n_frames):

		input = data[t,...,yfp]
		yfp_predictor.fit(X[~isnan(input)].reshape(-1,1),input[~isnan(input)].reshape(-1,1))
		yfp_predictions += [yfp_predictor.predict(x.reshape(-1,1))]

		input = data[t,...,cfp]
		cfp_predictor.fit(X[~isnan(input)].reshape(-1,1),input[~isnan(input)].reshape(-1,1))
		cfp_predictions += [cfp_predictor.predict(x.reshape(-1,1))]

	predictions = stack([yfp_predictions,cfp_predictions],axis=-1)
	return xx,tt,predictions,yfp_predictor,cfp_predictor


def get_args() :
	'''parse arguments from command line'''
	parser = ArgumentParser(description='creates bifrucation plot from crn file parameters')

	parser.add_argument('data_path', type=str, help='path to crn file')
	return vars(parser.parse_args())


def main(data_path) :

	cfp,yfp,_ = 0,1,2
	data = io.imread(data_path).astype(float); #data = log10(data+0.0001)
	xx,tt,predictions,yfp_steady_state,cfp_steady_state = fit(data)

	figure(figsize=(7,7))
	title('steady state fit', fontsize=16, y=1.02)
	plot(1.6*xx[0],data[-1,...,yfp].T,'.',color='#ffc000',alpha=0.1)
	plot(1.6*xx[0],data[-1,...,cfp].T,'.',color='#00b0f0',alpha=0.1)

	plot(1.6*xx[0],yfp_steady_state.predict(xx[0].reshape(-1,1)),color='#ffc000',linewidth=3)
	plot(1.6*xx[0],cfp_steady_state.predict(xx[0].reshape(-1,1)),color='#00b0f0',linewidth=3)

	xlabel('width, $x$ / cm',fontsize=16)
	ylabel('arbitrary fluoresence units',fontsize=16)
	show()

	figure(figsize=(7,7))

	contourf(1.6*xx,tt,predictions[...,0],cmap='yellow',alpha=1)
	contourf(1.6*xx,tt,predictions[...,1],cmap='cyan')

	difference = predictions[...,0].T-predictions[...,1].T
	difference[tt.T<8] = NaN
	contour(1.6*xx.T,tt.T,difference,levels=[0],colors=['black'],linewidths=[3])

	xlabel('width, $x$ / cm',fontsize=16)
	ylabel('time, $t$ / hours',fontsize=16)
	gca().invert_yaxis()
	show()

	boundary_index = argmin(abs(difference),axis=0)
	boundary_index = boundary_index[boundary_index>0]
	boundary = xx[0][boundary_index]
	print(max(boundary)-min(boundary))


# execute main program
if __name__ == '__main__' :
	args = get_args()
	main(**args)
