from __future__ import print_function
from PyDSTool import Vode_ODEsystem,args,ContClass,Generator
from PyDSTool.Toolbox.phaseplane import find_fixedpoints

from numpy import zeros
from matplotlib.pyplot import *

class Model(object) :

	def __init__(self, **kwargs) :
		
		self.system = args( name='system', **kwargs )
		self.state_keys = kwargs['varspecs'].keys()


	def initial_fixed_point(self) :
		print('finding initial fixed points...')

		self.odes = Vode_ODEsystem(self.system)
		steady_state = find_fixedpoints(self.odes, n=5, eps=1e-8)
		nS = len(steady_state)
		print('- found %d fixed point(s)'%nS)
		if nS > 0 :

			self.system.ics = steady_state[0]
			self.system.ttype = int
		
		else :
			raise Exception("cannot find initial fixed points")


	def get_cusp(self, params, maxstepsize=0.1, npoints=100):
		'''get cusp bifurcation diagrams via parameter continuation'''
		
		self.initial_fixed_point()
		self.bifurcations = ContClass(Generator.MapSystem(self.system))

		print('finding limit point...')
		self.bifurcations.newCurve(args(
			freepars = [params[0]], name='EQ1', type='EP-C', MaxNumPoints = npoints,
			LocBifPoints='LP', StepSize = 0.1, MaxStepSize=maxstepsize, StopAtPoints = ['B','LP']))

		print('- forwards')
		self.bifurcations['EQ1'].forward()

		if not self.bifurcations['EQ1'].getSpecialPoint('LP1') :
			self.bifurcations.delCurve('EQ1')
		
			self.bifurcations.newCurve(args(
				freepars = [params[0]], name='EQ1', type='EP-C', MaxNumPoints = npoints,
				LocBifPoints='LP', StepSize = 0.1, MaxStepSize=maxstepsize, StopAtPoints = ['B','LP'] ))

			print('- backwards')
			self.bifurcations['EQ1'].backward()
			
  
		if self.bifurcations['EQ1'].getSpecialPoint('LP1') :
		
			print('following limit curve...')
			self.bifurcations.newCurve(args(
				freepars = params, name='LC1back',type='LP-C', MaxNumPoints = npoints,
				initpoint='EQ1:LP1', LocBifPoints='CP', StepSize = 0.1, MaxStepSize=maxstepsize,
				StopAtPoints = ['B'] ))

			print('- backwards')
			self.bifurcations['LC1back'].backward()  
			
			#if not self.bifurcations['LC1'].getSpecialPoint('CP1') :
				#self.bifurcations.delCurve('LC1')
				
			self.bifurcations.newCurve(args(
				freepars = params, name='LC1for',type='LP-C', MaxNumPoints = npoints,
				initpoint='EQ1:LP1', LocBifPoints='CP', StepSize = 0.1, MaxStepSize=maxstepsize,
				StopAtPoints = ['B'] ))
		
			print('- forwards')
			self.bifurcations['LC1for'].forward()
				self.bifurcations['LC1'].forward()
				
			if not (self.bifurcations['LC1for'].getSpecialPoint('CP1')) and not (self.bifurcations['LC1back'].getSpecialPoint('CP1')):
				raise Exception('no cusp point found')

			print('done')

		else :
			 raise Exception('no limit points found!')