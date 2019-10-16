from __future__ import print_function
from PyDSTool import Vode_ODEsystem,args,ContClass,Generator
from numpy import zeros
from matplotlib.pyplot import *

class Model(object) :

	def __init__(self, **kwargs) :
		
		self.system = args( name='system', **kwargs )
		self.state_keys = kwargs['varspecs'].keys()

		init = zeros(len(self.state_keys))
		self.system.ics = dict(zip(self.state_keys,init))

		
	def integrate(self) :
		'''integrate system in given time domain'''

		self.odes = Vode_ODEsystem(self.system)
		trajectory = self.odes.compute('trajectory')
		trajectory = trajectory.sample()

		figure(figsize=(6,6))
		title('stability check', y=1.02, fontsize=16)
		plot(trajectory); xlabel('iteration, i',fontsize=16); ylabel('state values',fontsize=16)
		show()

		steady_state = trajectory.coordarray[:,-1]
		steady_state = dict(zip(self.state_keys,steady_state))
		self.system.ics = steady_state


	def get_cusp(self,*params) :
		'''get cusp bifurcation diagrams via parameter continuation'''

		self.system.ttype = int; del self.system.tdomain
		self.bifurcations = ContClass(Generator.MapSystem(self.system))
		
		print('finding limit point...')
		self.bifurcations.newCurve(args(
			freepars = [list(params)[0]], name='EQ1', type='EP-C', MaxNumPoints = 100,
			LocBifPoints='LP', StepSize = 0.1, MaxStepSize=0.1, StopAtPoints = ['B','LP'],
			TestTol=0.1 ))
		
		print(' - forwards')
		self.bifurcations['EQ1'].forward()
		
		if not self.bifurcations['EQ1'].getSpecialPoint('LP1') :
			self.bifurcations.delCurve('EQ1')
		
			self.bifurcations.newCurve(args(
				freepars = [list(params)[0]], name='EQ1', type='EP-C', MaxNumPoints = 100,
				LocBifPoints='LP', StepSize = 0.1, MaxStepSize=0.1, StopAtPoints = ['B','LP'] ))

			print(' - backwards')
			self.bifurcations['EQ1'].backward()
			
  
		if self.bifurcations['EQ1'].getSpecialPoint('LP1') :
		
			print('following limit curve...')
			self.bifurcations.newCurve(args(
				freepars = list(params), name='LC1',type='LP-C', MaxNumPoints = 100,
				initpoint='EQ1:LP1', LocBifPoints='CP', StepSize = 0.1, MaxStepSize=0.1,
				StopAtPoints = ['B'] ))

			print(' - forwards')
			self.bifurcations['LC1'].forward()  
			
			if not self.bifurcations['LC1'].getSpecialPoint('CP1') :
				self.bifurcations.delCurve('LC1')
				
				self.bifurcations.newCurve(args(
					freepars = list(params), name='LC1',type='LP-C', MaxNumPoints = 100,
					initpoint='EQ1:LP1', LocBifPoints='CP', StepSize = 0.1, MaxStepSize=0.1,
					StopAtPoints = ['B'] ))
			
				print(' - backwards')
				self.bifurcations['LC1'].backward()
				
				if not self.bifurcations['LC1'].getSpecialPoint('CP1') :
					raise Exception('no cusp point found')

			print('done')


		else :
			 raise Exception('no limit points found!')