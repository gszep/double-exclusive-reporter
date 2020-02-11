# kwargs for initialising Model object with double exclusive reporter
system_specifications = {

	# state and parameter domains
	'pdomain' : { 'c12': [-1.39,4.88], 'c6': [-1.39,4.88] },
	'xdomain' : {'luxR': [-5,5], 'lasR': [-5,5], 'lacI': [-5,5], 'tetR': [-5,5]},

	# rate functions
	'varspecs' : {
		'luxR': '( capacity*aR33*PTet(10**tetR) - 10**luxR*(growth+dR) )/(log(10)*10**luxR+0.0000001)',
		'lasR': '( capacity*aS175*PLac(10**lacI) - 10**lasR*(growth+dS) )/(log(10)*10**lasR+0.0000001)',
		'lacI': '( capacity*aL*P76(10**luxR,10**lasR) - 10**lacI*(growth+dL+iI*IPTG) )/(log(10)*10**lacI+0.0000001)',
		'tetR': '( capacity*aT*P81(10**luxR,10**lasR) - 10**tetR*(growth+dT+iA*ATC) )/(log(10)*10**tetR+0.0000001)',
	},

	# promoter activities
	'fnspecs' : {
		'boundLasR' : (['lasR'],'lasR**2 * ( (eS6*10**c6)**nS + (10**c12)**nS )'),
		'boundLuxR' : (['luxR'],'luxR**2 * ( (10**c6)**nR + (eR12*10**c12)**nR )'),
		'P76' : (['luxR','lasR'],'( eps76 + KGR_76*boundLuxR(luxR) + KGS_76*boundLasR(lasR) ) / ( 1.0 + KGR_76*boundLuxR(luxR) + KGS_76*boundLasR(lasR) )'),
		'P81' : (['luxR','lasR'],'( eps81 + KGR_81*boundLuxR(luxR) + KGS_81*boundLasR(lasR) ) / ( 1.0 + KGR_81*boundLuxR(luxR) + KGS_81*boundLasR(lasR) )'),
		'PLac' : (['lacI'],'1.0/(1.0+lacI**nL)'),
		'PTet' : (['tetR'],'1.0/(1.0+tetR**nT)')
	}
}

parameters = {

	# growth parameters
	'capacity' : 1.0,
	'growth' : 1.0,

	# signal affinity and crosstalk
	'KGR_76' : 10.0,
	'KGS_76' : 0.0001,
	'KGR_81' : 0.01,
	'KGS_81' : 10.0,
	'eS6' : 1e-6,
	'eR12' : 1e-8,

	# hill coefficients
	'nR' : 1.0,
	'nS' : 1.0,
	'nL' : 1.0,
	'nT' : 2.0,

	# activations
	'aR33' : 10.0,
	'aS175' : 3.0,
	'aL' : 10.0,
	'aT' : 10.0,

	'eps76' : 0.01,
	'eps81' : 0.005,
	'aYFP' : 1e4,
	'aCFP' : 1e4,

	# degredations
	'dR' : 50.0,
	'dS' : 50.0,
	'dYFP' : 0.1,
	'dCFP' : 0.4,
	'dL' : 1.0,
	'dT' : 1.0,

	# derepression
	'iA' : 0.2,
	'iI' : 200.0,

	# inputs
	#'c12': 4.88, 'c6' : 4.88,
	'c12': 3.0, 'c6' : 3.0,
	'ATC' : 0.0, 'IPTG' : 0.0,
	
	# relays
	'kC6': 0.0, 'kC12': 0.0,
	'dluxI': 0.00, 'dlasI': 0.0
}