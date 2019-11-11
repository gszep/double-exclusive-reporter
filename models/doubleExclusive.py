# kwargs for initialising Model object with double exclusive reporter
system_specifications = {

	# state and parameter domains
	'pdomain' : { 'c12': [-1.39,4.88], 'c6': [-1.39,4.88] },
	'xdomain' : {'luxR': [-5,5], 'lasR': [-5,5], 'lacI': [-5,5], 'tetR': [-5,5]},

	# rate functions
	'varspecs' : {
		'luxR': '( capacity*aR33*PTet(10**tetR) - 10**luxR*(growth+dR) )/(log(10)*10**luxR+0.01)',
		'lasR': '( capacity*aS175*PLac(10**lacI) - 10**lasR*(growth+dR) )/(log(10)*10**lasR+0.01)',
		'lacI': '( capacity*aL*P76(10**luxR,10**lasR) - 10**lacI*(growth+dL+iI*IPTG) )/(log(10)*10**lacI+0.01)',
		'tetR': '( capacity*aT*P81(10**luxR,10**lasR) - 10**tetR*(growth+dT+iA*ATC) )/(log(10)*10**tetR+0.01)',
	},

	# promoter activities
	'fnspecs' : {
		'boundLasR' : (['lasR'],'lasR**2 * ( (KS6*10**c6)**nS + (KS12*10**c12)**nS ) / (1.0 + KS6*10**c6 + KS12*10**c12 )**nS'),
		'boundLuxR' : (['luxR'],'luxR**2 * ( (KR6*10**c6)**nR + (KR12*10**c12)**nR ) / (1.0 + KR6*10**c6 + KR12*10**c12 )**nR'),
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
	'KGR_76' : 8.74464935835812, 
	'KGS_76' : 0.00040523437497702, 
	'KGR_81' : 0.0153475500621728,
	'KGS_81' : 10.9817728137119,
	'KR6' : 0.00014985256979569, 
	'KS6' : 4.22023069018813e-06, 
	'KR12' : 1.02380021889859e-08,
	'KS12' : 0.0246441407606443,

	# hill coefficients
	'nR' : 0.78145979866281,
	'nS' : 0.984847709819578,
	'nL' : 0.643742771597066,
	'nT' : 2.18076527621396,

	# activations
	'aR33' : 9.61357377397729,
	'aS175' : 2.92159578355242,
	'aL' : 65.2589987503168,
	'aT' : 60.6259602781465, 

	'eps76' : 0.0203965308312284,
	'eps81' : 0.00611932024058379, 
	'aYFP' : 33182.8359178623, 
	'aCFP' : 25061.89983953, 

	# degredations
	'dR' : 9.98791524173607, 
	'dYFP' : 0.146735252880641, 
	'dCFP' : 0.462404866514184, 
	'dL' : 0.647723328994184, 
	'dT' : 2.72162807214286, 

	# derepressors
	'iA' : 0.146707944645707, 
	'iI' : 129.229930836576,
	'ATC' : 0.0,
	'IPTG' : 0.0,

	# morphogens
	'c12': 4.88, 'c6' : 4.88,

	# relays
	'kC6': 0.0, 'kC12': 0.0,
	'dluxI': 0.00, 'dlasI': 0.0
}