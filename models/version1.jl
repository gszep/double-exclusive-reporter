P0 = Dict(
	# growth parameters
	"capacity" => 1.0,
	"growth" => 1.0,

	# signal affinity and crosstalk
	"KGR_76" => 8.74464935835812,
	"KGS_76" => 0.00040523437497702,
	"KGR_81" => 0.0153475500621728,
	"KGS_81" => 10.9817728137119,
	"KR6" => 0.00014985256979569,
	"KS6" => 4.22023069018813e-06,
	"KR12" => 1.02380021889859e-08,
	"KS12" => 0.0246441407606443,

	# hill coefficients
	"nR" => 0.78145979866281,
	"nS" => 0.984847709819578,
	"nL" => 0.643742771597066,
	"nT" => 2.18076527621396,

	# activations
	"aR33" => 9.61357377397729,
	"aS175" => 2.92159578355242,
	"aL" => 65.2589987503168,
	"aT" => 60.6259602781465,

	"e76" => 0.0203965308312284,
	"e81" => 0.00611932024058379,
	"aYFP" => 33182.8359178623,
	"aCFP" => 25061.89983953,

	# degredations
	"dR" => 9.98791524173607,
	"dS" => 10.0,
	"dYFP" => 0.146735252880641,
	"dCFP" => 0.462404866514184,
	"dL" => 0.647723328994184,
	"dT" => 2.72162807214286,

	# derepressors
	"iA" => 0.146707944645707,
	"iI" => 129.229930836576
	)

function F(x, p, c6, c12)
	luxR, lasR, lacI, tetR = x

	boundLasR(lasR) = lasR*lasR * ( (p["KS6"]*c6)^p["nS"] + (p["KS12"]*c12)^p["nS"] ) / (1.0 + p["KS6"]*c6 + p["KS12"]*c12 )^p["nS"]
	boundLuxR(luxR) = luxR*luxR * ( (p["KR6"]*c6)^p["nR"] + (p["KR12"]*c12)^p["nR"] ) / (1.0 + p["KR6"]*c6 + p["KR12"]*c12 )^p["nR"]
	P76(luxR,lasR) = ( p["e76"] + p["KGR_76"]*boundLuxR(luxR) + p["KGS_76"]*boundLasR(lasR) ) / ( 1.0 + p["KGR_76"]*boundLuxR(luxR) + p["KGS_76"]*boundLasR(lasR) )
	P81(luxR,lasR) = ( p["e81"] + p["KGR_81"]*boundLuxR(luxR) + p["KGS_81"]*boundLasR(lasR) ) / ( 1.0 + p["KGR_81"]*boundLuxR(luxR) + p["KGS_81"]*boundLasR(lasR) )
	PLac(lacI) = 1.0 / (1.0 + lacI^p["nL"])
	PTet(tetR) = 1.0 / (1.0 + tetR^p["nT"])

	dluxR = p["capacity"]*p["aR33"]*PTet(tetR) - luxR*(p["growth"]+p["dR"])
	dlasR = p["capacity"]*p["aS175"]*PLac(lacI) - lasR*(p["growth"]+p["dS"])
	dlacI = p["capacity"]*p["aL"]*P76(luxR,lasR) - lacI*(p["growth"]+p["dL"])
	dtetR = p["capacity"]*p["aT"]*P81(luxR,lasR) - tetR*(p["growth"]+p["dT"])

	return [dluxR, dlasR, dlacI, dtetR]
end

Bound(x, c6, c12, K6, K12, n) = x*x * ( (K6*10^c6)^n + (K12*10^c12)^n ) / (1.0 + K6*10^c6 + K12*10^c12 )^n
Prom(luxR,lasR,c6,c12,KR6,KR12,nR,KS6,KS12,nS,eps,KGR,KGS) = ( eps + KGR*Bound(luxR,c6, c12,  KR6, KR12, nR) + KGS*Bound(lasR, c6, c12, KS6, KS12, nS) ) / ( 1.0 + KGR*Bound(luxR, c6, c12, KR6, KR12, nR) + KGS*Bound(lasR, c6, c12, KS6, KS12, nS) )
IHill(x, n) = 1.0 / (1.0 + x^n)

function Flog(x, p, c6, c12; eps = 1e-7)
	luxR, lasR, lacI, tetR = x

	dluxR = ( p["aR33"]*IHill(10^tetR, p["nT"]) - 10^luxR*(p["growth"]+p["dR"]) ) / (log(10)*10^luxR + eps)
	dlasR = ( p["aS175"]*IHill(10^lacI, p["nL"]) - 10^lasR*(p["growth"]+p["dS"]) ) / (log(10)*10^lasR + eps)
	dlacI = ( p["aL"]*Prom(10^luxR,10^lasR,c6,c12,p["KR6"],p["KR12"],p["nR"],p["KS6"],p["KS12"],p["nS"],p["e76"],p["KGR_76"],p["KGS_76"]) - 10^lacI*(p["growth"]+p["dL"]) ) / (log(10)*10^lacI + eps)
	dtetR = ( p["aT"]*Prom(10^luxR,10^lasR,c6,c12,p["KR6"],p["KR12"],p["nR"],p["KS6"],p["KS12"],p["nS"],p["e81"],p["KGR_81"],p["KGS_81"]) - 10^tetR*(p["growth"]+p["dT"]) ) / (log(10)*10^tetR + eps)
	return [dluxR, dlasR, dlacI, dtetR]
end

dBound(x, c6, c12, K6, K12, n) = 2 * x * ( (K6*10^c6)^n + (K12*10^c12)^n ) / (1.0 + K6*10^c6 + K12*10^c12 )^n
dPromdluxR(luxR,lasR,c6,c12,KR6,KR12,nR,KS6,KS12,nS,e,KGR,KGS) = ((1.0-e)*KGR*dBound(luxR, c6, c12, KR6, KR12, nR)) / (1.0 + KGR*Bound(luxR, c6, c12, KR6, KR12, nR) + KGS*Bound(lasR, c6, c12, KS6, KS12, nS))^2
dPromdlasR(luxR,lasR,c6,c12,KR6,KR12,nR,KS6,KS12,nS,e,KGR,KGS) = ((1.0-e)*KGS*dBound(lasR, c6, c12, KS6, KS12, nS)) / (1.0 + KGR*Bound(luxR, c6, c12, KR6, KR12, nR) + KGS*Bound(lasR, c6, c12, KS6, KS12, nS))^2
dIHill(x, n) = -n*x^(n-1.0) / (1.0 + x^n)^2

function Jlog(x, p, c6, c12; eps = 1e-7)
	luxR, lasR, lacI, tetR = x
	n = 4
	J = zeros(4, 4)
	dP76dluxR = dPromdluxR(10^luxR, 10^lasR, c6, c12, p["KR6"], p["KR12"], p["nR"], p["KS6"], p["KS12"], p["nS"], p["e76"], p["KGR_76"], p["KGS_76"])
	dP76dlasR = dPromdlasR(10^luxR, 10^lasR, c6, c12, p["KR6"], p["KR12"], p["nR"], p["KS6"], p["KS12"], p["nS"], p["e76"], p["KGR_76"], p["KGS_76"])
	dP81dluxR = dPromdluxR(10^luxR, 10^lasR, c6, c12, p["KR6"], p["KR12"], p["nR"], p["KS6"], p["KS12"], p["nS"], p["e81"], p["KGR_81"], p["KGS_81"])
	dP81dlasR = dPromdlasR(10^luxR, 10^lasR, c6, c12, p["KR6"], p["KR12"], p["nR"], p["KS6"], p["KS12"], p["nS"], p["e81"], p["KGR_81"], p["KGS_81"])
	J[1, 1] = -((10^luxR*log(10)*((p["dR"] + p["growth"])*eps + p["aR33"]*log(10)*IHill(10^tetR, p["nT"]))) / (eps + 10^luxR*log(10))^2)
	J[1, 4] = (10^tetR*p["aR33"]*log(10)*dIHill(10^tetR, p["nT"])) / (eps + 10^luxR*log(10))
	J[2, 2] = -((10^lasR*log(10)*((p["dS"] + p["growth"])*eps + p["aS175"]*log(10)*IHill(10^lacI, p["nL"]))) / (eps + 10^lasR*log(10))^2)
	J[2, 3] = (10^lacI*p["aS175"]*log(10)*dIHill(10^lacI, p["nL"])) / (eps + 10^lasR*log(10))
	J[3, 1] = (10^luxR*p["aL"]*log(10)*dP76dluxR) / (eps + 10^lacI*log(10))
	J[3, 2] = (10^lasR*p["aL"]*log(10)*dP76dlasR) / (eps + 10^lacI*log(10))
	J[3, 3] = -((10^lacI*log(10)*((p["dL"] + p["growth"])*eps + p["aL"]*log(10)*Prom(10^luxR,10^lasR,c6,c12,p["KR6"],p["KR12"],p["nR"],p["KS6"],p["KS12"],p["nS"],p["e76"],p["KGR_76"],p["KGS_76"]))) / (eps + 10^lacI*log(10))^2)
	J[4, 1] = (10^luxR*p["aT"]*log(10)*dP81dluxR) / (eps + 10^tetR*log(10))
	J[4, 2] = (10^lasR*p["aT"]*log(10)*dP81dlasR) / (eps + 10^tetR*log(10))
	J[4, 4] = -((10^tetR*log(10)*((p["dT"] + p["growth"])*eps + p["aT"]*log(10)*Prom(10^luxR,10^lasR,c6,c12,p["KR6"],p["KR12"],p["nR"],p["KS6"],p["KS12"],p["nS"],p["e81"],p["KGR_81"],p["KGS_81"]))) / (eps + 10^tetR*log(10))^2)

	return J
end

# Test Jacobian
c6 = 5.
c12 = 5.
sol = [-1.,-1.,-1.,-1.]
J0 = Jlog(sol,P0,c6,c12)
function F(x)
	Flog(x, P0, c6, c12)
end
h = 1e-6
Japprox =
	hcat([ F(sol+[h,0.,0.,0.]) - F(sol)
	, F(sol+[0.,h,0.,0.]) - F(sol)
	, F(sol+[0.,0.,h,0.]) - F(sol)
	, F(sol+[0.,0.,0.,h]) - F(sol)
	]) / h
