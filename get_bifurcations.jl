include("FoldCurve.jl")

n_components = 4
parameters = (

    capacity = 1.2,
    growth = 1.0,
    K = 2.7609,

    # signal affinity and crosstalk
    KGR_76 = 8.74464935835812,
    KGS_76 = 0.00040523437497702,
    KGR_81 = 0.0153475500621728,
    KGS_81 = 10.9817728137119,
    KR6 = 0.00014985256979569,
    KS6 = 4.22023069018813E-06,
    KR12 = 1.02380021889859E-08,
    KS12 = 0.0246441407606443,

    # hill coefficients
    nR = 0.78145979866281,
    nS = 0.984847709819578,
    nL = 0.643742771597066,
    nT = 2.18076527621396,

    # activations
    aR33 = 9.61357377397729,
    aS175 = 2.92159578355242,
    aL = 65.2589987503168,
    aT = 60.6259602781465,

    e76 = 0.0203965308312284,
    e81 = 0.00611932024058379,
    aYFP = 33182.8359178623,
    aCFP = 25061.89983953,

    # degredations
    dR = 9.98791524173607,
    dYFP = 0.146735252880641,
    dCFP = 0.462404866514184,
    dL = 0.647723328994184,
    dT = 2.72162807214286,

    # derepressors
    iA = 0.146707944645707,
    iI = 129.229930836576,
    ATC = 0,
    IPTG = 0,

    # morphogens
	c12 = 4.88, c6 = -1.39
)


function rates( u, p )

    # unpacking parameters / variables
    capacity,growth,K,KGR_76,KGS_76,KGR_81,KGS_81,KR6,KS6,KR12,KS12,nR,nS,nL,nT,aR33,aS175,aL,aT,e76,e81,aYFP,aCFP,dR,dYFP,dCFP,dL,dT,iA,iI,ATC,IPTG,c12,c6 = p
	luxR,lasR,lacI,tetR = u
	f = similar(u)

	f[1] = ( capacity*aR33*PTet(10^tetR,nT) - 10^luxR*(growth+dR) )/(log(10)*10^luxR)
	f[2] = ( capacity*aS175*PLac(10^lacI,nL) - 10^lasR*(growth+dR) )/(log(10)*10^lasR)
	f[3] = ( capacity*aL*P76(10^luxR,10^lasR,c6,c12,nR,KR6,KR12,nS,KS6,KS12,e76,KGR_76,KGS_76) - 10^lacI*(growth+dL+iI*IPTG) )/(log(10)*10^lacI)
	f[4] = ( capacity*aT*P81(10^luxR,10^lasR,c6,c12,nR,KR6,KR12,nS,KS6,KS12,e81,KGR_81,KGS_81) - 10^tetR*(growth+dT+iA*ATC) )/(log(10)*10^tetR)

    return f
end

function Jacobian( u ,p )

    # unpacking parameters / variables
    capacity,growth,K,KGR_76,KGS_76,KGR_81,KGS_81,KR6,KS6,KR12,KS12,nR,nS,nL,nT,aR33,aS175,aL,aT,e76,e81,aYFP,aCFP,dR,dYFP,dCFP,dL,dT,iA,iI,ATC,IPTG,c12,c6 = p
	luxR,lasR,lacI,tetR = u
	J = zeros( length(u),  length(u))

	J[1,1] = ( capacity*aR33*PTet(10^tetR,nT) - 10^luxR*(growth+dR) )/(log(10)*10^luxR)
	J[2,2] = ( capacity*aS175*PLac(10^lacI,nL) - 10^lasR*(growth+dR) )/(log(10)*10^lasR)
	J[3,3]  = ( capacity*aL*P76(10^luxR,10^lasR,c6,c12,nR,KR6,KR12,nS,KS6,KS12,e76,KGR_76,KGS_76) - 10^lacI*(growth+dL+iI*IPTG) )/(log(10)*10^lacI)
	J[4,4] = ( capacity*aT*P81(10^luxR,10^lasR,c6,c12,nR,KR6,KR12,nS,KS6,KS12,e81,KGR_81,KGS_81) - 10^tetR*(growth+dT+iA*ATC) )/(log(10)*10^tetR)

    return J
end

boundLasR(lasR,c6,c12,nS,KS6,KS12) = lasR^2 * ( (KS6*10^c6)^nS + (KS12*10^c12)^nS ) / (1.0 + KS6*10^c6 + KS12*10^c12 )^nS
boundLuxR(luxR,c6,c12,nR,KR6,KR12) = luxR^2 * ( (KR6*10^c6)^nR + (KR12*10^c12)^nR ) / (1.0 + KR6*10^c6 + KR12*10^c12 )^nR
P76(luxR,lasR,c6,c12,nR,KR6,KR12,nS,KS6,KS12,e76,KGR_76,KGS_76) = ( e76 + KGR_76*boundLuxR(luxR,c6,c12,nR,KR6,KR12) + KGS_76*boundLasR(lasR,c6,c12,nS,KS6,KS12) ) / ( 1.0 + KGR_76*boundLuxR(luxR,c6,c12,nR,KR6,KR12) + KGS_76*boundLasR(lasR,c6,c12,nS,KS6,KS12) )
P81(luxR,lasR,c6,c12,nR,KR6,KR12,nS,KS6,KS12,e81,KGR_81,KGS_81) = ( e81 + KGR_81*boundLuxR(luxR,c6,c12,nR,KR6,KR12) + KGS_81*boundLasR(lasR,c6,c12,nS,KS6,KS12) ) / ( 1.0 + KGR_81*boundLuxR(luxR,c6,c12,nR,KR6,KR12) + KGS_81*boundLasR(lasR,c6,c12,nS,KS6,KS12) )
PLac(lacI,nL) = 1.0/(1.0+lacI^nL)
PTet(tetR,nT) = 1.0/(1.0+tetR^nT)

capacity,growth,K,KGR_76,KGS_76,KGR_81,KGS_81,KR6,KS6,KR12,KS12,nR,nS,nL,nT,aR33,aS175,aL,aT,e76,e81,aYFP,aCFP,dR,dYFP,dCFP,dL,dT,iA,iI,ATC,IPTG,c12,c6 = parameters
p = [capacity,growth,K,KGR_76,KGS_76,KGR_81,KGS_81,KR6,KS6,KR12,KS12,nR,nS,nL,nT,aR33,aS175,aL,aT,e76,e81,aYFP,aCFP,dR,dYFP,dCFP,dL,dT,iA,iI,ATC,IPTG,c12,c6]
curve = @time limit_curve(
	(u,c12,c6) -> rates(u,   [capacity,growth,K,KGR_76,KGS_76,KGR_81,KGS_81,KR6,KS6,KR12,KS12,nR,nS,nL,nT,aR33,aS175,aL,aT,e76,e81,aYFP,aCFP,dR,dYFP,dCFP,dL,dT,iA,iI,ATC,IPTG,c12,c6]),
	(u,c12,c6) -> Jacobian(u,[capacity,growth,K,KGR_76,KGS_76,KGR_81,KGS_81,KR6,KS6,KR12,KS12,nR,nS,nL,nT,aR33,aS175,aL,aT,e76,e81,aYFP,aCFP,dR,dYFP,dCFP,dL,dT,iA,iI,ATC,IPTG,c12,c6]),
	[0.0,0.0,0.0,0.0],2.0,2.0; maxSteps=200, pMin=-1.39, pMax=4.88 )

plot(curve.branch[2,:],curve.branch[1,:],
	label="fold curve",color="green")

rates([0.,0.,0.,0.],p)
