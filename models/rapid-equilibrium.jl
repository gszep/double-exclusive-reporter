using Parameters: @unpack
using LinearAlgebra
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
	"iI" => 129.229930836576,
	"atc" => 0.0,
	"iptg" => 0.0
	)

bound(x, c6, c12, K6, K12, n)  = x^2 * ( (K6*10^c6)^n + (K12*10^c12)^n ) / (1.0 + K6*10^c6 + K12*10^c12 )^n
∂bound(x, c6, c12, K6, K12, n) = 2*x * ( (K6*10^c6)^n + (K12*10^c12)^n ) / (1.0 + K6*10^c6 + K12*10^c12 )^n

hill(x, n)  = 1.0 / (1.0 + x^n)
∂hill(x, n) = -n*x^(n-1.0) / (1.0 + x^n)^2

function P76(c,c6,c12,p)
	@unpack KR6,KR12,KS6,KS12, nR,nS, e76,KGR_76,KGS_76 = p
	luxR, lasR, _, _ = 10 .^ c

	activation = KGR_76*bound(luxR,c6,c12,KR6,KR12,nR) + KGS_76*bound(lasR,c6,c12,KS6,KS12,nS)
	return ( e76 + activation ) / ( 1.0 + activation )
end

function ∂S76(c,c6,c12,p)
	@unpack KR6,KR12,KS6,KS12, nR,nS, e76,KGR_76,KGS_76 = p
	luxR, lasR, _, _ = 10 .^ c

	activation =   KGS_76* bound(lasR,c6,c12,KS6,KS12,nS) + KGR_76*bound(luxR,c6,c12,KR6,KR12,nR)
	return (1-e76)*KGS_76*∂bound(lasR,c6,c12,KS6,KS12,nS) / ( 1.0 + activation )^2
end

function ∂R76(c,c6,c12,p)
	@unpack KR6,KR12,KS6,KS12, nR,nS, e76,KGR_76,KGS_76 = p
	luxR, lasR, _, _ = 10 .^ c

	activation =   KGR_76* bound(luxR,c6,c12,KR6,KR12,nR) + KGS_76*bound(lasR,c6,c12,KS6,KS12,nS)
	return (1-e76)*KGR_76*∂bound(luxR,c6,c12,KR6,KR12,nR) / ( 1.0 + activation )^2
end

function P81(c,c6,c12,p)
	@unpack KR6,KR12,KS6,KS12, nR,nS, e81,KGR_81,KGS_81 = p
	luxR, lasR, _, _ = 10 .^ c

	activation = KGR_81*bound(luxR,c6,c12,KR6,KR12,nR) + KGS_81*bound(lasR,c6,c12,KS6,KS12,nS)
	return ( e81 + activation ) / ( 1.0 + activation )
end

function ∂S81(c,c6,c12,p)
	@unpack KR6,KR12,KS6,KS12, nR,nS, e81,KGR_81,KGS_81 = p
	luxR, lasR, _, _ = 10 .^ c

	activation =   KGS_81* bound(lasR,c6,c12,KS6,KS12,nS) + KGR_81*bound(luxR,c6,c12,KR6,KR12,nR)
	return (1-e81)*KGS_81*∂bound(lasR,c6,c12,KS6,KS12,nS) / ( 1.0 + activation )^2
end

function ∂R81(c,c6,c12,p)
	@unpack KR6,KR12,KS6,KS12, nR,nS, e81,KGR_81,KGS_81 = p
	luxR, lasR, _, _ = 10 .^ c

	activation =   KGR_81* bound(luxR,c6,c12,KR6,KR12,nR) + KGS_81*bound(lasR,c6,c12,KS6,KS12,nS)
	return (1-e81)*KGR_81*∂bound(luxR,c6,c12,KR6,KR12,nR) / ( 1.0 + activation )^2
end

function Flog(c, p, c6, c12)
	@unpack growth,capacity, aR33,aS175,aL,aT, nT,nL, dR,dS,dL,dT, iptg,atc, iI,iA = p
	luxR, lasR, lacI, tetR = 10 .^ c

	dluxR = ( capacity*aR33 * hill(tetR,nT) - (growth+dR) * luxR )
	dlasR = ( capacity*aS175 * hill(lacI,nL) - (growth+dS) * lasR )
	dlacI = ( capacity*aL * P76(c,c6,c12,p) - (growth+dL+iptg*iI) * lacI )
	dtetR = ( capacity*aT * P81(c,c6,c12,p) - (growth+dT+atc*iA) * tetR )
	return [dluxR, dlasR, dlacI, dtetR]
end

function Jlog(c, p, c6, c12)
	@unpack growth,capacity, aR33,aS175,aL,aT, nT,nL, dR,dS,dL,dT, iptg,atc, iI,iA = p
	luxR, lasR, lacI, tetR = 10 .^ c
	J = zeros(4,4)

	J[1, 1] = -(growth+dR)
	J[2, 2] = -(growth+dS)

	J[3, 3] = -(growth+dL+iptg*iI)
	J[4, 4] = -(growth+dT+atc*iA)

	J[1, 4] = capacity*aR33  * ∂hill(tetR,nT)
	J[2, 3] = capacity*aS175 * ∂hill(lacI,nL)

	J[3, 1] = capacity*aL * ∂R76(c,c6,c12,p)
	J[3, 2] = capacity*aL * ∂S76(c,c6,c12,p)

	J[4, 1] = capacity*aT * ∂R81(c,c6,c12,p)
	J[4, 2] = capacity*aT * ∂S81(c,c6,c12,p)

	return J * Diagonal(log(10) .* 10 .^ c)
end


function finiteDifferences(F, x::AbstractVector; δ = 1e-9)
	f = F(x)
	N = length(x)
	J = zeros(eltype(f), N, N)
	x1 = copy(x)
	for i=1:N
		x1[i] += δ
		J[:, i] .= (F(x1) .- F(x)) / δ
		x1[i] -= δ
	end
	return J
end

unit_test = fill(false,10000)
for i=1:10000
	sol,c6,c12 = randn(4),randn(),randn()

	unit_test[i] = all(isapprox.( Jlog(sol,P0,c6,c12),
		finiteDifferences( x->Flog(x,P0,c6,c12),sol),
		atol=1e-5*norm(Jlog(sol,P0,c6,c12))))
end
@assert all(unit_test)
