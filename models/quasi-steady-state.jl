using Parameters: @unpack
using LinearAlgebra
parameters = Dict(

	# growth parameters
	"capacity" => 1.0,
	"growth" => 1.0,

	# signal affinity and crosstalk
	"KGR_76" => 0.0342,
	"KGS_76" => 0.000165,
	"KGR_81" => 0.000742,
	"KGS_81" => 10.4,
	"eS6" => 1.01e-8,
	"eR12" => 4.2e-7,

	# hill coefficients
	"nR" => 0.618,
	"nS" => 0.324,
	"nL" => 0.58,
	"nT" => 1.65,

	# activations
	"aR33" => 16.7,
	"aS175" => 2.86,
	"aL" => 3480.0,
	"aT" => 67.4,

	"e76" => 0.0114,
	"e81" => 0.00267,

	# degredations
	"dR" => 24.6,
	"dS" => 8.58,
	"dL" => 1.0,
	"dT" => 1.0,

	# derepressors
	"iA" => 0.421,
	"iI" => 418.0,
	"atc" => 0.0,
	"iptg" => 0.0
)

bound(x, c6, c12, K6, K12, n)  = x^2 * ( (K6*10^c6)^n + (K12*10^c12)^n )
∂bound(x, c6, c12, K6, K12, n) = 2*x * ( (K6*10^c6)^n + (K12*10^c12)^n )

hill(x, n)  = 1.0 / (1.0 + x^n)
∂hill(x, n) = -n*x^(n-1.0) / (1.0 + x^n)^2

function P76(c,c6,c12,p)
	@unpack eS6,eR12, nR,nS, e76,KGR_76,KGS_76 = p
	luxR, lasR, _, _ = 10 .^ c;

	activation = KGR_76*bound(luxR,c6,c12,1,eR12,nR) + KGS_76*bound(lasR,c6,c12,eS6,1,nS)
	return ( e76 + activation ) / ( 1.0 + activation )
end

function ∂S76(c,c6,c12,p)
	@unpack eS6,eR12, nR,nS, e76,KGR_76,KGS_76 = p
	luxR, lasR, _, _ = 10 .^ c

	activation =   KGS_76* bound(lasR,c6,c12,eS6,1,nS) + KGR_76*bound(luxR,c6,c12,1,eR12,nR)
	return (1-e76)*KGS_76*∂bound(lasR,c6,c12,eS6,1,nS) / ( 1.0 + activation )^2
end

function ∂R76(c,c6,c12,p)
	@unpack eS6,eR12, nR,nS, e76,KGR_76,KGS_76 = p
	luxR, lasR, _, _ = 10 .^ c

	activation =   KGR_76* bound(luxR,c6,c12,1,eR12,nR) + KGS_76*bound(lasR,c6,c12,eS6,1,nS)
	return (1-e76)*KGR_76*∂bound(luxR,c6,c12,1,eR12,nR) / ( 1.0 + activation )^2
end

function P81(c,c6,c12,p)
	@unpack eS6,eR12, nR,nS, e81,KGR_81,KGS_81 = p
	luxR, lasR, _, _ = 10 .^ c

	activation = KGR_81*bound(luxR,c6,c12,1,eR12,nR) + KGS_81*bound(lasR,c6,c12,eS6,1,nS)
	return ( e81 + activation ) / ( 1.0 + activation )
end

function ∂S81(c,c6,c12,p)
	@unpack eS6,eR12, nR,nS, e81,KGR_81,KGS_81 = p
	luxR, lasR, _, _ = 10 .^ c

	activation =   KGS_81* bound(lasR,c6,c12,eS6,1,nS) + KGR_81*bound(luxR,c6,c12,1,eR12,nR)
	return (1-e81)*KGS_81*∂bound(lasR,c6,c12,eS6,1,nS) / ( 1.0 + activation )^2
end

function ∂R81(c,c6,c12,p)
	@unpack eS6,eR12, nR,nS, e81,KGR_81,KGS_81 = p
	luxR, lasR, _, _ = 10 .^ c

	activation =   KGR_81* bound(luxR,c6,c12,1,eR12,nR) + KGS_81*bound(lasR,c6,c12,eS6,1,nS)
	return (1-e81)*KGR_81*∂bound(luxR,c6,c12,1,eR12,nR) / ( 1.0 + activation )^2
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
	u,c6,c12 = randn(4),randn(),randn()

	unit_test[i] = all(isapprox.( Jlog(u,parameters,c6,c12),
		finiteDifferences( x->Flog(x,parameters,c6,c12),u),
		atol=1e-5*norm(Jlog(u,parameters,c6,c12))))
end
@assert all(unit_test)
