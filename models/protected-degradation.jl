using Parameters: @unpack
using LinearAlgebra

##############################  model parameters; change them here if you want to experiment
θ = Dict(

	# growth
	"ρ" => 0.002,
	"K" => 2.7,
	"r" => 1.0,

	# signal diffusion
	"D₆" => 0.0000018,
	"D₁₂" => 0.0000009,

	# dissociation constants
	"R⁻P76" => 0.0342,
	"S⁻P76" => 0.000165,
	"R⁻P81" => 0.000742,
	"S⁻P81" => 10.4,

	# crosstalk
	"S⁻C6" => 1.01e-8,
	"R⁻C12" => 4.2e-7,

	# leaky expression
	"ε₇₆" => 0.0114,
	"ε₈₁" => 0.00267,

	# hill coefficients
	"nᴿ" => 0.618,
	"nˢ" => 0.324,
	"nᴸ" => 0.58,
	"nᵀ" => 1.65,

	# synthesis rates
	"aᴿ" => 16.7,
	"aˢ" => 2.86,
	"aᴸ" => 2480.0,
	"aᵀ" => 67.4,
	"aᶜ" => 2.16E+04,
    "aʸ" => 1.62E+04,

	# degredations
	"dᴿ" => 24.6,
	"dˢ" => 8.58,
	"dᴸ" => 1.0,
	"dᵀ" => 1.0,
	"dᶜ" => 0.379,
    "dʸ" => 0.14,

	# derepressors
	"iᴬ" => 0.421,
	"iᴵ" => 418.0,
	"A" => 0.0,
	"I" => 0.0,

	#relay
	"k₆" => 300.0,
	"k₁₂" => 0.0, #186.113764100785,

	"dk₆" => 1.07,
    "dk₁₂" => 1.09
)

##############################  receiver binding
bound(x, cₓ, cₑ, ε, nₓ)  = x .* x .* ( abs.(cₓ) .^ nₓ .+ abs.(ε.*cₑ) .^ nₓ )
∂bound(x, cₓ, cₑ, ε, nₓ) = 2x .* ( abs.(cₓ) .^ nₓ .+ abs.(ε.*cₑ) .^ nₓ )

##############################  hill activations
hill(x, n)  = 1 ./ (1 .+ abs.(x).^n)
∂hill(x, n) = -n .* abs.(x).^(n.-1) ./ (1 .+ abs.(x).^n).^2

##############################  promoter activations
function P76( states, parameters )
	@unpack R⁻C12,S⁻C6, nˢ,nᴿ, R⁻P76,S⁻P76,ε₇₆ = parameters
	@unpack R,S, c₆,c₁₂ = states

	activation = R⁻P76 .* bound( R, c₆, c₁₂, R⁻C12, nᴿ ) .+ S⁻P76 .* bound( S, c₁₂, c₆, S⁻C6, nˢ )
	return  (ε₇₆ .+ activation) ./ (1 .+ activation)
end

function P81( states, parameters )
	@unpack R⁻C12,S⁻C6, nˢ,nᴿ, R⁻P81,S⁻P81,ε₈₁ = parameters
	@unpack R,S, c₆,c₁₂ = states

	activation = R⁻P81 .* bound( R, c₆, c₁₂, R⁻C12, nᴿ ) .+ S⁻P81 .* bound( S, c₁₂, c₆, S⁻C6, nˢ )
	return (ε₈₁ .+ activation) ./ (1 .+ activation)
end

##############################  partial derivatives used in jacobian
function ∂R76( states, parameters )
	@unpack R⁻C12,S⁻C6, nˢ,nᴿ, R⁻P76,S⁻P76,ε₇₆ = parameters
	@unpack R,S, c₆,c₁₂ = states

	activation = R⁻P76 .* bound( R, c₆, c₁₂, R⁻C12, nᴿ ) .+ S⁻P76 .* bound( S, c₁₂, c₆, S⁻C6, nˢ )
	return (1-ε₇₆) .* R⁻P76 .* ∂bound( R, c₆, c₁₂, R⁻C12, nᴿ ) ./ (1 .+ activation).^2
end

function ∂S76( states, parameters )
	@unpack R⁻C12,S⁻C6, nˢ,nᴿ, R⁻P76,S⁻P76,ε₇₆ = parameters
	@unpack R,S, c₆,c₁₂ = states

	activation = R⁻P76 .* bound( R, c₆, c₁₂, R⁻C12, nᴿ ) .+ S⁻P76 .* bound( S, c₁₂, c₆, S⁻C6, nˢ )
	return (1-ε₇₆) .* S⁻P76 .* ∂bound( S, c₁₂, c₆, S⁻C6, nˢ ) ./ (1 .+ activation).^2
end

function ∂R81( states, parameters )
	@unpack R⁻C12,S⁻C6, nˢ,nᴿ, R⁻P81,S⁻P81,ε₈₁ = parameters
	@unpack R,S, c₆,c₁₂ = states

	activation = R⁻P81 .* bound( R, c₆, c₁₂, R⁻C12, nᴿ ) .+ S⁻P81 .* bound( S, c₁₂, c₆, S⁻C6, nˢ )
	return (1-ε₈₁) .* R⁻P81 .* ∂bound( R, c₆, c₁₂, R⁻C12, nᴿ ) ./ (1 .+ activation).^2
end

function ∂S81( states, parameters )
	@unpack R⁻C12,S⁻C6, nˢ,nᴿ, R⁻P81,S⁻P81,ε₈₁ = parameters
	@unpack R,S, c₆,c₁₂ = states

	activation = R⁻P81 .* bound( R, c₆, c₁₂, R⁻C12, nᴿ ) .+ S⁻P81 .* bound( S, c₁₂, c₆, S⁻C6, nˢ )
	return (1-ε₈₁) .* S⁻P81 .* ∂bound( S, c₁₂, c₆, S⁻C6, nˢ ) ./ (1 .+ activation).^2
end

############################# rate functions used for spatial simulations
function rates( states::Dict, parameters::Dict, t::Float64 )
	@unpack r,K, D₆,D₁₂, aᴿ,aˢ,aᴸ,aᵀ,aᶜ,aʸ, nᵀ,nᴸ, dᴿ,dˢ,dᴸ,dᵀ,dᶜ,dʸ, I,A, iᴵ,iᴬ, k₆,k₁₂, dk₆,dk₁₂ = parameters
	@unpack R, S, L, T, c₆, c₁₂, CFP,YFP,ρ, luxI,lasI = states

	# growth
	γ = r .* (1 .- ρ./K )
	dρ = γ .* ρ

	# calculate reaction rates
	dR = aᴿ .* hill(T,nᵀ) .- (γ.+dᴿ) .* R
	dS = aˢ .* hill(L,nᴸ) .- (γ.+dˢ) .* S

	dL = aᴸ .* P76(states, parameters) - (γ.+dᴸ.+iᴵ*I) .* L
	dT = aᵀ .* P81(states, parameters) - (γ.+dᵀ.+iᴬ*A) .* T

	dCFP = aᶜ .* P76(states, parameters) - (γ.+dᶜ) .* CFP
	dYFP = aʸ .* P81(states, parameters) - (γ.+dʸ) .* YFP

	dluxI = P81(states, parameters) - (γ.+dk₆) .* luxI
	dlasI = P76(states, parameters) - (γ.+dk₁₂) .* lasI

	# calculate laplacian
	dc₆, dc₁₂ = zero(c₆), zero(c₁₂)
	if isa(c₆, Array)

		for i=2:length(c₆)-1
			dc₆[i] = D₆ * (c₆[i-1]-2*c₆[i]+c₆[i+1]) / Δx^2
			dc₁₂[i] = D₁₂ * (c₁₂[i-1]-2*c₁₂[i]+c₁₂[i+1]) / Δx^2
		end

		# with zero-flux boundaries
		dc₆[1] = D₆ * (c₆[2]-c₆[1]) / Δx^2
		dc₁₂[1] = D₁₂ * (c₁₂[2]-c₁₂[1]) / Δx^2

		dc₆[end] = D₆ * (c₆[end-1]-c₆[end]) / Δx^2
		dc₁₂[end] = D₁₂ * (c₁₂[end-1]-c₁₂[end]) / Δx^2

		# relay reaction
		dc₆ .+= ρ .* k₆ .* luxI
		dc₁₂ .+= ρ .* k₁₂ .* lasI
	end

	return Dict(
		"R" => dR, "S" => dS, "L" => dL,
		"T" => dT, "c₆" => dc₆, "c₁₂" => dc₁₂,
		"CFP" => dCFP, "YFP" => dYFP, "ρ" => dρ,
		"luxI" => dluxI, "lasI" => dlasI
	)
end

function rates( states::Array, parameters::Dict, t::Float64 )
	states = Dict(

		# growth
		"ρ" => states[:,1],

		# receivers
		"R" => states[:,2],
		"S" => states[:,3],

		# inhibitors
		"L" => states[:,4],
		"T" => states[:,5],

		# signals
		"c₆" => states[:,6],
		"c₁₂" => states[:,7],

		# responses
		"CFP" => states[:,8],
		"YFP" => states[:,9],

		# relay
		"luxI" => states[:,10],
		"lasI" => states[:,11]
	)

	@unpack ρ,R,S,L,T,c₆,c₁₂,CFP,YFP,luxI,lasI = rates(states,parameters,t)
	return [ρ R S L T c₆ c₁₂ CFP YFP luxI lasI]
end

############ jacobian and rate function in logspace used for bifurcation analysis
function rates(u::Array, p::NamedTuple{(:c₆,:c₁₂),Tuple{Float64,Float64}} ; parameters::Dict=θ)
	@unpack c₆,c₁₂ = p
	states = Dict(

		# cell density
		"ρ" => θ["ρ"],

		# receivers
		"R" => 10^u[1],
		"S" => 10^u[2],

		# inhibitors
		"L" => 10^u[3],
		"T" => 10^u[4],

		# signals
		"c₆" => 10^c₆,
		"c₁₂" => 10^c₁₂,

		# responses
		"CFP" => NaN,
		"YFP" => NaN,

		# relay
		"luxI" => NaN,
		"lasI" => NaN
	)

	t = NaN
	@unpack R,S,L,T = rates(states,parameters,t)
	return [R,S,L,T]
end

function jacobian( u::Array, p::NamedTuple{(:c₆,:c₁₂),Tuple{Float64,Float64}} ; parameters::Dict=θ)
	@unpack c₆,c₁₂ = p
	states = Dict(

		# cell density
		"ρ" => θ["ρ"],

		# receivers
		"R" => 10^u[1],
		"S" => 10^u[2],

		# inhibitors
		"L" => 10^u[3],
		"T" => 10^u[4],

		# signals
		"c₆" => 10^c₆,
		"c₁₂" => 10^c₁₂,

		# responses
		"CFP" => NaN,
		"YFP" => NaN,

		# relay
		"luxI" => NaN,
		"lasI" => NaN
	)

	@unpack r,K, D₆,D₁₂, aᴿ,aˢ,aᴸ,aᵀ, nᵀ,nᴸ, dᴿ,dˢ,dᴸ,dᵀ, I,A, iᴵ,iᴬ = parameters
	@unpack R, S, L, T, c₆, c₁₂, ρ = states

	γ = r .* (1 .- ρ./K )
	J = zeros(4,4)

	J[1, 1] = -(γ+dᴿ)
	J[2, 2] = -(γ+dˢ)

	J[3, 3] = -(γ+dᴸ+iᴵ*I)
	J[4, 4] = -(γ+dᵀ+iᴬ*A)

	J[1, 4] = aᴿ * ∂hill(T,nᵀ)
	J[2, 3] = aˢ * ∂hill(L,nᴸ)

	J[3, 1] = aᴸ * ∂R76(states, parameters)
	J[3, 2] = aᴸ * ∂S76(states, parameters)

	J[4, 1] = aᵀ * ∂R81(states, parameters)
	J[4, 2] = aᵀ * ∂S81(states, parameters)

	return J * Diagonal(log(10) .* [R, S, L, T] )
end
