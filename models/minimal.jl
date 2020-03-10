using Parameters: @unpack

# states
N,W = 101,0.016
space = range(0,W,length=N)

Δt = 0.007
Δx = step(space)

# model states and parameters
U = Dict(

	# receivers
	"R" => zeros(N),
	"S" => zeros(N),

	# inhibitors
	"L" => zeros(N),
	"T" => zeros(N),

	# signals
	"c₆" => zeros(N),
	"c₁₂" => zeros(N)

)
θ = Dict(

	# growth
	"γ" => 1.0,

	# signal diffusion
	"D₆" => 0.0000018,
	"D₁₂" => 0.0000009,

	# dissociation constants
	"Kᴿ" => 0.0012,
	"Kˢ" => 0.4,

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

	# degredations
	"dᴿ" => 24.6,
	"dˢ" => 8.58,
	"dᴸ" => 1.0,
	"dᵀ" => 1.0,

	# derepressors
	"iᴬ" => 0.421,
	"iᴵ" => 418.0,
	"A" => 0.0,
	"I" => 0.0
)

bound(x, cₓ, nₓ)  = x .^ 2 .* cₓ .^ nₓ
∂bound(x, cₓ, nₓ) = 2x .* cₓ .^ nₓ

hill(x, n)  = 1 ./ (1 .+ x.^n)
∂hill(x, n) = -n .* x.^(n.-1) ./ (1 .+ x.^n).^2

function P76( states = U, parameters = θ)
	@unpack nᴿ,Kᴿ = parameters
	@unpack R, c₆ = states

	activation = Kᴿ .* bound( R, c₆, nᴿ )
	return  activation ./ (1 .+ activation)
end

function ∂R76( states = U, parameters = θ)
	@unpack nᴿ,Kᴿ = parameters
	@unpack R, c₆ = states

	activation = Kᴿ .* bound( R, c₆, nᴿ )
	return Kᴿ .* ∂bound( R, c₆, nᴿ ) ./ (1 .+ activation).^2
end

function P81( states = U, parameters = θ)
	@unpack nˢ,Kˢ = parameters
	@unpack S, c₁₂ = states

	activation = Kˢ .* bound( S, c₁₂, nˢ )
	return activation ./ (1 .+ activation)
end

function ∂S81( states = U, parameters = θ)
	@unpack nˢ,Kˢ = parameters
	@unpack S, c₁₂ = states

	activation = Kˢ .* bound( S, c₁₂, nˢ )
	return Kˢ .* ∂bound( S, c₁₂, nˢ ) ./ (1 .+ activation).^2
end

function F( states = U, parameters = θ)
	@unpack γ,D₆,D₁₂, aᴿ,aˢ,aᴸ,aᵀ, nᵀ,nᴸ, dᴿ,dˢ,dᴸ,dᵀ, I,A, iᴵ,iᴬ = parameters
	@unpack R, S, L, T, c₆, c₁₂ = states

	# calculate reaction rates
	dR = aᴿ .* hill(T,nᵀ) .- (γ+dᴿ) .* R
	dS = aˢ .* hill(L,nᴸ) .- (γ+dˢ) .* S
	dL = aᴸ .* P76(states, parameters) - (γ+dᴸ+iᴵ*I) .* L
	dT = aᵀ .* P81(states, parameters) - (γ+dᵀ+iᴬ*A) .* T

	# calculate laplacian with zero-flux boundaries
	dc₆, dc₁₂ = zero(c₆), zero(c₁₂)
	if isa(c₆, Array)

		for i=2:N-1
			dc₆[i] = (c₆[i-1]-2*c₆[i]+c₆[i+1]) / Δx^2
			dc₁₂[i] = (c₁₂[i-1]-2*c₁₂[i]+c₁₂[i+1]) / Δx^2
		end

		# update states
		R .+= dR*Δt; S .+= dS*Δt;
		L .+= dL*Δt; T .+= dT*Δt;
		c₆ .+= dc₆*Δt; c₁₂ .+= dc₁₂*Δt
	end

	return Dict(
		"R" => dR, "S" => dS, "L" => dL,
		"T" => dT, "c₆" => dc₆, "c₁₂" => dc₁₂
	)
end
