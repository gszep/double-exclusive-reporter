using DifferentialEquations


function rates(u; parameters=θ)
	states = Dict(

		# receivers
		"R" => u[1],
		"S" => u[2],

		# inhibitors
		"L" => u[3],
		"T" => u[4],

		# signals
		"c₆" => u[5],
		"c₁₂" => u[6],

	)

	@unpack R,S,L,T,c₆,c₁₂ = F(states,parameters)
	return [R S L T c₆ c₁₂]
end


rates([[0.0] [0.1] [0.2] [3.] [4.] [1.]])

fff

function jacobian(u,c₆,c₁₂; parameters=θ)
	states = Dict(

		# receivers
		"R" => 10^u[1],
		"S" => 10^u[2],

		# inhibitors
		"L" => 10^u[3],
		"T" => 10^u[4],

		# signals
		"c₆" => 10^c₆,
		"c₁₂" => 10^c₁₂

	)

	@unpack γ,D₆,D₁₂, aᴿ,aˢ,aᴸ,aᵀ, nᵀ,nᴸ, dᴿ,dˢ,dᴸ,dᵀ, I,A, iᴵ,iᴬ = parameters
	@unpack R, S, L, T, c₆, c₁₂ = states
	J = zeros(4,4)

	J[1, 1] = -(γ+dᴿ)
	J[2, 2] = -(γ+dˢ)

	J[3, 3] = -(γ+dᴸ+iᴵ*I)
	J[4, 4] = -(γ+dᵀ+iᴬ*A)

	J[1, 4] = aᴿ * ∂hill(T,nᵀ)
	J[2, 3] = aˢ * ∂hill(L,nᴸ)

	J[3, 1] = aᴸ * ∂R76(states, parameters)
	J[4, 2] = aᵀ * ∂S81(states, parameters)

	return J * Diagonal(log(10) .* [R, S, L, T] )
end

function get_bifurcations(F, J;
	u=[-2.45,-1.13, 0.84, 1.37],
	c₆=0.0, c₁₂=5.0 )

	newtonOptions = Cont.NewtonPar(maxIter = 300, tol = 1e-10)
	parameters = Cont.ContinuationPar(
		dsmin=0.0001, ds=0.001, dsmax=0.001,
		pMin=-10.0, pMax=10.0, detect_fold=true,
		maxSteps=10000, newtonOptions = newtonOptions)
	fold_parameters = ContinuationPar(
		dsmin = 0.001, ds=-0.001, dsmax = 0.005,
		pMax=10.0, pMin=-10.0, maxSteps=5000,
		newtonOptions = Cont.NewtonPar(tol = 1e-4, maxIter = 300) )

	println("1. searching for initial equilibrium...")
	u, _, converged = Cont.newton(
		u -> F(u, c₆,c₁₂), u -> J(u, c₆,c₁₂),
		u, newtonOptions)

	if converged
		println("2. equilibrium curve along c6...")
		equilibrium_curve, _ = Cont.continuation(
			(u, p) -> F(u, p, c₁₂), (u, p) -> J(u, p, c₁₂),
			u, c₆, parameters, printsolution=x->x[1], plot=false)

		if length(equilibrium_curve.bifpoint) > 0
			printstyled("3. limit points ", color=:green)

			limit_curve, _, success = Cont.continuationFold(
				(u, α, β) -> F(u, α, β), (u, α, β) -> J(u, α, β),
				equilibrium_curve, 1, c₁₂, fold_parameters, plot=false)

			return u,limit_curve

		else
			println("no limit points\n")
			return u,nothing
		end
	else
		println("cannot find initial equilibrium\n")
		return [NaN,NaN,NaN,NaN],nothing
	end
end
