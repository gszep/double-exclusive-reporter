using PseudoArcLengthContinuation, LinearAlgebra, Plots
const Cont = PseudoArcLengthContinuation

function limit_curve( rates, Jacobian, u, β1,β2 ; kwargs...)
	"""computes limit curve using parameter continuation methods"""
	FOLD_INDEX = 1

	# integrate until steady state
	steady_state, _, stable = @time Cont.newton(
		u -> rates(u,β1,β2), u, Cont.NewtonPar() )

	#################################################################################
	if stable # find fold points
		branches, _, _ = Cont.continuation( (u,β1) -> rates(u,β1,β2), steady_state, β2, ContinuationPar(maxSteps=0))

		# forward/backward branches
		for σ ∈ [-1,1]
			branch, _, _ = Cont.continuation(
				(u,β1) -> rates(u,β1,β2),

				steady_state, β2,
				ContinuationPar(ds=σ*0.01; kwargs...))

			append!( branches.branch, branch.branch )
			append!( branches.bifpoint, branch.bifpoint )
		end

		#################################################################################
		if length(branches.bifpoint)>0 # continue on fold curve to find cusp
			β = β2*ones(length(branches.branch[1,:]))
			cusp, _, _ = Cont.continuationFold( rates, Jacobian, branches, FOLD_INDEX, β2, ContinuationPar(maxSteps=0))

			# forward/backward branches
			for σ ∈ [-1,1]
				branch, _, _ = Cont.continuationFold(
					rates, Jacobian,

					branches, FOLD_INDEX, β2,
					ContinuationPar(ds=σ*0.01;  kwargs...))

				append!( cusp.branch, branch.branch )
				append!( cusp.bifpoint, branch.bifpoint )

			end

		#################################################################################
		else
			printstyled(color=:red,"[error] no limit points found\n")
			return
		end
	else
		printstyled(color=:red,"[error] integration to steady state failed\n")
		return
	end
	return cusp
end

function rates(u,θ)
	f = similar(u)
	f[1] =  θ[2]-θ[1] + (θ[2]+θ[1])*u[1] - u[1]^3
	f[2] =  θ[2]-θ[1] + (θ[2]+θ[1])*u[2] - u[2]^3
	return f
end

function Jacobian(u,θ)
	J = zeros( length(u),  length(u))
	J[1,1] =  (θ[2]+θ[1]) - 3u[1]^2
	J[2,2] =  (θ[2]+θ[1]) - 3u[2]^2
	return J
end

θ = [1.0,-0.0]
curve = @time limit_curve(
	(u,θ1,θ2) -> rates(u,[θ1,θ2,θ[3:end]... ]),
	(u,θ1,θ2) -> Jacobian(u,[θ1,θ2,θ[3:end]... ]),
	[0.0,0.0],1.0,1.0; maxSteps=2000, pMin=-2.0, pMax=2.0 )

plot(curve.branch[2,:],curve.branch[1,:],
	label="fold curve",color="green")
