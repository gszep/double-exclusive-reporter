using PseudoArcLengthContinuation
const Cont = PseudoArcLengthContinuation

# algorithm for getting bistability region
function get_bifurcations(rates::Function, jacobian::Function;
	u::Vector{Float64}=[-2.45,-1.13, 0.84, 1.37],
	c₆::Float64=0.0, c₁₂::Float64=5.0 )

	# continuation parameters; if the proceedure fails try modifying these
	codim1_parameters = Cont.ContinuationPar( detect_fold=true,
		newtonOptions = Cont.NewtonPar(maxIter = 300, tol = 1e-10),
		dsmin=0.0001, ds=0.001, dsmax=0.001,
		pMin=-10.0, pMax=10.0, maxSteps=10000)

	codim2_parameters = Cont.ContinuationPar(
		newtonOptions = Cont.NewtonPar(tol = 1e-4, maxIter = 300),
		dsmin = 0.001, ds=-0.001, dsmax = 0.005,
		pMax=10.0, pMin=-10.0, maxSteps=5000)

	# search for initial equilibrium, u
	u, _, converged = Cont.newton(
		u -> rates(u, c₆,c₁₂), u -> jacobian(u, c₆,c₁₂),
		u, Cont.NewtonPar(maxIter = 300, tol = 1e-10) )

	if converged # continue equilibrium curve along c6 for fixed c12
		equilibrium_curve, _ = Cont.continuation(
			(u, p) -> rates(u, p, c₁₂), (u, p) -> jacobian(u, p, c₁₂),
			u, c₆, codim1_parameters, printsolution=x->x[1], plot=false)

		if length(equilibrium_curve.bifpoint) > 0
			# continue limit point in (c12,c6) region

			limit_curve, _, success = Cont.continuationFold(
				(u, α, β) -> rates(u, α, β), (u, α, β) -> jacobian(u, α, β),
				equilibrium_curve, 1, c₁₂, codim2_parameters, plot=false)
			return u,limit_curve

		else
			throw("no limit points\n")
		end
	else
		throw("cannot find initial equilibrium\n")
	end
end

# get bistability regions for each point in time
function get_bifurcations( solution::ODESolution, rates::Function, jacobian::Function;
	u::Vector{Float64}=[-2.45,-1.13, 0.84, 1.37], c₆::Float64=0.0, c₁₂::Float64=5.0,
	frames::Int=15 )

	# iterate through timepoints in solution
	bifurcations = []
	for i ∈ range(1,length(solution.u),length=frames)
		i = Int(floor(i))

		# pass cell density as parameter
		ρ, = eachcol(solution.u[i])
		θ["ρ"] = mean(ρ)

		try # perform continuation
			u,limit_curve = get_bifurcations(rates, jacobian; u = u, c₆=c₆, c₁₂=c₁₂)
			bistable_region = Shape(10 .^ limit_curve.branch[1,:], 10 .^ limit_curve.branch[2,:])
			push!(bifurcations,bistable_region)

		catch
			bistable_region = Shape([NaN],[NaN])
			push!(bifurcations,bistable_region)
		end
	end
	return bifurcations
end
