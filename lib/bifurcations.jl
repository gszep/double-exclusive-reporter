using PseudoArcLengthContinuation,Parameters,Setfield
const PALC = PseudoArcLengthContinuation
# algorithm for getting bistability region
function get_bifurcations(rates::Function, jacobian::Function;
	u::Vector{Float64}=[-2.27, -1.27, 1.15, 1.26],
	p::NamedTuple{(:c₆,:c₁₂),Tuple{Float64,Float64}}=(c₆=0.0,c₁₂=5.0) )

	# continuation parameters; if the proceedure fails try modifying these
	codim1_parameters = ContinuationPar( detectFold=true,
		newtonOptions = NewtonPar(maxIter = 300, tol = 1e-10),
		dsmin=0.001, ds=0.01, dsmax=0.01, a=1.0,
		pMin=-10.0, pMax=10.0, maxSteps=1000)

	codim2_parameters = ContinuationPar(
		newtonOptions = NewtonPar(maxIter = 300, tol = 1e-4),
		dsmin = 0.001, ds=-0.001, dsmax = 0.01, a=1.0,
		pMax=10.0, pMin=-10.0, maxSteps=2000)

	# search for initial equilibrium, u
	u, _, converged = newton( rates, jacobian, u, p,
		NewtonPar(maxIter = 300, tol = 1e-10, verbose = false) )

	if converged # continue equilibrium curve along c6 for fixed c12
		equilibrium_curve, _ = continuation( rates, jacobian, u, p, (@lens _.c₆),
			codim1_parameters; printSolution = (x, p) -> x[3], verbosity=0)

		if length(equilibrium_curve.foldpoint) > 0
			# continue limit point in (c12,c6) region

			limit_curve, _, success = continuation( rates, jacobian,
				equilibrium_curve, 1, p, (@lens _.c₆), (@lens _.c₁₂),
				codim2_parameters)

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
	u::Vector{Float64}=[-2.27, -1.27, 1.15, 1.26], p::NamedTuple{(:c₆,:c₁₂),Tuple{Float64,Float64}}=(c₆=0.0,c₁₂=5.0),
	frames::Int=15 )

	# iterate through timepoints in solution
	bifurcations = []
	for i ∈ range(1,length(solution.u),length=frames)
		i = Int(floor(i))

		# pass cell density as parameter
		ρ, = eachcol(solution.u[i])
		θ["ρ"] = mean(ρ)

		try # perform continuation
			u,limit_curve = get_bifurcations(rates, jacobian; u=u, p=p )
			bistable_region = Shape(10 .^ limit_curve.branch[1,:], 10 .^ limit_curve.branch[2,:])
			push!(bifurcations,bistable_region)

		catch
			bistable_region = Shape([NaN],[NaN])
			push!(bifurcations,bistable_region)
		end
	end
	return bifurcations
end
