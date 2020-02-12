using PseudoArcLengthContinuation, LinearAlgebra, Plots, Printf, Colors, NPZ
const Cont = PseudoArcLengthContinuation

module Parameters
	include("lib/parameters.jl")
end
module Version1
	include("models/version1.jl")
end
module Version2
	include("models/version2.jl")
end

function get_bifurcations(P, F, J;
		c6=5.0, c12=5.0, niters=1000)

	newtonOptions = Cont.NewtonPar(maxIter = 1000, tol = 1e-12)
	parameters = Cont.ContinuationPar(
		dsmin=0.001, ds=-0.01, dsmax=0.01,
		pMin=0.0, pMax=6.0, detect_fold=true,
		maxSteps=500, newtonOptions = newtonOptions)

	fold_parameters = ContinuationPar(
		dsmin = 0.001, ds=-0.001, dsmax = 0.005,
		pMax=5.1, pMin=0.0, maxSteps=1000,
		newtonOptions = Cont.NewtonPar(tol = 1e-5, maxIter = 300) )

	converged,i,x = false,0,fill(NaN,4)
	println("1. searching for initial equilibrium...")
	while !converged
		try
			x, _, converged = Cont.newton(
				x -> F(x, P, c6, c12), x -> J(x, P, c6, c12),
				randn(4), newtonOptions)
		catch
			converged = false
		end

		i += 1
		if i > niters
			printstyled("cannot find initial equilibrium", color=:red)
			return
		end
	end

	println("2. equilibrium curve along c6...")
	equilibrium_curve, _ = Cont.continuation(
		(x, p) -> F(x, P, p, c12), (x, p) -> J(x, P, p, c12),
		x, c6, parameters, printsolution=x->x[1], plot=false)

	if length(equilibrium_curve.bifpoint) > 0
		printstyled("3. limit points ", color=:green)

		limit_curve, _, success = Cont.continuationFold(
			(x, α, β) -> F(x, P, α, β), (x, α, β) -> J(x, P, α, β),
			equilibrium_curve, 1, c12, fold_parameters, plot=false)

		println("\n")
		return limit_curve
	else
		println("no limit points\n")
		return
	end
end

version = 1
if version == 1
	P0 = Main.Version1.P0
	Flog = Main.Version1.Flog
	Jlog = Main.Version1.Jlog

	seeds = range(0, stop=19) |> collect
	filter!(e->e∉[16],seeds)

elseif version == 2
	P0 = Main.Version2.P0
	Flog = Main.Version2.Flog
	Jlog = Main.Version2.Jlog

	seeds = range(0,19) |> collect
	filter!(e->e∉[0,6,11,16],seeds)
end

plot(xlims=(0,5), ylims=(0,5), label="",
	xlabel="C12",ylabel="C6", size=(500,500))

plot!([NaN],[NaN],color=:blue,legend=:topleft,label="iptg = 0.0")
plot!([NaN],[NaN],color=:green,legend=:topleft,label="iptg = 0.002")

for (i,seed) in enumerate(seeds)
	println("computing bifurcation for seed ", seed)

	P1 = merge(P0, Main.Parameters.parseSummary(@sprintf("inference_results/v%d/summary%d.txt", version, seed)))

	P1["iptg"] = 0.0
	limit_curve = get_bifurcations(P1, Flog, Jlog)
	if limit_curve != nothing
		plot!(limit_curve.branch[1,:], label="", alpha=0.5,
			limit_curve.branch[2,:], color=:blue ) |> display
	end

	P1["iptg"] = 0.002
	limit_curve = get_bifurcations(P1, Flog, Jlog)
	if limit_curve != nothing
		plot!(limit_curve.branch[1,:], label="", alpha=0.5,
			limit_curve.branch[2,:], color=:green ) |> display
	end

	# Write branch to file
	# npzwrite(@sprintf("bifurcation_results/branches%d_seed%d.npz", version, seed), transpose(outfoldco.branch[1:2,:]))
end
#npzwrite(@sprintf("bifurcation_results/cusps%d.npz", version), cusps)
