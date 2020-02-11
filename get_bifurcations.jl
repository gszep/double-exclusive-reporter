using PseudoArcLengthContinuation, LinearAlgebra, Plots, Printf, Colors, NPZ
const Cont = PseudoArcLengthContinuation

module Version1
	include("models/version1.jl")
end
module Version2
	include("models/version2.jl")
end
module Parameters
	include("lib/parameters.jl")
end

function evaluateCusp(P, F, J, x0; plot=false, c6=5.0, c12=5.0, ds=-0.1)
	newtonOptions = Cont.NewtonPar(maxIter = 200, tol = 1e-12)
	eq, hist, flag = Cont.newton(x -> F(x, P, c6, c12), x -> J(x, P, c6, c12), x0, newtonOptions)
		println("- Equilibrium found at ", eq)

	# Search in the C6 direction
	opts_br6 = Cont.ContinuationPar(dsmin=0.001, dsmax=0.01, ds=ds, pMin=0.0, pMax=6.0, detect_fold=true, maxSteps=250)
		opts_br6.newtonOptions = newtonOptions
	br6, u6 = Cont.continuation(
		(x, p) -> F(x, P, p, c12),
		(x, p) -> J(x, P, p, c12),
		eq, c6, opts_br6, printsolution = x -> x[1], plot=plot)
		#Cont.plotBranch(br6; xlabel="C6", ylabel="[luxR]", label="")

	# Compute the fold point more precisely
	indfold = 1
	optcontfold = ContinuationPar(dsmin = 0.001, dsmax = 0.005, ds=ds, pMax=5.1, pMin=0.0, maxSteps=1500)
		optcontfold.newtonOptions = Cont.NewtonPar(tol = 1e-6, maxIter = 300)
	outfoldco, hist, flag = Cont.continuationFold(
		(x, α, β) -> F(x, P, α, β),
		(x, α, β) -> J(x, P, α, β),
		br6, indfold, c12, optcontfold, plot=plot)
	return outfoldco
end

sol = fill(-1.0,4)
version = 1

if version == 1
	P0 = Main.Version1.P0
	Flog = Main.Version1.Flog
	Jlog = Main.Version1.Jlog
	seeds = [5, 6, 7, 9, 10, 11, 12, 17, 19]
elseif version == 2
	P0 = Main.Version2.P0
	Flog = Main.Version2.Flog
	Jlog = Main.Version2.Jlog
	seeds = [1, 3, 4, 8, 9, 10, 12, 19]
end

# Evaluate one with plotting
P = merge(P0, Main.Parameters.parseSummary(@sprintf("inference_results/v%d/summary%d.txt", version, seeds[2])))
evaluateCusp(P, Flog, Jlog, sol, plot=true)

# Evaluate all in a loop
plot(0,0); xlabel!("C12"); ylabel!("C6");
	xlims!((0.,5.));
	ylims!((0.,5.))
nseeds = length(seeds)
cusps = zeros(nseeds, 3)
for i in range(1, length = nseeds)
	seed = seeds[i]
	println("Computing bifurcation for seed ", seed)
	#i = seeds[2]
	P1 = merge(P0, Main.Parameters.parseSummary(@sprintf("inference_results/v%d/summary%d.txt", version, seed)))
	outfoldco = evaluateCusp(P1, Flog, Jlog, sol, ds=-0.001)
	#p = Cont.plotBranch!(outfoldco; label=@sprintf("Seed %d", i), legend=:topleft)
	bifpt = outfoldco.bifpoint[1]

	# Update cusps
	cusps[i,:] = [seed, bifpt.param, bifpt.printsol]

	# Write branch to file
	npzwrite(@sprintf("bifurcation_results/branches%d_seed%d.npz", version, seed), transpose(outfoldco.branch[1:2,:]))

	# Plot
	plot!(outfoldco.branch[1,:], outfoldco.branch[2,:], label=@sprintf("Seed %d", seed), legend=:topleft)
	p = scatter!([bifpt.param], [bifpt.printsol], color=:black, label="")
	display(p)
end
	npzwrite(@sprintf("bifurcation_results/cusps%d.npz", version), cusps)
