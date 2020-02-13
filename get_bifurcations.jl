using PseudoArcLengthContinuation, LinearAlgebra, Plots, Printf, Colors
const Cont = PseudoArcLengthContinuation

using DelimitedFiles,ColorSchemes
include("lib/parameters.jl")
pyplot()

# patches to make logscales work... damn it Julia
import Base: iterate,firstindex,lastindex,getindex
iterate(x::Surface{Array{Float64,2}}) = iterate(x.surf)
iterate(x::Surface{Array{Float64,2}}, i::Int64) = iterate(x.surf,i)
firstindex(x::Surface{Array{Float64,2}}) = firstindex(x.surf)
lastindex(x::Surface{Array{Float64,2}}) = lastindex(x.surf)
getindex(x::Surface{Array{Float64,2}}, i::UnitRange{Int64}) = getindex(x.surf,i)

# main bifurcation loop
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

module Version1
	include("models/version1.jl")
end

module Version2
	include("models/version2.jl")
end

version = 2
if version == 1
	P0 = Main.Version1.P0
	Flog = Main.Version1.Flog
	Jlog = Main.Version1.Jlog

	seeds = range(0,stop=19) |> collect
	filter!(e->e∉[16],seeds)

elseif version == 2
	P0 = Main.Version2.P0
	Flog = Main.Version2.Flog
	Jlog = Main.Version2.Jlog

	seeds = range(0,stop=19) |> collect
	filter!(e->e∉[0,6,11,16],seeds)
end

plot(xlims=(1.6,2.5e4), ylims=(1.6,2.5e4),
	label="", size=(400,400),
	xscale=:log, yscale=:log,
	xlabel="C12 (nM)",ylabel="C6 (nM)")

flowdata = readdlm("flowplot.dat",'\t')
flowdata = reshape(flowdata, (12,8,3))
heatmap!(unique(flowdata[:,:,1]),unique(flowdata[:,:,2]), colorbar=:none,
	flowdata[:,:,3], c=cgrad(ColorSchemes.Greys_5.colors))

plot!([NaN],[NaN],color=:red,legend=:topleft,label="limit curve estimates")
for (i,seed) in enumerate(seeds)
	println("computing bifurcation for seed ", seed)

	P1 = merge(P0,parseSummary(@sprintf("inference_results/v%d/summary%d.txt", version, seed)))

	P1["iptg"] = 0.0
	limit_curve = get_bifurcations(P1, Flog, Jlog)
	if limit_curve != nothing

		plot!(10 .^ limit_curve.branch[1,:], label="", alpha=0.5,
			10 .^ limit_curve.branch[2,:], color=:red )
	end
end

png("rapid-eqilibrium")

# P1["iptg"] = 0.002
# limit_curve = get_bifurcations(P1, Flog, Jlog)
# if limit_curve != nothing
#
# 	plot!(limit_curve.branch[1,:], label="", alpha=0.5,
# 		limit_curve.branch[2,:], color=:green ) |> display
# end

#plot!([NaN],[NaN],color=:green,legend=:topleft,label="iptg = 0.002")
