using PseudoArcLengthContinuation
const Cont = PseudoArcLengthContinuation
using LinearAlgebra: Diagonal

include("patches/plots.jl")

#using DelimitedFiles,ColorSchemes
#include("lib/parameters.jl")

function rates(u,c₆,c₁₂; parameters=θ)
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

	@unpack R,S,L,T = F(states,parameters)
	return [R,S,L,T]
end

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



include("models/minimal.jl")
	# yfp = reshape(readdlm("primed/dat2yfp.dat",'\t'), (12,8,3))[end:-1:1,end:-1:1,:]
	# cfp = reshape(readdlm("primed/dat2cfp.dat",'\t'), (12,8,3))[end:-1:1,end:-1:1,:]
	# heatmap!( unique(yfp[:,:,1]),unique(yfp[:,:,2]), yfp[:,:,3]'-cfp[:,:,3]',
	# 	c=cgrad([:cyan, :white, :yellow]), colorbar=true)


	x,limit_curve = get_bifurcations(rates, jacobian)
	if limit_curve != nothing

		plot(xlims=(0.25,5e4), ylims=(0.75,5e4), label="", size=(700,700), grid=false,
			xscale=:log, yscale=:log, xlabel="[C12] (nM)",ylabel="[C6] (nM)")
		plot!(10 .^ limit_curve.branch[1,:], legend=:topleft, label="",
			10 .^ limit_curve.branch[2,:], color=:red,linewidth=3) |> display
	end

	# plot!([25000.0  5000.0  1000.0  200.0  40.0  8.0  1.6], seriestype = :hline, label="")
	# plot!([25000.0  8333.0  2777.0  925.0  308.0 102.0 34.0  11.0 3.8 1.3 0.4], seriestype = :vline, label="")
	#

	# parameters["growth"] = 0.5
	# parameters["capacity"] = 1.0
	# x,limit_curve = get_bifurcations(parameters, Flog, Jlog)
	# if limit_curve != nothing
	#
	# 	plot!(10 .^ limit_curve.branch[1,:], label="0.5", legend=:topleft,
	# 		10 .^ limit_curve.branch[2,:], color=:blue,linewidth=3)
	# end

	# parameters["growth"] = 0.1
	# parameters["capacity"] = 1.0
	# x,limit_curve = get_bifurcations(parameters, Flog, Jlog)
	# if limit_curve != nothing
	#
	# 	plot!(10 .^ limit_curve.branch[1,:], label=L"\gamma_0 = 0.1", legend=:topleft,
	# 		10 .^ limit_curve.branch[2,:], color=:gold,linewidth=3)
	# end
	png("main-figure")





ddd






# version = 2
# if version == 1
# 	P0 = Main.Version1.P0
# 	Flog = Main.Version1.Flog
# 	Jlog = Main.Version1.Jlog
#
# 	seeds = range(0,stop=19) |> collect
# 	filter!(e->e∉[16],seeds)

for (i,seed) in enumerate([19])
	println("computing bifurcation for seed ", seed)

	# P1 = merge(P0,parseSummary(@sprintf("parameters/v%d/summary%d.txt", version, seed)))
	#
	# P1["iptg"] = 0.0
	# x,limit_curve = get_bifurcations(parameters, Flog, Jlog)
	# if limit_curve != nothing
	#
	# 	plot!(10 .^ limit_curve.branch[1,:], label="", alpha=1.0,
	# 		10 .^ limit_curve.branch[2,:], color=:red,linewidth=3)
	# end

	# P1["iptg"] = 0.001
	# limit_curve = get_bifurcations(P1, Flog, Jlog)
	# if limit_curve != nothing
	#
	# 	plot!(limit_curve.branch[1,:], label="", alpha=0.5,
	# 		limit_curve.branch[2,:], color=:green ) |> display
	# end
	#
	# plot!([NaN],[NaN],color=:green,legend=:topleft,label="iptg = 0.002")

end
png("main-figure")

P1["iptg"] = 0.0
plot(xlims=(1.6,2.5e4), ylims=(1.6,2.5e4),
	label="", size=(700,700),
	xscale=:log, yscale=:log,
	xlabel="C12 (nM)",ylabel="C6 (nM)")

	flowdata = readdlm("flowplot.dat",'\t')
	flowdata = reshape(flowdata, (12,8,3))
	heatmap!(unique(flowdata[:,:,1]),unique(flowdata[:,:,2]),colorbar=false,grid=:top,
		flowdata[:,:,3], c=cgrad(ColorSchemes.Greys_5.colors))

	for (i,seed) in enumerate([19])
		x = :none
		println("computing bifurcation for seed ", seed)
		P1 = merge(P0,parseSummary(@sprintf("inference_results/v%d/summary%d.txt", version, seed)))

		P1["atc"] = 0.0
		x,limit_curve = get_bifurcations(P1, Flog, Jlog; x=x)
		if limit_curve != nothing
			plot!(10 .^ limit_curve.branch[1,:], alpha=1.0,label="ATC = 0.0 mM",
				10 .^ limit_curve.branch[2,:], color=:red, linewidth=3 )
		end

		x,limit_curve = get_bifurcations(P1, Flog, Jlog; x=x, ds = 1.0)
		if limit_curve != nothing
			plot!(10 .^ limit_curve.branch[1,:], alpha=1.0,label="",
				10 .^ limit_curve.branch[2,:], color=:red, linewidth=3 )
		end

		for iptg=10 .^ (-1:0.05:0.5)
			P1["atc"] = iptg

			x,limit_curve = get_bifurcations(P1, Flog, Jlog; x=x)
			if limit_curve != nothing

				plot!(10 .^ limit_curve.branch[1,:], alpha=0.5,label="",
					10 .^ limit_curve.branch[2,:], color="gray", linewidth=2 )
			end

			x,limit_curve = get_bifurcations(P1, Flog, Jlog; x=x, ds = 1.0)
			if limit_curve != nothing

				plot!(10 .^ limit_curve.branch[1,:], alpha=0.5,label="",
					10 .^ limit_curve.branch[2,:], color="gray", linewidth=2 )
			end
		end

		P1["atc"] = 10^(0.5)
		x,limit_curve = get_bifurcations(P1, Flog, Jlog; x=x)
		if limit_curve != nothing
			plot!(10 .^ limit_curve.branch[1,:], alpha=1.0,label="ATC = 3.16 mM",
				10 .^ limit_curve.branch[2,:], color=:gold, linewidth=3 )
		end

		x,limit_curve = get_bifurcations(P1, Flog, Jlog; x=x, ds = 1.0)
		if limit_curve != nothing

			plot!(10 .^ limit_curve.branch[1,:], alpha=1.0,label="",
				10 .^ limit_curve.branch[2,:], color=:gold, linewidth=3 )
		end
	end

png("main-figure")

P1["iptg"]
P1["atc"] = 0.0
x,limit_curve = get_bifurcations(P1, Flog, Jlog; x=x)
if limit_curve != nothing
	plot!(10 .^ limit_curve.branch[1,:], alpha=1.0,label="ATC = 3.16 mM",
		10 .^ limit_curve.branch[2,:], color=:gold, linewidth=3 )
end

x,limit_curve = get_bifurcations(P1, Flog, Jlog; x=x, ds = 1.0)
if limit_curve != nothing

	plot!(10 .^ limit_curve.branch[1,:], alpha=1.0,label="",
		10 .^ limit_curve.branch[2,:], color=:gold, linewidth=3 )
end
png("main-figure")


plot(xlims=(1.6,2.5e4), ylims=(1.6,2.5e4),
	label="", size=(700,700),
	xscale=:log, yscale=:log,
	xlabel="[C12] (nM)",ylabel="[C6] (nM)")

	flowdata = readdlm("flowplot.dat",'\t')
	flowdata = reshape(flowdata, (12,8,3))
	# heatmap!(unique(flowdata[:,:,1]),unique(flowdata[:,:,2]),colorbar=false,grid=:top,
	# 	flowdata[:,:,3], c=cgrad(ColorSchemes.Greys_5.colors))

	plot!([NaN],[NaN],color=:red,legend=:topleft,label="Model",linewidth=3,grid=false)


	println("computing bifurcation for seed ", 19)

	P1 = merge(P0,parseSummary(@sprintf("inference_results/v%d/summary%d.txt", version, 19)))

	P1["iptg"] = 0.0
	x,limit_curve = get_bifurcations(P1, Flog, Jlog)
	if limit_curve != nothing

		plot!(10 .^ limit_curve.branch[1,:], label="", alpha=1.0,
			10 .^ limit_curve.branch[2,:], color=:red,linewidth=3)
	end
png("main-figure")
println(P1)

using NPZ
npzwrite("bifurcations.npz", transpose(limit_curve.branch[1:2,:]))
