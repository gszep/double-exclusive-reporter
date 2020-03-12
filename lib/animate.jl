include("patches/plots.jl")
using StatsBase

############################# generate animation
function generate_animation( solution::ODESolution, file_name::String; bifurcations=nothing,
	fps=15, frames=15, response_lim=(0,20), signal_lim=(0,3), c6_lim=(0.75,5e4), c12_lim=(0.25,5e4) )

	# iterate through time points
	frames,j = bifurcations == nothing ? frames : length(bifurcations), 1
	animation = @animate for i ∈ range(1,length(solution.u),length=frames)

		i = Int(floor(i))
		_,_,_,_,_,c₆,c₁₂,CFP,YFP = eachcol(solution.u[i])
		t = round(solution.t[i])

############################# real space plot
		real_space = plot([],[], ylim=response_lim, label="",
			xlabel=L"\mathrm{Space}, x\,/\,\mathrm{cm}", ylabel=L"\mathrm{Response}, X\mathrm{FP}\;/\;\mathrm{\mu M}",
			right_margin=30mm, top_margin=5mm, framestyle = :box)
		annotate!([(0.2, 1.0, text("t = $t h",24))])

		right_axis = twinx()
		plot!(right_axis, [],[], ylim=signal_lim, label="",
			ylabel=L"\mathrm{Signal}, \mathrm{C}X\,/\,\mathrm{\mu M}")

	    plot!(100*space, CFP/1e3, linewidth=5, ribbon=(CFP,zero(CFP)), color="#00b0f0", label="CFP", legend=:topleft)
		plot!(100*space, YFP/1e3, linewidth=5, ribbon=(YFP,zero(YFP)), color="#ffc000", label="YFP", legend=:topleft)

		plot!(right_axis,100*space, c₆/1e3,  linewidth=5, color="#000099", label="C6", legend=:topright)
		plot!(right_axis,100*space, c₁₂/1e3, linewidth=5, color="#ff6600", label="C12", legend=:topright)

############################# state space plot
		state_space = plot(xlims=c12_lim, ylims=c6_lim, label="", grid=false, xscale=:log, yscale=:log,
			xlabel=L"\mathrm{Signal}, \mathrm{C12}\,/\,\mathrm{nM}", ylabel=L"\mathrm{Signal}, \mathrm{C6}\,/\,\mathrm{nM}",
			right_margin=30mm, top_margin=5mm, framestyle = :box)

		C12,C6 = round(mean(c₁₂)),round(mean(c₆))
		annotate!([(60, 2.0, text("C12 = $C12 nM  C6 = $C6 nM",24))])

		mask = YFP .< CFP
		scatter!([C12],[C6], marker=:x, color=:black, markersize=15, label="Signal Spatial Average")
		scatter!(c₁₂[mask],c₆[mask], color="#00b0f0", markerstrokewidth=0, markersize=8, label="")
		scatter!(c₁₂[.~mask],c₆[.~mask], color="#ffc000", markerstrokewidth=0, markersize=8, label="")

############################# bistable region
		if bifurcations != nothing
			plot!(bifurcations[j],legend=:topleft, label="Bistable Region", color=:red, alpha=0.5)
			j += 1
		end

############################# final composite plot
		two_plots = @layout [a b]
		plot(real_space, state_space, layout = two_plots, size=(2100,1000),
			left_margin=25mm, bottom_margin=25mm, top_margin=10mm)
	end

	return gif(animation, file_name, fps = fps)
end
