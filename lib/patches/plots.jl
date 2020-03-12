using Plots, Plots.PlotMeasures, LaTeXStrings
font = Plots.font("Arial", 24)
pyplot(guidefont=font, xtickfont=font, ytickfont=font, legendfont=font)

# patches to make logscales work... damn it Julia
import Base: iterate,firstindex,lastindex,getindex
iterate(x::Surface{Array{Float64,2}}) = iterate(x.surf)
iterate(x::Surface{Array{Float64,2}}, i::Int64) = iterate(x.surf,i)
firstindex(x::Surface{Array{Float64,2}}) = firstindex(x.surf)
lastindex(x::Surface{Array{Float64,2}}) = lastindex(x.surf)
getindex(x::Surface{Array{Float64,2}}, i::UnitRange{Int64}) = getindex(x.surf,i)
