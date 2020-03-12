# double-exclusive-reporter
Code to produce figures for publication

## Modifying Parameters
The inferred model parameters are contained in `models/protected-degradation.jl`. Simply change or import new sets of parameters and re-run the simulation to see the effects

```
θ = Dict(

	# growth
	"ρ" => 0.002,
	"K" => 2.7,
	"r" => 1.0,
  
  ...
)
```

## Dependencies
[Download and install Julia 1.3.1](https://julialang.org) then run a julia session by typing `julia` in the terminal or command prompt. Install the dependencies - including the library `PseudoArcLengthContinuation.jl` written by R.Veltz - with the following commands
```julia
using Pkg
Pkg.add(["DifferentialEquations","Plots","LaTeXStrings","StatsBase","LinearAlgebra","Parameters"])
Pkg.add(PackageSpec(path="https://github.com/rveltz/PseudoArcLengthContinuation.jl"))
```
then clone this repo and run the main script
```bash
git clone https://github.com/gszep/double-exclusive-reporter.git
cd double-exclusive-reporter
julia main.jl
```
