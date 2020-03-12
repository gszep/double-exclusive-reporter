# Interpretation of morphogen gradients
# by a synthetic bistable circuit
Code to produce figures for publication

## Running the Simulation and Bifurcation Analysis
todo

## Modifying the Model and Parameters
If you would like to change the model you can modify `models/protected-degradation.jl`. Two methods need to be modified: the rate function `rates( states::Dict, parameters::Dict, t::Float64 )` and jacobian in logspace `jacobian( u::Array, c₆::Float64, c₁₂::Float64 ; parameters::Dict=θ)`. You are advised to use the testing utility `test/jacobian.jl` to ensure that the symbolic jacobian matches the finite difference approximation of the rates.

The inferred parameters are also contained in `models/protected-degradation.jl`. Simply change or import new sets of parameters to global variable `θ` and re-run the simulation to see the effects

```
θ = Dict(

	# growth
	"ρ" => 0.002,
	"K" => 2.7,
	"r" => 1.0,
  
  ...
)
```

## Installation and Dependencies
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
