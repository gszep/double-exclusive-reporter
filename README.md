# Interpretation of morphogen gradients by a bistable circuit

## Running the Simulation and Bifurcation Analysis
The parametrised differential equation model used to describe the genetic circuit is defined by the `rates` and `jacobian` functions in `models/protected-degradation.jl`. The parameters are contained in a dictionary `θ` which is then unpacked by symbol name in the two functions. Changing the model requires changing `θ`, `rates` and `jacobian`.

The initial conditions `u0` and other simulation parameters such as `t_final` can be changed in the simulation script `main.jl`. The `DifferentialEquations.jl` library is used to simulate the model. This is followed by bifurcation analysis for `15 frames` of the entire length of the simulation. The bifurcation algorithm parameters are contained in `get_bifurcations` in `lib/bifurcations.jl`

Running `main.jl` will produce an animation, the layout of which is defined in `lib/animate.jl`

![](_simulation.gif)

## Fluorescence Microscopy Data Processing
All microscopy data processing as described in the supplementary was done using `get_boundary.py` and the `bulk-processing.ipynb` notebook. The input data are large `TIFF` hyperstacks and are not uploaded here; but can be provided upon request
![](_kymographs.png)

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
Pkg.add(["DifferentialEquations","Plots","PyPlot","LaTeXStrings","StatsBase","LinearAlgebra","Parameters"])
Pkg.add(PackageSpec(path="https://github.com/rveltz/PseudoArcLengthContinuation.jl"))
```
then clone this repo and run the main script
```bash
git clone https://github.com/gszep/double-exclusive-reporter.git
cd double-exclusive-reporter
julia main.jl
```
