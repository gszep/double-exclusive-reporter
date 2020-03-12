# double-exclusive-reporter
Code to produce figures for publication

## Dependencies
[Download and install Julia 1.3.1](https://julialang.org) then run a julia session by typing `julia` in the terminal or command prompt. Install the dependencies with the following commands
```julia
using Pkg
Pkg.add(["DifferentialEquations","Plots","StatBase","LinearAlgebra"])
Pkg.add(PackageSpec(path="https://github.com/rveltz/PseudoArcLengthContinuation.jl"))
```
then clone this repo and run the main script
```bash
git clone https://github.com/gszep/double-exclusive-reporter.git
cd double-exclusive-reporter
julia main.jl
```
