using DifferentialEquations
include("lib/bifurcations.jl")
include("lib/animate.jl")

############################# import model with parameters
include("models/protected-degradation.jl")

############################# simulation grid
n_points,x_max,t_final = 101,0.016,24.0
	space = range(0,x_max,length=n_points)
	Δx, Δt = step(space), 0.007

############################# initial conditions
u0 = zeros(n_points,11)

	# cell density
	θ["ρ"] = 0.002
	u0[:,1] .= θ["ρ"]

	# signaling [nM]
	u0[space.-x_max/2 .> x_max/4, 6] .= 1000 # c6
	u0[space .< x_max/4, 7] .= 1000 # c12

	# derepressors
	θ["A"] = 0.0 # atc [ng/ml]
	θ["I"] = 0.0 # iptg [mM]

	#relay
	# θ["k₆"] = 300.0
	# θ["k₁₂"] = 186.113764100785

############################# simulate
problem = ODEProblem( rates, u0, (0.0,t_final), θ )
	solution = solve( problem, Euler(), dt=Δt )

############################# bifurcation analysis
bifurcations = get_bifurcations( solution, rates, jacobian, frames=15 )

#################### generate animation simulate
generate_animation(solution, "output.gif", bifurcations = bifurcations, fps=5, frames=15)
