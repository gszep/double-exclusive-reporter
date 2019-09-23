using PseudoArcLengthContinuation, LinearAlgebra, Plots
const Cont = PseudoArcLengthContinuation

function rates(u, β1 = 0, β2 = 0)
	f = similar(u)
	f[1] =  β1+β2 + (β1-β2)*u[1] + u[1]^3
	return f
end

function Jacobian(u, β1 = 0, β2 = 0)
	J = zeros( length(u),  length(u))
	J[1,1] = (β1-β2) + 3*u[1]^2
	return J
end

β1,β2 = -0.5,0.75
initial_condition = [0.0]

# integrate until steady state
steady_state, trajectory, success = @time Cont.newton(
	u -> rates(u,β1,β2),
	initial_condition,
	Cont.NewtonPar() )

if success

	# find fold points
	branch, _, _ = @time Cont.continuation(
		(u,β1) -> rates(u,β1,β2),
		steady_state, β1,
		Cont.ContinuationPar(pMax = 2, pMin = -2))

	# locate first fold precisely
	fold_index = 1
	# fold_points, _, success = @time Cont.newtonFold(
	# 	(u,β1) -> rates(u,β1,β2), (u,β1) -> Jacobian(u,β1,β2),
	# 	branch, fold_index, NewtonPar())

	Cont.plotBranch(branch)

	# continue on fold curve to find cusp
	cusp, _, _ = @time Cont.continuationFold(
		(u,β1,β2) -> rates(u,β1,β2), (u,β1,β2) -> Jacobian(u,β1,β2),
		FoldPoint(branch,fold_index), β2, branch.bifpoint[fold_index][end-1],
		ContinuationPar(pMax = 2, pMin = -2))

	Cont.plotBranch(cusp)

else
	raise("[error] integration to steady state failed")
end
