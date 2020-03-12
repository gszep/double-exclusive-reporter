using Test


function finiteDifferences(F, x::AbstractVector; δ = 1e-9)
	f = F(x)
	N = length(x)
	J = zeros(eltype(f), N, N)
	x1 = copy(x)
	for i=1:N
		x1[i] += δ
		J[:, i] .= (F(x1) .- F(x)) / δ
		x1[i] -= δ
	end
	return J
end


function test_jacobian(rates, jacobian, n=10000)
	unit_test = fill(false,n)

	for i=1:n

		u,c6,c12 = randn(4),randn(),randn()
		approximation = finiteDifferences( x->rates(x,c6,c12), u )
		tolerance = 1e-5 * norm(jacobian(u,c6,c12))

		unit_test[i] = all(
				isapprox.( jacobian(u,c6,c12), approximation,
			atol=tolerance)
		)
	end

	@assert all(unit_test)
end
