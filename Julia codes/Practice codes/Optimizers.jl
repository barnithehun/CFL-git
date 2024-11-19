using Optim
using Optim: converged, maximum, maximizer, minimizer, iterations #some extra functions

result = optimize(x -> -(x^2), -2.0, 1.0)

converged(result) || error("Failed to converge in $(iterations(result)) iterations")
xmin = result.minimizer


f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
x_iv = [0.0, 0.0]


results = optimize(f, x_iv, LBFGS())


results = optimize(f, x_iv, f_tol = 0.0001) # i.e. optimize(f, x_iv, NelderMead())
results.iterations
typeof(result)