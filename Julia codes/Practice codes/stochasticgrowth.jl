# Codes that are related to the stochastic growth model lectures
# https://julia.quantecon.org/dynamic_programming/optgrowth.html#equation-fcbell20-optgrowth
# The model itself is just an example, the emphasis should be on the recoursice solution techniques

using LinearAlgebra, Statistics, LaTeXStrings, Plots, Interpolations, NLsolve, Optim, Random
using Optim: maximum, maximizer

# Bellman operator - just 1 application of it, the tolerance is for the maximization function in it
function T(w; p, tol = 1e-10)
    # it would make a name tuple and unpack them only within the function - this is good practice supposedly
    (; beta, u, f, Xi, y) = p # unpack parameters
    w_func = LinearInterpolation(y, w)   # this creates  the value function 

    Tw = similar(w)    # similar() defines a vector that has the same dimensions and types as the original
    sigma = similar(w)

    for (i, y_val) in enumerate(y)
        # solve maximization for each point in y, using y itself as initial condition.
        results = maximize(c -> u(c; p) +   # current return - the u() fc should be given as a parameter of the fc
                                beta * mean(w_func.(f(y_val - c; p) .* Xi)), # expected value of future return

                                # the mean fucntion is needed bc. the stochastic part is computed by Monte-Car
                           tol, y_val)
        Tw[i] = maximum(results)
        sigma[i] = maximizer(results)
    end
    return (; w = Tw, sigma) # returns named tuple of results
end

### First check an analytical result of the model  -
    # Note that an analytical solution is only possible with these specific function forms

Random.seed!(42) # for reproducible results
u(c; p) = log(c) # utility
f(k; p) = k^p.alpha # deterministic part of production function

# this just creates paramter values in a tupple
function OptimalGrowthModel(; alpha = 0.4, beta = 0.96, mu = 0.0, s = 0.1,
                            u = u, f = f, # defaults defined above
                            y = range(1e-5, 4.0, length = 200), # grid on y
                            Xi = exp.(mu .+ s * randn(250)))
    return (; alpha, beta, mu, s, u, f, y, Xi)
end # named tuples defaults

# True value and policy function
function v_star(y; p)
    (; alpha, mu, beta) = p
    c1 = log(1 - alpha * beta) / (1 - beta)
    c2 = (mu + alpha * log(alpha * beta)) / (1 - alpha)
    c3 = 1 / (1 - beta)
    c4 = 1 / (1 - alpha * beta)
    return c1 + c2 * (c3 - c4) + c4 * log(y)
end
c_star(y; p) = (1 - p.alpha * p.beta) * y

### Now see how close to this is the numerical solution of the model
p = OptimalGrowthModel() # use all default parameters from named tuple
w_star = v_star.(p.y; p)  # evaluate closed form value along grid

# here you need to apply the Bellman only once, because you know you alraedy start with the optimal value*

w = T(w_star; p).w # evaluate operator, access Tw results

# almost perfect, apart from a small numerical errors
plt = plot(ylim = (-35, -24))
plot!(plt, p.y, w, linewidth = 2, alpha = 0.6, label = L"T(v^*)")
plot!(plt, p.y, w_star, linewidth = 2, alpha = 0.6, label = L"v^*")
plot!(plt, legend = :bottomright)


### Now see how much iteration it takes if the starting guess is not perfect
w = 5 * log.(p.y)  # An initial condition -- fairly arbitrary
n = 35

plot(xlim = (extrema(p.y)), ylim = (-50, 10))

lb = "initial condition"
plt = plot(p.y, w, color = :black, linewidth = 2, alpha = 0.8, label = lb)

for i in 1:n
    w = T(w; p).w
    plot!(p.y, w, color = RGBA(i / n, 0, 1 - i / n, 0.8), linewidth = 2,
          alpha = 0.6,
          label = "")
end

lb = "true value function"
plot!(plt, p.y, v_star.(p.y; p), color = :black, linewidth = 2, alpha = 0.8,
      label = lb)
plot!(plt, legend = :bottomright)


function solve_optgrowth(initial_w; p, iterations = 500, m = 3, show_trace = false)
    results = fixedpoint(w -> T(w; p).w, initial_w; iterations, m, show_trace) # Anderson iteration
    v_star = results.zero
    sigma = T(results.zero; p).sigma
    return (; v_star, sigma, results)
end

initial_w = 5 * log.(p.y)
sol = solve_optgrowth(initial_w; p)
v_star_approx = sol.v_star
println("Converged in $(sol.results.iterations) to an ||residuals||_∞ of $(sol.results.residual_norm)")

plt = plot(ylim = (-35, -24))
plot!(plt, p.y, v_star_approx, linewidth = 2, alpha = 0.6,
      label = "approximate value function")
plot!(plt, p.y, v_star.(p.y; p), linewidth = 2, alpha = 0.6,
      label = "true value function")
plot!(plt, legend = :bottomright)


