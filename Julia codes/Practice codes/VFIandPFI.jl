# This is a Julia code on the basic structure on VFI and PFI
# The example model is a simple growth model with recursive formulation: 
#  y = k^α, k_0 is given
#  Bellman: V(k) = max_k'  ln(f(k)-k') + βV(k')
# No uncertainty, one state var: k, and one decisions var k'

using Roots, Optim, Interpolations, Plots

### Problem paramaters ###

# paramaters
alpha    = 0.65
beta      = 0.95
# setting up the grid
grid_max  = 2       # upper bound of capital grid
n         = 200     # number of grid points
N_iter    = 3000    # number of iterations
kgrid     = 1e-2:(grid_max-1e-2)/(n-1):grid_max  # equispaced grid
f(x) = x^alpha      # defines the production function f(k)
tol = 1e-9

### Analytical solution ###
# this is only used to check the validity of the results
ab        = alpha * beta
c1        = (log(1 - ab) + log(ab) * ab / (1 - ab)) / (1 - beta) 
c2        = alpha / (1 - ab)

v_star(k) = c1 .+ c2 .* log.(k)  # Max value, given k (?)
k_star(k) = ab * k.^alpha        # Optimal k' given k
c_star(k) = (1-ab) * k.^alpha    # Implies consumption
ufun(x) = log.(x)                # Utility function

# Function for the bellman operator
# Takes a grid of state variables and computes the next iterate of the value function
# Basically for some V0 it loops over every possible k to find the one that gives the largest value
#         - this will be V1
#         - the k that maximizes V1 is the policy
function bellman_operator(grid,v0) # out contraction mapping
    
    v1  = zeros(n)         # next guess
    pol = zeros(Int,n)     # policy function
    w   = zeros(n)         # temporary vector 

    # loop over current states
    # current capital
    for (i,k) in enumerate(grid)  # for each i - as in for each value for the state variable
       # due to enumerate() , i refers to the index and k refers to the value

        # loop over all possible kprime choices
        for (iprime,kprime) in enumerate(grid) # makes a vector w[] with values for each k'
            if f(k) - kprime < 0   #check for negative consumption
                w[iprime] = -Inf
            else
                w[iprime] = ufun(f(k) - kprime) + beta * v0[iprime]
            end
        end
        # for the given each 
        v1[i], pol[i] = findmax(w)     # stores Value and policy (index of optimal choice)
    end  # you need to find the best k' for each i that is why you need two loops
    return (v1,pol)   # return both value and policy function

end
# now you have the contraction mapping, applying it over and over again will yield the optimal value function

function VFI()
    v_init = zeros(n)     # initial guess - this VFI will autmatically use a V0 = 0 starting guess
    for iter in 1:N_iter
        v_next = bellman_operator(kgrid,v_init)  # returns a tuple: (v1,pol)
        # thats why you need v_next[1] below, because you are just checking the errors
        
        # check convergence
        if maximum(abs,v_init.-v_next[1]) < tol
            return v_next

            verrors = maximum(abs,v_next[1].-v_star(kgrid))         # just comparing it to the ideal solution
            perrors = maximum(abs,kgrid[v_next[2]].-k_star(kgrid))  # just comparing it to the ideal solution
            println("Found solution after $iter iterations")
            println("maximal value function error = $verrors")
            println("maximal policy function error = $perrors")
        elseif iter==N_iter
            warn("No solution found after $iter iterations")
            return v_next
        end
        v_init = v_next[1]  # update guess 
    end
end 


function plotVFI()
    v = VFI()
    p = Any[]
    
    # value and policy functions
    push!(p,plot(kgrid,v[1],
            lab="V",
            ylim=(-50,-30),legend=:bottomright),
            plot(kgrid,kgrid[v[2]],
            lab="policy",legend=:bottomright))
    
    # errors of both
    push!(p,plot(kgrid,v[1].-v_star(kgrid),
        lab="V error",legend=:bottomright),
        plot(kgrid,kgrid[v[2]].-k_star(kgrid),
        lab="policy error",legend=:bottomright))

    plot(p...,layout=grid(2,2) )
    
end

plotVFI()


### PFI ###

function policy_iter(grid,c0,u_prime,f_prime)
    
    c1  = zeros(length(grid))     # next guess
    pol_fun = extrapolate(interpolate((collect(grid),), c0, Gridded(Linear()) ) , Interpolations.Flat())
    
    
    # loop over current states
    # of current capital
    for (i,k) in enumerate(grid)
        objective(c) = u_prime(c) - beta * u_prime(pol_fun(f(k)-c)) * f_prime(f(k)-c)
        c1[i] = fzero(objective, 1e-10, f(k)-1e-10) 
    end
    return c1
end

uprime(x) = 1.0 ./ x
fprime(x) = alpha * x.^(alpha-1)



function PFI()
    c_init = kgrid
    for iter in 1:N_iter
        c_next = policy_iter(kgrid,c_init,uprime,fprime)  
        # check convergence
        if maximum(abs,c_init.-c_next) < tol
            perrors =  maximum(abs,c_next.-c_star(kgrid))
            println("PFI:")
            println("Found solution after $iter iterations")
            println("max policy function error = $perrors")
            return c_next
        elseif iter==N_iter
            warn("No solution found after $iter iterations")
            return c_next
        end
        c_init = c_next  # update guess 
    end
end
function plotPFI()
    v = PFI()
    plot(kgrid,[v v.-c_star(kgrid)],
            lab=["policy" "error"],
            legend=:bottomright,
            layout = 2)
end

plotPFI()
