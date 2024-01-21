# Example solution of a discrete DP

# consumption-SAVING model
# Household stats with stock of good s, (state variable)
# Stores a (action variable)
# Consumes s-a,
# Has a sochastic income {0, ..., U}
# Stock evolves according to  s' = a + U
# discount factor is beta

using BenchmarkTools, LaTeXStrings, LinearAlgebra, Plots, QuantEcon, Statistics, SparseArrays, BenchmarkTools, Plots, QuantEcon

# this is a way to make parameters
# calling SimpleOG() makes a named tuple with the parameter values
SimpleOG(; B = 10, M = 5, alpha = 0.5, beta = 0.9) = (; B, M, alpha, beta)


function transition_matrices(g)
    (; B, M, alpha, beta) = g       # will define the parameters withing the function - in the spirit of avoiding global variables
    u(c) = c^alpha
    n = B + M + 1
    m = M + 1

    R = zeros(n, m)
    Q = zeros(n, m, n)

    for a in 0:M
        Q[:, a + 1, (a:(a + B)) .+ 1] .= 1 / (B + 1)
        for s in 0:(B + M)
            R[s + 1, a + 1] = (a <= s ? u(s - a) : -Inf)
        end
    end

    return (; Q, R)
end

# same thing but verbose  - this is way more intuitive and readable
function verbose_matrices(g)
    (;B, M, alpha, beta) = g
    u(c) = c^alpha

    #Matrix dimensions. The +1 is due to the 0 state.
    n = B + M + 1
    m = M + 1

    R = fill(-Inf, n, m) #Start assuming nothing is feasible
    Q = zeros(n,m,n) #Assume 0 by default

    #Create the R matrix
    #Note: indexing into matrix complicated since Julia starts indexing at 1 instead of 0
    #but the state s and choice a can be 0
    for a in 0:M
         for s in 0:(B + M)
            if a <= s #i.e. if feasible
                R[s + 1, a + 1] = u(s - a)
            end
        end
    end

    #Create the Q multi-array
    for s in 0:(B+M) #For each state
        for a in 0:M #For each action
            for sp in 0:(B+M) #For each state next period
                if( sp >= a && sp <= a + B) # The support of all realizations
                    Q[s + 1, a + 1, sp + 1] = 1 / (B + 1) # Same prob of all
                end
            end
            @assert sum(Q[s + 1, a + 1, :]) ≈ 1 #Optional check that matrix is stochastic
         end
    end
return (;Q,R)
end

g = SimpleOG()
Q, R = verbose_matrices(g)

ddp = DiscreteDP(R, Q, g.beta);
# this function just takes in the arrays R and Q and the discount paramter beta
# then it outputs it in a single array (they call it a DiscreteDP object) that can be handed to the solver

# solve is not a generic function it is part of the QuantEcon package just as DiscreteDP
@elapsed results = solve(ddp, VFI)
results.num_iter


fieldnames(typeof(results))

# Doing the same thing with state-action pairs,
B = 10
M = 5
alpha = 0.5
beta = 0.9
u(c) = c^alpha
n = B + M + 1
m = M + 1

s_indices = Int64[]
a_indices = Int64[]
Q = zeros(0, n)
R = zeros(0)

b = 1 / (B + 1)

# If I understand correctly,
#  - this approach stacks the 3 dimensional Q into a 2 dimensional 
#  - stacks R matrix into a column vector
#  - these will be sparse matrices but if you DiscreteDP() as below, it will handle it  
for s in 0:(M + B)
    for a in 0:min(M, s)
        s_indices = [s_indices; s + 1]
        a_indices = [a_indices; a + 1]
        q = zeros(1, n)
        q[(a + 1):((a + B) + 1)] .= b
        Q = [Q; q]
        R = [R; u(s - a)]
    end
end


ddp = DiscreteDP(R, Q, beta, s_indices, a_indices);
@elapsed results = solve(ddp, PFI) # look like it is even faster
