# Example solution of a discrete DP

# consumption-SAVING model
# Household starts with stock of good s, (state variable)
# Stores a (action variable)
# Consumes s-a,
# Has a sochastic income {0, ..., U}
# Stock evolves according to  s' = a + U
# discount factor is beta

using BenchmarkTools, LaTeXStrings, LinearAlgebra, Plots, QuantEcon, Statistics, SparseArrays

# this is a way to make parameters
# calling SimpleOG() makes a named tuple with the parameter values
SimpleOG(; B = 3, M = 2, alpha = 0.5, beta = 0.9) 

function transition_matrices(g)
    (; B, M, alpha, beta) = g       # define the parameters within the function -in the spirit of avoiding global variables
    u(c) = c^alpha
    n = B + M + 10
    m = M + 1

    R = zeros(n, m)
    Q = zeros(n, m, n)

    for a in 0:M
        Q[:, a + 1, (a:(a + B)) .+ 1] .= 1 / (B + 1)   # here the state variable will be reset every period by the action, so you can define entire columns for the Q matrix! 
        for s in 0:(B + M)
            R[s + 1, a + 1] = (a <= s ? u(s - a) : -Inf)
        end
    end

    return (; Q, R)
end

g = SimpleOG();
Q, R = transition_matrices(g);

# Lets check
Q[:,:,1]
# This means that if:
#    - storage (action) = 0, s' = 0 with a 0.25 prob, for any s = 0:5
#    - storage (action) = 1:2, s' = 0 with a 0 prob, for any s = 0:5
#    - this makes sense because you cannot have a negative shock to your storage so it wont be less...

# Now see: 
Q[:,:,2]
# This means that if:
#    - storage (action) = 0:1, s' = 2 with a 0.25 prob, for any s = 0:5
#    - storage (action) = 2, s' = 0 with a 0 prob, for any s = 0:5
#    - this makes sense because you cannot have a negative shock to your storage so it wont be less...

# Now see: 
Q[:,:,6]
# This means that next period state variable can only be 5 if you store the maximum possible of 2
# Then it has a 25% chance that you will have B=3 and so s' = 5


# same thing but verbose  - this is way more intuitive and readable
function verbose_matrices(g)
    (; B, M, alpha, beta) = g
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
    for s in 0:(B+M) # For each state
        for a in 0:M # For each action
            for sp in 0:(B+M) # For each state next period
                if( sp >= a && sp <= a + B) # The support of all realizations
                    Q[s + 1, a + 1, sp + 1] = 1 / (B + 1) # Same prob of all
                end
            end
            @assert sum(Q[s + 1, a + 1, :]) ≈ 1 #Optional check that matrix is stochastic
         end
    end
return (;Q,R)
end

# produces the same result
g = SimpleOG()
Q, R = verbose_matrices(g)

ddp = DiscreteDP(R, Q, g.beta);
# this function just takes in the arrays R and Q and the discount paramter beta
# then it outputs it in a single array (they call it a DiscreteDP object) that can be handed to the solver

# solve is not a generic function it is part of the QuantEcon package just as DiscreteDP
@elapsed results = solve(ddp, PFI)
results.num_iter

# all the 
fieldnames(typeof(results))

# mc stands for markov chain
# this is the transition probability given the optimal policy!!
results.mc
# it gives the dynamics of the state when the agent follows the optimal policy

# now, it is simple to find the corresponding stationary distribution
stationary_distributions(results.mc)[1]


########## Doing the same thing with state-action pairs #############
B = 3
M = 2
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

# You just define the feasible state-action pairs
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

# This desribes all feasible(!) state action pairs
[s_indices a_indices]
# Corresponding return functions
[s_indices a_indices R] # notice how there is no -INF in the return function, that is because we only consider the feasible combinations
# Corresponding transition probabilities
[s_indices a_indices Q]
#  Q[1, :] says that if the state is 0 and the actions is 0  - you have a 25% chance that the state for each future state 0-3


ddp = DiscreteDP(R, Q, beta, s_indices, a_indices);
@elapsed results = solve(ddp, PFI) # look like it is faster


