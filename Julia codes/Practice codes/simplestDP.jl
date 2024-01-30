# This is the simplest possible application of the discrete DP
# The problem is discussed at: https://quanteconpy.readthedocs.io/en/latest/markov/ddp.html#solution-algorithms
using BenchmarkTools, LaTeXStrings, LinearAlgebra, Plots, QuantEcon, Statistics, SparseArrays

#= 
    Set of states S = {0, 1}
    Set of actions A = {0, 1}
    Set of feasible state-action pairs SA = {(0, 0), (0, 1), (1, 0)}
    Rewards 
        r(s, a): r(0, 0) = 5, r(0, 1) =10, r(1, 0) = -1
    Transition probabilities q(s_next|s, a):
        q(0|0, 0) = 0.5, q(1|0, 0) = 0.5, q(0|0, 1) = 0, q(1|0, 1) = 1, q(0|1, 0) = 0, q(1|1, 0) = 1
    Discount factor 0.95 
=#

# Standard approach
beta = 0.95
R = [5 10;-1 -float(Inf)] 
Q = cat([0.5 0; 0 0.5], [0.5 1; 1 0.5]; dims=3)

ddp = DiscreteDP(R, Q, beta)
results = solve(ddp, VFI)

# State - action formulation
# In this case, you only have to cover the feasible cases
# Note indicies must be positive integers!!!! - this is not the state and action but their corresponding indicies!
s_indices = [1, 1, 2]  # State indices
a_indices = [1, 2, 1]  # Action indices

R = [5, 10, -1]
Q = [0.5 0.5; 0 1; 0 1]
beta = 0.95

[s_indices a_indices]
# Corresponding return functions
[s_indices a_indices R] # notice how there is no -INF in the return function, that is because we only consider the feasible combinations
# Corresponding transition probabilities
[s_indices a_indices Q]

R = convert(Vector{Float64}, R)
Q = convert(Matrix{Float64}, Q)

ddp = DiscreteDP(R, Q, beta, s_indices, a_indices)
results = solve(ddp, VFI)

