# Chase Abram
# June 11, 2023

# Basic Hopenhayn (1992) model of firm dynamics with entry and exit
# Continuous time, based on notes by Moll and MATLAB by Cioffi

# Load packages
using LinearAlgebra

using SparseArrays

using Plots
default(linewidth = 3.0,  
legendfontsize = 12,
size = (800,600))

using LaTeXStrings

using NLsolve

using Parameters

##
# Model struct
@with_kw mutable struct HopenhaynModel
    
    # Grid points
    nz = 1000

    # Goods demand
    epsilon = 0.5
    pf = arg -> arg^(-epsilon)
    # Initial goods quantity
    Q = 1.0
    # Intial goods price
    p = 0.5

    # Labor supply
    phi = 0.5
    wf = arg -> arg^(phi)
    # Initial labor quantity
    N = 1.0
    # Initial wage
    w = 1.0

    # Production tech
    alpha = 0.5
    f = (arg1, arg2) -> arg1*arg2^alpha

    # operating cost
    c_f = 0.05

    # Scrap Value
    vartheta = 0.0
    vartheta_vec = vartheta .* ones(nz)

    # discount rate
    rho = 0.05

    # z values
    z_lb = 0.0
    z_ub = 1.0
    z = LinRange(z_lb, z_ub, nz)
    dz = diff(z)

    # employment at z
    n = zeros(nz)

    # output at z 
    y = zeros(nz)

    # profits at z
    prof = zeros(nz)

    # Value function
    v = zeros(nz)
    # Index of entry/exit cutoff
    vartheta_ind = -1

    # Equilibrium distribution 
    g = ones(nz)

    # Entry distribution
    h_lb = (z_ub - z_lb) * 0.7
    beta = 0.0
    h = entry_dist(z, h_lb, beta)

    # Entry cost 
    c_e = 0.6

    # Entry rate
    m = 1.0
    mbar = 0.1
    eta = 1000.0

    # Diffusion process terms
    mu = -0.01
    sigma = 0.01
    A = diffusion_op(z,mu,sigma)
    Atilde = deepcopy(A)

    # Objects for LCP problem
    B = rho * I - A
    x = zeros(nz)
    q = ones(nz)

    # Numerical parameters
    w_update = 0.2
    p_update = 0.001
    maxit = 1000
    tol = 1e-5
end

##

# Build diffusion process operator
function diffusion_op(z, mu, sigma)
    
    nz = length(z)
    dz = diff(z)

    # Initialize (diffusion is sparse)
    A = spzeros(nz,nz)

    # Diffusion proportional to TFP
    sig = sigma .* z

    # Reflect at bottom
    A[1,2] = 1/2*sig[1]^2/dz[1]^2 + max(mu,0)/dz[1]
    A[1,1] = -A[1,2]

    for i in 2:(nz-1)
        A[i,i-1] = 1/2*sig[i]^2/dz[i]^2 - min(mu,0)/dz[i]
        A[i,i+1] = 1/2*sig[i]^2/dz[i]^2 + max(mu,0)/dz[i]
        A[i,i] = -(A[i,i-1] + A[i,i+1])
    end

    # Reflect at top
    A[nz,nz-1] = 1/2*sig[nz]^2/dz[end]^2 - min(mu,0)/dz[end]
    A[nz,nz] = -A[nz,nz-1]

    return A
end

# Trapezoidal integration
function trap_sum(a, dx)
    return sum((a[1:end-1] + a[2:end])/2 .* dx)
end

# Entry distibution
function entry_dist(z, lb, beta)
    
    # Shifted exponential
    # beta = 0 => Uniform
    h = exp.(-beta .* z) .* (z .>= lb)
    h ./= trap_sum(h, diff(z))

    return h
end

# Solve firm problem, given prices
function firm_solve(model)
    @unpack_HopenhaynModel model

    # Labor demand
    n .= (alpha .* z .* p ./ w) .^(1/(1-alpha))
    # Output supply
    y .= z .* n .^ alpha
    # Profits
    prof .= p .* y .- w .* n .- c_f

    @pack_HopenhaynModel! model
end

# Solve HJBVI via LCP solver
function value_solve(model)
    @unpack_HopenhaynModel model

    # Set up LCP
    q .= -prof + B * vartheta_vec

    function F!(Func, x)
        Func .= B*x + q
    end

    lb = zeros(nz)
    ub = Inf .* ones(nz)
    
    # Solve LCP (NLsolve is just one choice of solver)
    soln = mcpsolve(F!, lb, ub, v .- vartheta_vec)
    @assert converged(soln) == true
    
    # Implied value function
    v .= soln.zero .+ vartheta_vec
    @pack_HopenhaynModel! model
end

# Solve KF equation
function KF_solve(model)
    @unpack_HopenhaynModel model

    # Expected value of entry (trapezoid sum)
    v_e = trap_sum(v .* h, dz)

    # Entry rate (some elasticity for numerical stability)
    m = mbar * exp(eta * (v_e - c_e))

    # find exit index
    vartheta_ind = sum(v .<= vartheta_vec)

    # Implied flows given entry/exit behavior
    # Follow drift process
    Atilde .= A
    # Except in exit region no flows (so absorbing at boundary)
    Atilde[:,1:vartheta_ind-1] .= 0
    for i in 1:(vartheta_ind - 1)
        Atilde[i, i] = 1
    end

    # Solve KF
    g .= -Atilde' \ (m .* h)
    # For numerical stability, keep dist >= 0
    g .= max.(g,0)

    @pack_HopenhaynModel! model
end

# Compute aggregates
function aggs(model)
    @unpack_HopenhaynModel model

    # The bounds add stability
    N = max(min(trap_sum(n .* g, dz), 1e4), 1e-4)
    Q = max(min(trap_sum(y .* g, dz), 1e4), 1e-4)

    @pack_HopenhaynModel! model
end

# Solve equilibrium
function solve_eq(model)
    @unpack_HopenhaynModel model

    # Initialize wage loop
    it_w = 1
    err_w = Inf
    # Until labor market clears or maxit
    while it_w < maxit && err_w > tol
        
        # Initialize price loop
        it_p = 1
        err_p = Inf

        # Until goods market clears or maxit
        while it_p < maxit && err_p > tol

            @pack_HopenhaynModel! model
            # Solve (static) firm problems
            firm_solve(model)

            # Solve HJBVI
            value_solve(model)

            # Solve KF
            KF_solve(model)

            # Implied aggregates
            aggs(model)
            @unpack_HopenhaynModel model

            # Goods market residual
            err_p = abs(pf(Q) - p)
            # Update goods price
            p = p_update * pf(Q) + (1-p_update) * p

            if (it_p + 1) % 100 == 0
                println("it_w: ", it_w, ", err_w: ", err_w, ", it_p: ", it_p, ", pf: ", pf(Q), ", err_p: ", err_p)
            end

            it_p += 1
        end

        # Labor market residual
        err_w = abs(wf(N) - w)
        # Update wage
        w = w_update * wf(N) + (1-w_update) * w
        println("it_w: ", it_w, ", err_w: ", err_w, ", it_p: ", it_p, ", pf: ", pf(Q), ", err_p: ", err_p)
        it_w += 1
    end

    @pack_HopenhaynModel! model
end

##
# Lower entry cost
model0 = HopenhaynModel(c_e = 0.6)
solve_eq(model0)

# Higher entry cost
model1 = HopenhaynModel(c_e = 0.7, w_update = 0.5, p_update = 1e-4)
@elapsed solve_eq(model1)

##

# Value functions
pltv = plot(model0.z, model0.v, label = L"v_{c_e = 0.6}(z)", xlabel = L"z")
plot!(pltv, model1.z, model1.v, label = L"v_{c_e = 0.7}(z)", xlabel = L"z")
plot!(pltv, model0.z, model0.vartheta_vec, label = L"\vartheta")
display(pltv)

# TFP distributions
pltg = plot(model0.z, model0.g, label = L"g_{c_e = 0.6}(z)", xlabel = L"z")
plot!(pltg, model1.z, model1.g, label = L"g_{c_e = 0.7}(z)", xlabel = L"z")
# Entry dist
plot!(pltg, model0.z, model0.h, label = L"h(z)", xlabel = L"z")
display(pltg)

# Employment distributions
pltn = plot(model0.n, model0.g, label = L"g_{c_e = 0.6}(n)", xlabel = L"n")
plot!(pltn, model1.n, model1.g, label = L"g_{c_e = 0.7}(n)", xlabel = L"n")
display(pltn)

# Output distributions
plty = plot(model0.y, model0.g, label = L"g_{c_e = 0.6}(y)", xlabel = L"y")
plot!(plty, model1.y, model1.g, label = L"g_{c_e = 0.7}(y)", xlabel = L"y")
display(plty)

