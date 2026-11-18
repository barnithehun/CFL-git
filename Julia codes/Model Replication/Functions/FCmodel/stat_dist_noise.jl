####### STATIONARY DISTRIBUTION ########
# zeta_R: fixed reorganization cost for the liquidation decision; defaults to the global zeta_Rl
function stat_dist(SumPol, Fmat, f0; zeta_R = zeta_Rl)

    # Exiting firms + defaulting firms
    n = size(SumPol,1)

    # "liq" vector: firms that are liquidated in default
    phi_a = parameters().phi_a
    phi_c_hh = parameters().phi_c_hh
    delta = parameters().delta

    # liq decision works with policies works bc. of stationary distribution
    # with valuation noise (parameters().sigma_v > 0, Baseline_noise.jl) the liquidation probability of FirmOptim's fn_Gam
    sigma_v = parameters().sigma_v
    L = phi_a .* (1-delta) .* SumPol[1:n-1, 3]; R = phi_c_hh .* SumPol[1:n-1, end]
    liq = sigma_v == 0 ? Int.( L .>= (R .- zeta_R) ) :
          [r <= 0 ? 1.0 : l + zeta_R <= 0 ? 0.0 : cdf(Normal(), log((l + zeta_R) / r) / sig) for (l, r, sig) in
           zip(L, R, (haskey(parameters(), :sv_slope) ? sigma_at.(SumPol[1:n-1, 3]) : fill(sigma_v, n-1)))]

    # a default state exits the market, it means that its need no updating after Fmat is changed to reflect Q
    xpol = [SumPol[1:n-1,6] .+ SumPol[1:n-1,7] .* liq ; 1]
    Ident = Matrix(I,n,n)

    xpol_mat = Ident - Diagonal(xpol)   # I - diag(X)
    f0 = xpol_mat*f0                    # (I - diag(X))f0
    Mmat = Fmat*xpol_mat                # M = F(I - diag(X))

    # unscaled stationary distribution
    mu_0 = inv(Ident - Mmat)*f0         # inv(I-M)*f0

    # ok, bc exit and default implies k = 0, n = 0
    Nd = transpose(SumPol[1:n,10])*mu_0
    Ns = 10000

    # This is given that labour supply is one (10000) inelastically
    m = Ns/Nd
    mu = m.*mu_0

    return ( mu, m , xpol )

    println(zeta_Rl)

end