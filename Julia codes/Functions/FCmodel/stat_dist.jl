####### STATIONARY DISTRIBUTION ########
function stat_dist(SumPol, Fmat, f0)

    # Exiting firms + defaulting firms
    n = size(SumPol,1)
    # here it states that a firm in default state exits the market, it means that it need no updating after Fmat is changed to reflect Q
    xpol = [SumPol[1:n-1,6] + SumPol[1:n-1,7] ; 1] 
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

end