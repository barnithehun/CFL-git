############ General equilibrium: solving for wage ##############
function FindZeta(wage, phi_c, zeta_Rl , zeta_Rs)  
    
    beta = parameters().beta
    SumPol, e_chain, Fmat = FirmOptim(wage, phi_c, zeta_Rl, zeta_Rs_alt)
    c_e, f0 = EntryValue(SumPol, e_chain) ;
    mu, m, xpol = stat_dist(SumPol, Fmat, f0);

    liqprob = sumSSsme(SumPol,Fmat,f0)[10,1]

    return(liqprob)

end