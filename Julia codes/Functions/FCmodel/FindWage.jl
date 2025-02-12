############ General equilibrium: solving for wage ##############
function FindWage(wage, phi_c, zeta_Rl, zeta_Rs;)  

    beta = parameters().beta
    SumPol, e_chain, _ = FirmOptim(wage, phi_c, zeta_Rl, zeta_Rs)

    # entrant productivities set to be equal to the stationary distribution of e_chain    
    e_entry  = reduce(+,stationary_distributions(e_chain))

    # entrant X distribution - x_e = 0  in every case 
    x_vals = unique(SumPol[:, 1])
    zero_index = findall(x -> x == 0.0, x_vals)
    x_entry = zeros(length(x_vals))
    x_entry[zero_index .+ 0] .= 1 # if x = 0 prob = 1 

    # (x,e) are independent, the joint of the two distribution is their product
    xe_entry = [kron(e_entry, x_entry); 0] # also the f0 vector

    # map the entry probabilities to values
    Ve = transpose(xe_entry) * (SumPol[:,end])*beta

    return ( Ve )    

end