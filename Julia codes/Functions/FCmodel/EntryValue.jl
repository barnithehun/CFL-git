####### ENTRANTS #######
function EntryValue(SumPol, e_chain)  

    # entrant ln(e) distribution - equal to the stationary distribution of e_chain    
    e_entry  = reduce(+, stationary_distributions(e_chain))

    # entrant X distribution - x_e = 0  in every case 
    x_vals = unique(SumPol[:, 1])
    zero_index = findall(x -> x == 0.0, x_vals)
    x_entry = zeros(length(x_vals))
    x_entry[zero_index .+ 0] .= 1 # if x = 0 prob = 1 

    # (x,e) are independent, the joint of the two distribution is their product
    xe_entry = [kron(e_entry, x_entry); 0] # also the f0 vector

    # map the entry probabilities to values
    beta = parameters().beta
    Ve = transpose(xe_entry) * (SumPol[:,end])*beta

    return ( Ve, xe_entry )    

end
