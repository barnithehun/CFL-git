function ErrorFunc1(initvec)

    ipc = initvec[1]
    ipdef_exo = initvec[2]
    izeta_Rl = initvec[3]
    iphi_c_hh  = initvec[4]
    iphi_c  = initvec[5]

    # Check if any value in initvec is less than or equal to zero
    if any(initvec .<= 0)
        return ( Inf )
    end

    SumPol, e_chain, Fmat = FirmOptim_Ext1(1, ipc, ipdef_exo, izeta_Rl, iphi_c_hh, iphi_c)
    c_e, f0 = EntryValue(SumPol, e_chain) 
    SumRes = sumSS(SumPol,Fmat,f0)   

    ResultValues = SumRes[[7, 8, 10, 12, 13],1] |> Vector
    TargetValues = [ 0.49, 5.3, 0.48, 0.77, 0.43 ]
    ModelError = sum( abs.((ResultValues .- TargetValues)) ./ TargetValues )

    return ( ModelError )
        
end


function ErrorFunc_CF(initvec)

    ipc = initvec[1]
    ipdef_exo = initvec[2]
    izeta_Rl = initvec[3]

    # Check if any value in initvec is less than or equal to zero
    if any(initvec .<= 0)
        return ( Inf )
    end

    SumPol, e_chain, Fmat = FirmOptim_Ext_CF(1, ipc, ipdef_exo, izeta_Rl)
    c_e, f0 = EntryValue(SumPol, e_chain) 
    SumRes = sumSS(SumPol,Fmat,f0)   

    ResultValues = SumRes[[7, 8, 10],1] |> Vector
    TargetValues = [ 0.49, 5.3, 0.48 ]
    ModelError = sum( abs.((ResultValues .- TargetValues)) ./ TargetValues )

    return ( ModelError )
        
end


function ErrorFunc_AB(initvec)

    ipc = initvec[1]
    ipdef_exo = initvec[2]
    izeta_Rl = initvec[3]

    # Check if any value in initvec is less than or equal to zero
    if any(initvec .<= 0)
        return ( Inf )
    end

    phi_c = 0
    SumPol, e_chain, Fmat = FirmOptim_Ext_AB(1, phi_c, ipc, ipdef_exo, izeta_Rl)
    c_e, f0 = EntryValue(SumPol, e_chain) 
    SumRes = sumSS(SumPol,Fmat,f0)   

    ResultValues = SumRes[[7, 8, 10],1] |> Vector
    TargetValues = [ 0.49, 5.3, 0.48 ]
    ModelError = sum( abs.((ResultValues .- TargetValues)) ./ TargetValues )

    return ( ModelError )
        
end


