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


function ErrorFunc2(initvec)

    ipc = initvec[1]
    ipdef_exo_l = initvec[2]
    ipdef_exo_s = initvec[3]
    izeta_Rl = initvec[4]
    izeta_Rs = initvec[5]
    iphi_a = initvec[6]
    iphi_c  = initvec[7]
    iFcut  = initvec[8]

    # Check if any value in initvec is less than or equal to zero
    if any(initvec .<= 0)
        return ( Inf )
    end

    SumPol, e_chain, Fmat = FirmOptim_Ext2(1, ipc, ipdef_exo_l, ipdef_exo_s, izeta_Rl, izeta_Rs, iphi_a, iphi_c, iFcut)
    c_e, f0 = EntryValue(SumPol, e_chain) 
    SumRes = sumSSsme(SumPol,Fmat,f0)   

    ResultValues = [ SumRes[[7,8,10,12],1] ; SumRes[[7,8,10,12],2]] |> Vector
    TargetValues = [ 0.36 , 5.54 , 0.702, 0.32, 0.68, 4.67, 0.322, 0.52 ]
    ModelError = sum( abs.((ResultValues .- TargetValues)) ./ TargetValues )
    # ModelError = sum( abs.(ResultValues .- TargetValues))

    return ( ModelError )
        
end