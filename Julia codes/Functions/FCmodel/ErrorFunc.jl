function ErrorFunc(initvec)

    include("C:/Users/szjud/OneDrive/Asztali gép/EBCs/CFL-git/Julia codes/Functions/FCmodel/FirmOptim.jl")

    ipc = initvec[1]
    ipdef_exo = initvec[2]
    idelta = initvec[3]
    izeta_R = initvec[4]
    izeta_L = initvec[5]

    # Check if any value in initvec is less than or equal to zero
    if any(initvec .<= 0)
        return ( Inf )
    end

    SumPol, e_chain, Fmat = FirmOptim(1, ipc, ipdef_exo, idelta, izeta_R, izeta_L, phi_c = 0.8)
    c_e, f0 = EntryValue(SumPol, e_chain) 
    SumRes = sumSS(SumPol,Fmat,f0)   

    ResultValues = SumRes[[7,8,9,12,13],1] |> Vector
    TargetValues = [ 0.51 , 4.9 , 2.9, 0.767 , 0.45 ]
    ModelError = sum( abs.((ResultValues .- TargetValues)) ./ TargetValues )

    return ( ModelError )
        
end