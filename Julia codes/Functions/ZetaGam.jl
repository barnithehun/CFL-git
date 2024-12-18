function ZetaGam(zeta_R_vec::AbstractVector)

    wage = 1
    phi_c = 0.4073866012405985;
    zeta_Rl = 8459.339189508773;

    
    function gridsize()
        # grid sizes - x, k, b should be even numbers!!
        x_size::Int = 46
        e_size::Int = 23
        k_size::Int = 32
        b_size::Int = 32
        return (x_size = x_size, e_size = e_size, k_size = k_size, b_size = b_size)
    end
    n,_, _,_,_, _, _,_, _,_ = GridMake()

    SumPolZeta::Array{Float64, 3} = zeros(n, 18, length(zeta_R_vec))
    for (ind, zeta_RR) in enumerate(zeta_R_vec)

        println("Iteration: $ind, Value: $zeta_RR")

        SumPol, _, _ = FirmOptim(wage, phi_c, zeta_RR, zeta_Rl)
        
        SumPolZeta[:,:,ind] = SumPol


    end

    return  SumPolZeta
    
end 
