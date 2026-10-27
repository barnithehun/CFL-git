#### UPDATE THIS TO THE LATEST VER IF YOU USE IT ####

function GridMake()

        fn_L(k,e) =  (nu*e*k^alpha / wage)^(1/(1-nu))
        fn_Y(k, e) =  e*k^alpha*fn_L(k,e)^nu
        fn_Pi(k, e) = (1-nu)*fn_Y(k,e)-pc
        fn_X(k,b,e) =  fn_Pi(k, e) + (1-delta) * k - b
        fn_D(next_k, next_b, x, q) =  x - next_k + q * next_b
        # calling parameters
        rho_e, sigma_e, nul_e, alpha, nu, pc, beta, delta, pdef_exo, discount, phi_a, tauchen_sd, kappa, zeta_R, zeta_L, tau_vec = parameters()

        # calling grid size
        x_size, e_size, k_size, b_size = gridsize()
        
        # Setting the state-space
        # productivity process
        e_chain = tauchen(e_size, rho_e, sigma_e, (1-rho_e)*nul_e, tauchen_sd)
        e_vals = exp.(e_chain.state_values) 
        # adding exogeneous default shocks
        e_ptrans = e_chain.p

        # Log-grids
        k_grid = [0;exp.(range(log(10), log(10*10^5), k_size-1))]
        b_grid = [0;exp.(range(log(10), log(10*10^5), b_size-1))]  # no savings

        x_low = fn_X(k_grid[1],b_grid[end], e_vals[1])
        x_high = fn_X(k_grid[end], b_grid[1], e_vals[end])
        x_grid_low = sort(-exp.(range(log(10), log(-x_low), ceil(div(x_size, 3))))[2:end], rev = false) # set the negative part of the x-grid size to 1/3rd
        x_grid_high = exp.(range(log(10), log(x_high), (x_size - length(x_grid_low)-1))) 
        x_grid = [ x_grid_low; 0; x_grid_high ]

        # Define Q(n,m,n) matrix
        n = x_size * e_size  # all possible states
        m = k_size * b_size  # all possible actions

        # total number of possible states, +1 state for being quit 
        s_vals = [gridmake(x_grid, e_vals) zeros(n)]           
        s_vals = [s_vals; [0 0 1]]

        s_i_vals = [gridmake(1:x_size, 1:e_size) zeros(n)] 
        s_i_vals = Int.([s_i_vals; [div(x_size,2) 1 1]])  # productivity after default does not matter

        a_vals = [gridmake(k_grid, b_grid) zeros(m,2) ]
        a_vals = [a_vals; [0 0 1 0]; [0 0 0 1]]

        a_i_vals = [gridmake(1:k_size, 1:b_size) zeros(m,2)]
        a_i_vals = Int.([a_i_vals; [1 div(b_size,2) 1 0]; [1 div(b_size,2) 0 1]])

        # adjusting the gridsize
        n = n+1
        m = m+2


        return(n,m, e_vals,e_chain,e_ptrans, x_grid, s_vals,s_i_vals, a_vals,a_i_vals)
    
end