function ZetaGamPlot(SumPolZeta,zeta_R_vec, e_i, x_i)
    
x_size = gridsize().x_size    
xe_i_1 = x_i + (e_i-1)*x_size   
xe_i_2 = x_i + (e_i-2)*x_size   
xe_i_3 = x_i + (e_i-3)*x_size   
x_val = Int(round(SumPol[xe_i_1, 1],digits = 0))
e_val_1 = round(SumPol[xe_i_1, 2],digits = 1)
e_val_2 = round(SumPol[xe_i_2, 2],digits = 1)
e_val_3 = round(SumPol[xe_i_3, 2],digits = 1)


p1 = plot(zeta_R_vec ./ 1000, SumPolZeta[xe_i_1, 4, :], 
title="Debt", titlefont=font(14), xlabel="Zeta_R (thousands)", ylabel="Debt", 
legend=:topright, linecolor=:blue, linewidth=2, label="Prod = $e_val_1")

# Add xe_i_2 and xe_i_3 series to p1
plot!(zeta_R_vec ./ 1000, SumPolZeta[xe_i_2, 4, :], linecolor=:red, linewidth=2, label="Prod = $e_val_2")
plot!(zeta_R_vec ./ 1000, SumPolZeta[xe_i_3, 4, :], linecolor=:green, linewidth=2, label="Prod = $e_val_3")

p2 = plot(zeta_R_vec ./ 1000, SumPolZeta[xe_i_1, 9, :], 
    title="Inverse Interest Rate", titlefont=font(14), xlabel="Zeta_R (thousands)", ylabel="Inverse Interest rate rate (q)", 
    legend=:topright, linecolor=:blue, linewidth=2, label="Prod = $e_val_1")

# Add xe_i_2 and xe_i_3 series to p2
plot!(zeta_R_vec ./ 1000, SumPolZeta[xe_i_2, 9, :], linecolor=:red, linewidth=2, label="Prod = $e_val_2")
plot!(zeta_R_vec ./ 1000, SumPolZeta[xe_i_3, 9, :], linecolor=:green, linewidth=2, label="Prod = $e_val_3")

p3 = plot(zeta_R_vec ./ 1000, SumPolZeta[xe_i_1, 14, :], 
    title="Liquidation Probability", titlefont=font(14), xlabel="Zeta_R (thousands)", ylabel="Liquidation Probability", 
    legend=:bottomright, linecolor=:blue, linewidth=2, label="Prod = $e_val_1")

# Add xe_i_2 and xe_i_3 series to p3
plot!(zeta_R_vec ./ 1000, SumPolZeta[xe_i_2, 14, :], linecolor=:red, linewidth=2, label="Prod = $e_val_2")
plot!(zeta_R_vec ./ 1000, SumPolZeta[xe_i_3, 14, :], linecolor=:green, linewidth=2, label="Prod = $e_val_3")

p4 = plot(zeta_R_vec ./ 1000, SumPolZeta[xe_i_1, 18, :], 
    title="Firm Value", titlefont=font(14), xlabel="Zeta_R (thousands)", ylabel="FirmValue", 
    legend=:bottomright, linecolor=:blue, linewidth=2, label="Prod = $e_val_1")

# Add xe_i_2 and xe_i_3 series to p4
plot!(zeta_R_vec ./ 1000, SumPolZeta[xe_i_2, 18, :], linecolor=:red, linewidth=2, label="Prod = $e_val_2")
plot!(zeta_R_vec ./ 1000, SumPolZeta[xe_i_3, 18, :], linecolor=:green, linewidth=2, label="Prod = $e_val_3")

# Combine the plots into a 2x2 layout
plot(p1, p3, p2, p4, layout=(2, 2), size=(1000, 800))

end