
#=
########################## Improvements For Small firms
c_e_baseline = 537.5330125606192

zeta_Rl_vec = 0:500:3000
zeta_Rs_vec = 0:500:3000

prodimp = zeros(length(zeta_Rl_vec), length(zeta_Rs_vec))
liqprob_SME = zeros(length(zeta_Rl_vec), length(zeta_Rs_vec))
liqprob_LE = zeros(length(zeta_Rl_vec), length(zeta_Rs_vec))
d2c_SME = zeros(length(zeta_Rl_vec), length(zeta_Rs_vec))
d2c_LE = zeros(length(zeta_Rl_vec), length(zeta_Rs_vec))
CFrel_SME = zeros(length(zeta_Rl_vec), length(zeta_Rs_vec))
CFrel_LE = zeros(length(zeta_Rl_vec), length(zeta_Rs_vec))

for (ind_l, zeta_Rl) in enumerate(zeta_Rl_vec)
    for (ind_s, zeta_Rs) in enumerate(zeta_Rl_vec)

        
        SumPol, e_chain, Fmat = FirmOptim(wage, phi_c, zeta_Rs, zeta_Rl)
        c_e, f0 = EntryValue(SumPol, e_chain) ;
        mu, m, xpol = stat_dist(SumPol, Fmat, f0);

        prodimp[ind_l, ind_s] = c_e / c_e_baseline
        liqprob_SME[ind_l, ind_s] = sumSSsme(SumPol,Fmat,f0)[10,1]
        liqprob_LE[ind_l, ind_s] = sumSSsme(SumPol,Fmat,f0)[10,2]
        d2c_SME[ind_l, ind_s] = sumSSsme(SumPol,Fmat,f0)[7,1]
        d2c_LE[ind_l, ind_s] = sumSSsme(SumPol,Fmat,f0)[7,2]
        CFrel_SME[ind_l, ind_s] = sumSSsme(SumPol,Fmat,f0)[12,1]
        CFrel_LE[ind_l, ind_s] = sumSSsme(SumPol,Fmat,f0)[12,2]

    end
end

# Define the red point coordinates
x_target = 2590
y_target = 2590

ps = heatmap(zeta_Rs_vec, zeta_Rl_vec, (prodimp .-1) ./ 2,
    xlabel="Reorg. Cost - Large Firms", ylabel="Reorg. Cost - Small Firms",
    color=:viridis, title="Productivity Improvements", size=(600,400))
scatter!(ps, [x_target], [y_target], color=:red, markersize=5, label=nothing)
savefig("productivity_heatmap.png")

# Create individual heatmaps and add the red point without legend
p1 = heatmap(zeta_Rs_vec, zeta_Rl_vec, liqprob_SME,
    xlabel="Reorg. Cost - Large Firms", ylabel="Reorg. Cost - Small Firms",
    color=:viridis, title="Liquidation Probability - SMEs")
scatter!(p1, [x_target], [y_target], color=:red, markersize=5, label=nothing)

p2 = heatmap(zeta_Rs_vec, zeta_Rl_vec, liqprob_LE,
    xlabel="Reorg. Cost - Large Firms", ylabel="Reorg. Cost - Small Firms",
    color=:viridis, title="Liquidation Probability - LEs")
scatter!(p2, [x_target], [y_target], color=:red, markersize=5, label=nothing)

p3 = heatmap(zeta_Rs_vec, zeta_Rl_vec, d2c_SME,
    xlabel="Reorg. Cost - Large Firms", ylabel="Reorg. Cost - Small Firms",
    color=:viridis, title="Debt-to-Collateral - SMEs")
scatter!(p3, [x_target], [y_target], color=:red, markersize=5, label=nothing)

p4 = heatmap(zeta_Rs_vec, zeta_Rl_vec, d2c_LE,
    xlabel="Reorg. Cost - Large Firms", ylabel="Reorg. Cost - Small Firms",
    color=:viridis, title="Debt-to-Collateral - LEs")
scatter!(p4, [x_target], [y_target], color=:red, markersize=5, label=nothing)

p5 = heatmap(zeta_Rs_vec, zeta_Rl_vec, CFrel_SME,
    xlabel="Reorg. Cost - Large Firms", ylabel="Reorg. Cost - Small Firms",
    color=:viridis, title="CFL Reliance - SMEs")
scatter!(p5, [x_target], [y_target], color=:red, markersize=5, label=nothing)

p6 = heatmap(zeta_Rs_vec, zeta_Rl_vec, CFrel_LE,
    xlabel="Reorg. Cost - Large Firms", ylabel="Reorg. Cost - Small Firms",
    color=:viridis, title="CFL Reliance - LEs")
scatter!(p6, [x_target], [y_target], color=:red, markersize=5, label=nothing)

# Plot combined figure
plot(p1, p5, p2, p6, layout=(2,2), size=(1100,700), margin=5mm)
savefig("combined_heatmaps.png")




=#

############ The effects of reform across firm sizes
using Plots

function plotPolHeatmaps(SumPol, SumPol0, mu, muflag)
    # Filtering rows based on conditions
    keepvec = (SumPol[:, 1] .>= 0) .&& (SumPol[:, 1] .<=  50000) .&&  (SumPol[:, 2] .>= 3)
    
    x = SumPol0[keepvec, 1]
    e = SumPol0[keepvec, 2]
    
    SumPol = SumPol[keepvec, :]
    SumPol0 = SumPol0[keepvec, :]
    mu = mu[keepvec]
    
    # Compute required variables
    b2a = (SumPol0[:, 4] ./ (SumPol0[:, 3] .+ SumPol0[:, 1])) .- (SumPol[:, 4] ./ (SumPol[:, 3] .+ SumPol[:, 1]))
    b2a .= ifelse.(b2a .< -0.3 .|| b2a .> 0.3 , NaN, b2a)
    
    CFL = SumPol0[:, 17] .- SumPol[:, 17]
    gam = SumPol0[:, 14] .- SumPol[:, 14]
    q = ((1 ./SumPol0[:, 9]) .- (1 ./SumPol[:, 9])) .* 100
    q .= ifelse.(q .< -4 .|| q .> 4, NaN, q)
    
    # Mask values where mu is too small
    if muflag == 1
        mask = mu .<= 0.000001
        b2a[mask] .= NaN
        CFL[mask] .= NaN
        gam[mask] .= NaN
        q[mask] .= NaN
    end
        
    function plot_heatmap(x, e, values, title)
        x_labels = unique(string.(round.(Int, x)))
        e_labels = unique(string.(round.(e, digits=1)))
    
        resmat = fill(NaN, length(e_labels), length(x_labels))
    
        for e_i in 1:length(e_labels)
            for x_i in 1:length(x_labels)
                ind = (e_i - 1) * length(x_labels) + x_i 
                resmat[e_i, x_i] = values[ind]
            end
        end
    
        heatmap(x_labels, e_labels, resmat, color=:viridis, 
        xlabel="Cash on Hand", ylabel="Productivity", title=title, 
        legend=true, xrotation=45, xlabelfontsize=12, 
        ylabelfontsize=12, tickfontsize=9,  titlefontsize=13)

        
    end
    
    plt1 = plot_heatmap(x, e, gam, "A: Liquidation Probability")
    plt2 = plot_heatmap(x, e, CFL, "B: CFL reliance")
    plt3 = plot_heatmap(x, e, b2a, "C: Debt to Collateral")
    plt4 = plot_heatmap(x, e, q, "D: Interest Rate")
    
    plot(plt1, plt2, plt3, plt4, layout=(2, 2), size=(900, 600), margin=5mm)
end


plotPolHeatmaps(SumPol, SumPol0, mu, 1)




########################## Shifting variable costs 

wage = 1;
zeta_Rl =  1295       # reform fixed costs
c_e_baseline = 537.53 # entry value in the baseline case

shift_vec = 0:0.02:0.2

prodimp = zeros(length(shift_vec))
liqprob = zeros(length(shift_vec))
d2c = zeros(length(shift_vec))
CFrel = zeros(length(shift_vec))
intrate = zeros(length(shift_vec))

liqprob_SME = zeros(length(shift_vec))
liqprob_LE = zeros(length(shift_vec))
d2c_SME = zeros(length(shift_vec))
d2c_LE = zeros(length(shift_vec))
CFrel_SME = zeros(length(shift_vec))
CFrel_LE = zeros(length(shift_vec))

for (ind, shift) in enumerate(shift_vec)

        
    SumPol, e_chain, Fmat = FirmOptim(wage, zeta_Rl, shift)
    c_e, f0 = EntryValue(SumPol, e_chain) ;
    mu, m, xpol = stat_dist(SumPol, Fmat, f0);

    prodimp[ind] = (c_e / c_e_baseline) / 2
    liqprob[ind] = sumSS(SumPol,Fmat,f0)[10,1]
    d2c[ind] = sumSS(SumPol,Fmat,f0)[7,1]
    CFrel[ind] = sumSS(SumPol,Fmat,f0)[13,1]
    intrate[ind] = sumSS(SumPol,Fmat,f0)[8,1]

    liqprob_SME[ind] = sumSSsme(SumPol,Fmat,f0)[10,1]
    liqprob_LE[ind] = sumSSsme(SumPol,Fmat,f0)[10,2]
    d2c_SME[ind] = sumSSsme(SumPol,Fmat,f0)[7,1]
    d2c_LE[ind] = sumSSsme(SumPol,Fmat,f0)[7,2]
    CFrel_SME[ind] = sumSSsme(SumPol,Fmat,f0)[12,1]
    CFrel_LE[ind] = sumSSsme(SumPol,Fmat,f0)[12,2]

end

p1 = bar(shift_vec, (prodimp .* 2 .- 1) ./ 2, title="ProdImp * 2", xlabel="Shift Vec", legend=false)
p2 = bar(shift_vec, liqprob, title="LiqProb", xlabel="Shift Vec", legend=false)
p3 = bar(shift_vec, d2c, title="D2C", xlabel="Shift Vec", legend=false)
p4 = bar(shift_vec, CFrel, title="CFRel", xlabel="Shift Vec", legend=false)
p5 = bar(shift_vec, intrate, title="IntRate", xlabel="Shift Vec", legend=false)
p6 = bar(shift_vec, liqprob_SME, title="LiqProb_SME", xlabel="Shift Vec", legend=false)
p7 = bar(shift_vec, liqprob_LE, title="LiqProb_LE", xlabel="Shift Vec", legend=false)
p8 = bar(shift_vec, d2c_SME, title="D2C_SME", xlabel="Shift Vec", legend=false)
p9 = bar(shift_vec, d2c_LE, title="D2C_LE", xlabel="Shift Vec", legend=false)
p10 = bar(shift_vec, CFrel_SME, title="CFRel_SME", xlabel="Shift Vec", legend=false)
p11 = bar(shift_vec, CFrel_LE, title="CFRel_LE", xlabel="Shift Vec", legend=false)

# Arrange plots in a grid (4 rows, 3 columns)
plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, layout=(4, 3), size=(1200, 900))