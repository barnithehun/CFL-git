wage = 1
phi_c = 0.4073866012405985;
zeta_Rs = 3882.0297136399195

# Cases are added manually now - if you change the number of cases, adjsut the parmater 'run'
nrun = 3

# ABL case
SumPol, e_chain, Fmat =  FirmOptim(wage, 0, zeta_Rs, zeta_Rl)
SumPolBig = fill(NaN,  size(SumPol, 1) ,  size(SumPol, 2) , nrun);
FmatBig = fill(NaN,  size(Fmat, 1) ,  size(Fmat, 2) , nrun);

SumPolBig[:,:,1] = SumPol
FmatBig[:,:,1] = Fmat

# Baseline case
SumPol, e_chain, Fmat =  FirmOptim(wage, phi_c, zeta_Rs, zeta_Rl)
SumPolBig[:,:,2] = SumPol
FmatBig[:,:,2] = Fmat

# Zeta = half case 
FirmOptim(wage, phi_c, zeta_Rs/4, zeta_Rl/4)
SumPolBig[:,:,3] = SumPol
FmatBig[:,:,3] = Fmat

# Zeta = quarter case
#SumPol, e_chain, Fmat =  FirmOptim(wage, phi_c, zeta_Rs, zeta_Rl)
#SumPolBig[:,:,4] = SumPol
#FmatBig[:,:,4] = Fmat

# Zeta = tenth case
#SumPol, e_chain, Fmat = FirmOptim(wage, phi_c, zeta_Rs, zeta_Rl)
#SumPolBig[:,:,5] = SumPol
# FmatBig[:,:,5] = Fmat


# simulates a random draw from Fvec
function next_si(Fvec)
    
    r = rand()  # uniform distribution between 0-1
    cumulative_prob = 0.0
    
    for (i, prob) in enumerate(Fvec)
        cumulative_prob += prob
        if r < cumulative_prob
            return i
        end
    end 

    return length(Fvec)  
    
end



# Setting up the estimation
simn_length = 100000
simt_length = 10
e_i = 19


x_size, _, _, _ = gridsize()
x_i = findall(x -> x == 0.0, unique(SumPol[:, 1]))[1]

Bsimmat = fill(NaN, simn_length, simt_length, nrun)
Qsimmat = fill(NaN, simn_length, simt_length, nrun)
Vsimmat = fill(NaN, simn_length, simt_length, nrun)
Lsimmat = fill(NaN, simn_length, simt_length, nrun)

for run in 1:nrun  
    for simn in 1:simn_length 
        
        s_i = x_i + (e_i-1)*x_size
    
        # Loop for 0
        for simt in 1:simt_length 
    
            xpol = SumPolBig[s_i, 6, run] + SumPolBig[s_i, 7, run] 
    
            if xpol != 1

                            
                # B
                if simt == 1
                    Bsimmat[simn, simt, run] = SumPolBig[s_i, 4, run]
                else
                    # Fmat is the transposed(!) transition matrix   
                    Fvec = FmatBig[: , s_i, run]
                    s_i = next_si(Fvec)
                    Bsimmat[simn, simt, run] = SumPolBig[s_i, 4, run]
                end
                
                # Q
                if simt == 1
                    Qsimmat[simn, simt, run] = SumPolBig[s_i, 9, run]
                else
                    # Fmat is the transposed(!) transition matrix   
                    Fvec = FmatBig[: , s_i, run]
                    s_i = next_si(Fvec)
                    Qsimmat[simn, simt, run] = SumPolBig[s_i, 9, run]
                end
                

                # Val
                if simt == 1
                    Vsimmat[simn, simt, run] = SumPolBig[s_i, end, run]
                else
                    # Fmat is the transposed(!) transition matrix   
                    Fvec = FmatBig[: , s_i, run]
                    s_i = next_si(Fvec)
                    Vsimmat[simn, simt, run] = SumPolBig[s_i, end, run]
                end
                

                # Uncond. Pliq
                if simt == 1
                    Lsimmat[simn, simt, run] = SumPolBig[s_i, 14, run]
                else
                    # Fmat is the transposed(!) transition matrix   
                    Fvec = FmatBig[: , s_i, run]
                    s_i = next_si(Fvec)
                    Lsimmat[simn, simt, run] = SumPolBig[s_i, 14, run]
                end
    
            else
                break
            end
        end
    end
end

# 4 stands for the  for variables I am studying
Bmeanv = zeros(simt_length, 1, nrun) 
Qmeanv = zeros(simt_length, 1, nrun) 
Vmeanv = zeros(simt_length, 1, nrun) 
Lmeanv = zeros(simt_length, 1, nrun) 

for run in 1:nrun 
    for simt in 1:simt_length 

        Bmeanv[simt, 1, run] = mean(filter(!isnan, Bsimmat[:, simt, run]))
        Qmeanv[simt, 1, run] = mean(filter(!isnan, Qsimmat[:, simt, run]))
        Vmeanv[simt, 1, run] = mean(filter(!isnan, Vsimmat[:, simt, run]))
        Lmeanv[simt, 1, run] = mean(filter(!isnan, Lsimmat[:, simt, run]))

    end
end


x_axis = 1:simt_length
function create_plot(matrix, title_label)
    p = plot(title=title_label, xlabel="simt", ylabel="Mean Value", legend=false)
    for run in 1:nrun
        plot!(x_axis, [mean(filter(!isnan, matrix[:, t, run])) for t in x_axis], label="", lw=1, alpha=0.5)
    end
    return p
end

# Generate individual plots
p1 = create_plot(Bsimmat, "Bmeanv")
p2 = create_plot(Qsimmat, "Qmeanv")
p3 = create_plot(Vsimmat, "Vmeanv")
p4 = create_plot(Lsimmat, "Lmeanv")

# Combine plots in a grid layout
plot(p1, p2, p3, p4, layout=(2,2), size=(1200, 800))


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






########################## Improvements: shifting costs
c_e_baseline = 537.5330125606192

zeta_Rs_vec = 200:200:2600

b_prodimp = zeros( length(zeta_Rs_vec));
b_liqprob_SME = zeros( length(zeta_Rs_vec));
b_liqprob_LE = zeros( length(zeta_Rs_vec));
b_d2c_SME = zeros( length(zeta_Rs_vec));
b_d2c_LE = zeros( length(zeta_Rs_vec));
b_intrate_SME = zeros( length(zeta_Rs_vec));
b_intrate_LE = zeros( length(zeta_Rs_vec));
b_CFrel_SME = zeros( length(zeta_Rs_vec));
b_CFrel_LE = zeros( length(zeta_Rs_vec));

for (ind_s, zeta_Rl) in enumerate(zeta_Rs_vec)

    zeta_Rs = 0

    SumPol, e_chain, Fmat = FirmOptim(wage, phi_c, zeta_Rl, zeta_Rs)
    c_e, f0 = EntryValue(SumPol, e_chain) ;
    mu, m, xpol = stat_dist(SumPol, Fmat, f0);

    b_prodimp[ind_s] = c_e / c_e_baseline
    b_liqprob_SME[ind_s] = sumSSsme(SumPol,Fmat,f0)[10,1]
    b_liqprob_LE[ind_s] = sumSSsme(SumPol,Fmat,f0)[10,2]
    b_d2c_SME[ind_s] = sumSSsme(SumPol,Fmat,f0)[7,1]
    b_d2c_LE[ind_s] = sumSSsme(SumPol,Fmat,f0)[7,2]
    b_intrate_SME[ind_s] = sumSSsme(SumPol,Fmat,f0)[8,1]
    b_intrate_LE[ind_s] = sumSSsme(SumPol,Fmat,f0)[8,2]
    b_CFrel_SME[ind_s] = sumSSsme(SumPol,Fmat,f0)[12,1]
    b_CFrel_LE[ind_s] = sumSSsme(SumPol,Fmat,f0)[12,2]

end


function create_long_df(b_var, var, var_name, ind_s)
    return DataFrame(
        ind_s = repeat(ind_s, 2),  # Repeat indices for both baseline & counterfactual
        variable = var_name,
        type = vcat(fill("Shifting costs", length(ind_s)), fill("Reducing Costs", length(ind_s))),
        value = vcat(var, b_var)  # Stack the values
    )
end

indices = collect(1:length(prodimp))

df_prodimp = create_long_df(b_prodimp, prodimp, "prodimp", indices)
df_liqprob_SME = create_long_df(b_liqprob_SME, liqprob_SME, "liqprob_SME", indices)
df_liqprob_LE = create_long_df(b_liqprob_LE, liqprob_LE, "liqprob_LE", indices)
df_d2c_SME = create_long_df(b_d2c_SME, d2c_SME, "d2c_SME", indices)
df_d2c_LE = create_long_df(b_d2c_LE, d2c_LE, "d2c_LE", indices)
df_intrate_SME = create_long_df(b_intrate_SME, intrate_SME, "intrate_SME", indices)
df_intrate_LE = create_long_df(b_intrate_LE, intrate_LE, "intrate_LE", indices)
df_CFrel_SME = create_long_df(b_CFrel_SME, CFrel_SME, "CFrel_SME", indices)
df_CFrel_LE = create_long_df(b_CFrel_LE, CFrel_LE, "CFrel_LE", indices)

zeta_Rs = [zeta_Rs_vec; zeta_Rs_vec] 
# Define font size properties
title_fontsize = 20
xlabel_fontsize = 16
ylabel_fontsize = 16
legend_fontsize = 18
xtick_fontsize = 14  # Set for x-axis ticks
ytick_fontsize = 14  # Set for y-axis ticks

# Now apply to the plots
p1 = groupedbar(zeta_Rs, (df_prodimp.value .- 1) ./ 2, group=df_prodimp.type, 
                xlabel="Reorg. Cost - Small Firms", ylabel="Productivity relative to Baseline", 
                title="Productivity Improvements", palette=colors, 
                ylims=(minimum((df_prodimp.value .- 1) ./ 2), maximum((df_prodimp.value .- 1) ./ 2)), 
                size=(600,600), legendfontsize=legend_fontsize, titlefontsize=title_fontsize, 
                xlabelfontsize=xlabel_fontsize, ylabelfontsize=ylabel_fontsize,
                xtickfontsize=xtick_fontsize, ytickfontsize=ytick_fontsize)
Plots.savefig("shiftingCosts1.png")

p2 = groupedbar(zeta_Rs, df_liqprob_SME.value, group=df_liqprob_SME.type, 
                xlabel="Reorg. Cost - Small Firms", ylabel="Liquidation Probability", 
                title="Liquidation Probability - SMEs", palette=colors, 
                ylims=(minimum(df_liqprob_SME.value), maximum(df_liqprob_SME.value)),
                titlefontsize=title_fontsize, xlabelfontsize=xlabel_fontsize, 
                ylabelfontsize=ylabel_fontsize, xtickfontsize=xtick_fontsize, 
                ytickfontsize=ytick_fontsize)

p3 = groupedbar(zeta_Rs, df_liqprob_LE.value, group=df_liqprob_LE.type, 
                xlabel="Reorg. Cost - Small Firms", ylabel="Liquidation Probability", 
                title="Liquidation Probability - LEs", palette=colors, 
                ylims=(minimum(df_liqprob_LE.value), maximum(df_liqprob_LE.value)),
                titlefontsize=title_fontsize, xlabelfontsize=xlabel_fontsize, 
                ylabelfontsize=ylabel_fontsize, xtickfontsize=xtick_fontsize, 
                ytickfontsize=ytick_fontsize)

p4 = groupedbar(zeta_Rs, df_d2c_SME.value, group=df_d2c_SME.type, 
                xlabel="Reorg. Cost - Small Firms", ylabel="Debt to Collateral", 
                title="Debt-to-Collateral - SMEs", palette=colors, 
                ylims=(minimum(df_d2c_SME.value), maximum(df_d2c_SME.value)),
                titlefontsize=title_fontsize, xlabelfontsize=xlabel_fontsize, 
                ylabelfontsize=ylabel_fontsize, xtickfontsize=xtick_fontsize, 
                ytickfontsize=ytick_fontsize)

p5 = groupedbar(zeta_Rs, df_d2c_LE.value, group=df_d2c_LE.type, 
                xlabel="Reorg. Cost - Small Firms", ylabel="Debt to Collateral", 
                title="Debt-to-Collateral - LEs", palette=colors, 
                ylims=(minimum(df_d2c_LE.value), maximum(df_d2c_LE.value)),
                titlefontsize=title_fontsize, xlabelfontsize=xlabel_fontsize, 
                ylabelfontsize=ylabel_fontsize, xtickfontsize=xtick_fontsize, 
                ytickfontsize=ytick_fontsize)

p6 = groupedbar(zeta_Rs, df_intrate_SME.value, group=df_intrate_SME.type, 
                xlabel="Reorg. Cost - Small Firms", ylabel="Interest rate", 
                title="Interest Rate - SMEs", palette=colors, 
                ylims=(minimum(df_intrate_SME.value), maximum(df_intrate_SME.value)),
                titlefontsize=title_fontsize, xlabelfontsize=xlabel_fontsize, 
                ylabelfontsize=ylabel_fontsize, xtickfontsize=xtick_fontsize, 
                ytickfontsize=ytick_fontsize)

p7 = groupedbar(zeta_Rs, df_intrate_LE.value, group=df_intrate_LE.type, 
                xlabel="Reorg. Cost - Small Firms", ylabel="Interest rate", 
                title="Interest Rate - LEs", palette=colors, 
                ylims=(minimum(df_intrate_LE.value), maximum(df_intrate_LE.value)),
                titlefontsize=title_fontsize, xlabelfontsize=xlabel_fontsize, 
                ylabelfontsize=ylabel_fontsize, xtickfontsize=xtick_fontsize, 
                ytickfontsize=ytick_fontsize)

p8 = groupedbar(zeta_Rs, df_CFrel_SME.value, group=df_CFrel_SME.type, 
                xlabel="Reorg. Cost - Small Firms", ylabel="CF Reliance", 
                title="CFL Reliance - SMEs", palette=colors, 
                ylims=(minimum(df_CFrel_SME.value), maximum(df_CFrel_SME.value)),
                titlefontsize=title_fontsize, xlabelfontsize=xlabel_fontsize, 
                ylabelfontsize=ylabel_fontsize, xtickfontsize=xtick_fontsize, 
                ytickfontsize=ytick_fontsize)

p9 = groupedbar(zeta_Rs, df_CFrel_LE.value, group=df_CFrel_LE.type, 
                xlabel="Reorg. Cost - Small Firms", ylabel="CF Reliance", 
                title="CFL Reliance - LEs", palette=colors, 
                ylims=(minimum(df_CFrel_LE.value), maximum(df_CFrel_LE.value)),
                titlefontsize=title_fontsize, xlabelfontsize=xlabel_fontsize, 
                ylabelfontsize=ylabel_fontsize, xtickfontsize=xtick_fontsize, 
                ytickfontsize=ytick_fontsize)



plot(p2, p3, p8, p9, layout=(2, 2), size=(900,800), margin=5mm)
Plots.savefig("shiftingCosts2.png")


plot(p4, p5, p7, p7, layout=(2, 2), size=(900,800), margin=5mm)
Plots.savefig("shiftingCosts3.png")

