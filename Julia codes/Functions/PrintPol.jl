# Print policies
function PrintPol()    
    mu, _, _ = stat_dist(SumPol, Fmat, f0)
    column_names = [:x, :e, :k, :b, :next_x, :exit, :def, :pdef, :q, :l, :y, :Pi, :d, :gam, :Pi_liq, :Pi_reo, :tau, :val, :SSpercent]
    SumPoll = hcat(SumPol[1:end, :],round.(mu ./ sum(mu) .* 100, digits=4))
    SumPol_df = DataFrame(SumPoll, column_names)
    time_string = "$(Dates.day(now()))_$(Dates.hour(now()))_$(Dates.minute(now()))_$(Dates.second(now()))"
    XLSX.writetable("$time_string.xlsx", sheetname="sheet1", SumPol_df)
end