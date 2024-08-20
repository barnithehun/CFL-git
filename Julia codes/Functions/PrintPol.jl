# Print policies
function PrintPol(SumPol, mu, SumPol0, mu0)    


    column_names = [:x, :e, :k, :b, :next_x, :exit, :def, :pdef, :q, :l, :y, :Pi, :d, :gam, :Pi_liq, :Pi_reo, :tau, :val, :SSpercent]

    
    SumPoll = hcat(SumPol, round.(mu ./ sum(mu) .* 100, digits=4))
    SumPol_df = DataFrame(SumPoll, column_names)
    time_string = "CFL_$(Dates.day(now()))_$(Dates.hour(now()))_$(Dates.minute(now()))"
    XLSX.writetable("$time_string.xlsx", sheetname="sheet1", SumPol_df)

    #=
    SumPoll = hcat(SumPol0, round.(mu0 ./ sum(mu) .* 100, digits=4))
    SumPol_df = DataFrame(SumPoll, column_names)
    time_string = "ABL_$(Dates.day(now()))_$(Dates.hour(now()))_$(Dates.minute(now()))"
    XLSX.writetable("$time_string.xlsx", sheetname="sheet1", SumPol_df)

    SumPoll = hcat(SumPol, round.(mu ./ sum(mu) .* 100, digits=4)) .- hcat(SumPol0, round.(mu0 ./ sum(mu) .* 100, digits=4))
    SumPol_df = DataFrame(SumPoll, column_names)
    time_string = "DIFF_$(Dates.day(now()))_$(Dates.hour(now()))_$(Dates.minute(now()))"
    XLSX.writetable("$time_string.xlsx", sheetname="sheet1", SumPol_df)
    =#

end