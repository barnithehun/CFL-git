############ E18 and E19: two aggregate feedbacks on the zeta/3 reform ############
# E18 fire sales: the resale value of capital falls when more capital is liquidated,
#       phi_a = phi_a0 * (liquidated capital / its baseline level)^(-eta).
# E19 court congestion: the cost of a reorganization rises with the number of cases,
#       zeta = (statutory zeta) * (reorganizations / their baseline level)^(eta).
# Both are 1 at the baseline, so the baseline is unchanged, and the reform is a fixed point: solve the GE, update the
# price or the cost from the flows, repeat (damped) until it settles. The elasticities are not calibrated: a range.
include(joinpath(@__DIR__, "common.jl"))
dir = outdir("E18")
B = base_eq(); ob = outcomes(B); zeta, phi0 = B.params.zeta_Rl, B.params.phi_a
o3 = outcomes(reform_eq(3))

# fixed point of x = f(outcomes at x), with damping; x0 is the value without the feedback
function fixed_point(solve, f, x0; tol = 2e-4, maxit = 8, damp = 0.5)
    x, o = x0, nothing
    for it in 1:maxit
        o = solve(x)
        x_new = f(o)
        println("  iteration $it: x = ", round(x, digits = 5), " -> implied ", round(x_new, digits = 5))
        abs(x_new - x) <= tol * abs(x0) && return (x = x, o = o, it = it, converged = true)
        x = damp * x + (1 - damp) * x_new
    end
    return (x = x, o = o, it = maxit, converged = false)
end

etas = [0.25, 0.5, 1.0]
fire = cached("fire_sale") do
    map(etas) do eta
        println("fire sales, eta = $eta")
        fixed_point(phi -> outcomes(ge((zeta_Rl = zeta / 3, phi_a = phi, phi_af = phi); bracket = (0.8, 1.3))),
                    o -> phi0 * (o["liqK"] / ob["liqK"])^(-eta), phi0)
    end
end
cong = cached("congestion") do
    map(etas) do eta
        println("congestion, eta = $eta")
        fixed_point(z -> outcomes(ge((zeta_Rl = z,); bracket = (0.8, 1.3))),
                    o -> (zeta / 3) * (o["reorganizations"] / ob["reorganizations"])^eta, zeta / 3)
    end
end

cols = [SWEEP_COLS[1:5]; "liq. capital %"; "reorganizations %"; "phi_a"; "zeta / base"; "iterations"; "converged"]
row(o, phi, z, it, conv) = [sweep_row(o, ob)[1:5]; pct(o["liqK"], ob["liqK"]); pct(o["reorganizations"], ob["reorganizations"]); phi; z / zeta; it; conv]
save_table(dir, "feedbacks.csv", cols,
           [["zeta/3, no feedback" => row(o3, phi0, zeta / 3, 0, 1)];
            ["fire sales, eta $eta" => row(r.o, r.x, zeta / 3, r.it, r.converged) for (eta, r) in zip(etas, fire)];
            ["congestion, eta $eta" => row(r.o, phi0, r.x, r.it, r.converged) for (eta, r) in zip(etas, cong)]])

plt = small_plot(title = "Y/L gain of zeta/3 with feedbacks, %", xlabel = "elasticity eta")
plot!(plt, [0.0; etas], [pct(o3["wage"], ob["wage"]); [pct(r.o["wage"], ob["wage"]) for r in fire]], marker = :circle, label = "fire sales")
plot!(plt, [0.0; etas], [pct(o3["wage"], ob["wage"]); [pct(r.o["wage"], ob["wage"]) for r in cong]], marker = :square, label = "court congestion")
savefig(plt, joinpath(dir, "feedbacks.png"))
