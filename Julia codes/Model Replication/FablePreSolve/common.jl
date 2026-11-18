############ FablePreSolve: common code for the extension files E01-E22 ############
# Loads the model (Model_ext.jl: Baseline_noise46.jl plus two switches that are off by default) and Functions/Results.jl,
# and defines what the extensions share:
#   base_eq(), reform_eq(d)   the baseline and the zeta/d reform at its GE wage (cached in out/cache)
#   events(E)                 what happens next period to every continuing firm (production, repayment, default,
#                             liquidation or reorganization), exactly as FirmOptim's Fmat has it
#   kernel, statevecs         the same events as transition matrices and as per-firm expected flows by state
#   welfare(E)                main5's C = Y - I - Psi, from the events
#   path(E, d0)               the perfect-foresight path of the firm distribution from arrivals d0
#   small helpers for tables and plots
# Each E file starts with include("common.jl") and writes its tables and figures to out/<name>/.

ENV["GKSwstype"] = "100"     # GR without a display
using Serialization, Printf

const FPS_DIR = @__DIR__
const MODEL_DIR = dirname(FPS_DIR)
const OUT_DIR = joinpath(FPS_DIR, "out")
let code = read(joinpath(FPS_DIR, "Model_ext.jl"), String)
    # loaded under a file name in Model Replication, so its cd(@__DIR__) and its includes of Functions/ resolve there
    fake = joinpath(MODEL_DIR, "Model_ext.jl")
    task_local_storage(() -> include_string(Main, code, fake), :SOURCE_PATH, fake)
end
include(joinpath(MODEL_DIR, "Functions", "Results.jl"))

outdir(name) = mkpath(joinpath(OUT_DIR, name))

# cached(f, name): f()'s result, saved in out/cache/<name>.jls (delete the folder after changing Model_ext.jl)
function cached(f, name)
    file = joinpath(OUT_DIR, "cache", name * ".jls")
    isfile(file) && return deserialize(file)
    r = f()
    mkpath(dirname(file)); serialize(file, r)
    return r
end
base_eq() = cached(() -> (reset_params!(); solve_eq()), "base")
# the zeta/d reform for all firms, at the wage where free entry holds with the baseline entry cost
reform_eq(d) = cached(() -> (B = base_eq(); ge_eq((zeta_Rl = B.params.zeta_Rl / d,), B.c_e)), "reform_zeta_over_$d")
ge(vals; kw...) = ge_eq(vals, base_eq().c_e; kw...)

continuing(E) = [E.SumPol[s, 6] == 0 && E.SumPol[s, 7] == 0 for s in 1:size(E.SumPol, 1) - 1]
# mass of firms that continue, by state (the last state, default, is left out everywhere below)
cont_mass(E) = vec(E.mu)[1:end-1] .* continuing(E)
# entrants that stay, per unit of potential entrants
entrant_dist(E) = ((1 .- vec(E.xpol)) .* vec(E.f0))[1:end-1]

# events(E): one row per (continuing state s, next productivity e', x' grid point), with the probability of the row
# for a firm in s (prob), the default probability at that grid point (pdef: chosen there, or the exogenous shock), the
# liquidation probability in default (liq), and where a reorganized firm lands (r1 with weight rw, r2 with the rest).
# rev_shock: an unexpected loss of that share of revenue next period, with wages paid in full (E22).
function events(E; rev_shock = 0.0)
    S, p, w = E.SumPol, E.params, E.wage
    n = size(S, 1) - 1
    x_size, e_size = E.gridsize.x_size, E.gridsize.e_size
    alpha, nu, delta, pc = p.alpha, p.nu, p.delta, p.pc
    x_grid = S[1:x_size, 1]                         # x varies fastest in the state index
    e_vals = exp.(E.e_chain.state_values)
    P = E.e_chain.p
    V = S[:, 18]
    wedge = hasproperty(p, :liq_wedge) ? p.liq_wedge : 0.0

    fn_L(k, e) = (nu * e * k^alpha / w)^(1 / (1 - nu))
    fn_Y(k, e) = e * k^alpha * fn_L(k, e)^nu
    fn_Gam(val, k, zeta_R) = p.sigma_v == 0 ? Float64(p.phi_a * (1 - delta) * k >= p.phi_c_hh * val * exp(wedge) - zeta_R) :
                             p.phi_c_hh * val <= 0 ? 1.0 : (p.phi_a * (1 - delta) * k + zeta_R) <= 0 ? 0.0 :
                             cdf(Normal(), (log((p.phi_a * (1 - delta) * k + zeta_R) / (p.phi_c_hh * val)) - wedge) / p.sigma_v)
    function V_at(x, e_i)                           # V at x, interpolated as in FirmOptim's q update
        base = (e_i - 1) * x_size
        j = searchsortedlast(x_grid, x)
        1 <= j < x_size || return V[(j <= 0 ? 1 : x_size) + base]
        t = (x - x_grid[j]) / (x_grid[j+1] - x_grid[j])
        (1 - t) * V[j + base] + t * V[j + 1 + base]
    end

    c = (s = Int[], en = Int[], land = Int[], prob = Float64[], pdef = Float64[], chosen = Float64[], liq = Float64[],
         k = Float64[], b = Float64[], y = Float64[], lab = Float64[], profit = Float64[], xnext = Float64[],
         vreo = Float64[], zeta = Float64[], theta = Float64[], large = Bool[], tau = Float64[],
         r1 = Int[], rw = Float64[], r2 = Int[])
    cont = continuing(E)
    for s in 1:n
        cont[s] || continue
        k, b = S[s, 3], S[s, 4]
        large = (p.size_by_debt == 1 ? b : k) >= p.Fcut
        zeta_R = large ? p.zeta_Rl : p.zeta_Rs
        theta = large ? p.phi_c : p.phi_c_s
        pdef_exo = large ? p.pdef_exo_l : p.pdef_exo_s
        e_i = (s - 1) ÷ x_size + 1
        for e_n in 1:e_size
            pt = P[e_i, e_n]
            pt == 0 && continue
            base = (e_n - 1) * x_size
            y = fn_Y(k, e_vals[e_n])
            profit = (1 - nu) * y - rev_shock * y - pc
            x_next = profit + (1 - delta) * k - b   # cash on hand after repaying b'
            val_reo = V_at(x_next + b, e_n)         # reorganization: a fresh start at x' + b'
            liq = fn_Gam(val_reo, k, zeta_R)
            r_close, rw, r_far = landing(x_next + b, x_grid)
            x_close, cw, x_far = landing(x_next, x_grid)
            for (x_i, wgt) in ((x_close, cw), (x_far, 1 - cw))
                (x_i == 0 || wgt == 0) && continue
                chosen = S[x_i + base, 7]
                push!(c.s, s); push!(c.en, e_n); push!(c.land, x_i + base); push!(c.prob, pt * wgt)
                push!(c.pdef, chosen + (1 - chosen) * pdef_exo); push!(c.chosen, chosen); push!(c.liq, liq)
                push!(c.k, k); push!(c.b, b); push!(c.y, y); push!(c.lab, fn_L(k, e_vals[e_n])); push!(c.profit, profit)
                push!(c.xnext, x_next); push!(c.vreo, val_reo); push!(c.zeta, zeta_R); push!(c.theta, theta)
                push!(c.large, large); push!(c.tau, S[s, 17])
                push!(c.r1, r_close + base); push!(c.rw, rw); push!(c.r2, r_far == 0 ? 0 : r_far + base)
            end
        end
    end
    return (; c..., n = n, x_grid = x_grid, e_vals = e_vals)
end

# what the lender recovers per unit of face value in a liquidation and in a reorganization, row by row (fn_Tau_Q's payoffs)
function recoveries(E, ev)
    p = E.params
    Pl = max.(0, p.phi_af .* (1 - p.delta) .* ev.k .- p.zeta_L)
    liq = min.(ev.b, (1 .- ev.tau) .* Pl .+ ev.tau .* p.kappa .* Pl)
    reo = min.(ev.b, (1 .- ev.tau) .* Pl .+ ev.tau .* ev.theta .* ev.vreo)
    return (liq = liq, reo = reo)
end

# kernel(ev): T[s', s] = probability that a firm continuing in s arrives in s' next period (repaying, or reorganized)
function kernel(ev)
    n = ev.n
    reo = ev.prob .* ev.pdef .* (1 .- ev.liq)
    two = ev.r2 .!= 0
    Trep = sparse(ev.land, ev.s, ev.prob .* (1 .- ev.pdef), n, n)
    Treo = sparse(ev.r1, ev.s, reo .* ev.rw, n, n) + sparse(ev.r2[two], ev.s[two], reo[two] .* (1 .- ev.rw[two]), n, n)
    return (T = Trep + Treo, Trep = Trep, Treo = Treo)
end

# by_state(ev, v): sum over a state's rows of prob * v, i.e. the expected value of v next period for a firm in that state
by_state(ev, v) = (out = zeros(ev.n); for i in eachindex(ev.s); out[ev.s[i]] += ev.prob[i] * v[i]; end; out)

# statevecs(E, ev): expected next-period flows of one firm continuing in each state
function statevecs(E, ev)
    d, l, r = ev.pdef, ev.pdef .* ev.liq, ev.pdef .* (1 .- ev.liq)
    S = E.SumPol
    return (Y = by_state(ev, ev.y), L = by_state(ev, ev.lab), def = by_state(ev, d), liq = by_state(ev, l), reo = by_state(ev, r),
            liqK = by_state(ev, l .* ev.k), liqPi = by_state(ev, l .* ev.profit), reoZeta = by_state(ev, r .* ev.zeta),
            reoV = by_state(ev, r .* ev.vreo), K = S[1:ev.n, 3] .* continuing(E), B = S[1:ev.n, 4] .* continuing(E),
            exit = Float64.(S[1:ev.n, 6] .!= 0))
end

# resources(E, sv, c, c_next, m_next): what is left for consumption in the period in which the firms c (continuing mass by
# state, chosen last period) produce, when c_next continue into the following period and m_next potential entrants enter.
#   C = Y - fixed operating costs - entry costs - investment - profits and capital lost in liquidation - reorganization costs
# In a steady state c_next = c and investment is delta K plus the capital destroyed in liquidations.
function resources(E, sv, c, c_next, m_next; sv_next = sv)
    p = E.params
    r = Dict{String, Float64}()
    r["Y"] = sv.Y' * c
    K, K_next, liqK = sv.K' * c, sv_next.K' * c_next, sv.liqK' * c   # sv_next: the policies of the firms that continue
    r["I: net investment"] = K_next - (1 - p.delta) * K
    r["Psi: fixed operating costs"] = p.pc * sum(c)
    r["Psi: entry costs"] = m_next * E.c_e
    r["Psi: capital lost in liquidation"] = (1 - p.phi_a) * (1 - p.delta) * liqK
    r["Psi: profits lost in liquidation"] = sv.liqPi' * c
    r["Psi: reorganization, fixed (zeta)"] = sv.reoZeta' * c
    r["Psi: reorganization, variable"] = (1 - p.phi_c_hh) * (sv.reoV' * c)
    r["Psi"] = sum(v for (k, v) in r if startswith(k, "Psi: "))
    r["C"] = r["Y"] - r["I: net investment"] - r["Psi"]
    r["firms"] = sum(c); r["K"] = K; r["L"] = sv.L' * c
    r["defaults"] = sv.def' * c; r["liquidations"] = sv.liq' * c; r["reorganizations"] = sv.reo' * c
    return r
end

# welfare(E): the steady state's resources, scaled so that labor demand (production timed with next period's
# productivity, as in the law of motion) is 10,000. "C" charges every potential entrant c_e, as the code's free entry
# does; "C, main5 entry" charges only the entrants that stay, as main5's Psi.
function welfare(E; ev = events(E), sv = statevecs(E, ev))
    c = cont_mass(E)
    sc = 10000 / (sv.L' * c)
    r = resources(E, sv, sc .* c, sc .* c, sc * E.m)
    stay = sc * E.m * sum(entrant_dist(E))
    r["Psi: entry costs, main5"] = stay * E.c_e
    r["C, main5 entry"] = r["C"] + r["Psi: entry costs"] - r["Psi: entry costs, main5"]
    r["C, reorganization costs as transfers"] = r["C"] + r["Psi: reorganization, fixed (zeta)"] + r["Psi: reorganization, variable"]
    r["scale"] = sc; r["potential entrants"] = sc * E.m; r["entrants that stay"] = stay
    r["exits"] = sc * (sv.exit' * vec(E.mu)[1:end-1]); r["wage"] = E.wage
    r["law of motion gap"] = maximum(abs.(kernel(ev).T * c .+ E.m .* entrant_dist(E) .- vec(E.mu)[1:end-1])) / sum(c)
    return r
end

# path(E, d0; T): the firm distribution from the arrivals d0 (mass by state, before exit) in period 0, under E's
# policies and wage. Each period potential entrants m_t enter so that labor demand next period is Ns; while m_t > 0 the
# wage that free entry pins down stays E's wage, so this is the perfect-foresight path (values and loan prices are
# those of E's steady state). Returns per period: the continuing mass c_t, m_t and the arrivals d_t.
function path(E, d0; T = 60, Ns = 10000.0, ev = events(E), sv = statevecs(E, ev), K = kernel(ev))
    cont, ent = continuing(E), entrant_dist(E)
    d = copy(d0)
    cs, ms, ds = Vector{Float64}[], Float64[], Vector{Float64}[]
    corners = 0
    for t in 0:T
        c_inc = d .* cont
        m = (Ns - sv.L' * c_inc) / (sv.L' * ent)
        m < 0 && (m = 0.0; corners += 1)        # incumbents alone demand more labor than there is: the wage would have to rise
        c = c_inc .+ m .* ent
        push!(cs, c); push!(ms, m); push!(ds, copy(d))
        d = K.T * c
    end
    return (c = cs, m = ms, d = ds, sv = sv, corners = corners)
end

# path_static(E, d0, m; T): an approximate path, used where path() fails (entry at the corner). The firms follow E's
# steady-state policies, m potential entrants enter every period (E's steady-state number), and each period the wage
# clears the labor market given the capital in place: labor demand at wage w is L(w_E) (w_E / w)^(1/(1-nu)), so
# w_t = w_E (L_t / Ns)^(1-nu) and output is w_t Ns / nu. Firms do not foresee these wage deviations.
function path_static(E, d0, m; T = 60, Ns = 10000.0, ev = events(E), sv = statevecs(E, ev), K = kernel(ev))
    cont, ent = continuing(E), entrant_dist(E)
    d = copy(d0)
    cs, wages = Vector{Float64}[], Float64[]
    for t in 0:T
        c = d .* cont .+ m .* ent
        push!(cs, c); push!(wages, E.wage * ((sv.L' * c) / Ns)^(1 - E.params.nu))   # the wage when the firms c produce (next period)
        d = K.T * c
    end
    return (c = cs, wage = wages, sv = sv)
end

# liq_prob: FirmOptim's liquidation probability in default, for vectors of reorganization values and capital, at cost zeta
function liq_prob(E, vreo, k, zeta)
    p = E.params
    wedge = hasproperty(p, :liq_wedge) ? p.liq_wedge : 0.0
    L = p.phi_a .* (1 - p.delta) .* k
    return [p.sigma_v == 0 ? Float64(l >= p.phi_c_hh * v * exp(wedge) - z) :
            p.phi_c_hh * v <= 0 ? 1.0 : l + z <= 0 ? 0.0 : cdf(Normal(), (log((l + z) / (p.phi_c_hh * v)) - wedge) / p.sigma_v)
            for (l, v, z) in zip(L, vreo, zeta isa Number ? Iterators.repeated(zeta) : zeta)]
end

# follow(E, K, sv, d0; T): a group of firms arriving with the distribution d0, followed for T years under E's policies:
# the share still operating at the start of each year, and the cumulative shares that defaulted and that were liquidated
function follow(E, K, sv, d0; T = 5)
    cont, n0, d = continuing(E), sum(d0), copy(d0)
    alive, cdef, cliq = [1.0], [0.0], [0.0]
    for t in 1:T
        c = d .* cont
        push!(cdef, cdef[end] + (sv.def' * c) / n0); push!(cliq, cliq[end] + (sv.liq' * c) / n0)
        d = K.T * c
        push!(alive, sum(d) / n0)
    end
    return (alive = alive, cdef = cdef, cliq = cliq)
end

# state_stats(E, wts): averages of the firms' choices over the states, with the masses wts (continuing firms).
# Debt to capital, the rate and CF reliance are over borrowers, with sumSS's filters
function state_stats(E, wts)
    S, n = E.SumPol, length(wts)
    k, b, q, tau, e = S[1:n, 3], S[1:n, 4], S[1:n, 9], S[1:n, 17], S[1:n, 2]
    bor = b .> 0
    d2c = b ./ max.(k, 1e-12)
    okd = bor .& (k .> 0) .& (q .>= 0.5) .& (d2c .<= E.params.d2c_cut)
    okq = bor .& (q .>= 0.5)
    tot = sum(wts)
    return (mass = tot, k = sum(wts .* k) / tot, lab = sum(wts .* S[1:n, 10]) / tot, e = sum(wts .* e) / tot,
            borrowing = sum(wts[bor]) / tot, tau = wmean(tau, wts, bor), d2c = wmean(d2c, wts, okd),
            rate = wmean(100 .* (1 ./ max.(q, 1e-12) .- 1), wts, okq), debt = sum(wts .* b) / tot)
end

# outcomes(E): the headline numbers of an equilibrium, for the sweeps. "liq share" is the share of borrowers' defaults
# that end in liquidation (a flow); "liq prob" is sumSS's firm average of the liquidation probability.
function outcomes(E; ev = events(E))
    sv = statevecs(E, ev)
    w = welfare(E; ev = ev, sv = sv)
    rw = cont_mass(E)[ev.s] .* ev.prob
    bor = ev.b .> 0
    d = rw .* ev.pdef
    return Dict("wage" => E.wage, "C" => w["C"], "C main5" => w["C, main5 entry"], "firms" => w["firms"], "K" => w["K"],
                "default rate" => 100 * w["defaults"] / w["firms"], "liq per 100" => 100 * w["liquidations"] / w["firms"],
                "reo per 100" => 100 * w["reorganizations"] / w["firms"], "liq share" => sum(d[bor] .* ev.liq[bor]) / sum(d[bor]),
                "liq prob" => E.ss["avg_liqprob", "values"], "CF reliance" => E.ss["CF reliance", "values"],
                "CF share" => E.ss["CF share", "values"], "rate" => E.ss["avg intrate", "values"], "d2c" => E.ss["avg_b2a", "values"],
                "liq small" => E.sme["liqprob", "SME Values"], "liq large" => E.sme["liqprob", "LE Values"],
                "rate small" => E.sme["intrate", "SME Values"], "rate large" => E.sme["intrate", "LE Values"],
                "CFrel small" => E.sme["CF reliance", "SME Values"], "CFrel large" => E.sme["CF reliance", "LE Values"],
                "liqK" => w["scale"] * (sv.liqK' * cont_mass(E)), "reorganizations" => w["reorganizations"],
                "reo zeta" => w["Psi: reorganization, fixed (zeta)"])
end
# the standard rows of a sweep table: changes from the baseline outcomes ob
sweep_row(o, ob) = [pct(o["wage"], ob["wage"]), pct(o["C"], ob["C"]), pct(o["C main5"], ob["C main5"]), pct(o["firms"], ob["firms"]),
                    100 * o["liq share"], o["CF reliance"], o["CF share"], o["rate small"], o["default rate"]]
const SWEEP_COLS = ["Y/L %", "C %", "C main5 %", "firms %", "liq share", "CF reliance", "CF share", "rate small", "default rate"]

# weighted mean of v with weights w over the entries where keep is true
wmean(v, w, keep = trues(length(v))) = sum(w[keep] .* v[keep]) / sum(w[keep])

# tables: a vector of (label => values) rows and column names, printed and saved as CSV
function save_table(dir, file, cols, rows; digits = 3)
    df = DataFrame(; item = [string(k) for (k, v) in rows])
    for (j, cname) in enumerate(cols)
        df[!, cname] = [Float64(v[j]) for (k, v) in rows]
    end
    CSV.write(joinpath(dir, file), df)
    println("\n", file)
    println(rpad("", 46), join(lpad.(cols, 13)))
    for (k, v) in rows
        println(rpad(k, 46), join([lpad(round(x, digits = digits), 13) for x in v]))
    end
    return df
end
pct(a, b) = 100 * (a / b - 1)
small_plot(args...; kw...) = plot(args...; size = (560, 360), dpi = 110, legend = :best, titlefontsize = 10, guidefontsize = 9,
                                  tickfontsize = 8, legendfontsize = 8, left_margin = 3Measures.mm, bottom_margin = 3Measures.mm, kw...)
