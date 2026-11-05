####### GRID EDGES ########
# check_grid: warns when an equilibrium runs up against the edge of a grid. E is an equilibrium from solve_eq, or any
# NamedTuple with SumPol, mu and gridsize. solve_eq and model_fit call it after the stationary distribution.
#  - e: the lowest point stands for every productivity below it (tauchen_lo in parameters()), so no firm should
#       continue there; and firms should not pile up at the highest point (twice as many firms there as one point
#       below: the top point holds the whole upper tail, so it is often about as large as its neighbour)
#  - k', b': no state should choose the top of the k or b grid (top, 10^5 in FirmOptim), whatever its mass
#  - x: no firms at the lowest or highest x point
# Shares of firms below tol are ignored. Returns the numbers behind the checks.
function check_grid(E; tol = 1e-4, top = 10.0^5)
    S = E.SumPol
    x_size, e_size, _, _ = E.gridsize
    n = size(S, 1) - 1                                  # the last state is default
    w = E.mu[1:n] ./ sum(E.mu[1:n])                     # share of firms in each (x, e) state; x runs fastest
    stays = (S[1:n, 6] .== 0) .& (S[1:n, 7] .== 0)       # neither exits nor defaults
    e_mass = vec(sum(reshape(w, x_size, e_size), dims = 1))
    x_mass = vec(sum(reshape(w, x_size, e_size), dims = 2))

    c = (e_low = sum(w[1:x_size] .* stays[1:x_size]), e_top = e_mass[end], e_below_top = e_mass[end-1],
         k_top = count(stays .& (S[1:n, 3] .>= 0.999 * top)), b_top = count(stays .& (S[1:n, 4] .>= 0.999 * top)),
         x_low = x_mass[1], x_high = x_mass[end])

    pct(v) = round(100 * v, sigdigits = 2)
    c.e_low > tol && @warn "Grid edge: $(pct(c.e_low))% of firms continue at the lowest productivity point; lower tauchen_lo"
    c.e_top > max(2 * c.e_below_top, tol) &&
        @warn "Grid edge: firms pile up at the highest productivity point ($(pct(c.e_top))% of firms, $(pct(c.e_below_top))% one point below); raise tauchen_sd"
    c.k_top > 0 && @warn "Grid edge: $(c.k_top) states choose the top of the k grid ($top)"
    c.b_top > 0 && @warn "Grid edge: $(c.b_top) states choose the top of the b grid ($top)"
    max(c.x_low, c.x_high) > tol && @warn "Grid edge: firms at the ends of the x grid (lowest point $(pct(c.x_low))%, highest $(pct(c.x_high))%)"
    return c
end
