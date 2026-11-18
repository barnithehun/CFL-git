# screen_noise_rank.py <label> [top]: ranks the points of screens/noise_<label>_*_w*.csv (screen_noise.jl) for every tag in
# the files. Per tag: the calibration objective (five aggregates + the small minus large gaps in liquidation and CF reliance,
# each |model - data| / |data|), whether small and large CF reliance are within 0.10 of the data, and liquidation within
# 0.12 / 0.10. Prints the best points per tag and writes screens/rank_<label>.csv.
import csv, glob, math, os, sys

here = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.normpath(os.path.join(here, "..", "..", "..", "Data Replication"))
label = sys.argv[1]
top = int(sys.argv[2]) if len(sys.argv) > 2 else 8

rows = []
for f in sorted(glob.glob(os.path.join(here, "screens", f"noise_{label}_[0-9]*_w*.csv"))):
    rows += list(csv.DictReader(open(f, encoding="utf-8")))
if not rows:
    sys.exit("no rows")
tags = sorted({k[len("s_liq_"):] for k in rows[0] if k.startswith("s_liq_")})

def targets(tag):
    r = {x["moment"]: x for x in csv.DictReader(open(os.path.join(data_dir, f"calib_targets_{tag}.csv"), encoding="utf-8"))}
    g = lambda m, c: float(r[m][c])
    return dict(d2c=g("Debt to Collateral", "weighted"), rate=g("Interest rate", "weighted"), liq=g("Liquidation prob.", "weighted"),
                CFrel=g("CF-reliance", "weighted"), CFsh=g("CF-based to total debt", "weighted"),
                s_liq=g("Liquidation prob.", "small"), l_liq=g("Liquidation prob.", "large"),
                s_CFrel=g("CF-reliance", "small"), l_CFrel=g("CF-reliance", "large"),
                s_rate=g("Interest rate", "small"), l_rate=g("Interest rate", "large"))

F = lambda x: float(x) if x not in ("", "NaN", "nan") else math.nan
out = []
for tag in tags:
    d = targets(tag)
    for r in rows:
        m = {k: F(r[k]) for k in ("d2c", "rate", "liq", "CFrel", "CFsh")}
        s = {k: F(r[f"{k}_{tag}"]) for k in ("s_liq", "l_liq", "s_CFrel", "l_CFrel", "s_rate", "l_rate")}
        agg = sum(abs(m[k] - d[k]) / abs(d[k]) for k in m)
        gl = abs((s["s_liq"] - s["l_liq"]) - (d["s_liq"] - d["l_liq"])) / abs(d["s_liq"] - d["l_liq"])
        gc = abs((s["s_CFrel"] - s["l_CFrel"]) - (d["s_CFrel"] - d["l_CFrel"])) / abs(d["s_CFrel"] - d["l_CFrel"])
        gr = abs((s["s_rate"] - s["l_rate"]) - (d["s_rate"] - d["l_rate"])) / abs(d["s_rate"] - d["l_rate"])
        cf_ok = abs(s["s_CFrel"] - d["s_CFrel"]) <= 0.10 and abs(s["l_CFrel"] - d["l_CFrel"]) <= 0.10
        liq_ok = abs(s["s_liq"] - d["s_liq"]) <= 0.12 and abs(s["l_liq"] - d["l_liq"]) <= 0.10
        if math.isnan(agg + gl + gc):
            continue
        o = dict(tag=tag, point=r["point"], obj7=agg + gl + gc, obj8=agg + gl + gc + gr, agg=agg, cf_ok=int(cf_ok), liq_ok=int(liq_ok),
                 cf_miss=max(abs(s["s_CFrel"] - d["s_CFrel"]), abs(s["l_CFrel"] - d["l_CFrel"])))
        o.update({k: r[k] for k in ("pc", "pdef_exo_l", "zeta_Rl", "phi_c_hh", "phi_c", "sigma_v") + (("sv_slope",) if "sv_slope" in r else ())})
        o.update(s); o.update(m)
        out.append(o)

with open(os.path.join(here, "screens", f"rank_{label}.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)

print(f"{len(rows)} points, tags {', '.join(tags)}")
for tag in tags:
    o = sorted([x for x in out if x["tag"] == tag], key=lambda x: x["obj7"])
    n_cf = sum(x["cf_ok"] for x in o); n_both = sum(x["cf_ok"] and x["liq_ok"] for x in o)
    best_cf = min(x["cf_miss"] for x in o)
    d = targets(tag)
    print(f"\n== {tag}: data liq {d['s_liq']:.2f}/{d['l_liq']:.2f}, CF rel {d['s_CFrel']:.2f}/{d['l_CFrel']:.2f};"
          f" points with CF s/l ok {n_cf}, CF and liq ok {n_both}; smallest CF miss {best_cf:.3f}")
    for x in o[:top]:
        sv = f" sl {float(x['sv_slope']):+.2f}" if "sv_slope" in x else ""
        print(f"  pt {x['point']:>4} obj7 {x['obj7']:.3f} agg {x['agg']:.3f} liq {x['s_liq']:.2f}/{x['l_liq']:.2f} CF {x['s_CFrel']:.2f}/{x['l_CFrel']:.2f}"
              f" ok {x['cf_ok']}{x['liq_ok']} | pc {float(x['pc']):.0f} Px {float(x['pdef_exo_l']):.4f} z {float(x['zeta_Rl']):.0f}"
              f" hh {float(x['phi_c_hh']):.2f} th {float(x['phi_c']):.2f} sv {float(x['sigma_v']):.2f}{sv}")
    ok = sorted([x for x in o if x["cf_ok"]], key=lambda x: x["obj7"])[:4]
    for x in ok:
        print(f"  CF-ok pt {x['point']:>4} obj7 {x['obj7']:.3f} agg {x['agg']:.3f} liq {x['s_liq']:.2f}/{x['l_liq']:.2f} CF {x['s_CFrel']:.2f}/{x['l_CFrel']:.2f}"
              f" | hh {float(x['phi_c_hh']):.2f} th {float(x['phi_c']):.2f} sv {float(x['sigma_v']):.2f}")
