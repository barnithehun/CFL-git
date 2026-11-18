# summary.py [wave ...]: compact view of runs/*/DIAG.csv (all waves if none given); also writes tracker.csv
import csv, glob, os, sys

here = os.path.dirname(os.path.abspath(__file__))
rows = []
for f in sorted(glob.glob(os.path.join(here, "runs", "*", "DIAG.csv"))):
    with open(f, newline="", encoding="utf-8") as fh:
        r = list(csv.DictReader(fh))
        if r:
            rows.append(r[0])
cols = []
for r in rows:
    for c in r:
        if c not in cols:
            cols.append(c)
with open(os.path.join(here, "tracker.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=cols)
    w.writeheader()
    for r in rows:
        w.writerow(r)

waves = sys.argv[1:]
if waves:
    rows = [r for r in rows if r.get("wave") in waves]

def g(r, k, fmt="{:.2f}"):
    v = r.get(k, "")
    if v in ("", None):
        return "-"
    if v in ("true", "false"):
        return "Y" if v == "true" else "."
    try:
        return fmt.format(float(v))
    except ValueError:
        return v

view = [("id", "id", "{}"), ("st", "status", "{}"), ("err", "error", "{:.3f}"),
        ("liq s/l", None, None), ("rate s/l", None, None), ("d2c s/l", None, None), ("CFrel s/l", None, None),
        ("sz6", "err_size6", "{:.2f}"), ("YL", "r_YL", "{:.2f}"), ("mass", "r_mass", "{:.1f}"),
        ("rS b>r", None, None), ("rL b>r", None, None), ("M", "dc_M", "{:.2f}"), ("KY", "dc_KY", "{:.2f}"),
        ("SBt/u/c", None, None), ("tau", "tau_interior", "{:.2f}"), ("score", "score", "{}")]
print("  ".join(h for h, _, _ in view) + "   label")
for r in rows:
    out = []
    for h, k, fmt in view:
        if k:
            out.append(g(r, k, fmt))
        elif h == "liq s/l":
            out.append(g(r, "s_liq") + "/" + g(r, "l_liq"))
        elif h == "rate s/l":
            out.append(g(r, "s_rate", "{:.1f}") + "/" + g(r, "l_rate", "{:.1f}"))
        elif h == "d2c s/l":
            out.append(g(r, "s_d2c") + "/" + g(r, "l_d2c"))
        elif h == "CFrel s/l":
            out.append(g(r, "s_CFrel") + "/" + g(r, "l_CFrel"))
        elif h == "rS b>r":
            out.append(g(r, "s_rate", "{:.1f}") + ">" + g(r, "r_s_rate", "{:.1f}"))
        elif h == "rL b>r":
            out.append(g(r, "l_rate", "{:.1f}") + ">" + g(r, "r_l_rate", "{:.1f}"))
        elif h == "SBt/u/c":
            out.append(g(r, "sb_t_YL") + "/" + g(r, "sb_u_YL") + "/" + g(r, "sb_c_YL"))
    print("  ".join(out) + "   " + r.get("label", ""))
# data rows of the target sets in use
seen = set()
for r in rows:
    t = r.get("targets")
    if t and t not in seen:
        seen.add(t)
        print(f"data {t}: liq {g(r,'ds_liq')}/{g(r,'dl_liq')}  rate {g(r,'ds_rate','{:.1f}')}/{g(r,'dl_rate','{:.1f}')}  d2c {g(r,'ds_d2c')}/{g(r,'dl_d2c')}  CFrel {g(r,'ds_CFrel')}/{g(r,'dl_CFrel')}  share {g(r,'share_data')}")
