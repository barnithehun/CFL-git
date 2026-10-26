# Predicted liquidation probability and spreads of CF- vs AB-based debt

> **Update (2026-10-14):** this report documents the full specification search, now archived in `liqprob_specsearch.R`.
> `liqprob.R` has been rewritten around the FJC liquidation rate, with outputs in `liqprob_out/fjc/`. It includes the decisive flexible-size test.
> **The small-firm FJC result does not survive asset-bucket fixed effects** (CF 0.07 / 0.20 / 0.49, n.s.) **or firm fixed effects.**
> It is identified from the asset-size profile and cannot be separated from a nonlinear size effect, so treat the recommendation in section 7 as superseded.

Script: `R/Replication/liqprob_specsearch.R` · Outputs: `R/Replication/liqprob_out/` · Run time ≈ 12.5 min (≈ 7 min without the bootstrap)

## 1. Bottom line

1. **In the full sample there is no robust result.** Across ~100 identified phat constructions, adding phat to the main regression gives a significantly
   positive CF coefficient in about half of the specifications. In most of those, the AB coefficient is just as large or larger.
   - The original V3 phat (`phat_stata.dta`) and its re-estimated version both work the wrong way: AB spreads load on phat 3x more than CF spreads.
   - The full-sample specifications that do give "CF > AB" usually break a main result. Typically pledgeability becomes significant for CF debt.
2. **In the smaller firms the first-best pattern appears**: phat raises CF spreads and has no effect on AB spreads, and the main regression stays intact.
   The cleanest version uses an external liquidation rate. This is the Ch.7 share of business bankruptcies in the FJC bankruptcy data
   (about 200k cases), by asset × liability bucket. It was fixed in advance, before I looked at results.
   - It holds in all three samples, in the bottom half and the bottom tercile of firms, and in the within-firm design.
   - The CF − AB difference is significant at 5% in every one of those cells.
   - It remains after controlling for predicted default probability (CF − AB p = 0.047 / 0.004 / 0.001; section 4).
3. **Phats estimated from Compustat defaults are fragile.** The probability models are fit on only 300–1,000 default events, and the
   linear model extrapolates far outside the [0, 1] range for the larger firms in the sample. Once the first-stage noise is bootstrapped, only the
   pre-default (4-quarter-lagged) probit keeps a borderline significant small-firm CF − AB gap (one-sided p ≈ 0.03–0.07).
4. **Industry level does not rescue the result.** Across SIC2 industries the sign is *wrong*: industries with a higher Ch.7 share have *lower* CF spreads.
   Within an industry over time there is weak support: CF spreads rise with the industry's predicted liquidation probability in the pre-2021
   samples (p = 0.02–0.06), while AB spreads don't move.
5. **Default risk is a confound you need to address.** Predicted *default* probability itself raises CF spreads and not AB spreads. So "CF spreads respond to
   distress risk, AB spreads don't" is very robust. "CF spreads respond to *liquidation-given-default* risk" is only robust for smaller firms.

**What I would put in the paper:** the FJC size-cell liquidation rate (`idb_liq_final_AL`) in the main spread regression, estimated
for smaller firms (or as a phat × small-firm interaction), with the within-firm version and the default-probability control as robustness checks. Section 7 has the details.

## 2. What was built

### Regression sample (unchanged from `results.R`)
- New issuances from `R/Replication/dlevel_full.csv`. Spec `f_spread`: lebitda, lassets, pledge, levrg, lage, lemp, Nconts, CFshare + sic_group, rating, quarter FE.
- The baseline reproduces `table_outcome.tex` exactly: N = 22,820 (CF 13,772 / AB 9,048).
- Three samples:
  - full 2010–2024: N = 22,820
  - year < 2021: N = 18,363 (the paper / do-file sample)
  - year < 2019: N = 13,850
- Inference:
  - firm-clustered SEs throughout
  - HC1 (the paper's `, robust`) also stored
  - industry-clustered SEs added for industry-level phats
- Note: with firm clustering, even the baseline CFshare coefficient for CF debt (−0.24) is not significant (p = 0.12).

### phat constructions (118 variants, all standardised: coefficients are percentage points of spread per 1 sd of phat)

| Family | What it is | # |
|---|---|---|
| Compustat default model | Pr(Ch.7 \| default) estimated on `V3/liquidations.dta` and predicted on the regression sample. Varies along 4 dimensions (listed below). | 92 |
| FJC IDB size-cell rate | Ch.7 share of business Ch.7+Ch.11 cases 2005–2024 by asset bucket, liability bucket, asset×liability, and each × filing year. Two outcomes: final chapter (not-confirmed Ch.11 counted as liquidation, as in `IDB_easst.R`) or chapter at filing. Cells are shrunk toward coarser rates (empirical Bayes). | 12 |
| Industry | SIC2 leave-own-firm-out Ch.7 share of Compustat defaults, and SIC2 (×year) means of firm-level phat over the Compustat universe | 5 |
| Firm mean | Time-invariant firm average of phat (the V3 `mphat` analogue) | 3 |
| Original | `phat1` and `mphat` from `V3/phat_stata.dta` and the V4 copy | 4 |
| Contrast / placebo | Unconditional Pr(Ch.7 within 8q); Pr(any default within 8q) | 2 |

The four dimensions of the Compustat default models:
- **Event sample:**
  - all flagged quarters (the V3 replication; N = 996, 18% Ch.7)
  - first flag per firm (N = 727)
  - covariates 4 quarters before the first flag (N = 678)
- **Outcome:** the V3 flags, or a strict version where Ch.11 = footnote TL or deletion code 02.
- **Covariates:**
  - S1: size
  - S2: size + EBITDA + pledge (= V3)
  - S3: S2 + leverage + employment
  - S4: *only* variables excluded from the spread regression (market-to-book, PP&E share of collateral, cash)
  - S5: S2 + S4
  - S6: S5 + SG&A
- **Link:** LPM (primary), probit, logit.

The re-estimated V3 probit reproduces the V3 coefficients: lassets −0.155, lebitda −0.090, pledge −0.330.

**Identification check.** For every spec the script computes the R² of phat on the controls and FE.
- An LPM built only from variables already in the spread regression (S1–S3) is an exact linear combination of the controls. Its coefficient is not identified, and those 18 variants are reported as "not identified" rather than estimated.
- In an earlier draft I clipped LPM predictions to [0, 1], which created a kink that looked identified. The script now uses unclipped values.
- Where phat is identified, the controls typically explain 85–96% of its variance. The coefficient comes from the remaining nonlinear or excluded variation.

### Designs (each run on every phat × sample)
- **D1:** main regression + phat. All coefficients and FE are interacted with CF, so the point estimates equal the CF/AB split of `results.R` and β_CF − β_AB gets a proper test.
- **D1pd:** D1 + predicted default probability (logit on Compustat, 303k firm-quarters).
- **D1x:** D1pd + phat × default probability, testing whether liquidation risk matters in proportion to default risk.
- **D3:** horse race. phat replaces its own inputs (lassets, lebitda, pledge).
- **D4:** within firm (firm FE).
- **D5:** smaller half by assets, larger half, bottom tercile.
- **D5pd:** smaller half + default probability control.
- **I1/I2:** industry cross-section and industry × year panel (section 5).
- **D6:** SBRA triple difference.

**"Main results intact".** Each spec is checked against the baseline of the same sample and design, using the paper's own HC1 inference. The four key coefficients are CFshare and pledge in the CF and AB columns. A coefficient counts as intact if it keeps both its sign and its significance status at 10%.

## 3. Robustness across all specifications

Share of identified specifications (firm-clustered SEs) in which:
- **CF+**: β_CF > 0 at 5%
- **first best**: CF+ and β_AB insignificant at 10%
- **CF>AB**: β_CF − β_AB > 0 at 5%
- **intact**: main results unchanged

| Design | Sample | CF+ | first best | CF>AB | intact | first best & intact |
|---|---|---|---|---|---|---|
| D1 all firms | 2010–24 | 50% | 16% | 11% | 33% | 3% |
| D1 all firms | <2021 | 51% | 9% | 13% | 34% | 1% |
| D1 all firms | <2019 | 54% | 10% | 7% | 33% | 2% |
| D5 smaller half | 2010–24 | 24% | 23% | 19% | 59% | 20% |
| D5 smaller half | <2021 | 21% | 8% | 11% | 61% | 7% |
| D5 smaller half | <2019 | 28% | 26% | 13% | 59% | 25% |
| D5 bottom tercile | 2010–24 | 16% | 15% | 13% | 78% | 14% |
| D5pd smaller half + PD control | 2010–24 | 21% | 14% | 13% | 11%* | 11% |
| D5pd smaller half + PD control | <2021 | 18% | 7% | 13% | 2%* | 2% |
| D4 within firm | 2010–24 | 12% | 3% | 2% | 41% | 1% |

\* Low because the PD control itself changes the small-firm pledge coefficient for CF debt (see section 4a).

How to read this:
- **All firms.** phat is "priced" in both CF and AB debt. In the median spec the AB coefficient (0.28) is slightly larger than the CF one (0.24).
- **Smaller firms.** The AB coefficient collapses to about 0: median 0.09 in the smaller half, −0.02 in the bottom tercile. Whether the CF coefficient is positive then depends on which variables drive phat (next point).
- **The main split is the covariates of the Compustat model.**
  - Models that include size, profitability and market-to-book (S5) or pre-default leverage (lagged S3) give CF ≫ AB for small firms.
  - Models driven by cash alone (S4) give a *negative* CF coefficient, because cash raises Ch.7 probability but also lowers CF spreads directly.
  - The V3 specification (S2) loads on AB debt.

See `liqprob_out/spec_curve.png`: in the right-hand panels the red AB dots sit on zero, while the blue CF line rises above zero for roughly the top 40–50% of specifications.

## 4. The constructions in detail

### (a) Fixed in advance (headline). Firm-clustered SE; *p* in brackets.

**FJC IDB rate, assets × liabilities** (`idb_liq_final_AL`; mean 0.14, sd 0.20). The first stage uses about 200k cases, so the generated-regressor problem is negligible.

| | CF | AB | CF − AB | Main results intact |
|---|---|---|---|---|
| All firms 2010–24 | 0.52 [<0.001] | 0.41 [<0.001] | 0.11 [0.47] | yes |
| All firms <2021 | 0.66 [<0.001] | 0.48 [<0.001] | 0.18 [0.27] | yes |
| **Smaller half 2010–24** | **0.50 [0.049]** | **−0.13 [0.44]** | **0.62 [0.03]** | yes |
| **Smaller half <2021** | **0.74 [0.005]** | **−0.17 [0.32]** | **0.91 [0.003]** | yes |
| **Smaller half <2019** | **0.98 [0.001]** | **−0.21 [0.34]** | **1.19 [0.001]** | yes |
| Bottom tercile 2010–24 / <2021 / <2019 | 0.44 [0.11] / 0.63 [0.03] / 0.83 [0.01] | −0.22 / −0.30 / −0.37 (n.s.) | diff p = 0.04 / 0.007 / 0.003 | yes |
| Within firm (firm FE) 2010–24 / <2021 / <2019 | 0.64 [0.006] / 0.84 [0.003] / 0.88 [0.01] | 0.06 / −0.03 / 0.01 (n.s.) | diff p = 0.04 / 0.01 / 0.04 | pledge-CF loses 10% significance (−0.51 → −0.40) |
| Smaller half + default probability control 2010–24 / <2021 / <2019 | 0.44 [0.09] / 0.69 [0.01] / 0.91 [0.002] | −0.16 / −0.22 / −0.28 (n.s.) | diff p = 0.047 / 0.004 / 0.001 | see note* |

\* The default probability control itself absorbs the small-firm pledge coefficient for CF debt (0.66*** → 0.10), whether or not phat is included. So the "intact" flag fails because of the control, not because of phat.

- **Magnitude.** In the smaller half (2010–24), a 10 percentage point higher liquidation rate raises CF spreads by about 25 bp; in the pre-2019 sample, about 50 bp.
- **Other IDB cells.**
  - The final-outcome cells (A, AL, AY, L, ALY) behave similarly.
  - The *filing*-chapter cells do not: in small firms their CF coefficient is negative or zero.
  - So the result relies on *eventual* liquidation, which includes Ch.11 cases that end up dismissed or converted. Economically that is the right object, but it is a choice you should state.

**Primary LPM** (`fs_v3_first_S5_lpm`: first default per firm, core + market-to-book/PP&E share/cash, 337 events).
- **All firms:** CF 0.94 [<0.001], AB 0.28 [0.13], diff p = 0.014. But pledge becomes significant for CF debt (0.90***), so the main result is *not* intact.
- **Smaller half:** CF 0.99 / 1.23 / 1.61, AB 0.19 / 0.34 / 0.21 (n.s.) for 2010–24 / <2021 / <2019. Intact.
- **The bootstrap removes it.** Resampling default events and firms (200 reps) gives a smaller-half CF interval of [−0.80, 2.08] and P(diff ≤ 0) = 0.10–0.22.
- **Extrapolation.** 45% of the LPM fitted values for the regression sample lie outside [0, 1]: the regression firms are much larger than the defaulted Compustat firms.
- I would not use it.

**V3 replication** (`fs_v3_all_S2_probit`, = `phat_stata.dta`).
- **All firms:** CF 0.24 [0.08], AB 0.79 [<0.001], diff −0.55 [0.016].
- **Smaller half:** CF ≈ 0, AB 0.4–0.55.
- Wrong direction everywhere, and the bootstrap confirms it (P(diff ≤ 0) = 0.92–0.98).
- The original `phat_stata.dta` files give the same picture (AB 0.8–1.2 vs CF 0.3–0.4).

### (b) Found by the specification search (treat as exploratory)

**Pre-default probit** (`fs_v3_lag4_S3_probit`: covariates 4 quarters before the default, incl. leverage and employment).
- **All firms:** CF 0.75–0.94***, AB 0.18–0.41 (n.s.), diff p = 0.03–0.17. Pledge for CF debt becomes significant, so not intact.
- **Smaller half:** CF 1.09–1.20***, AB −0.01–0.11 (n.s.), diff p = 0.003–0.015, intact.
- **Bootstrap:** smaller-half CF − AB 0.75–0.86, 95% CI about [0, 1.9], one-sided p = 0.03 / 0.07 / 0.03. This is the only Compustat-based phat that survives the bootstrap, and only marginally.

**IDB liabilities × year, final outcome** (`idb_liq_final_LY`).
- **All firms:** first best and intact in 2010–24 (CF 0.31 [0.001], AB 0.06 [0.46], diff p = 0.04) and in <2021 (diff p = 0.01). Not in <2019.

### (c) Liquidation risk × default risk (D1x)
- If liquidation risk matters in proportion to default risk, phat × PD should load on CF spreads.
- **IDB family:** the interaction is positive and significant for CF in 100% of variants and samples. For AB it is significant in at most 25%, and the CF − AB difference is significant in 50–83%.
- **Compustat-model family:** CF-positive in about half.
- This is the theory-consistent pattern, but the main effect of phat in these models still loads on AB debt for all firms.

### (d) Default probability placebo
- Predicted default probability (Pr(default within 8q)), entered in place of phat:
  - all firms: CF 0.16 [0.005] vs AB −0.21 [0.03]
  - smaller half: CF 0.92–1.07 vs AB −0.15 to −0.29
- "CF debt prices distress risk, AB debt doesn't" is robust.
- To credit the effect to liquidation risk specifically you need the D1pd/D5pd controls (included above). The phat × PD interaction (D1x) is the natural way to state it.

## 5. Industry level

**I1, SIC2 cross-section.** Spreads are residualised on rating and quarter FE plus the controls, and averaged by industry (≈ 59 industries).
- Industries with a higher Ch.7 share among defaults have **lower** CF spreads: −0.13 to −0.14 per sd, p ≤ 0.01. AB spreads show no relationship.
- The industry mean of phat gives the same sign, insignificant.
- Industry composition (e.g., tech and services have high Ch.7 shares and also many CF borrowers with high market-to-book) dominates any liquidation channel.

**I2, SIC2 × year panel with industry and year FE.**
- Using the industry-year mean of the S5 phat:
  - CF spreads rise: 0.54 [0.02] (<2019), 0.39 [0.06] (<2021), 0.21 [0.19] (full).
  - AB spreads don't: −0.07 to 0.01.
  - CF − AB gap p = 0.05 / 0.07 / 0.13.
- The S4-based measure gives a large, significant gap, but it is driven by AB spreads *falling*.
- This is weak, supportive evidence. It is about 700 industry-years with noisy industry phats.

## 6. Other findings

- **SBRA triple difference** (small = CapIQ debt ≤ $7.5m, post = 2020Q2+).
  - Spreads of small firms fall after 2020 for both loan types: −1.46 [<0.001] for AB, with a CF-specific extra −0.24 [0.62].
  - No differential CF effect. The IDB confirms that the Ch.7 share of small-business cases fell (83% → 75%), but in the spread data this does not show up specifically in CF debt.
- **Main regression within size groups.** In the smaller half of firms the paper's main results already look different:
  - pledge is *positive and significant* for CF debt (0.64)
  - CFshare is positive (insignificant) for CF debt

  The "intact" check for size groups therefore compares with the size-group baseline. If the paper moves towards small firms, that table changes too.
- **Not identified.** 18 LPM variants per design are exact linear combinations of the controls. In the within-firm design 27–36 more are nearly collinear (R² > 0.999). These are reported as "not identified" in `spec_grid.csv`.

## 7. Recommendation

1. **Main specification.** Use the FJC asset × liability liquidation rate. It is external, has a huge N, and has no generated-regressor problem. Include it either as `phat × small` interacted with CF, or by estimating the spread regression on the smaller half of firms.
   Report the CF and AB columns with the CF − AB test (p ≤ 0.03 in all three samples).
2. **Robustness table:**
   - bottom tercile
   - firm FE
   - + predicted default probability
   - the phat × PD interaction
   - the pre-default Compustat probit with bootstrapped SEs
3. **Be transparent about the scope.** The effect is a small-firm effect. In the full sample both debt types load on phat. This fits the paper's own story: liquidation is a realistic outcome only for small firms.
4. **Drop the V3 phat.** It loads on AB debt in every sample.

## 8. Caveats

- **Specification search.** About 2,700 regressions were run (5,830 rows in `spec_grid.csv` across SE types). Only the three a-priori headline phats and the pre-specified size split should be treated as confirmatory. The lag-4 probit and IDB-LY results are exploratory.
- **Functional-form identification.** IDB phat is a step function of assets and liabilities, and its CF effect is identified by the 4–14% of its variance not explained by linear size, leverage and the FE. The within-firm design helps, but not completely.
- **Small default sample.** The Compustat flags were reconstructed from deletion codes (Ch.7 = deletion reason 03). There are only 177–188 Ch.7 events, and they are public firms, much larger than the paper's small firms.
- **Clustering.** With firm-clustered SEs the baseline CFshare coefficient for CF debt is already insignificant (p = 0.12). The "intact" check uses HC1 to match the paper.

## 9. Files

| File | Content |
|---|---|
| `liqprob_specsearch.R` | The whole analysis (archived; `liqprob.R` is the FJC-focused version). `LIQPROB_OUT`, `LIQPROB_BOOT`, `LIQPROB_NBOOT`, `LIQPROB_HEADLINE` env variables redirect output / switch the bootstrap |
| `liqprob_out/spec_grid.csv` | One row per phat × design × sample × SE type: β_CF, β_AB, diff, SEs, the main coefficients, intact flags, identification R² |
| `liqprob_out/spec_summary_by_design.csv`, `spec_summary_by_family.csv` | The shares in section 3 |
| `liqprob_out/spec_curve.png` | Specification curves (all firms vs smaller half, three samples) |
| `liqprob_out/headline_{all,small}_{full,lt2021,lt2019}.{tex,md}` | Regression tables: baseline + 3 headline phats, CF and AB columns |
| `liqprob_out/first_stage.tex`, `first_stage_samples.csv` | First-stage models and default-event samples |
| `liqprob_out/phat_variants.csv` | All 118 phats: definition, first-stage N and R², coverage, share outside [0, 1] |
| `liqprob_out/bootstrap_summary.csv`, `bootstrap_draws.csv` | Two-sample cluster bootstrap (200 reps) |
| `liqprob_out/industry_results.csv`, `sbra_results.csv` | Sections 5 and 6 |
| `liqprob_out/run_log.txt` | Log of the final run |
