# Data Replication

This folder holds the R code for the paper's empirical facts and for the calibration targets of the model.

**Rule:** all data-replication work goes in this folder. That includes new scripts, intermediate `.RDS` files, tables and plots. Model work goes in `Julia codes/Model Replication/`.

Existing scripts break this rule in a few places, and those are left as they are:
- Raw inputs are read from outside the repo and are read-only:
  - `EBCs/Data/R/V4/raw_dat/` holds the WRDS downloads.
  - `EBCs/Data/R/IDB dataset/` holds the FJC IDB snapshots.
  - `EBCs/Data/R/V4/exrate.csv` and `connected_firms.csv`.
- `DefaultData.R` and `results.R` save some paper figures straight to `Latex codes/Plots/`.

## Scripts

Run them in this order. `results.R` is an analysis script and comes last.

| Script | Does | Main outputs |
|---|---|---|
| `data_fetch.R` | Downloads Compustat and Capital IQ, 2010–2024, from WRDS | Raw CSVs in `EBCs/Data/R/V4/raw_dat/` |
| `data_cleaning.R` | Cleans and links Compustat and Capital IQ (sections 1–7) | `comps.RDS`, `capiq.RDS`, `dlevel_cleaned.RDS`, `flevel.RDS`, `comps_fjc.RDS`, `flevel_full.RDS`, `dlevel_full.RDS` |
| `DefaultData.R` | Builds business bankruptcy cases from the FJC IDB | `idb_cases.RDS`; `liqprob2.png` and `Liqprob_OT.png` (saved to `Latex codes/Plots/`) |
| `firm_weights.R` | Counts US corporations by asset size from IRS SOI, downloading files into `soi_raw/` if they are missing | `soi_weights.RDS` |
| `calib_targets.R` | Computes the small- and large-firm moments and the weighted calibration targets. The weights come from the SOI, the IDB or Compustat-CapIQ (`weight_source`); with Compustat-CapIQ only the liquidation probability is reweighted, by post-stratifying IDB cases to the Compustat-CapIQ size distribution | None saved; the tables are printed to the console. These are the targets for step 3 in `Julia codes/Model Replication/` |
| `results.R` | Summary statistics, credit-spread regressions and firm-level results | `spreads.png`, `spreads_age.png`, `mixplot.png`; `CFshare_firms.png` (saved to `Latex codes/Plots/`) |

## Requirements

- **R packages:** tidyverse, tictoc, lubridate, RPostgres, RPostgreSQL, readxl, fixest, estimatr, modelsummary.
- **Data access:** a WRDS account for Compustat and Capital IQ, used by `data_fetch.R`. The FJC IDB and IRS SOI data are public.
- **Not in git:** `.RDS`, `.csv`, `.xlsx` and `.zip` files are gitignored. Regenerate them by running the scripts in the order above.
