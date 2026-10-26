# CALIBRATION TARGETS: DATA MOMENTS FOR SMALL AND LARGE FIRMS AND THEIR WEIGHTED AVERAGE - 2010-2024
# Small firms are underrepresented in Compustat-CapIQ. The moments are calculated separately for small and large firms
# (Compustat-CapIQ; liquidation probability from the IDB) and mapped onto a firm size distribution:
#     target = w * small + (1 - w) * large,    w = share of small firms among non-micro firms
# w comes from the US population of corporations (SOI, firm_weights.R) or from the IDB. CF-based to total debt is not
# weighted: it is the aggregate CF-based debt over all debt in the sample.
#
# INPUTS (all in CFL-git/Data Replication/):
#   dlevel_cleaned.RDS - contract level data before the firm level filters (data_cleaning.R, section 3)
#   comps_fjc.RDS      - cleaned Compustat (data_cleaning.R, section 5)
#   comps.RDS          - raw Compustat in USD, only for total liabilities ltq (data_cleaning.R, section 1)
#   idb_cases.RDS      - bankruptcy cases, one row per case (DefaultData.R, section 2)
#   soi_weights.RDS    - number of US corporations by asset class and year (firm_weights.R)
# OUTPUTS: none - the tables are printed in the console
#
# Sections 1 and 2 repeat sections 4 and 6 of data_cleaning.R, so that capital leases and the interest rate
# weighing can be switched. The "old" rate scheme follows R/V4/dlevel_v4.R.
#
# Structure:
#   0: Settings
#   1: Flevel, collapsing contract level data   (data_cleaning.R, section 4)
#   2: Flevel Full, connecting Compustat          (data_cleaning.R, section 6)
#   3: IDB, small and large firms
#   4: Compustat-CapIQ, small and large firms
#   5: Results: shares of small firms, calibration targets

library(tidyverse)
library(tictoc)

rm(list = ls())

setwd("C:/Users/szjud/OneDrive/Asztali gép/EBCs/CFL-git/Data Replication")
source("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V1//s_tab.R")


# 0: Settings ---------------------------------------------------------------------------------------------------

end_year         <- 2024     # sample: 2010Q1 - end_year Q4 (Compustat-CapIQ quarters and IDB filing years)

# Firm - debt level
# Defaults = the settings closest to the paper's targeted and untargeted moments (Tables tab:targets and tab:untargeted)
cap_lease_stock  <- FALSE    # capital leases in the debt stocks (value, CF/AB volumes, CFshare, d2c and the ±30% match filter)
cap_lease_rates  <- FALSE    # capital leases in the firm level interest rates (data_cleaning.R: excluded)

# Firm level aggregation of interest rates
rate_scheme      <- "new"    # "new": data_cleaning.R - weighted.mean(), contracts without a rate are dropped
                             # "old": R/V4/dlevel_v4.R - rate averaged over the contract's life, then sum(value*rate)/sum(value):
                             #        contracts without a rate stay in the denominator, which gives lower rates

# Small and large firms
size_var         <- "assets" # firm size: "assets" - Compustat total assets (atq)    <-> IDB estimated assets (EASST)
                             #            "ltq"    - Compustat total liabilities (ltq) <-> IDB estimated liabilities (ELBLTS)
                             #            "value"  - CapIQ total debt                  <-> IDB estimated liabilities (ELBLTS)
cutoff_usd       <- 1e8      # small = IDB size bin up to the cutoff (USD); the cutoff has to be an upper bound of an IDB bin:
                             #   5e4  ($0 to $50k)       1e5  ($50k to $100k)    5e5  ($100k to $500k)   1e6  ($500k to $1m)
                             #   1e7  ($1 to $10m)       5e7  ($10 to $50m)      1e8  ($50 to $100m)     5e8  ($100 to $500m)
                             #   1e9  ($500 to $1000m)   5e9  ($1b + (pre 2015)) 1e10 ($1-$10b)          5e10 ($10- $50b)

micro_usd        <- 1e6      # micro firms (size up to micro_usd) are dropped in Compustat-CapIQ, IDB and SOI; 0 = keep them
                             # SOI: 0, 5e5 or 1e6 (smallest SOI class: $1 under $500k); IDB: 0 or an IDB bin bound

# Weights
weight_source    <- "SOI"    # "SOI": US corporations (firm_weights.R) - needs size_var = "assets" and a cutoff that is also an
                             #        SOI class bound: 1e6, 1e7, 5e7, 1e8, 5e8
                             # "IDB": share of small cases in the IDB

# IDB sample
nonlead          <- "exclude"   # "include" / "exclude": cases that are members of a joint or consolidated filing
                                # (JNTLEAD or CNSLLEAD points to the lead case; it is NA for the lead case itself)

idb_bins <- c(5e4, 1e5, 5e5, 1e6, 1e7, 5e7, 1e8, 5e8, 1e9, 5e9, 1e10, 5e10)
stopifnot(cutoff_usd %in% idb_bins, micro_usd %in% c(0, idb_bins), micro_usd < cutoff_usd,
          size_var %in% c("assets", "ltq", "value"),
          rate_scheme %in% c("new", "old"),
          weight_source %in% c("SOI", "IDB"),
          nonlead %in% c("include", "exclude"))
if (weight_source == "SOI") stopifnot(size_var == "assets", cutoff_usd %in% c(1e6, 1e7, 5e7, 1e8, 5e8), micro_usd %in% c(0, 5e5, 1e6))

tic("calib_targets.R")


# 1: Flevel, collapsing contract level data ------------------------
# data_cleaning.R, section 4

dlevel <- readRDS("dlevel_cleaned.RDS") %>% filter(year >= 2010, year <= end_year)

# A: debt stocks by type

# capital leases in the debt stocks or not
dlevel_stock <- dlevel
if (!cap_lease_stock) dlevel_stock <- dlevel_stock %>% filter(debttext != "Capital leases")

flevel <- dlevel_stock %>%
  group_by(companyid, companyname, gvkey, yearfounded, country, CF, year, quarter) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = CF, values_from = value, values_fill = 0) %>%
  rename(AB_val = `0`, CF_val = `1`) %>%
  mutate(CFshare = if_else((CF_val + AB_val) > 0, CF_val / (CF_val + AB_val), NA_real_),
         dateper = (year - 2010) * 4 + quarter)

# B: weighted rates, spreads and maturities

# capital leases in the rates or not - data_cleaning.R drops them: accounting changes at the end of 2018 distort the data
dlevel_misc <- dlevel
if (!cap_lease_rates) dlevel_misc <- dlevel_misc %>% filter(debttext != "Capital leases")

if (rate_scheme == "new") {
  # data_cleaning.R: contracts without a rate in the given quarter are dropped from the weighted mean
  flevel_intrate <- dlevel_misc %>%
    mutate(w = coalesce(value, 0)) %>%
    mutate( dateper = (year - 2010) * 4 + quarter) %>%
    group_by(companyid, gvkey, dateper) %>%
    summarise(maturity = weighted.mean(matumax, w, na.rm = TRUE),
              intrate = weighted.mean(intrate, w, na.rm = TRUE),
              spread = weighted.mean(spread, w, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(across(where(is.numeric), ~ replace(.x, is.nan(.x), NA)))
}

if (rate_scheme == "old") {
  # R/V4/dlevel_v4.R: the rate of a contract is its average over the contract's life, spreads filled only where missing
  flevel_intrate <- dlevel_misc %>% group_by(id) %>%
    mutate(intrate = mean(intrate, na.rm = T),
           spread = ifelse(is.na(spread),  mean(spread, na.rm = T) , spread)) %>% ungroup() %>%
    mutate( dateper = (year - 2010) * 4 + quarter) %>%
    group_by(companyid, gvkey, dateper) %>%
    summarise(
      maturity = sum(value * matumax, na.rm = TRUE) / sum(value, na.rm = TRUE),
      intrate = sum(value * intrate, na.rm = TRUE) / sum(value, na.rm = TRUE),
      spread = sum(value * spread, na.rm = TRUE) / sum(value, na.rm = TRUE),
      .groups = "drop")
}

# Join: one row per firm-quarter
flevel <- flevel %>%
  left_join(flevel_intrate, by = c("companyid", "gvkey", "dateper"))

rm(flevel_intrate, dlevel_misc, dlevel_stock, dlevel)


# 2: Flevel Full, Connecting Compustat to Flevel -------------------------------------------------------
# data_cleaning.R, section 6

comps <- readRDS("comps_fjc.RDS")

# total liabilities (ltq) is not kept in comps_fjc.RDS - taken from the raw Compustat data (already in USD)
comps_ltq <- readRDS("comps.RDS") %>%
  group_by(gvkey, year, quarter) %>% filter(row_number() == 1) %>% ungroup() %>%   # gvkey duplicates, as in data_cleaning.R
  select(c("gvkey", "year", "quarter", "ltq"))

# joining flevel and comps
flevel <- flevel %>% left_join(comps, by = c("gvkey", "year", "quarter")) %>%
  left_join(comps_ltq, by = c("gvkey", "year", "quarter"))
rm(comps, comps_ltq)

# !!! dropping all observations that are not in CAPIQ !!!
flevel <-flevel %>% filter(!is.na(datadate))

# New variables that are calculable only with CAPIQ
flevel <- flevel %>%   rename(compsdebt = tot_debt) %>% mutate(
  emp = if_else(emp == 0, NaN , emp),
  sec_share = sec_debt/compsdebt,
  age = year - yearfounded + 1,
  value = (AB_val+CF_val)/1000000,
  d2a = value / assets,
  consist_flevel = value/compsdebt)

flevel <- flevel %>%  select(-c("compsdebt", "yearfounded", "datadate", "sic"))


# Truncating outliers - calculated ratios can take extreme values
truncate_p <- function(x, p) {

  lower_limit <- quantile(x, probs = p, na.rm = T)
  upper_limit <- quantile(x, probs = 1-p , na.rm = T)

  x[x < lower_limit] <- NA
  x[x > upper_limit] <- NA

  return(x)
}


# Only accepting observations where compustat and capiq reports similar data
sum(0.7 < flevel$consist_flevel & 1.3 > flevel$consist_flevel, na.rm = T) / sum(!is.na(flevel$consist_flevel))
# filter if the match quality - compustat and capiq - is poor
flevel <- flevel %>% filter(0.7 < consist_flevel & 1.3 > consist_flevel,
                            companyname!= "Berkshire Hathaway Inc.",
                            sic_group != ("Nonclassifiable"))

# Filtering unlikely observations
flevel <- flevel %>%
  mutate(dateper = (year - 2010)*4 + quarter,
         levrg = ifelse(levrg < 0 | levrg > 1, NA, levrg),
         pledge = ifelse(pledge < 0 | pledge > 1, NA, pledge),
         intcov = ifelse(intcov == -Inf | intcov > 500, NA, intcov),
         liq = ifelse(liq < 0, NA, liq),
         inv_rate = ifelse(inv_rate > 20 | inv_rate < - 20, NA, inv_rate),
         inv = ifelse(inv_rate > 20 | inv_rate < -20, NA, inv),
         age = ifelse(age < 0, NA, age))

# Further variables for the multivariate analysis
flevel <- flevel %>% group_by(companyid) %>% arrange(companyid, dateper) %>%
  mutate(ebit2a = ebitda / assets ,
         d2c = value / collat,
         d2c = ifelse(d2c < 0 | d2c > 1.5, NA, d2c),
         AB2a =  AB_val / 1000000 / assets,
         CF2a = CF_val / 1000000 / assets)


flevel <- flevel %>% ungroup() %>%
  mutate(d2a = truncate_p(d2a, 0.02),
         prod = truncate_p(prod, 0.02),
         AB2a = truncate_p(AB2a, 0.02),
         CF2a = truncate_p(CF2a, 0.02))


# 3: IDB, small and large firms ---------------------------------------------------------------------------------
# idb_cases.RDS: one row per case (DefaultData.R, section 2); chapter = final status, unconfirmed Ch.11 cases count as Ch.7

idb <- readRDS("idb_cases.RDS") %>% filter(year(ORGFLDT) >= 2010, year(ORGFLDT) <= end_year)

# size bin of the case: estimated assets or estimated liabilities (bin upper bound, USD)
idb <- idb %>% mutate(size = if (size_var == "assets") EASSTnum else ELBLTSnum) %>% filter(!is.na(size))

# micro firms
idb <- idb %>% filter(size > micro_usd)

# members of joint / consolidated filings
s_tab(nonlead_case, idb %>% mutate(nonlead_case = !is.na(JNTLEAD) | !is.na(CNSLLEAD)))
if (nonlead == "exclude") idb <- idb %>% filter(is.na(JNTLEAD) & is.na(CNSLLEAD))

idb <- idb %>% mutate(small = ifelse(size <= cutoff_usd, 1, 0))
s_tab(small, idb)


# 4: Compustat-CapIQ, small and large firms ----------------------------------------------------------------------
# assets, ltq and value are in million USD; small = up to the cutoff, as in the IDB bins

flevel <- flevel %>% mutate(size = .data[[size_var]]) %>%
  filter(!is.na(size), size * 1e6 > micro_usd) %>%
  mutate(small = ifelse(size * 1e6 <= cutoff_usd, 1, 0))
s_tab(small, flevel)


# 5: Results ------------------------------------------------------------------------------------------------------

cat("\nSettings: years 2010 -", end_year, "| capital leases in stocks:", cap_lease_stock, "| in rates:", cap_lease_rates,
    "| rate scheme:", rate_scheme, "| size:", size_var, "| cutoff:", format(cutoff_usd, big.mark = ",", scientific = FALSE),
    "USD | micro up to:", format(micro_usd, big.mark = ",", scientific = FALSE), "USD | non-lead cases (IDB):", nonlead,
    "| weights:", weight_source, "\n\n")

# SOI: non-financial corporations above the micro cut, all tax years pooled (firm_weights.R)
soi <- readRDS("soi_weights.RDS") %>% filter(lo >= micro_usd) %>%
  mutate(n = n_nf, small = ifelse(hi <= cutoff_usd, 1, 0))

# Table 1: share of small and large firms in Compustat-CapIQ, IDB and SOI
# Compustat-CapIQ firms are classified by their average size over the sample
firms <- flevel %>% group_by(companyid) %>% summarise(size = mean(size, na.rm = T)) %>%
  mutate(small = ifelse(size * 1e6 <= cutoff_usd, 1, 0))

shares <- tibble(
  data   = c("Compustat-CapIQ (firm-quarters)", "Compustat-CapIQ (firms)", "IDB (cases)", "SOI non-financial (returns)"),
  N      = c(nrow(flevel), nrow(firms), nrow(idb), sum(soi$n)),
  small  = c(mean(flevel$small), mean(firms$small), mean(idb$small), sum(soi$n * soi$small) / sum(soi$n)),
  large  = 1 - small)
print(shares)

# weight of small firms
w <- if (weight_source == "SOI") last(shares$small) else mean(idb$small)

# Table 2: calibration targets for small and large firms and their weighted average
moments <- function(f, d) {
  tibble(
    moment = c("Debt to Collateral", "Interest rate", "Liquidation prob.", "CF-reliance", "CF-based to total debt"),
    value  = c(mean(f$d2c, na.rm = T),
               mean(f$intrate, na.rm = T),
               mean(d$chapter == 7),
               mean(f$CFshare, na.rm = T),
               sum(f$CF_val, na.rm = T) / sum(f$CF_val + f$AB_val, na.rm = T)))
}

targets <- moments(filter(flevel, small == 1), filter(idb, small == 1)) %>% rename(small = value) %>%
  left_join(moments(filter(flevel, small == 0), filter(idb, small == 0)) %>% rename(large = value), by = "moment") %>%
  mutate(weighted = w * small + (1 - w) * large) %>%
  left_join(moments(flevel, idb) %>% rename(unweighted = value), by = "moment") %>%   # Compustat-CapIQ (IDB) sample averages
  mutate(weighted = ifelse(moment == "CF-based to total debt", unweighted, weighted))   # aggregate ratio, not weighted

cat("\n", weight_source, "weight of small firms: w =", round(w, 3),
    " \n Comps-CapIQ share of small firms:", round(mean(flevel$small), 3), "(firm-quarters)\n")
print(targets %>% mutate(across(where(is.numeric), ~ round(.x, 3))))

toc()
