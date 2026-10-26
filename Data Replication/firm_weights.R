# FIRM SIZE WEIGHTS: NUMBER OF US CORPORATIONS BY SIZE OF TOTAL ASSETS - IRS STATISTICS OF INCOME (SOI), 2010-2022
# Population weights for the calibration targets (calib_targets.R): Compustat-CapIQ underrepresents small firms.
# SOI counts returns of active corporations (C and S corps, REITs, RICs); a consolidated group counts once.
#
# INPUTS:  soi_raw/ (downloaded in section 1 if missing), idb_cases.RDS, flevel_full.RDS, comps_fjc.RDS (comparison plots only)
# OUTPUT:  soi_weights.RDS - number and share of corporations by asset class and year, all industries and non-financial
#
# DOWNLOAD LINKS (each file is listed in soi_files, section 1):
#   Complete Report, Table 2.1 (2014-2022): https://www.irs.gov/statistics/soi-tax-stats-corporation-income-tax-returns-complete-report-publication-16
#   Source Book flat files (2004-2013):      https://www.irs.gov/statistics/soi-tax-stats-corporation-source-book-data-file
#   Explanation of terms:                    https://www.irs.gov/statistics/soi-tax-stats-section-5-explanation-of-terms
#
# Non-financial firms: industry detail exists only in the Source Book (2010-2013). For 2014-2022 the 2010-2013 industry mix
# within each asset class is carried forward: non-financial count = all industries count * (1 - excluded share of the class).
#
# Structure:
#   0: Settings
#   1: SOI, downloading and reading
#   2: SOI, non-financial firms and the weights table
#   3: Comparison: SOI vs IDB vs Compustat

library(tidyverse)
library(readxl)

rm(list = ls())

setwd("C:/Users/szjud/OneDrive/Asztali gép/EBCs/CFL-git/Data Replication")

my_theme <- theme(plot.title = element_text(face = "bold", hjust = 0.5, size =14),
                  panel.spacing = unit(2, "lines"), strip.background=element_rect(fill="white"),
                  axis.ticks.length=unit(0, "cm"),
                  axis.text.y = element_text(size = 12, angle = 0, vjust = 0.5, hjust=0.5),
                  axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
                  axis.title.x=element_text(size = 12), axis.title.y=element_text(size = 12),
                  panel.background = element_rect(fill = "white"),
                  panel.grid.major = element_line(linewidth = 0.2, linetype = 'solid', colour = "grey80"),
                  panel.grid.minor = element_line(linewidth = 0.15, linetype = 'solid', colour = "white"),
                  panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.5))


# 0: Settings ---------------------------------------------------------------------------------------------------

soi_path     <- "soi_raw/"
drop_sectors <- c(22, 52, 53, 55)  # NAICS: utilities, finance and insurance, real estate, holding companies
                                   # (Compustat sample drops SIC 4900-4999 and 6000-6799, data_cleaning.R)
years        <- 2010:2024          # comparison plots: IDB filing years, Compustat years
nonlead      <- "exclude"          # comparison plots: IDB member cases of joint / consolidated filings


# 1: SOI, downloading and reading --------------------------------------------------------------------------------

soi_files <- tribble(
  ~year, ~file,             ~url,
  2010,  "10sbfltfile.zip", "https://www.irs.gov/pub/irs-soi/10sbfltfile.zip",   # Source Book, industry x asset class
  2011,  "11sbfltfile.zip", "https://www.irs.gov/pub/irs-soi/11sbfltfile.zip",
  2012,  "12sbfltfile.zip", "https://www.irs.gov/pub/irs-soi/12sbfltfile.zip",
  2013,  "13sbfltfile.zip", "https://www.irs.gov/pub/irs-soi/13sbfltfile.zip",
  2014,  "14co21ccr.xlsx",  "https://www.irs.gov/pub/irs-soi/14co21ccr.xlsx",    # Table 2.1, all industries by asset class
  2015,  "15co21ccr.xlsx",  "https://www.irs.gov/pub/irs-soi/15co21ccr.xlsx",
  2016,  "16co21ccr.xlsx",  "https://www.irs.gov/pub/irs-soi/16co21ccr.xlsx",
  2017,  "17co21ccr.xlsx",  "https://www.irs.gov/pub/irs-soi/17co21ccr.xlsx",
  2018,  "18co21ccr.xlsx",  "https://www.irs.gov/pub/irs-soi/18co21ccr.xlsx",
  2019,  "19co21ccr.xlsx",  "https://www.irs.gov/pub/irs-soi/19co21ccr.xlsx",
  2020,  "20co21ccr.xlsx",  "https://www.irs.gov/pub/irs-soi/20co21ccr.xlsx",
  2021,  "21co21ccr.xlsx",  "https://www.irs.gov/pub/irs-soi/21co21ccr.xlsx",
  2022,  "22co21ccr.xlsx",  "https://www.irs.gov/pub/irs-soi/22co21ccr.xlsx")

dir.create(soi_path, showWarnings = FALSE)
for (i in seq_len(nrow(soi_files))) {
  f <- paste0(soi_path, soi_files$file[i])
  if (!file.exists(f)) download.file(soi_files$url[i], f, mode = "wb")
  sb_dir <- paste0(soi_path, "sb", soi_files$year[i] %% 100)
  if (grepl("zip$", f) & !dir.exists(sb_dir)) unzip(f, exdir = sb_dir)
}

# asset classes: AC = Source Book code, same order as the columns of Table 2.1; lo / hi in USD
soi_classes <- tibble(
  AC    = 3:13,
  class = c("$1 under $500k", "$500k under $1m", "$1 under $5m", "$5 under $10m", "$10 under $25m", "$25 under $50m",
            "$50 under $100m", "$100 under $250m", "$250 under $500m", "$500m under $2.5bn", "$2.5bn or more"),
  lo    = c(0, 5e5, 1e6, 5e6, 1e7, 2.5e7, 5e7, 1e8, 2.5e8, 5e8, 2.5e9),
  hi    = c(5e5, 1e6, 5e6, 1e7, 2.5e7, 5e7, 1e8, 2.5e8, 5e8, 2.5e9, Inf))

# Source Book 2010-2013: WT_CT = number of returns; INDY_CD 1 = all industries, 2 digits = sectors (41 = 42 + 44 + 46, dropped)
soi_sb <- map_dfr(2010:2013, function(y) {
  read_csv(sprintf("%ssb%d/%dsb1.csv", soi_path, y %% 100, y), show_col_types = FALSE) %>%
    transmute(year = y, sector = INDY_CD, AC, returns = WT_CT)
}) %>% filter(AC %in% soi_classes$AC, sector == 1 | (sector < 100 & sector != 41))   # zero assets (AC 2) and total (AC 1) out

# Table 2.1 2014-2022: row "Number of returns", columns 4-14 = AC 3-13 (col 2 total, col 3 zero assets)
soi_t21 <- map_dfr(2014:2022, function(y) {
  x <- read_excel(sprintf("%s%dco21ccr.xlsx", soi_path, y %% 100), col_names = FALSE, col_types = "text", .name_repair = "minimal")
  n <- unlist(x[str_squish(x[[1]]) %in% "Number of returns", 4:14])
  tibble(year = y, sector = 1, AC = 3:13, returns = as.numeric(n))
})


# 2: SOI, non-financial firms and the weights table ---------------------------------------------------------------

# share of returns in the dropped sectors, by asset class (2010-2013 pooled)
excl <- soi_sb %>% filter(sector != 1) %>%
  group_by(AC) %>% summarise(excl_sh = sum(returns[sector %in% drop_sectors], na.rm = T) / sum(returns, na.rm = T))

# non-financial counts: actual sector sums in 2010-2013, all industries * (1 - excluded share) in 2014-2022
soi_nf <- soi_sb %>% filter(sector != 1, !(sector %in% drop_sectors)) %>%
  group_by(year, AC) %>% summarise(n_nf = sum(returns, na.rm = T), .groups = "drop")

soi_weights <- bind_rows(soi_sb %>% filter(sector == 1), soi_t21) %>%
  select(year, AC, n_all = returns) %>%
  left_join(excl, by = "AC") %>%
  left_join(soi_nf, by = c("year", "AC")) %>%
  mutate(nf_imputed = is.na(n_nf),
         n_nf = ifelse(nf_imputed, n_all * (1 - excl_sh), n_nf)) %>%
  left_join(soi_classes, by = "AC") %>%
  group_by(year) %>% mutate(sh_all = n_all / sum(n_all), sh_nf = n_nf / sum(n_nf)) %>% ungroup() %>%
  select(year, AC, class, lo, hi, n_all, n_nf, sh_all, sh_nf, nf_imputed)

saveRDS(soi_weights, "soi_weights.RDS")

# excluded share by class, and the share of firms by class (non-financial, 2010-2022 pooled)
soi_weights %>% group_by(AC, class) %>% summarise(n_nf = sum(n_nf), .groups = "drop") %>%
  mutate(sh_nf = n_nf / sum(n_nf)) %>% left_join(excl, by = "AC") %>%
  mutate(across(c(sh_nf, excl_sh), ~ round(.x, 4))) %>% print()


# 3: Comparison: SOI vs IDB vs Compustat ----------------------------------------------------------------------------

# number of corporations by asset class over time
p1 <- soi_weights %>% pivot_longer(c(n_all, n_nf), names_to = "sample", values_to = "n") %>%
  mutate(sample = recode(sample, n_all = "All industries", n_nf = "Non-financial")) %>%
  ggplot(aes(x = year, y = n, color = fct_reorder(class, AC), linetype = sample)) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = seq(2010, 2022, by = 2)) +
  scale_y_log10(labels = scales::comma) +
  labs(x = "", y = "Number of returns (log scale)", color = "Total assets", linetype = "",
       title = "SOI: active corporations by size of total assets (non-financial 2014-2022: imputed)") +
  my_theme
print(p1)

# common asset classes: bounds that exist in both SOI and IDB
common_edges  <- c(0, 5e5, 1e6, 1e7, 5e7, 1e8, 5e8, Inf)
common_labels <- c("< $500k", "$500k - 1m", "$1 - 10m", "$10 - 50m", "$50 - 100m", "$100 - 500m", "> $500m")

idb <- readRDS("idb_cases.RDS") %>% filter(year(ORGFLDT) %in% years, !is.na(EASSTnum))
if (nonlead == "exclude") idb <- idb %>% filter(is.na(JNTLEAD) & is.na(CNSLLEAD))

# Compustat: one observation per firm-year (last quarter), assets in USD
comps <- readRDS("comps_fjc.RDS") %>% filter(year %in% years, assets > 0, sic_group != "Nonclassifiable") %>%
  group_by(gvkey, year) %>% filter(quarter == max(quarter)) %>% ungroup() %>%
  transmute(sample = "Compustat (all non-fin.)", assets = assets * 1e6)
flevel <- readRDS("flevel_full.RDS") %>% filter(year %in% years, assets > 0) %>%
  group_by(companyid, year) %>% filter(quarter == max(quarter)) %>% ungroup() %>%
  transmute(sample = "Compustat-CapIQ sample", assets = assets * 1e6)

dist_assets <- bind_rows(
  soi_weights %>% transmute(sample = "SOI all industries", upper = hi, n = n_all),
  soi_weights %>% transmute(sample = "SOI non-financial", upper = hi, n = n_nf),
  idb %>% transmute(sample = "IDB", upper = EASSTnum, n = 1),                                # IDB bins: upper bound
  bind_rows(comps, flevel) %>% transmute(sample, upper = assets * (1 + 1e-9), n = 1)) %>%     # "under" the bound, as in SOI
  mutate(size_class = cut(upper, common_edges, labels = common_labels, right = TRUE)) %>%
  group_by(sample, size_class) %>% summarise(n = sum(n), .groups = "drop") %>%
  group_by(sample) %>% mutate(share = n / sum(n)) %>% ungroup()

dist_assets %>% select(-n) %>% pivot_wider(names_from = sample, values_from = share) %>%
  mutate(across(where(is.double), ~ round(.x, 3))) %>% print(width = Inf)

p2 <- ggplot(dist_assets, aes(x = size_class, y = share, fill = sample)) +
  geom_col(position = "dodge", alpha = 0.8) +
  labs(x = "Total assets", y = "Share of firms", fill = "", title = "Firm size distribution by total assets") +
  my_theme + theme(legend.position = "bottom")
print(p2)

p3 <- ggplot(bind_rows(comps, flevel), aes(x = assets, fill = sample)) +
  geom_histogram(bins = 60, alpha = 0.5, position = "identity") +
  geom_vline(xintercept = common_edges[2:7], linetype = "dashed", color = "grey40") +   # common asset classes
  scale_x_log10(labels = scales::label_dollar(scale_cut = scales::cut_short_scale())) +
  scale_fill_manual(values = c("#377EB8", "#FF5555")) +
  labs(x = "Total assets (log scale)", y = "Firm-years", fill = "", title = "Compustat: total assets (dashed: common asset classes)") +
  my_theme + theme(legend.position = "bottom")
print(p3)
