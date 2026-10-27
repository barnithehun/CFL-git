# DEFAULT DATA: US BANKRUPTCY CASES OF BUSINESS DEBTORS (FJC INTEGRATED DATABASE)
# source: https://www.fjc.gov/research/idb
#     https://www.fjc.gov/research/idb/bankruptcy-cases-filed-terminated-and-pending-fy-2008-present
# INPUTS: raw FJC IDB snapshots in R/IDB dataset/ (DB10-DB25 .RDS) - read only
# OUTPUTS: idb_cases.RDS - cleaned bankruptcy cases, one row per case (used by liqprob.R and liqprob_specsearch.R)
#          Liqprob_OT.png - share of Ch.7 filings over time by liabilities (written to the paper's Plots folder)
#
# Structure:
#   1: Collection - assembling the business debtors from the raw IDB snapshots
#   2: Cleaning   - one case per snapshot, final chapter, unconfirmed Ch.11 cases, size bins
#   3: Asset based analysis       (firm size = estimated total assets, EASST)    - one row per case
#   4: Liabilities based analysis (firm size = estimated total liabilities, ELBLTS) - SBRA eligibility is defined on debt
#
# Firm size is only reported in bins (EASST / ELBLTS, A-L). Section 3 is the basis of the liquidation share figure (liqprob2.png)
# and the Ch.11 duration figure (daysin11.png) in the paper; section 4 gives the small/large split at 10m USD liabilities
# and the Ch.7 filing shares around the SBRA.

library(tidyverse)
library(tictoc)

rm(list = ls())

setwd("C://Users//szjud//OneDrive//Asztali gép//EBCs//CFL-git//Data Replication")
source("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V1//s_tab.R")

idb_path  <- "C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//IDB dataset//"   # raw FJC IDB files (read only)
plot_path <- "C://Users//szjud//OneDrive//Asztali gép//EBCs//CFL-git//Latex codes//Plots//"

# sample period of the paper (filing years) - used in the figures over time (sections 4c, 4d)
start_year <- 2010
end_year   <- 2024

my_theme <- theme(plot.title = element_text(face = "bold", hjust = 0.5, size =14),
                  panel.spacing = unit(2, "lines"), strip.background=element_rect(fill="white"),
                  axis.ticks.length=unit(0, "cm"),
                  axis.text.y = element_text(size = 12, angle = 0, vjust = 0.5, hjust=0.5),
                  axis.text.x = element_text(size = 12, angle = 75, vjust = 1, hjust=1),   # size bins do not fit horizontally
                  axis.title.x=element_text(size = 12), axis.title.y=element_text(size = 12),
                  panel.background = element_rect(fill = "white"),
                  panel.grid.major = element_line(size = 0.2, linetype = 'solid', colour = "grey80"),
                  panel.grid.minor = element_line(size = 0.15, linetype = 'solid', colour = "white"),
                  panel.border = element_rect(color = "grey50", fill = NA, size = 0.5))


# 1: Collection, assembling the full dataset --------------------------------------------------------------------

# things were originally downloaded in "sas7bdat" format but it took too much disk space
# I formatted them to RDS format.. ser original files at the FJC link above
# DB14to23 overlaps with DB23: 2023 snapshots are taken from DB23

#read in like this due to memory constraints
tic("Reading raw IDB files")
dat1 <- bind_rows(
  readRDS(paste0(idb_path, "DB10.RDS")),
  readRDS(paste0(idb_path, "DB11.RDS")),
  readRDS(paste0(idb_path, "DB12.RDS")),
  readRDS(paste0(idb_path, "DB13.RDS"))) %>% filter(DBTRTYP == "u" & NTRDBT == "b")  # I only keep the business debts by corporations!
dat2 <- readRDS(paste0(idb_path, "DB14to23.RDS")) %>% filter(year(SNAPSHOT) < 2023) %>% filter(DBTRTYP == "u" & NTRDBT == "b")
dat3 <- bind_rows(
  readRDS(paste0(idb_path, "DB23.rds")),
  readRDS(paste0(idb_path, "DB24.rds")),
  readRDS(paste0(idb_path, "DB25.rds"))) %>% filter(DBTRTYP == "u" & NTRDBT == "b")
toc()

# wide set of variables
idb_raw <- bind_rows(dat1, dat2, dat3) %>% select(c("SNAPSHOT", "SNAPFILE", "SNAPPEND", "SNAPCLOS", "CASEKEY", "ORGFLDT", "ORGFLCHP", "CRNTCHP", "ORGFLSUBCHP",
                                                   "EASST", "ELBLTS", "ASSTCASE", "SMLLBUS", "TOTASSTS", "SECURED", "UNSECPR", "UNSECNPR", "CLOSEDT", "CLCHPT",
                                                   "D1FDSP", "C11DVDND", "C11FTRPAY", "CNSLLEAD", "JNTLEAD", "ORIGIN" ))
rm(dat1, dat2, dat3)

# 2: Cleaning ---------------------------------------------------------
# starts from idb_raw (section 1), so the section can be re-run on its own without re-reading the raw files

## Firm Size
# TOTASSTS and EASST are not consistent! get rid of TOTASSTS an co.
# ggplot(idb_raw, aes(x = EASST, y = TOTASSTS)) +  geom_point() +  labs(x = "Category", y = "Value")
master <- idb_raw %>% select(-c("TOTASSTS", "SECURED", "UNSECPR", "UNSECNPR"))

# Clean up EASST
master <- master %>% filter(!grepl("^\\d+$|^$", EASST))

# GOAL: currently I have a case for each date-case - double check this. I need a single observation for each case.
# If it has been closed, use the closetype and values from the close date
# If it has not been closed use the currenttype

## Duplications
master <- master %>% group_by(SNAPSHOT, CASEKEY) %>%
  mutate(dup = ifelse(row_number() > 1, row_number(), 0))
s_tab(dup,master) # no duplications
master <- master %>%  select(-c("dup"))

# For each case, flag the last time it was observed: the final status of the case
# (filing date, filing chapter and size bins are practically the same in every snapshot of a case)
master <- master %>% arrange(SNAPSHOT) %>% group_by(CASEKEY) %>%
  mutate(last_obs = ifelse(row_number() == max(row_number()), 1, 0)) %>% ungroup()

# how long the case was open
master <- master %>% mutate(case_length = as.numeric(CLOSEDT - ORGFLDT))

# Only CH7 and CH11
master <- master %>% filter(!((CRNTCHP == 7 & CLCHPT == 11) | (CRNTCHP == 11 & CLCHPT == 7))) # weird observations in general

# If case is closed then CLCHPT if not closed then CRNTCHP
master <- master %>% mutate(chapter = ifelse(SNAPCLOS == 1, CLCHPT, CRNTCHP)) %>%
  filter((chapter == 7 | chapter == 11) & EASST != "M" & (ORGFLCHP == 7 | ORGFLCHP == 11))

# Consider cases that are probably not confirmed
s_tab( D1FDSP, master) # final disposition of debtor 1 (D1FDSP), all cases
s_tab( D1FDSP, master %>% filter(SNAPCLOS == 1 & CLCHPT == '11')) # final disposition of debtor 1, closed Ch.11 cases
master <- master %>%
  mutate(
    notconf = case_when(
      D1FDSP %in% c("A","B","C","D") ~ 0L,      # A–D
      D1FDSP %in% LETTERS[5:26]      ~ 1L,      # E–Z
      TRUE                           ~ NA_integer_))  # everything else (NA, blanks, weird codes) )

master <- master %>%  mutate(notconf = ifelse(is.na(notconf), 0, notconf),
                             chapter = ifelse(chapter == 11 & notconf == 1, 7, chapter))


## Size bins: EASST and ELBLTS use the same codes (cases filed from Dec 2007)
# value = upper bound of the bin (needed for interpolation later), label = used for EASST and on the figures' axes
size_bins <- tibble(
  code  = LETTERS[1:12],
  value = c(5e4, 1e5, 5e5, 1e6, 1e7, 5e7, 1e8, 5e8, 1e9, 5e9, 1e10, 5e10),
  label = c("$0 to $50k", "$50k to $100k", "$100k to $500k", "$500k to $1m", "$1 to $10m", "$10 to $50m",
            "$50 to $100m", "$100 to $500m", "$500 to $1000m", "$1b + (pre 2015)", "$1-$10b", "$10- $50b"))

# numeric bin values: liabilities (ELBLTSnum) and assets (EASSTnum); other codes (e.g. M, blank) get NA
master <- master %>% mutate(ELBLTSnum = size_bins$value[match(ELBLTS, size_bins$code)],
                            EASSTnum  = size_bins$value[match(EASST, size_bins$code)])

# axis labels for the size bins, keyed by the bin value (EASSTnum / ELBLTSnum) - used in the figures of sections 3 and 4
size_labels <- setNames(size_bins$label, as.character(size_bins$value))

# EASST gets its labels (ELBLTS stays in letters)
master <- master %>% mutate(EASST = size_bins$label[match(EASST, size_bins$code)])

## Final dataset: one row per case (last snapshot) - used in sections 3 and 4, saved for liqprob.R and liqprob_specsearch.R
# type_convert() gives the same column types as the old write.csv / read_csv round trip
idb_cases <- master %>% filter(last_obs == 1) %>% type_convert()


saveRDS(idb_cases, "idb_cases.RDS")
rm(master)


# 3a: Asset based analysis, liquidation share across firm sizes -------------------------------------------------
master <- idb_cases

# How often do firms convert from liquidation to reorg (closed cases, one row per case)
# conv_rates() prints the share of closed Ch.7 filings that end in Ch.11 and the other way around:
#   court record: chapter at closing (CLCHPT) - conversions ordered by the court
#   final status: chapter - also counts unconfirmed Ch.11 cases (notconf == 1) as Ch.7, as in the analysis below
conv_rates <- function(d) {
  for (from in c(7, 11)) {
    to <- ifelse(from == 7, 11, 7)
    x <- filter(d, ORGFLCHP == from)
    cat(sprintf("Filed as Ch.%d, ended as Ch.%d: %.2f%% (court record), %.2f%% (final status) - of %d closed Ch.%d filings\n",
                from, to, 100 * mean(x$CLCHPT == to), 100 * mean(x$chapter == to), nrow(x), from))
  }
}

master_closed <- filter(master, SNAPCLOS == 1)
conv_rates(master_closed) # Liquidaiton filins almost never get overturned into reorganization

# master <- master %>% filter(is.na(JNTLEAD), is.na(CNSLLEAD) , is.na(ORIGIN))

# summarized dataset for liquidation probability - size classification based on ASSETS (EASST)
sum_master <- master %>% group_by(EASST, chapter) %>%
  summarise(count = n(), totlib = mean(EASSTnum)) %>% spread(chapter, count) %>% mutate(pr_liq = `7` / (`7`+`11`)) %>% arrange(totlib)


# Panel A: number of Ch.7 and Ch.11 cases by asset bin
p1 <- ggplot(master %>% filter(year(SNAPSHOT) >= start_year), aes(x = as.factor(EASSTnum), fill = as.factor(chapter))) +
  geom_bar(position = "dodge", stat = "count", alpha = 0.8) +
  labs(x = "Estimated total assets", y = "Number of cases (thousands)", fill = "Chapter:") +
  ggtitle("Chapter 7 and Chapter 11 Filings") +
  scale_y_continuous(expand = c(0.02, 0), labels = scales::comma_format(scale = 1e-3)) +
  scale_fill_manual(values = c("#FF5555", "#377EB8")) +
  scale_x_discrete(labels = function(b) size_labels[b]) +
  my_theme +
  theme(legend.position = c(0.8, 0.9),  # Adjust the x and y coordinates as needed
        legend.text = element_text(size = 10),
        legend.key.size = unit(1.1, "lines"))


# Panel B: share of Ch.7 cases by asset bin (liqprob2.png in the paper)
p2 <- ggplot(sum_master, aes(x = as.factor(totlib), y = pr_liq, group = 1)) +
  geom_line(color = "#377EB8", linewidth = 1.4, alpha = 1) +
  labs(x = "Estimated Total Assets", y = "Chapter 7 to Chapter 11 filings") +
  ggtitle("Share of Liqudations across Firms Sizes") +
  scale_y_continuous(expand = c(0.01,0), limits = c(0,1)) +
  scale_x_discrete(labels = function(b) size_labels[b]) +
  my_theme

cowplot::plot_grid(p1, p2, ncol = 2)
#ggsave(paste0(plot_path, "liqprob2.png"), dpi = 400, height = 6, width = 10)


# 3b: Calibration Targets - Based on assets  ----------------------------------

# Calibration targets
master <- idb_cases
# master <- idb_cases %>% filter(is.na(JNTLEAD), is.na(CNSLLEAD) , is.na(ORIGIN))

# number of closed cases (one row per case)
sum(!is.na(master$CLCHPT))

# Smaller than 10M in assets
summaster <- master %>% filter(EASSTnum > 100000) %>%  mutate(Less10m = ifelse(EASSTnum <= 10^7, 1, 0))
s_tab(Less10m, summaster)

# Liquidation share by size (one row per case)
summaster %>% group_by(Less10m, chapter) %>% summarise(count = n()) %>% spread(chapter, count) %>% mutate(pr_liq = `7` / (`7`+`11`))


# 3c: Asset based analysis, length of Ch.11 cases (reorg length) ----------------------------------------------------------
master <- idb_cases

# Collapse the first two asset bins into one: "$0 to $100k"
sum_master <- master %>%
  filter(chapter == 11) %>%
  mutate(
    totlib = ifelse(EASSTnum <= 1e5, 1e5, EASSTnum),
    EASST  = ifelse(EASSTnum <= 1e5, "$0 to $100k", EASST)
  ) %>%
  group_by(EASST) %>%
  summarise(
    count = n(),
    case_length = mean(case_length, na.rm = TRUE),
    length_scaled = mean(case_length / totlib, na.rm = TRUE) * 1e6,
    totlib = first(totlib),
    .groups = "drop"
  ) %>%
  arrange(totlib)

# Plot 1: case_length by asset bin
p1 <- ggplot(sum_master,
             aes(x = factor(totlib, levels = sum_master$totlib),
                 y = as.numeric(case_length),
                 group = 1)) +
  geom_line(color = "#377EB8", linewidth = 1.4) +
  labs(x = "Estimated Total Assets",
       y = "Days in Chapter 11") +
  scale_x_discrete(labels = sum_master$EASST) +   # sum_master is sorted by totlib, as the x axis
  my_theme

# Plot 2: length_scaled by asset bin
p2 <- ggplot(sum_master,
             aes(x = factor(totlib, levels = sum_master$totlib),
                 y = as.numeric(length_scaled),
                 group = 1)) +
  geom_line(color = "#377EB8", linewidth = 1.4) +
  labs(x = "Estimated Total Assets",
       y = "Days / Millions of Assets") +
  scale_x_discrete(labels = sum_master$EASST) +   # sum_master is sorted by totlib, as the x axis
  my_theme

cowplot::plot_grid(p1, p2, ncol = 2)


# 4a: Liabilities based analysis, liquidation share across firm sizes ----------------------------------------------------
master <- idb_cases

# summarized dataset for liquidation probability - size classification based on LIABILITIES (ELBLTS)
sum_master <- master %>% group_by(ELBLTS, chapter) %>%
  summarise(count = n(), totlib = mean(ELBLTSnum)) %>% spread(chapter, count) %>% mutate(pr_liq = `7` / (`7`+`11`)) %>% arrange(totlib)


p1 <- ggplot(master %>% filter(year(SNAPSHOT) >= start_year), aes(x = as.factor(ELBLTSnum), fill = as.factor(chapter))) +
  geom_bar(position = "dodge", stat = "count", alpha = 0.8) +
  labs(x = "Estimated total liabilities", y = "Number of cases (thousands)", fill = "Chapter:") +
  ggtitle("Chapter 7 and Chapter 11 Filings") +
  scale_y_continuous(expand = c(0.02, 0), labels = scales::comma_format(scale = 1e-3)) +
  scale_fill_manual(values = c("#FF5555", "#377EB8")) +
  scale_x_discrete(labels = function(b) size_labels[b]) +
  my_theme +
  theme(legend.position = c(0.8, 0.9),  # Adjust the x and y coordinates as needed
        legend.text = element_text(size = 10),
        legend.key.size = unit(1.1, "lines"))


p2 <- ggplot(sum_master, aes(x = as.factor(totlib), y = pr_liq, group = 1)) +
  geom_line(color = "#FF5555", linewidth = 2, alpha = 1) +
  labs(x = "Estimated Total Liabilities", y = "Chapter 7 to Chapter 11 filings") +
  ggtitle("Share of Liquidations") +
  scale_y_continuous(expand = c(0.01,0), limits = c(0,1)) +
  scale_x_discrete(labels = function(b) size_labels[b]) +
  my_theme

cowplot::plot_grid(p1, p2, ncol = 2)
# it would be nice too see these things over time as well

# 4b: Calibration targets based on Liabilities  --------------------------------------------------------

master <- idb_cases
# master <- idb_cases %>% filter(is.na(JNTLEAD), is.na(CNSLLEAD) , is.na(ORIGIN))

# number of closed cases (one row per case)
sum(!is.na(master$CLCHPT))

# Smaller than 10M in assets
summaster <- master %>% filter(ELBLTSnum > 100000) %>%  mutate(Less10m = ifelse(ELBLTSnum <= 10^7, 1, 0))
# summaster <- master %>% filter(EASSTnum > 100000) %>%  mutate(Less10m = ifelse(ELBLTSnum <= 10^7, 1, 0))# I think this was originally used for the calibration targets in the paper!
s_tab(Less10m, summaster)

# Liquidation share by size (one row per case)
summaster %>% group_by(Less10m, chapter) %>% summarise(count = n()) %>% spread(chapter, count) %>% mutate(pr_liq = `7` / (`7`+`11`))


# 4c: Liabilities based analysis, liquidation probability over time (SBRA) ------------------------------------------------
# status at filing: filing chapter (ORGFLCHP) by filing year
master <- idb_cases %>%
  mutate(Less50m = case_when(
    ELBLTSnum <= 1e7                       ~ 1L,  # ≤ 10m
    ELBLTSnum > 1e7 & ELBLTSnum <= 5e7   ~ 2L,  # 10m - 50m
    TRUE                                  ~ 3L   # > 50m
  )) # %>% filter(is.na(JNTLEAD), is.na(CNSLLEAD) , is.na(ORIGIN))

sum_master <- master %>% group_by(year(ORGFLDT), Less50m, ORGFLCHP) %>%
  summarise(count = n()) %>%
  spread(ORGFLCHP, count) %>% mutate(pr_liq = `7` / (`7`+`11`)) %>%
  rename("year" = "year(ORGFLDT)", "SBRA" = "Less50m") %>% filter(year >= start_year, year <= end_year)


year_labels <- seq(start_year, end_year, by = 2)  # every 2nd year
year_breaks <- year_labels    #

ggplot(sum_master, aes(x =year , y = pr_liq, color = as.factor(SBRA))) +
  my_theme +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 0.6),                      # 0% to 60%
    breaks = seq(0, 0.6, by = 0.1)           # ticks every 10%
  ) +
  scale_x_continuous(
    breaks = year_breaks,
    labels = year_labels) +
  scale_color_manual(
    name   = "Liabilities M USD",
    values = c("1" = "#377EB8", "2" = "#FF5555", "3" = "green4"),
    labels = c("1" = "< 10M USD ", "2" = "10–50M USD", "3" = "50+M USD")) +
  labs(x = "", y = "Share of Chapter 7 Filings", title = "Share of Liquidations" ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),   # rotate x-axis ticks
    legend.position = "bottom",                         # move legend
    legend.text = element_text(size = 13),              # bigger legend text
    legend.title = element_text(size = 13),             # bigger legend title
    plot.title = element_text(size = 20))

 ggsave(paste0(plot_path, "Liqprob_OT.png"), dpi = 400, height = 5, width = 6)


# 4d: Liabilities based analysis, filings by chapter and Subchapter V (liabilities <= 10m) ------------------------------------------------
# one row per case, so each case is counted once (the Subchapter V flag is most complete in the last snapshot)
sum_master <- idb_cases %>% # filter(is.na(JNTLEAD), is.na(CNSLLEAD) , is.na(ORIGIN)) %>%
  filter(year(ORGFLDT) >= start_year, year(ORGFLDT) <= end_year, ELBLTSnum <= 1e7) %>%
  mutate(
    ORGFLCHP = ifelse(!is.na(ORGFLSUBCHP) & ORGFLSUBCHP == "V", 5, ORGFLCHP)
  ) %>%
  group_by(year = year(ORGFLDT), ORGFLCHP) %>%
  summarise(count = n(), .groups = "drop") %>%
  spread(ORGFLCHP, count, fill = 0) %>%    # fill missing chapter counts with 0
  print(n = 100) 

# Reshape back to long for plotting (stacked bars need long format)
sum_master_long <- sum_master %>%
  pivot_longer(-year, names_to = "Chapter", values_to = "Count")

# Make Chapter a factor so "5" shows up nicely
sum_master_long <- sum_master_long %>%
  mutate(Chapter = factor(Chapter, levels = c("7","11","5")))

# Plot
ggplot(sum_master_long, aes(x = factor(year), y = Count, fill = Chapter)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = c("7" = "#377EB8", "11" = "#FF5555", "5" = "green4"),
    labels = c("7" = "Chapter 7", "11" = "Chapter 11", "5" = "Subchapter V")
  ) +
  labs(x = "Filing Year", y = "Number of Cases", fill = "Filing Chapter") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14))

