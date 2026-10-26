# RESULTS: SUMMARY STATISTICS, CREDIT SPREAD REGRESSIONS AND FIRM LEVEL RESULTS - 2010-2024
# INPUTS ARE THE OUTPUTS OF data_cleaning.R: dlevel_full.RDS AND flevel_full.RDS

library(tidyverse)
library(lubridate)
library(estimatr)      # lm_robust: OLS with Stata-type robust standard errors
library(modelsummary)  # regression tables (replaces outreg2)
library(fixest)        # fixed effects regressions (replaces xtreg, fe)

rm(list = ls())

setwd("C:/Users/szjud/OneDrive/Asztali gép/EBCs/CFL-git/Data Replication")
source("C:/Users/szjud/OneDrive/Asztali gép/EBCs/Data/R/V1/s_tab.R")

plot_path <- "C://Users//szjud//OneDrive//Asztali gép//EBCs//CFL-git//Latex codes//Plots//"

# sample period of the paper: 2010Q1 - 2024Q4
end_year <- 2024

my_theme <- theme(plot.title = element_text(face = "bold", hjust = 0.5, size =20),
                  panel.spacing = unit(2, "lines"), strip.background=element_rect(fill="white"),
                  axis.ticks.length=unit(0, "cm"),
                  axis.text.y = element_text(size = 16, angle = 0, vjust = 0.5, hjust=0.5),
                  axis.text.x = element_text(size = 16, angle = 0, vjust = 0.5, hjust=0.5),
                  axis.title.x=element_text(size = 18), axis.title.y=element_text(size = 18),
                  panel.background = element_rect(fill = "white"),
                  panel.grid.major = element_line(size = 0.2, linetype = 'solid', colour = "grey80"),
                  panel.grid.minor = element_line(size = 0.15, linetype = 'solid', colour = "white"),
                  panel.border = element_rect(color = "grey50", fill = NA, size = 0.5))

# log that returns NA for non-positive values - same behaviour as log() in Stata
slog <- function(x) if_else(x > 0, log(x), NA_real_)


# 1: Dlevel, summary statistics ---------------------------------------
dlevel <- readRDS("dlevel_full.RDS") %>% filter(year <= end_year)

# share of CF-based debt by volume
dlevel %>%
  filter(!is.na(CF)) %>%
  summarise(CF_val = sum(value[CF == 1], na.rm = TRUE),
            AB_val = sum(value[CF == 0], na.rm = TRUE),
            CFshare_count = mean(CF == 1)) %>%
  mutate(CFshare_volume = CF_val / (CF_val + AB_val))

# number of debt contracts
unique(dlevel$id) %>% length()
# number of firms
unique(dlevel$companyid) %>% length()
# number of firms that hold only a single debt contract (footnote 27)
s_tab(Ncontr, dlevel %>% group_by(companyid, year, quarter) %>% summarise(Ncontr = n()))

# table A5 - debt level summary statistics
# NOTE: only CF-based contracts, as in the original code
summary_stats <- dlevel %>% ungroup() %>% mutate(value = value / 10^6) %>% filter(CF == 1) %>%
  select(c("value", "intrate", "matumax", "spread")) %>%
  summarise_all(list(
    mean = ~mean(., na.rm = TRUE),
    p10 = ~quantile(., 0.1, na.rm = TRUE),
    p25 = ~quantile(., 0.25, na.rm = TRUE),
    median = ~median(., na.rm = TRUE),
    p75 = ~quantile(., 0.75, na.rm = TRUE),
    p90 = ~quantile(., 0.9, na.rm = TRUE)
  )) %>%
  gather(statistic, value) %>%
  separate(statistic, into = c("variable", "statistic"), sep = "_") %>%
  spread(statistic, value) %>%
  xtable::xtable(caption = "Summary Statistics") %>%
  print(type = "latex", include.rownames = TRUE, booktabs = TRUE, sanitize.text.function = function(x) x)


# _Misc:: plotting conditional means (not in the paper) ---------------------------------
ggplot(data = dlevel %>% filter(inflow == 1 & lassets > 5)) +
  geom_point(aes(x = lassets, y = spread, color = as.factor(CF)), alpha = 0.04, size = 1) +
  geom_smooth(aes(x = lassets, y = spread, color = as.factor(CF)), method = "lm", se = FALSE, size = 1) +
  scale_color_manual(name = "Loan-Type", values = c("0" = "firebrick", "1" = "blue4"), labels = c("AB debt", "CF debt")) +
  labs(title = "Credit Spreads Over Firm Assets", x = "Log of Assets", y = "Credit Spread") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = c(0.98, 0.98), legend.justification = c("right", "top"), legend.key.size = unit(0.8, 'cm'), legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)) +
  my_theme

ggsave("spreads.png", dpi = 150, height = 6, width = 10)

ggplot(data = dlevel %>% filter(inflow == 1 & age < 100 & lassets > 5)) +
  geom_point(aes(x = age, y = spread, color = as.factor(CF)), alpha = 0.04, size = 1) +
  geom_smooth(aes(x = age, y = spread, color = as.factor(CF)), method = "lm", se = FALSE, size = 1) +
  scale_color_manual(name = "Loan-Type", values = c("0" = "firebrick", "1" = "blue4"), labels = c("AB debt", "CF debt")) +
  labs(title = "Credit Spreads Over Firm Age", x = "Firm Age by Years", y = "Credit Spread") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = c(0.98, 0.98), legend.justification = c("right", "top"), legend.key.size = unit(0.8, 'cm'), legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)) +
  my_theme

ggsave("spreads_age.png", dpi = 150, height = 6, width = 10)


# 2: Dlevel, determinants of credit spreads (replicates dlevel_regress.do) ---------------------------------------

# number of contracts per firm-quarter - counted before keeping only inflows, as in Stata
dreg <- dlevel %>%
  group_by(companyid, yq) %>%
  mutate(Nconts = sum(AB %in% c(0, 1))) %>% 
  ungroup() 

# only looking at current inflows (new issuances)
dreg <- dreg %>% filter(inflow == 1) %>% 
  mutate(rating = relevel(factor(rating), ref = "A+")) 

# factors (= encode in Stata) and log transformations
dreg <- dreg %>% 
  mutate(sic_group = factor(sic_group),
         rating = factor(rating),
         dateper = factor(dateper),
         lemp = slog(emp),
         lage = slog(age),
         lebitda = sign(ebitda) * log(abs(ebitda) + 1))   # log-modulus transformation

# MAIN RESULTS - TABLE 1 AND A7: WHAT EXPLAINS THE CREDIT SPREAD FOR NEW ISSUANCES
f_spread <- spread ~ lebitda + lassets + pledge + levrg + lage + lemp + Nconts + CFshare + sic_group + rating + dateper

# se_type = "stata" gives the same robust standard errors as ', robust' in Stata
reg_all <- lm_robust(f_spread, data = dreg, se_type = "stata")
# Only CF debt contracts
reg_cf <- lm_robust(f_spread, data = droplevels(filter(dreg, CF == 1)), se_type = "stata")
# Only AB debt contracts
reg_ab <- lm_robust(f_spread, data = droplevels(filter(dreg, CF == 0)), se_type = "stata")

modelsummary(list("All contracts" = reg_all, "CF contracts" = reg_cf, "AB contracts" = reg_ab),
             output = "table_outcome.tex",
             coef_omit = "dateper",
             stars = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
             gof_map = c("nobs", "r.squared"))


# ROBUSTNESS - SAME PICTURE WITH INTEREST RATES AND MATURITY AS CONTROL (not in the paper)
f_intrate <- intrate ~ lebitda + pledge + lassets + levrg + lage + lemp + matumax + sic_group + rating + dateper

rob_all <- lm_robust(f_intrate, data = dreg, se_type = "stata")
# Only CF debt contracts - non-robust SE, as in the original code
rob_cf <- lm_robust(f_intrate, data = droplevels(filter(dreg, CF == 1)), se_type = "classical")
# Only AB debt contracts - non-robust SE, as in the original code
rob_ab <- lm_robust(f_intrate, data = droplevels(filter(dreg, CF == 0)), se_type = "classical")

modelsummary(list("All contracts" = rob_all, "CF contracts" = rob_cf, "AB contracts" = rob_ab),
             coef_omit = "dateper",
             stars = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
             gof_map = c("nobs", "r.squared"))


rm(list = setdiff(ls(), c("slog", "s_tab", "my_theme", "plot_path", "end_year")))

# 3: Flevel, SBRA eligible and SME classifications ------------------------------------------------------------
flevel <- readRDS("flevel_full.RDS") %>% filter(year <= end_year)

flevel %>%
  ungroup() %>%
  summarise(CF_val = sum(CF_val, na.rm = TRUE),
            AB_val = sum(AB_val, na.rm = TRUE)) %>%
  mutate(CFshare_volume = CF_val / (CF_val + AB_val))


# SME classification (Not used in the paper)
flevel <- flevel %>% ungroup() %>%
  mutate(
    SME = case_when(
      (sic_group == "Manufacturing" | sic_group == "Construction") & emp <= 0.5 ~ 1,
      (sic_group == "Mining" |  sic_group == "Agriculture, Forestry, & Fishing") & emp <= 0.5  ~ 1,
      sic_group == "Retail Trade" & emp <= 0.5  ~ 1,
      sic_group == "Services" & emp <= 0.5  ~ 1,
      sic_group == "Transportation & P.U" & emp <= 0.5  ~ 1,
      sic_group == "Wholesale Trade" & emp <= 0.5 ~ 1,
      TRUE ~ 0 # Default for exceeding caps or unknown sic_groups
    ))

# These are the firms that can be classified as SME - with a pretty generous classification
s_tab(SME, flevel)
flevel$assets[flevel$SME == 0]  %>%  summary()


# 4: Flevel, tables -----------------------------------------------------------

# table A4 - firm level summary statistics
summary_stats <- flevel %>% ungroup() %>% filter(year >= 2010) %>%
  select(c("assets", "revenue", "emp", "eq", "inv_rate",  "prod", "age", "value", "d2c", "d2a", "levrg", "pledge", "liq", "maturity", "intrate", "spread", "CFshare")) %>%
  rename(invrate = inv_rate) %>%
  mutate(pledge = pledge * 100) %>%
  summarise_all(list(
    mean = ~mean(., na.rm = TRUE),
    p10 = ~quantile(., 0.1, na.rm = TRUE),
    p25 = ~quantile(., 0.25, na.rm = TRUE),
    median = ~median(., na.rm = TRUE),
    p75 = ~quantile(., 0.75, na.rm = TRUE),
    p90 = ~quantile(., 0.9, na.rm = TRUE)
  )) %>%
  gather(statistic, value) %>%
  separate(statistic, into = c("variable", "statistic"), sep = "_") %>%
  spread(statistic, value) %>%
  xtable::xtable(caption = "Summary Statistics") %>%
  print(type = "latex", include.rownames = TRUE, booktabs = TRUE, sanitize.text.function = function(x) x)


# Calibration targets 


# SBRA classification (used in the paper - will need to change later)
# This will be used to weigh data moments for small firms and non-small firms
flevel <- flevel %>%  mutate(SBRA = ifelse(liab < 10, 1,0))
s_tab(SBRA, flevel)
# the value of assets for these firms
flevel$assets[flevel$SBRA == 1]  %>%  summary()


# Target moments for calibration
summary_stats <- flevel %>% ungroup() %>% filter(SBRA == 1 & year >= 2010) %>%
  select(c("d2a", "d2c", "spread", "pledge", "intrate", "CFshare")) %>%
  summarise_all(list(
    mean = ~mean(., na.rm = TRUE),
    p10 = ~quantile(., 0.1, na.rm = TRUE),
    p25 = ~quantile(., 0.25, na.rm = TRUE),
    median = ~median(., na.rm = TRUE),
    p75 = ~quantile(., 0.75, na.rm = TRUE),
    p90 = ~quantile(., 0.9, na.rm = TRUE)
  )) %>%
  gather(statistic, value) %>%
  separate(statistic, into = c("variable", "statistic"), sep = "_") %>%
  spread(statistic, value) %>%
  xtable::xtable(caption = "Summary Statistics") %>%
  print(type = "latex", include.rownames = TRUE, booktabs = TRUE, sanitize.text.function = function(x) x)

# Share of ABL and CFL firms: 0 = only AB, 1 = both, 2 = only CF (footnote 27)
flevel <- flevel %>%  group_by(companyid) %>%
  mutate(CFcat = ifelse(CFshare == 0, 0,
                        ifelse(CFshare > 0 & CFshare < 1, 1, 2)))
s_tab(CFcat, flevel)[["percent"]]

# share of firms, assets and debt by CFcat
flevel %>% group_by(CFcat) %>%
  summarise(numshare = n(),
            totass = sum(assets, na.rm = T),
            totloan = sum(value, na.rm = T)) %>%
  mutate(numshare = numshare / sum(numshare) * 100,
         totass = totass / sum(totass) * 100,
         totloan = totloan / sum(totloan) * 100) %>%
  xtable::xtable( caption = "Summary Statistics") %>%
  print(type = "latex", include.rownames = TRUE, booktabs = TRUE,
        sanitize.text.function = function(x) x)

# share of CF-based debt by volume
CFshare <- flevel %>% ungroup() %>%  filter(year >= 2010) %>%
  mutate(ABdebt = (1-CFshare)*value ,
         CFdebt = CFshare*value) %>%
  summarise(CFdebt = sum(CFdebt, na.rm = T), ABdebt = sum(ABdebt, na.rm = T)) %>%
  mutate(CFshare = CFdebt / (CFdebt + ABdebt))
CFshare


# 5: Flevel, figure 5 - share of CF-based debt across firm classes -------------------------------------------------------------------------
# NOTE: classes are based on total debt from CapIQ (value, millions USD), not on liabilities
sumdat <- flevel %>% ungroup() %>%
  mutate(SBRA = case_when(
    value < 10 ~ 1,               # SBRA group 1
    value >= 10 & value < 50 ~ 2, # SBRA group 2
    TRUE ~ 3                      # SBRA group 3 (rest)
  )) %>%
  group_by(year, SBRA) %>%
  summarize(
    AB_val = sum(AB_val, na.rm = TRUE),
    CF_val = sum(CF_val, na.rm = TRUE),
    CFrel  = mean(CFshare, na.rm = TRUE)
  ) %>%
  mutate(CFshare = CF_val / (AB_val + CF_val))

# label only by years (one tick per two years)
year_labels <- seq(2010, 2024, by = 2)

ggplot(sumdat, aes(x = year, y = CFshare, color = as.factor(SBRA))) +
  my_theme +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  scale_x_continuous(
    breaks = year_labels,
    labels = year_labels) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(0.1, 0.9, 0.2),   # 10%, 30%, 50%, 70%, 90%
    limits = c(0.1, 0.9)           # optional: force axis range 10–90%
  ) +
  scale_color_manual(
    name   = "Liabilities M USD",
    values = c("1" = "#377EB8", "2" = "#FF5555", "3" = "green4"),
    labels = c("1" = "< 10M USD ", "2" = "10–50M USD", "3" = "50+M USD")) +
  labs(x = "", y = "CF-based Debt Share", title = "Reliance on CF-based debt" ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    plot.title = element_text(size = 20))

ggsave(paste0(plot_path, "CFshare_firms.png"), dpi = 400, height = 5, width = 6)


# _Misc:: CF-share and pledgeability (not in the paper) ----------------------------------------------

# Define the intervals
intervals <- c(0, 0.2, 0.4, 0.6, 0.8)

# Define a function to run regression for each interval
run_regression <- function(lower) {

  data <- flevel %>% ungroup() %>%
    filter(CFshare >= lower, CFshare < lower + 0.2,
           !is.na(spread) & !is.infinite(spread),
           !is.na(ebit2a) & !is.infinite(ebit2a),
           !is.na(pledge) & !is.infinite(pledge),
           !is.na(lassets) & !is.infinite(lassets),
           !is.na(levrg) & !is.infinite(levrg),
           !is.na(age) & !is.infinite(age),
           !is.na(emp) & !is.infinite(emp),
           !is.na(sic_group),
           !is.na(rating),
           !is.na(dateper))

  res <- lm(spread ~ ebit2a + pledge + lassets + levrg + age + emp + factor(sic_group) + factor(rating) + dateper, data)  %>% summary()

  pledge <- cbind(res[["coefficients"]][3,1], res[["coefficients"]][3,2] * 1.96)
  assets <- cbind(res[["coefficients"]][4,1], res[["coefficients"]][4,2] * 1.96)

  rbind(pledge, assets) %>%  return()
}

dat <- t(sapply(intervals, function(lower) run_regression(lower)))
dat <- rbind(dat[, 1:2],dat[, 3:4]) %>% cbind(intervals) %>%  cbind(c(rep("Pledgeability", 5), rep("Log of Assets", 5))) %>%
  cbind(rep(c("0 - 0.2", "0.2 - 0.4", "0.4 - 0.6", "0.6 - 0.8", "0.8 - 1"), 2)) %>% as.data.frame() %>%
  mutate(V1 = as.numeric(V1), V2 = as.numeric(V2), intervals = as.numeric(intervals))
colnames(dat) <- c("Estimate", "CI", "Interval", "Var", "INTV")

ggplot(dat, aes(x = INTV, y = Estimate , colour = Var)) +
  geom_hline(yintercept = 0, size = 1, color = "gray", linetype = "dashed") +
  geom_errorbar(aes(ymin = Estimate - CI, ymax = Estimate + CI),
                size = 1.1, width = 0.4, alpha = 0.8, position = position_dodge(width = 0.2)) +
  geom_point(shape = 15, size = 2, alpha = 1, position = position_dodge(width = 0.2)) +
  geom_line(size = 0.8, alpha = 0.8, position = position_dodge(width = 0.2)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  xlab("CFL reliance") +
  ylab("Coefficient Estimate") +
  scale_color_brewer(palette = "Set1") + my_theme +
  ggtitle("Coefficient Estimates: Assets and Pledgeability") +
  theme(plot.title = element_text(size = 10), legend.title = element_blank(),
        legend.text = element_text(size = 8), legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"))

ggsave("mixplot.png", dpi = 400, height = 6, width = 10)
