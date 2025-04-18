library(tidyverse)
library(tictoc)
library(RPostgreSQL)
library(RPostgres)
library(tictoc)
library(lubridate)


rm(list = ls())

setwd("C:/Users/szjud/OneDrive/Asztali gép/EBCs/Data/R/V3")
source("C:/Users/szjud/OneDrive/Asztali gép/EBCs/Data/R/V1/s_tab.R")

my_theme <- theme(plot.title = element_text(face = "bold", hjust = 0.5, size =18),
                  panel.spacing = unit(2, "lines"), strip.background=element_rect(fill="white"),
                  axis.ticks.length=unit(0, "cm"),
                  axis.text.y = element_text(size = 10, angle = 0, vjust = 0.5, hjust=0.5),
                  axis.text.x = element_text(size = 10, angle = 0, vjust = 0.5, hjust=0.5),
                  axis.title.x=element_text(size = 12), axis.title.y=element_text(size = 12),
                  panel.background = element_rect(fill = "white"),
                  panel.grid.major = element_line(linewidth = 0.2, linetype = 'solid', colour = "grey80"),
                  panel.grid.minor = element_line(linewidth = 0.15, linetype = 'solid', colour = "grey90"),
                  panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.5))



# Connecting and compustat cleaning ---------------------------------------
  flevel <- read_csv("flevel.csv")
  comps <- read_csv("comps.csv")
  
  
  # drop of the financial companies
  comps <- comps %>%
    mutate(fincomp = ifelse(sic >= 6000 & sic <= 6799, 1, 0),                # by SIC code
           fincomp = ifelse(sic >= 9000 & indfmt == "FS", 1, fincomp)) %>%   # by Compustat classific
    filter(fincomp != 1) %>% 
    filter(!(sic > 4900 & sic < 4999)) # these are utility companies
  
  # gvkey duplicates
  comps <- comps %>% 
    group_by(gvkey, year, quarter) %>%
    mutate(dup = ifelse(row_number() > 1, row_number(), 0)) %>%  
    filter(dup < 1) %>% select(-c("fincomp", "indfmt"))
  
  # Convert 'sic' to string
  comps <- comps %>%
    mutate(sic = as.character(sic),
           aux = substr(sic, 1, 2),
           aux = as.numeric(aux))
  
  # Define the lookup table with 'aux2' and corresponding 'sic_group'
  lookup_table <- data.frame(
    aux = c(1:9, 10:14, 15:17, 20:39, 40:49, 50:51, 52:59, 60:67, 70:89, 90:98, 99),
    sic_group = c(
      rep("Agriculture, Forestry, & Fishing", 9), rep("Mining", 5),  rep("Construction", 3), rep("Manufacturing", 20),
      rep("Transportation & P.U", 10), rep("Wholesale Trade", 2), rep("Retail Trade", 8), rep("Finance, Insurance, & Real Estate", 8),
      rep("Services", 20), rep("Public Admin", 9),"Nonclassifiable"))
  
  # Perform the join operation to assign 'sic_group'
  comps <- comps %>%
    left_join(lookup_table, by = "aux")
  rm(list = "lookup_table")
  
  # create derived variables: 
  comps <- comps %>%  
    mutate(atq = if_else(atq == 0, NaN , atq), 
           tot_debt = dlcq + dlttq,
           net_debt = tot_debt - chq,
           levrg = tot_debt/atq,
           collat = ppentq + invtq + rectq + chq,
           pledge = collat/atq,
           intcov = oibdpq /xintq,
           eq = atq - ltq,
           liq = chq / atq,
           prod = ((revtq - cogsq)/1000) / ((ppegtq/1000)^(1/3) + emp^(2/3))) %>%  
    group_by(gvkey, year) %>%  
    mutate(sppeq = ifelse(quarter == 1, sppey, sppey - lag(sppey)),
           capxq = ifelse(quarter == 1, capxy, capxy - lag(capxy)),
           capxq = ifelse(capxq < 0, NA, capxq)) %>%
    group_by(gvkey) %>% 
    mutate(inv = capxq - sppeq,
           inv_rate =  100*(capxq - sppeq)/ lag(atq),
           lassets = log10(atq*1000000)) 
  # the correct definition would be using 'ppegtq' instead of assets but that has to many missing values 
  
  
    comps <- comps %>% rename(ebitda = oibdpq, assets = atq, liab = lctq, sec_debt = dm, rating = spcsrc, revenue = revtq) %>%  
    select(c("datadate", "revenue","sic","rating","year","quarter","emp","sic_group", "ebitda", "sec_debt",
             "tot_debt","net_debt","levrg","collat","pledge","intcov","eq",  "liq", "assets", "liab", "inv", "inv_rate", "prod", "lassets"))
  
  write_csv(comps, "comps_fjc.csv")
  
  
  # joining flevel and comps
  comciq <- flevel %>% left_join(comps, by = c("gvkey", "year", "quarter"))
  rm(list = c("comps", "flevel"))
  
  # dropping all observations that are not in CAPIQ
  comciq <-comciq %>% filter(!is.na(datadate))
  
  # new variables that could be generated after connecting CAPIQ
  comciq <- comciq %>%   rename(compsdebt = tot_debt) %>% mutate(
    emp = if_else(emp == 0, NaN , emp),
    sec_share = sec_debt/compsdebt,
    age = year - yearfounded + 1,
    value = (AB_val+CF_val)/1000000,
    d2a = value / assets,
    consist_comciq = value/compsdebt)
  
    comciq <- comciq %>%  select(-c("compsdebt", "yearfounded", "datadate", "sic"))
  
  # Cleaning and trimming unlikely observations ----------------------------------------------------------------
  truncate_p <- function(x, p) {
    
    lower_limit <- quantile(x, probs = p, na.rm = T)
    upper_limit <- quantile(x, probs = 1-p , na.rm = T)
    
    x[x < lower_limit] <- NA
    x[x > upper_limit] <- NA
    
    return(x)
  }
  
  
# Only accepting observations where compustat and capiq reports similar data
sum(0.7 < comciq$consist_comciq & 1.3 > comciq$consist_comciq, na.rm = T) / sum(!is.na(comciq$consist_comciq))
# filter if the match quality - compustat and capiq - is poor
comciq <- comciq %>% filter(0.7 < consist_comciq & 1.3 > consist_comciq,
                            companyname!= "Berkshire Hathaway Inc.",
                            sic_group != ("Nonclassifiable"))

# filtering unlikely observations 
comciq <- comciq %>% 
  mutate(dateper = (year - 2010)*4 + quarter,
         levrg = ifelse(levrg < 0 | levrg > 1, NA, levrg),
         pledge = ifelse(pledge < 0 | pledge > 1, NA, pledge),
         intcov = ifelse(intcov == -Inf | intcov > 500, NA, intcov),
         liq = ifelse(liq < 0, NA, liq),
         inv_rate = ifelse(inv_rate > 20 | inv_rate < - 20, NA, inv_rate),
         inv = ifelse(inv_rate > 20 | inv_rate < -20, NA, inv),
         age = ifelse(age < 0, NA, age)) 

# Further variables for the multivariate analysis
comciq <- comciq %>% group_by(companyid) %>% arrange(companyid, dateper) %>% 
  mutate(ebit2a = ebitda / assets ,
         d2c = value / collat,
         d2c = ifelse(d2c < 0 | d2c > 2.5, NA, d2c),
         AB2a =  AB_val / 1000000 / assets,
         CF2a = CF_val / 1000000 / assets)


comciq <- comciq %>% ungroup() %>%
  mutate(ebit2a = truncate_p(ebit2a, 0.01),
#          d2c = truncate_p(d2c, 0.1),
         d2a = truncate_p(d2a, 0.02),
         prod = truncate_p(prod, 0.02),
         AB2a = truncate_p(AB2a, 0.02),
         CF2a = truncate_p(CF2a, 0.02))


# dropping Canada 
comciq <- comciq %>% filter(country == "United States") %>% select(-c("country"))


# I leave the multivariate analysis for STATA
write_csv(comciq, "comciq_stata.csv")
# haven::write_dta(comciq, "comciq_stata.dta", version = 13)


# Adding interest rate data -----------------------------------------------
flevel_intrate <- read_csv("flevel_intrate.csv") %>% select(-c("gvkey"))
comciq <- read_csv("comciq_stata.csv") %>% select(-c("gvkey"))

comciq <- comciq %>% left_join(flevel_intrate, by=c("companyid", "dateper"))

# haven::write_dta(comciq, "flevel_full.dta", version = 13)
# write_csv(comciq, "flevel_full.csv")

# duplications
comciq <- comciq %>% group_by(companyid, dateper) %>%
  mutate(dup = ifelse(row_number() > 1, row_number(), 0)) 
s_tab(dup, comciq)

# Define classification using dplyr and case_when
comciq <- comciq %>% ungroup() %>% 
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
s_tab(SME, comciq)
comciq$assets[comciq$SME == 0]  %>%  summary()


# These are the firms eligible for SBRA
comciq <- comciq %>%  mutate(SBRA = ifelse(liab < 10, 1,0))
s_tab(SBRA, comciq)
# the value of assets for these firms
comciq$assets[comciq$SBRA == 1]  %>%  summary()


# classifying firms by debt-value
summary(comciq$value)

comciq <- comciq %>%
  mutate(category = case_when(
    value <= 0.05 ~ "A",           # 50,000 / 1,000,000
    value <= 0.1 ~ "B",            # 100,000 / 1,000,000
    value <= 0.5 ~ "C",            # 500,000 / 1,000,000
    value <= 1 ~ "D",              # 1,000,000 / 1,000,000
    value <= 10 ~ "E",             # 10,000,000 / 1,000,000
    value <= 50 ~ "F",             # 50,000,000 / 1,000,000
    value <= 100 ~ "G",            # 100,000,000 / 1,000,000
    value <= 500 ~ "H",            # 500,000,000 / 1,000,000
    value <= 1000 ~ "I",           # 1,000,000,000 / 1,000,000
    value <= 5000 ~ "J",           # 5,000,000,000 / 1,000,000
    value <= 10000 ~ "K",          # 10,000,000,000 / 1,000,000
    value <= 50000 ~ "L",          # 50,000,000,000 / 1,000,000
    TRUE ~ NA_character_           # Assign NA for values outside the defined intervals
  ))


s_tab(category, comciq)

# liquidation probability you should calibrate towards
pr_liq <- c(0.665, 0.810, 0.749, 0.606, 0.441, 0.310, 0.188, 0.112, 0.0643, 0.0120, 0.00917, 0.00748)
percent <- c(0.015, 0.01, 0.042, 0.03, 0.17, 0.149, 0.057, 0.186, 0.096, 0.171, 0.037, 0.033)
# Weighted average calculation
weighted_avg_pr_liq <- sum(pr_liq * percent) / sum(percent)


# Tables -----------------------------------------------------------
# # share of CF held by the top x% of companies 
# binnum = 100
# comciq <- comciq %>%  mutate(bin_assets = statar::xtile(assets, binnum))
# aux <- filter(comciq,!is.na(bin_assets)) %>%  group_by(bin_assets) %>%  summarize(sumCFL = sum(CF_val)) %>% arrange(bin_assets)
# print(paste0("Share of the top ", 100/binnum , " percent: ", aux$sumCFL[binnum] / sum(aux$sumCFL) * 100 ))

# table 1 - summary statistics 
summary_stats <- comciq %>% ungroup() %>% filter(year >= 2010) %>% 
  select(c("assets", "revenue", "emp", "eq", "inv_rate",  "prod", "age", "value", "d2c", "d2a", "levrg", "pledge", "liq", "matur", "intrate", "spread", "CFshare")) %>%
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


# table 2 - targets
summary_stats <- comciq %>% ungroup() %>% filter(SBRA == 0 & year >= 2010) %>% 
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

# Share of ABL and CFL firms
comciq <- comciq %>%  group_by(companyid) %>% 
  mutate(CFcat = ifelse(CFshare == 0, 0, 
                        ifelse(CFshare > 0 & CFshare < 1, 1, 2)))
s_tab(CFcat, comciq)[["percent"]]

# table 2 - the share of firms who have multiple debt contracts
comciq %>% group_by(CFcat) %>%
  summarise(numshare = n(), 
            totass = sum(assets, na.rm = T),
            totloan = sum(value, na.rm = T)) %>% 
  mutate(numshare = numshare / sum(numshare) * 100,
         totass = totass / sum(totass) * 100,
         totloan = totloan / sum(totloan) * 100) %>%  
  xtable::xtable( caption = "Summary Statistics") %>%
  print(type = "latex", include.rownames = TRUE, booktabs = TRUE,
        sanitize.text.function = function(x) x)


dlevel <- read_csv("dlevel_imputed.csv")
# number of firms that hold only a single debt contract
s_tab(Ncontr, dlevel %>% group_by(companyid, year, quarter) %>% summarise(Ncontr = n()))

# share of CFbased debt
CFshare <- comciq %>% group_by(SBRA) %>%  filter(year >= 2010) %>% 
  mutate(ABdebt = (1-CFshare)*value , 
         CFdebt = CFshare*value) %>% 
  summarise(CFdebt = sum(CFdebt), ABdebt = sum(ABdebt)) %>% 
  mutate(CFshare = CFdebt / (CFdebt + ABdebt))
CFshare

# CF-share and pledgeability ----------------------------------------------

# Define the intervals
intervals <- c(0, 0.2, 0.4, 0.6, 0.8)

# Define a function to run regression for each interval
run_regression <- function(lower) {

data <- comciq %>% mutate(lage = log10(age),
                  lemp = log10(emp),
                  lebitda = sign(ebitda) * log(abs(ebitda) + 1)) %>% 
    filter(CFshare >= lower, CFshare < lower + 0.2,
           !is.na(spread) & !is.infinite(spread),
           !is.na(lebitda) & !is.infinite(lebitda),
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
  labs(color = "Transition:") +
  theme(legend.text = element_text(size = 5), legend.title = element_text(size = 7),
        legend.position = "bottom", legend.box = "horizontal") +
  xlab("CFL reliance") +
  ylab("Coefficient Estimate") + 
  scale_color_brewer(palette = "Set1") + my_theme + 
  #scale_color_manual(values = wes_palette("Darjeeling2"))
  ggtitle("Coefficient Estimates: Assets and Pledgeability") + 
  theme(plot.title = element_text(size = 10), legend.title = element_blank(),
        legend.text = element_text(size = 8), legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"))


cowplot::plot_grid(p1, p2, ncol = 2)
ggsave("mixplot.png", dpi = 400, height = 6, width = 10)
# it would be nice too see these things over time as well


