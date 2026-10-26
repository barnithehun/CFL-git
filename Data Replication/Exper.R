library(tidyverse)
library(tictoc)
library(RPostgreSQL)
library(RPostgres)
library(tictoc)
library(lubridate)
library(did)

rm(list = ls())

setwd("C:/Users/szjud/OneDrive/Asztali gép/EBCs/Data/R/V4")
source("C:/Users/szjud/OneDrive/Asztali gép/EBCs/Data/R/V1/s_tab.R")

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

# Problem #1 
# Inflows and outflows are always going to be tricky
# A lot of cases, for example with revolving credit lines, CapitalIQ does not log a new debt ID when the firm accesses the new credit, yet it is an inflow of debt
# The inflow outflow analysis, especially when it comes to credit spreads of new contracts is still valid, but when you want to analyse the volume of inflows this gets into a bigger problem
# FOR NOW: I am disregarding this problem 

# Problem #2
# For the Calloway and Santa Anna type of DID estimator you need a balanced panel
# This is especially problematic when you have quarterly data because almost no firms get CF inflows every quarter
# Therefore only a few variables work, such as CFshare
# With yearly data this is less of a problem, but this bias is still there

# ---------------------------------------
# Quarterly data
flevel <- read_csv("flevel_full.csv") %>%  mutate(treated  = ifelse(value <= 7.5, 1, 0))  %>% 
  filter(dateper %in% seq(1,41) | dateper %in% seq(42,62)) # %>% filter(value <= 2 |  value >= 10)
s_tab(treated, flevel)

# this is the Santa'Anna treatment timing variable
flevel <- flevel %>%
  mutate(G = ifelse(treated == 1, 41, Inf))  # 41 is the first treatment period


att <- att_gt(
  yname = "CFshare",         # your outcome variable name
  tname = "dateper",         # time variable
  idname = "companyid",      # firm identifier
  gname = "G",               # treatment timing variable
  data = flevel,
  est_method = "dr",         # doubly robust (can change to "ipw" or "reg
  panel = TRUE)              # panel data

ggdid(att)


# total share of CF based inflows in the economy per quarter
sumdat <- flevel %>% 
  mutate(dateper = year) %>% 
  group_by(dateper) %>%
  # mutate(size_q = ifelse(value < 3*10^6, 1,2)) %>% 
  mutate(size_q = case_when(value < 7.5 ~ 1, value >= 7.5 ~ 2, TRUE ~ NA_integer_ )) %>% 
  group_by(dateper, size_q) %>%
  summarise(CFshare_inflow = mean(CFshare_inflow, na.rm = T)) # %>% filter(!is.na(size_q))


ggplot(data = sumdat) + 
  geom_line(aes(x = dateper, y = CFshare_inflow, color = as.factor(size_q)), linewidth = 1, alpha = 0.5) +
  ggtitle("") + 
  ylab("Share of Debt") + xlab("Year") + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  my_theme


flevel <- flevel %>% group_by(dateper) %>%  mutate(size_q = ntile(assets, 2)) %>%  ungroup()
ggplot(data = flevel) + 
  geom_point(aes(x = dateper, y = CFshare_inflow, colour = as.factor(size_q)), linewidth = 1, alpha = 0.5) +
  geom_smooth(aes(x = dateper, y = CFshare_inflow, colour = as.factor(size_q), method = "gam"), linewidth = 1, alpha = 0.5) +
  ggtitle("") + 
  ylab("Share of Debt") + xlab("Year") + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  my_theme


# -----------------------------------------------------------------------------------------
# Yearly data
flevel <- read_csv("flevel_full.csv") 

# collapse further to yearly aggregates acrss firms
yr_flevel <- flevel %>% filter(dateper <= 60) %>%  group_by(year) %>% 
  mutate(size_q = ntile(assets, 3)) %>%
  mutate(size_q = case_when(value < 7.5 ~ 1, value >= 7.5 & value < 50 ~ 2 , value >= 50 ~ 3, TRUE ~ NA_integer_ ))  %>% group_by(year, size_q) %>%
  summarise(AB_val = sum(AB_val, na.rm = T), 
            CF_val = sum(CF_val, na.rm = T),
            CF_inflow = sum(CF_inflow, na.rm = T),
            AB_inflow = sum(AB_inflow, na.rm = T),
            spread_CF        = mean(spread_CF, na.rm = TRUE),
            spread_AB        = mean(spread_AB, na.rm = TRUE),
            spread_CF_inflow = mean(spread_CF_inflow, na.rm = TRUE),
            spread_AB_inflow = mean(spread_AB_inflow, na.rm = TRUE))  %>% 
  
  mutate(CFshare = if_else((CF_val + AB_val) > 0, CF_val / (CF_val + AB_val), NA_real_),
         CFshare_inflow = if_else((CF_inflow + AB_inflow) > 0, CF_inflow / (CF_inflow + AB_inflow), NA_real_)) %>% ungroup() %>% 
  select(-c("CF_val", "AB_val")) # these are only good for calculating shares


# Plotting stuff
ggplot(data = yr_flevel) + 
  geom_line(aes(x = year, y = CFshare, color = as.factor(size_q)), linewidth = 1, alpha = 0.5) +
  ggtitle("") + 
  ylab("Share of Debt") + xlab("Year") + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  my_theme

ggplot(data = yr_flevel) + 
  geom_line(aes(x = year, y = CFshare_inflow, color = as.factor(size_q)), linewidth = 1, alpha = 0.5) +
  ggtitle("") + 
  ylab("Share of Debt") + xlab("Year") + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  my_theme



# -----------------------------------------------------------------------------------------
# Yearly data - firm level
flevel <- read_csv("flevel_full.csv") 

yr_flevel <- flevel %>% filter(dateper != 63) %>%  group_by(companyid, year) %>% 
  summarise(value  = max(value, na.rm = T),
            AB_val = sum(AB_val, na.rm = T), 
            CF_val = sum(CF_val, na.rm = T),
            CF_inflow = sum(CF_inflow, na.rm = T),
            AB_inflow = sum(AB_inflow, na.rm = T),
            spread_CF        = mean(spread_CF, na.rm = TRUE),
            spread_AB        = mean(spread_AB, na.rm = TRUE),
            spread_CF_inflow = mean(spread_CF_inflow, na.rm = TRUE),
            spread_AB_inflow = mean(spread_AB_inflow, na.rm = TRUE),
            
            assets = mean(assets, na.rm = TRUE),
            emp = mean(emp, na.rm = TRUE),
            ebitda = sum(ebitda, na.rm = TRUE),
            age = mean(age, na.rm = TRUE),
            levrg = mean(levrg, na.rm = T)) %>% 
  
  mutate(CFshare = if_else((CF_val + AB_val) > 0, CF_val / (CF_val + AB_val), NA_real_),
         CFshare_inflow = if_else((CF_inflow + AB_inflow) > 0, CF_inflow / (CF_inflow + AB_inflow), NA_real_)) %>% ungroup() %>% 
  select(-c("CF_val", "AB_val")) %>%  # these are only good for calculating shares
  
  mutate(lasset    = log10(assets),
         lemp     = log10(emp),
         lage     = log10(age),
         lebitda  = sign(ebitda) * log10(abs(ebitda) + 1))  %>%  group_by(year) %>% 
  # mutate(size_q = case_when(value < 7.5 ~ 1, value >= 7.5 ~ 2, TRUE ~ NA_integer_ )) 
  mutate(size_q = ntile(assets, 2)) 


# Plotting stuff
ggplot(data = yr_flevel) + 
  geom_smooth(aes(x = year, y = CFshare, color = as.factor(size_q)), linewidth = 1, alpha = 0.5) +
  ggtitle("") + 
  ylab("Share of Debt") + xlab("Year") + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  my_theme

ggplot(data = yr_flevel) + 
  geom_smooth(aes(x = year, y = CFshare_inflow, color = as.factor(size_q)), linewidth = 1, alpha = 0.5) +
  ggtitle("") + 
  ylab("Share of Debt") + xlab("Year") + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  my_theme


ggplot(data = yr_flevel) + 
  geom_smooth(aes(x = year, y = spread_CF_inflow, color = as.factor(size_q)), linewidth = 1, alpha = 0.5) +
  ggtitle("") + 
  ylab("Share of Debt") + xlab("Year") + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  my_theme


# Santa Anna DID
estimdat <- yr_flevel %>%  filter(year %in% seq(2012,2020) | year %in% seq(2021,2024)) %>%  filter(value < 1000000) %>% 
  mutate(treated  = ifelse(value <= 5, 1, 0), 
         G = ifelse(treated == 1, 2020, Inf))

s_tab(treated, estimdat)

att <- att_gt(
  yname = "CFshare_inflow",         # your outcome variable name
  tname = "year",         # time variable
  idname = "companyid",      # firm identifier
  gname = "G",               # treatment timing variable
  data = estimdat,
  est_method = "dr",
  panel = TRUE)

p1 <- ggdid(att) + my_theme + ggtitle("CF-share of New Inflows")
 

att <- att_gt(
  yname = "CFshare",         # your outcome variable name
  tname = "year",         # time variable
  idname = "companyid",      # firm identifier
  gname = "G",               # treatment timing variable
  data = estimdat, 
  est_method = "dr",
  panel = TRUE)

p2 <- ggdid(att) + my_theme + ggtitle("CF-share of Outstanding Debt")

cowplot::plot_grid(p1, p2, ncol = 2)
ggsave("C://Users//szjud//OneDrive//Asztali gép//EBCs//CFL-git//Latex codes//Plots//DID.png", dpi = 400, height = 5, width = 9)


# TWFE
library(fixest)

flevel_twfe <- yr_flevel %>% filter(year %in% seq(2010, 2025), value < 11.5 & value > 4) %>% mutate(treated = ifelse(value <= 7.5, 1, 0))

# Run TWFE regression
# Example with outcome spread_CF_inflow and firm and year fixed effects
twfe_model <- feols(
  CFshare_inflow ~ treated * (year >= 2020) + lasset + lebitda + levrg | as.factor(companyid) + as.factor(year),
  data = flevel_twfe)

summary(twfe_model)


# RDD ---------------------------------------------------------------------
# firm level
flevel <- read_csv("flevel_full.csv") %>% filter(dateper %in% seq(49,58)) %>% 
  mutate(
    lasset    = log10(assets),
    lemp     = log10(emp),
    lage     = log10(age),
    lvalue   = log10(value),
    lebitda  = sign(ebitda) * log10(abs(ebitda) + 1))

# Run sharp RDD
RDestimate(spread_CF_inflow ~ value, bw = 3.5 , data = flevel, cutpoint = 7.5) %>% summary()
RDestimate(spread_CF_inflow ~ value | lasset + levrg + lebitda , bw = 3.5 , data = flevel, cutpoint = 7.5) %>% summary()


# debt level 
dlevel <- read_csv("dlevel_full.csv") %>% filter(dateper %in% seq(53,58)) %>% 
  mutate(
    lasset    = log10(assets),
    lemp     = log10(emp),
    lage     = log10(age),
    lvalue   = log10(value),
    lebitda  = sign(ebitda) * log10(abs(ebitda) + 1))

# Run sharp RDD
RDestimate(spread ~ value | lebitda +lasset+ pledge+ levrg, bw = 3.5*10^6, data = dlevel %>%  filter(CF == 1 & inflow == 1), cutpoint = 7.5*10^6) %>% summary()
RDestimate(spread ~ value, bw = 3.5*10^6, data = dlevel %>%  filter(CF == 0 & inflow == 1), cutpoint = 7.5*10^6) %>% summary()


# fuzzyRDD is not feasibe in this setup



#  Predicting default probability -----------------------------------------
flevel <- read_csv("flevel_full.csv")
IDB <- readRDS("C://Users//szjud//OneDrive//Asztali gép//EBCs//CFL-git//Data Replication//idb_cases.RDS")   # from IDB_easst.R
linking <- read_csv("FJC_link.csv")
cik <- read_csv("cik.csv") %>% group_by(gvkey, cik) %>%  summarise()


test <- left_join(flevel, cik, by = "gvkey") %>% 
  left_join(linking, by = "gvkey") %>% 
  left_join(IDB, by = c("CASEKEY"))

