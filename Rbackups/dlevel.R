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
                  axis.title.x=element_text(size = 10), axis.title.y=element_text(size = 10),
                  panel.background = element_rect(fill = "white"),
                  panel.grid.major = element_line(size = 0.2, linetype = 'solid', colour = "grey80"),
                  panel.grid.minor = element_line(size = 0.15, linetype = 'solid', colour = "grey90"),
                  panel.border = element_rect(color = "grey50", fill = NA, size = 0.5))

# data import  ---------------------------------------------------------
dlevel <- read_csv("dlevel_imputed.csv") 
comciq <- read_csv("comciq_stata.csv")

# Fixing inflow and outflows - 
# if a firm appears first in capiq it should not be counted
comp_id <- dlevel %>% group_by(companyid, year, quarter) %>%  summarise() %>%
  group_by(companyid) %>% 
  mutate(rn = row_number(),
         firstobs_firm = ifelse(rn == 1 , 1 , 0),
         lasttobs_firm = ifelse(rn == max(rn), 1 , 0))


dlevel <- dlevel %>%  left_join(comp_id, by = (c("companyid", "year", "quarter")))

dlevel <- dlevel %>% group_by(id) %>% arrange(companyid, id, year, quarter) %>% 
  mutate(rn = row_number(),
         yq = str_c(as.character(year), "Q", as.character(quarter)),
         inflow = ifelse(row_number() == 1 & yq != "2010Q1" & firstobs_firm != 1 , 1 ,0),
         outflow = ifelse(row_number() == max(rn) & yq != "2023Q3" & lasttobs_firm != 1, 1,0)) %>% 
  select(-c("rn"))

s_tab(inflow, dlevel)
s_tab(outflow, dlevel)

# finding maturity
dlevel <- dlevel %>% 
  mutate(dateper = (year - 2000)*4 + quarter,
         matuper = (year(maturity) - 2000)*4 + quarter(maturity)) %>% group_by(id) %>% 
  mutate(matumax = max(matuper-dateper),
         matumax = ifelse(matumax < 4, NA, matumax),
         matumax = ifelse(matumax < 0, 0, matumax)) %>% ungroup() %>% filter(year >= 2010) %>% select(-c("dateper", "matuper"))

# importing benchmark rates
benchmark <- left_join(read_csv("C:/Users/szjud/OneDrive/Asztali gép/EBCs/Data/R/V3/raw_dat/treasuries.csv"),
                       read_csv("C:/Users/szjud/OneDrive/Asztali gép/EBCs/Data/R/V3/raw_dat/tbills.csv"), by = c("DATE")) %>% rename("date" = "DATE") %>% 
  mutate(year = year(date), quarter = quarter(date)) %>% select(-c("date")) %>% group_by(year, quarter) %>% 
  summarise_all(mean, na.rm = TRUE)

bonds_and_notes <- c("Bonds and notes", "Commercial paper", "Debentures", "Mortgage bonds", "Mortgage notes", "Notes payable")
bank_loans <- c("Bank loans", "Capital leases", "General borrowings", "Mortgage loans", "Other borrowings", "Revolving credit",  "Term loan")

# defining the benchmark for each debt contract
dlevel <- dlevel %>% ungroup() %>% 
  left_join(benchmark, by = c("year", "quarter")) %>%
  mutate(benchmark = NA, benchmarktype = NA) %>%
  mutate(benchmark = ifelse(debttext %in% bank_loans & matumax == 0, TB4WK, benchmark),
         benchmark = ifelse(debttext %in% bank_loans & matumax == 1, TB3MS, benchmark),
         benchmark = ifelse(debttext %in% bank_loans & matumax == 2, TB6MS, benchmark),
         benchmark = ifelse(debttext %in% bank_loans & matumax == 3, TB1YR, benchmark),
         benchmark = ifelse(debttext %in% bank_loans & matumax >= 4 & matumax <= 6, GS1, benchmark),
         benchmark = ifelse(debttext %in% bank_loans & matumax >= 7 & matumax <= 12, GS2, benchmark),
         benchmark = ifelse(debttext %in% bank_loans & matumax >= 13 & matumax <= 32, GS5, benchmark),
         benchmark = ifelse(debttext %in% bank_loans & matumax >= 33, GS10, benchmark),
         benchmark = ifelse(debttext %in% bonds_and_notes & matumax <= 6, GS1, benchmark),
         benchmark = ifelse(debttext %in% bonds_and_notes & matumax >= 7 & matumax <= 12, GS2, benchmark),
         benchmark = ifelse(debttext %in% bonds_and_notes & matumax >= 13 & matumax <= 32, GS5, benchmark),
         benchmark = ifelse(debttext %in% bonds_and_notes & matumax >= 33 & matumax <= 60, GS10, benchmark),
         benchmark = ifelse(debttext %in% bonds_and_notes & matumax >= 61 & matumax <= 100, GS20, benchmark),
         benchmark = ifelse(debttext %in% bonds_and_notes & matumax >= 101, GS30, benchmark),
         
         benchmarktype = ifelse(debttext %in% bank_loans & matumax == 0, "TB4WK", benchmarktype),
         benchmarktype = ifelse(debttext %in% bank_loans & matumax == 1, "TB3MS", benchmarktype),
         benchmarktype = ifelse(debttext %in% bank_loans & matumax == 2, "TB6MS", benchmarktype),
         benchmarktype = ifelse(debttext %in% bank_loans & matumax == 3, "TB1YR", benchmarktype),
         benchmarktype = ifelse(debttext %in% bank_loans & matumax >= 4 & matumax <= 6, "GS1", benchmarktype),
         benchmarktype = ifelse(debttext %in% bank_loans & matumax >= 7 & matumax <= 12, "GS2", benchmarktype),
         benchmarktype = ifelse(debttext %in% bank_loans & matumax >= 13 & matumax <= 32, "GS5", benchmarktype),
         benchmarktype = ifelse(debttext %in% bank_loans & matumax >= 33, "GS10", benchmarktype),
         benchmarktype = ifelse(debttext %in% bonds_and_notes & matumax <= 6, "GS1", benchmarktype),
         benchmarktype = ifelse(debttext %in% bonds_and_notes & matumax >= 7 & matumax <= 12, "GS2", benchmarktype),
         benchmarktype = ifelse(debttext %in% bonds_and_notes & matumax >= 13 & matumax <= 32, "GS5", benchmarktype),
         benchmarktype = ifelse(debttext %in% bonds_and_notes & matumax >= 33 & matumax <= 60, "GS10", benchmarktype),
         benchmarktype = ifelse(debttext %in% bonds_and_notes & matumax >= 61 & matumax <= 100, "GS20", benchmarktype),
         benchmarktype = ifelse(debttext %in% bonds_and_notes & matumax >= 101, "GS30", benchmarktype)) %>%
  select(-c("GS10", "GS1", "GS2", "GS5", "GS30", "GS20", "TB3MS", "TB1YR", "TB4WK", "TB6MS", "intbenchmark"))

# dropping outliers
dlevel <- dlevel %>% mutate(spread = intrate - benchmark) %>% mutate(spread = ifelse(spread > 15 | spread < 0, NA, spread),
                                                                     intrate = ifelse(intrate > 30 | intrate < 0, NA, intrate))

# dropping Canada
dlevel <- dlevel %>% filter(country == "United States")

# dropping extra variables
dlevel <- dlevel %>% select(-c("companyname", "gvkey", "country", "unittypeid", "maturity", "atq", "d2a", "imputed", "yearfounded"))
comciq <- comciq %>% select(-c("AB_val", "CF_val", "sec_debt", "net_debt", "sec_share", "consist_comciq", "d2a")) %>% rename("tot_debt" = "value")

# joining firm level 
dlevel <- left_join(dlevel, comciq, by = c("companyid", "year", "quarter")) 
dlevel <- dlevel %>% filter(!is.na(sic_group))


dlevel <- dlevel %>% group_by(id) %>%
  mutate(intrate = mean(intrate, na.rm = T),
         spread = ifelse(is.na(spread),  mean(spread, na.rm = T) , spread))

# picking up IDB bankruptcy data
IDB <- read_csv("sum_master.csv") %>%  select("EASST", "pr_liq") %>% filter(EASST != "$1b + (pre 2015)") %>% 
  mutate(
    min_ass = c(0, 1e9, 1e6, 1e10, 1e7, 1e8, 1e5, 5e7, 5e8, 5e5, 5e4),
    max_ass = c(5e4, 1e10, 1e7, Inf, 5e7, 5e8, 5e5, 1e8, 1e9, 1e6, 1e5))

dlevel <- dlevel %>% mutate(assets = assets * 1000000) %>%
  rowwise() %>%
  mutate(
    EASST = IDB$EASST[which(assets >= IDB$min_ass & assets < IDB$max_ass)]
  ) %>%
  ungroup()

s_tab(EASST, dlevel)


# # how often do interest rates change within debt contract? -- NEVER
dlevel <- dlevel %>% group_by(id) %>%
  mutate(intmean = mean(spread, na.rm = T),
         intmax = max(spread, na.rm = T),
         neq = ifelse(spread != intmean, 1,0))
s_tab(neq, dlevel)

dlevel <- dlevel %>% left_join(IDB %>% select(c("EASST", "pr_liq")), by = "EASST")

write_csv(dlevel, "dlevel_full.csv")
haven::write_dta(dlevel, "dlevel_full.dta", version = 13)

# Firm-level interest rates  ----------------------------------------------
dlevel <- read_csv("dlevel_full.csv") %>% arrange(companyid, dateper) 

flevel_intrate <- dlevel  %>% group_by(companyid, gvkey, dateper) %>% 
  mutate(dweight = value / sum(value)) %>% 
  summarise(matur = sum(dweight*matumax, na.rm = T),
            intrate = sum(dweight*intrate, na.rm = T),
            spread = sum(dweight*spread, na.rm = T),
            pr_liq = mean(pr_liq))

write_csv(flevel_intrate, "flevel_intrate.csv")

ggplot(data = dlevel %>% filter(inflow == 1)) + 
  geom_point(aes(x = lassets, y = spread), color = "pink", alpha = 0.3) +
  geom_smooth(aes(x = lassets, y = spread),  method = "lm", fill = "pink", color = "firebrick", se = TRUE) +
  ggtitle("") + 
  ylab("Credit Spread") + xlab("Log of Assets") + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + my_theme


# Summary statistics -----------------------------------------------------
# number of debt contracts
unique(dlevel$id) %>% length()
# number of firms 
unique(dlevel$companyid) %>% length()

# table 2 - summary statistics 
summary_stats <- dlevel %>% ungroup() %>% mutate(value = value / 10^6) %>% 
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


# Multivariate analysis ---------------------------------------------------
dlevel <- read_csv("dlevel_full.csv") %>% 
  mutate(lage = log10(age),
         lemp = log10(emp))

dlevel <- dlevel %>% filter(inflow == 1)

dlevel <- dlevel %>% 
  group_by(companyid, dateper) %>%
  mutate(dup = ifelse(row_number() > 1, row_number(), 0)) 

dlevel <- dlevel %>%
  mutate(sic_group = factor(sic_group),
         rating = factor(rating),
         dateper = factor(dateper),
         companyid = factor(companyid))

lm(spread ~ ebit2a + pledge + lassets  + levrg + age + emp + 
     + sic_group + rating + dateper, data = dlevel) %>% summary()

lm(spread ~ ebit2a + pledge + lassets + levrg + age + emp + 
               + sic_group + rating + dateper, data = filter(dlevel, CF == 1)) %>% summary()

lm(spread ~ ebit2a + pledge + lassets + levrg + age + emp + 
               + sic_group + rating + dateper, data = filter(dlevel, CF == 0)) %>% summary


# # FIRM LSDV REGRESSIONS - they are a bit weird 
# lm(spread ~ ebit2a + pledge + lassets + levrg + age + emp + CF +
#      companyid, data = dlevel) %>% summary()
# 
# lm(spread ~ ebit2a + pledge + lassets + levrg + age + emp +
#      companyid, data = filter(dlevel, CF == 1)) %>% summary()
# 
# lm(spread ~ ebit2a + pledge + lassets + levrg + age + emp +
#      companyid, data = filter(dlevel, CF == 0)) %>% summary



#  2-variate analysis - conditional means ---------------------------------
dlevel <- read.csv("C:/Users/szjud/OneDrive/Asztali gép/EBCs/Data/R/V3//dlevel_full.csv")

ggplot(data = dlevel %>% filter(inflow == 1)) + 
  geom_point(aes(x = lassets, y = spread), color = "pink", alpha = 0.3) +
  geom_smooth(aes(x = lassets, y = spread),  method = "lm", fill = "pink", color = "firebrick", se = TRUE) +
  ggtitle("") + 
  ylab("Credit Spread") + xlab("Log of Assets") + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + my_theme


ggplot(data = dlevel %>% filter(inflow == 1 & AB == 1)) + 
  geom_point(aes(x = lassets, y = spread), color = "pink", alpha = 0.4) +
  geom_smooth(aes(x = lassets, y = spread), fill = "pink", color = "firebrick", se = TRUE) +
  ggtitle("") + 
  ylab("Credit Spread") + xlab("Log of Assets") + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + my_theme


ggplot(data = dlevel %>% filter(inflow == 1 & AB == 0)) + 
  geom_point(aes(x = lassets, y = spread), color = "lightblue", alpha = 0.4) +
  geom_smooth(aes(x = lassets, y = spread), method = "lm", fill = "lightblue", color = "blue4", se = TRUE) +
  ggtitle("") + 
  ylab("Credit Spread") + xlab("Log of Assets") + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + my_theme
