# COMPUSTAT THAT CAPIQ PANEL! DATA - 2010- 2024
# I RUN INTO MEMORY ISSUES, SO A LOT OF COMPUTATION IS SOLVED WITH FOR LOOPS

library(tidyverse)
library(tictoc)
library(RPostgreSQL)
library(RPostgres)
library(tictoc)

rm(list = ls())

setwd("C:/Users/szjud/OneDrive/Asztali gép/EBCs/CFL-git/Data Replication")
source("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V1//s_tab.R")

# 1: CapitalIQ cleaning and duplication control ---------------------------------------
exrate <- read_csv("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V4//exrate.csv")

# ONLY READ THIS IF IT IS THE 2ND TIME RUNNING THE CODE
# matchf <- read_csv("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V4//connected_firms.csv") 

start_year <- 2000
end_year <- 2024
year_increment <- 2

for (i in seq(start_year, end_year, year_increment)) {
  start_date <- i
  end_date <- i + year_increment - 1
  
  comps <- read_csv( paste0("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V4//raw_dat//", "comps_", start_date, "_", end_date, ".csv"))
  capiq <- read_csv( paste0("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V4//raw_dat//", "capiq_", start_date, "_", end_date, ".csv"))
  
  # ONLY READ THIS IF IT IS THE 2ND TIME RUNNING THE CODE
  # Keeps firms that can be matched to compustat later - this cuts comp time significantly, but only works if you have already matched the firms to compustat
  # capiq <- filter(capiq, companyid %in% unique(matchf$value)) 
  
  # 1) Multiplying up with unittype
  capiq <- capiq %>%
    mutate(dataitemvalue = case_when(
      unittypeid == 0 ~ dataitemvalue * 1,
      unittypeid == 1 ~ dataitemvalue * 1000,
      unittypeid == 2 ~ dataitemvalue * 1000000)) 
  
  # 2) Bringing everything to USD
  capiq <- capiq %>%
    mutate(year = year(periodenddate),
           quarter = quarter(periodenddate)) %>% 
    rename(currencyid = issuedcurrencyid) %>% 
    left_join(exrate,
              by = c("year", "quarter", "currencyid")) %>% 
    mutate(avg_priceclose = ifelse(currencyid == 30, NA, avg_priceclose), # this is just a bad observation
           dataitemvalue = dataitemvalue / avg_priceclose) 
  
  # 3) Renaming
  capiq <- capiq %>%
    rename(id = componentid,
           date = periodenddate,
           description = capitalstructuredescription,
           subtypeid = capitalstructuresubtypeid,
           debttext = descriptiontext,
           value = dataitemvalue,
           currency = currencyid,
           maturity = maturityhigh,
           intbenchmark = interestratebenchmarktypeid,
           intrate = interestratehighvalue,
           secured = securedtypeid) %>% 
    select(-c("isocountry3", "countryid", "currency", "avg_priceclose"))
  
  # 4) Duplications 
  
  # By all variables - need to break up the sample otherwise i run into memory constraints
  capiq_aux <- capiq %>% mutate(dupall = 1) %>% filter(FALSE)
  for (i in unique(capiq$year)) {
    for (j in unique(capiq$quarter)) {
      
      capiq_aux2 <- filter(capiq, year == i, quarter == j)
      capiq_aux2$dupall <- select(capiq_aux2, -c("filingdate")) %>% duplicated %>% as.numeric()
      
      capiq_aux <- rbind(capiq_aux, capiq_aux2)
    }
  }
  
  capiq <- capiq_aux
  rm(list = c("capiq_aux", "capiq_aux2"))
  
  # 4a) For these observations, the only difference is the filingdate, 
  capiq <- filter(capiq, dupall == 0) %>% select(-c("dupall"))
  
  # 4b) facilities are always registered twice
  capiq <- capiq %>%
    filter(!(grepl("facility", debttext, ignore.case = TRUE)))
  
  # There are still a few duplication left due to subsidiaries, but I need Compustat data to root these out
  
  # 4c) handling subsidiaries
  
  # add assets from compustat, to be able to drop subsidiearies and find outliers
  capiq <- left_join(capiq, comps %>% mutate(year = year(datadate), quarter = quarter(datadate)) %>% select(c("atq", "gvkey", "year", "quarter")), by = c("gvkey", "year", "quarter"))
  
  # Managing duplications: 2nd round -- gvkeys sometimes change for the same firm
  # this is due to subsidiaries: drop duplicates with the lowest assets keep highest consolidation level
  
  capiq <- capiq %>% 
    group_by(companyid, year, quarter) %>% 
    mutate(aux = length(unique(gvkey))) %>% 
    arrange(companyid, id, desc(atq)) %>%  
    group_by(id, companyid, year, quarter) %>%
    mutate(dup2 = ifelse(aux == 2, row_number(), 0)) %>%
    ungroup() %>% filter(dup2 != 2)
  
  # dropping remaing dulicates these are usually just the same debt filed twice in a quarter
  capiq <- capiq %>% group_by(id, year, quarter) %>%
    mutate(dup = ifelse(row_number() > 1, row_number(), 0)) %>% filter(dup < 1)
  
  write_csv(capiq, paste0("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V4//raw_dat//", "capiq2_", start_date, "_", end_date, ".csv"))
  
}

# merging the raw datasets for capiq into one dataset - same for compustat
capiq0 <- read_csv("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V4//raw_dat//capiq2_2010_2011.csv")  %>% filter(FALSE)
comps0 <- read_csv("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V4//raw_dat//comps_2010_2011.csv")  %>% filter(FALSE)
for (i in seq(start_year, end_year, year_increment)) {
  start_date <- i
  end_date <- i + year_increment - 1
  
  capiq <- read_csv(paste0("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V4//raw_dat//", "capiq2_", start_date, "_", end_date, ".csv"))
  comps <- read_csv( paste0("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V4//raw_dat//", "comps_", start_date, "_", end_date, ".csv"))
  
  
  comps0 <- rbind(comps0, comps)
  capiq0 <- rbind(capiq0, capiq)
}
capiq <- capiq0
comps <- comps0  %>% rename(quarter = qtr)
dlevel <- capiq
capiq %>% group_by(year, quarter) %>%  summarise(aaa= sum(value, na.rm = T)) %>% print(n = 105)

# COMPUSTAT: some firms report in Canadian dollars, bring all monetary values to USD
monetary_vars <- c("ppegtq", "sppey", "dpq", "dlttq", "dlcq", "dvpq", "cdvcy", "pdvcy", "prstkcy", "cheq", "saley", "actq", "atq", "lctq", "ltq", "oibdpq", "xintq", "invtq", "rectq", "ppentq", "ibq", "cogsq", "xsgaq", "dvy", "dpactq", "chq", "capxy", "revtq", "dm")

cad_rate <- exrate %>%
  filter(currencyid == 27) %>%
  select(year, quarter, cad_rate = avg_priceclose)

comps <- comps %>%
  filter(curcdq %in% c("USD", "CAD")) %>%
  left_join(cad_rate, by = c("year", "quarter")) %>%
  mutate(across(all_of(monetary_vars), ~ if_else(curcdq == "CAD", .x / cad_rate, .x))) %>%
  select(-cad_rate)

saveRDS(comps, "comps.RDS")
saveRDS(capiq, "capiq.RDS")

rm(list = setdiff(ls(), c("comps", "dlevel", "s_tab")))


# 2: CapitalIQ, classification of Debt contracts ---------------------------------------
dlevel <- readRDS("capiq.RDS")

# First round of classification
dlevel <- dlevel %>%
  mutate(mortgage_flag = ifelse(grepl("Mortgage", debttext), 1, 0),
         mortgage_flag = ifelse(!(grepl("Mortgage", debttext)), 0, mortgage_flag),
         mortgage_flag = ifelse(grepl("mortgage", tolower(description)), 1, mortgage_flag),
         
         caplease_flag = ifelse(grepl("Capital leases", debttext), 1, 0),
         caplease_flag = ifelse(!(grepl("Capital leases", debttext)), 0, caplease_flag),
         caplease_flag = ifelse(grepl("capital", tolower(description)) & grepl("lease", tolower(description)), 1, caplease_flag),
         
         bond_flag = ifelse(grepl("Bonds and notes", debttext), 1, 0),
         bond_flag = ifelse(!(grepl("Bonds and notes", debttext)), 0, bond_flag),
         bond_flag = ifelse(grepl("bond", tolower(description)) & !grepl("mortgage", tolower(description)), 1, bond_flag),
         
         debent_flag = ifelse(grepl("Debentures", debttext), 1, 0),
         debent_flag = ifelse(!(grepl("Debentures", debttext)), 0, debent_flag),
         debent_flag = ifelse(grepl("debenture", tolower(description)), 1, debent_flag),
         
         secured = ifelse(secured == 4, 2, secured),
         secured = ifelse(secured > 4, 3, secured),
         unsec_flag = ifelse(secured == 3, 1, 0),
         unsec_flag = ifelse(secured != 3, 0, unsec_flag),
         
         ABkey = 0,
         ABkey = ifelse(tolower(description) %in% c("borrowing base", "real estate", "plant", "property", "collateral", "building", "machine"), 1, ABkey),
         ABkey = ifelse(tolower(description) %in% c("inventory", "working capital", "automobile", "vehicle", "aircraft", "asset-based", "asset based"), 1, ABkey),
         ABkey = ifelse(tolower(description) %in% c("sba", "reserve-based", "project finance", "construction", "finance company"), 1, ABkey),
         CFkey = 0,
         CFkey = ifelse(tolower(description) %in% c("substantially all", "lien", "term facility", "term loan", "syndicated", "tranche"), 1, CFkey),
         CFkey = ifelse(tolower(description) %in% c("acquisition line", "bridge loan", "senior", "subordinated", "notes"), 1, CFkey))

# Classification based on flags and assets types
# Asset-based loans
dlevel <- dlevel %>%
  mutate(AB = ifelse(mortgage_flag == 1 | caplease_flag == 1 | ABkey == 1, 1, 0)) %>%
  mutate(AB = ifelse(AB != 1, 0, AB))

# Cash-flow based loans
dlevel <- dlevel %>%
  mutate(CF = ifelse((bond_flag == 1 | debent_flag == 1 | unsec_flag == 1 | CFkey == 1) & mortgage_flag != 1, 1, 0)) %>%
  mutate(CF = ifelse(CF != 1, 0, CF))

# # Second round of classification
dlevel <- dlevel %>%
  mutate(CF = ifelse(subtypeid == 1, 1, CF),                      # Commercial Papers
         AB = ifelse(subtypeid == 2 & secured != 3, 1, AB),       # Revolving credit
         AB = ifelse(subtypeid == 3 & secured != 3, 1, AB),       # Term loans
         CF = ifelse(subtypeid == 4 & mortgage_flag != 1, 1, CF), # Bonds and Notes
         AB = ifelse(subtypeid == 7 & secured != 3, 1, AB))       # Other borrowings

# Look at overlaps
table(dlevel$CF, dlevel$AB)

# Fix overlaps - If both, classify as asset based (conservative classification as in Lian and Ma, 2021)
dlevel <- dlevel %>%
  mutate(CF = ifelse(CF == 1 & AB == 1, 0, CF))

# dropping debttypes I am not interested in (less than 0.1% off all contracts)
debt_types <- c("Trust Preferred Securities", "Bank overdraft", "Bills payable", "FEDERAL RESERVE BANK BORROWINGS", "FHLB borrowings",
                "Federal Funds Purchased", "Letter of Credit Outstanding", "Securities loaned", "Securities sold under agreement to repurchase")
dlevel <- dlevel %>% filter(!debttext %in% debt_types )

# dropping classification flags
dlevel <- dlevel %>%
  select(-c("mortgage_flag", "caplease_flag", "debent_flag", "unsec_flag", "ABkey", "bond_flag", "CFkey", "aux", "dup", "dup2", "filingdate", "zipcode", "region"))


# Finding inflows and outflows
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
         inflow = ifelse(row_number() == 1  & firstobs_firm != 1 , 1 ,0),
         outflow = ifelse(row_number() == max(rn) & lasttobs_firm != 1, 1,0)) %>% 
  select(-c("rn")) %>% ungroup()

# 3: Dlevel, Calculating Credit Spreads and Outliers ---------------------------------------

# Summary stats do not make much sense without the pre-filtering
# # The share of CF based credit over the entire sample
# sum_dlevel <- function(dat){
#   
#   CFsize <- dlevel %>% mutate(ones = 1) %>% group_by(CF) %>%
#     summarise(value = sum(value, na.rm = T),
#               number = sum(ones)) 
#   print(paste0("The share of credit identitfied as CF - by value: ",  CFsize[2,2]/(CFsize[1,2]+CFsize[2,2]))) %>% return()
#   print(paste0("The share of credit identitfied as CF - by number: ",  CFsize[2,3]/(CFsize[1,3]+CFsize[2,3]))) %>% return()
#   
#   CFsize <- dlevel %>% mutate(AB_val = value * AB,
#                               CF_val = value * CF) %>% group_by(year, quarter) %>% 
#     summarise(AB_val = sum(AB_val, na.rm = T),
#               CF_val = sum(CF_val, na.rm = T),
#               nobs = max(row_number())) %>% 
#     mutate(CFshare =  CF_val / (CF_val + AB_val)) %>% print(n = 100) %>% return()
# }
# 
# 
# sum_dlevel(dlevel) %>% print(n = 104)

# Longitudinal checks
dlevel <- dlevel %>%
  filter(!(value == 0 | is.na(value))) %>%  # dropping zero values and NAs, about 18% of observations 
  arrange(id, year, quarter) %>%
  mutate(aux_value = (value - lag(value)) / lag(value), .by = id) %>%
  filter(aux_value < 2 | is.na(aux_value)) %>%  #  Change in value too large, drops around 2.3% of observations
  select(-aux_value)

# Unrealistic Leverage Values
dlevel <- dlevel %>%
  filter(country %in% c("United States", "Canada")) %>%
  mutate(d2a = value / atq / 1000000,
         value = if_else(d2a > 3, NA_real_, value, missing = value))


# calculating maturity
dlevel <- dlevel %>% 
  mutate(dateper = (year - 2000)*4 + quarter,
         matuper = (year(maturity) - 2000)*4 + quarter(maturity)) %>% group_by(id) %>% 
  mutate(matumax = max(matuper-dateper),
         matumax = ifelse(matumax < 4, NA, matumax),
         matumax = ifelse(matumax < 0, 0, matumax)) %>% ungroup() %>% filter(year >= 2010) %>% select(-c("dateper", "matuper"))


# importing benchmark rates - new data so you need some renaming
benchmark <- read_csv("C:/Users/szjud/OneDrive/Asztali gép/EBCs/Data/R/V4/raw_dat/benchmark.csv") %>% filter(year(date) >= 2010) %>%
  mutate(year = year(date), quarter = quarter(date)) %>% select(-c("date")) %>% group_by(year, quarter) %>% 
  summarise_all(mean, na.rm = TRUE) %>% rename_with(toupper, -c(year, quarter))


bonds_and_notes <- c("Bonds and notes", "Commercial paper", "Debentures", "Mortgage bonds", "Mortgage notes", "Notes payable")
bank_loans <- c("Bank loans", "Capital leases", "General borrowings", "Mortgage loans", "Other borrowings", "Revolving credit",  "Term loan")

benchmark_long <- benchmark %>%
  pivot_longer(-c(year, quarter), names_to = "benchmarktype", values_to = "benchmark")

dlevel <- dlevel %>%
  ungroup() %>%
  mutate(benchmarktype = case_when(
    debttext %in% bank_loans & matumax == 0 ~ "TB4WK",
    debttext %in% bank_loans & matumax == 1 ~ "TB3MS",
    debttext %in% bank_loans & matumax == 2 ~ "TB6MS",
    debttext %in% bank_loans & matumax == 3 ~ "TB1YR",
    debttext %in% bank_loans & matumax <= 6 ~ "GS1",
    debttext %in% bank_loans & matumax <= 12 ~ "GS2",
    debttext %in% bank_loans & matumax <= 32 ~ "GS5",
    debttext %in% bank_loans & matumax > 32 ~ "GS10",
    debttext %in% bonds_and_notes & matumax <= 6 ~ "GS1",
    debttext %in% bonds_and_notes & matumax <= 12 ~ "GS2",
    debttext %in% bonds_and_notes & matumax <= 32 ~ "GS5",
    debttext %in% bonds_and_notes & matumax <= 60 ~ "GS10",
    debttext %in% bonds_and_notes & matumax <= 100 ~ "GS20",
    debttext %in% bonds_and_notes & matumax > 100 ~ "GS30"
  )) %>%
  left_join(benchmark_long, by = c("year", "quarter", "benchmarktype")) %>%
  select(-intbenchmark)

# dropping unrealistic values
dlevel <- dlevel %>% mutate(spread = intrate - benchmark) %>%
  mutate(spread = ifelse(spread > 15 | spread < 0, NA, spread),
         intrate = ifelse(intrate > 30 | intrate < 0, NA, intrate))

# getting rid of a few stray variables 
dlevel <- dlevel %>%
  select(-c("benchmarktype", "benchmark", 
            "firstobs_firm", "lasttobs_firm"))

saveRDS(dlevel, "dlevel_cleaned.RDS")


# 4: Flevel, collapsing contract level data ------------------------

# A: debt stocks by type
dlevel <- readRDS("dlevel_cleaned.RDS") %>% filter(year >= 2010)

flevel <- dlevel %>%
  group_by(companyid, companyname, gvkey, yearfounded, country, CF, year, quarter) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = CF, values_from = value, values_fill = 0) %>%
  rename(AB_val = `0`, CF_val = `1`) %>%
  mutate(CFshare = if_else((CF_val + AB_val) > 0, CF_val / (CF_val + AB_val), NA_real_),
         dateper = (year - 2010) * 4 + quarter)

# B: weighted rates, spreads, maturities and new issuance

# drop capital leases - accounting changes at the end of 2018 distort the data
dlevel_misc <- dlevel %>%
  filter(debttext != "Capital leases")

flevel_intrate <- dlevel_misc %>%
  mutate(w = coalesce(value, 0),
         is_cf = CF %in% 1,
         is_ab = CF %in% 0,
         is_new = inflow %in% 1) %>%
  mutate( dateper = (year - 2010) * 4 + quarter) %>% 
  group_by(companyid, gvkey, dateper) %>%
  summarise(maturity = weighted.mean(matumax, w, na.rm = TRUE),
            intrate = weighted.mean(intrate, w, na.rm = TRUE),
            spread = weighted.mean(spread, w, na.rm = TRUE),
            CF_inflow = sum(w[is_cf & is_new]),
            AB_inflow = sum(w[is_ab & is_new]),
            CFshare_inflow = if_else(CF_inflow + AB_inflow > 0, CF_inflow / (CF_inflow + AB_inflow), NA_real_),
            spread_CF = weighted.mean(spread[is_cf], w[is_cf], na.rm = TRUE),
            spread_AB = weighted.mean(spread[is_ab], w[is_ab], na.rm = TRUE),
            spread_CF_inflow = weighted.mean(spread[is_cf & is_new], w[is_cf & is_new], na.rm = TRUE),
            spread_AB_inflow = weighted.mean(spread[is_ab & is_new], w[is_ab & is_new], na.rm = TRUE),
            .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ replace(.x, is.nan(.x), NA)))

# Join: one row per firm-quarter
flevel <- flevel %>%
  left_join(flevel_intrate, by = c("companyid", "gvkey", "dateper"))

rm(flevel_intrate, dlevel_misc)

saveRDS(flevel, "flevel.RDS")

# 5: Compustat, Cleaning and connecting on to CapIq ---------------------------------------
comps <- readRDS("comps.RDS")

# Keep non-financial companies only
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

# Create derived variables
comps <- comps %>%  
  mutate(atq = if_else(atq == 0, NaN , atq), 
         tot_debt = dlcq + dlttq,
         net_debt = tot_debt - chq,
         levrg = tot_debt/atq,
         collat = ppentq + invtq + rectq,
         pledge = collat/atq,
         intcov = oibdpq /xintq,
         eq = atq - ltq,
         liq = chq / atq,
         prod = ((revtq - cogsq)/1000) / ((ppegtq/1000)^(1/3) * emp^(2/3))) %>%  
  ungroup() %>% arrange(gvkey, datadate) %>% 
  group_by(gvkey, year) %>%  
  mutate(sppeq = ifelse(quarter == 1, sppey, sppey - lag(sppey)),
         capxq = ifelse(quarter == 1, capxy, capxy - lag(capxy)),
         capxq = ifelse(capxq < 0, NA, capxq)) %>%
  group_by(gvkey) %>% 
  mutate(inv = capxq - sppeq,
         inv_rate =  100*(capxq - sppeq)/ lag(atq),
         lassets = log10(atq*1000000)) 
# the correct definition would be using 'ppegtq' instead of assets but that has to many missing values 

# keep only the relevant variables
comps <- comps %>% rename(ebitda = oibdpq, assets = atq, liab = lctq, sec_debt = dm, rating = spcsrc, revenue = revtq, capit = ppentq) %>%  
  select(c("datadate", "revenue","sic","rating","year","quarter","emp","sic_group", "ebitda", "sec_debt",
           "tot_debt","net_debt","levrg", "capit", "collat","pledge","intcov","eq",  "liq", "assets", "liab", "inv", "inv_rate", "prod", "lassets"))

saveRDS(comps, "comps_fjc.RDS")


# 6: Flevel Full, Connecting Compustat to Flevel -------------------------------------------------------

# joining flevel and comps
flevel <- readRDS("flevel.RDS") %>% left_join(comps, by = c("gvkey", "year", "quarter"))


# !!! dropping all observations that are not in CAPIQ !!!
flevel <-flevel %>% filter(!is.na(datadate)) 
# save connecting firms so that you do not have to do the whole download and connecting 
write.csv(as_tibble(unique(flevel$companyid)), "connected_firms.csv")


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


# I leave the multivariate analysis for STATA
saveRDS(flevel, "flevel_full.RDS")

# 7: Dlevel full, Connecting Joined data to Dlevel ----------------------------------------
# I connect to dlevel here to carry over the data cleaning done on the flevel

# dropping extra variables
dlevel <- dlevel %>%
  select(-c("companyname", "gvkey", "country", "unittypeid", "atq", "d2a", "yearfounded"))

# keeping only the firm-level variables needed for the Stata regressions
flevel <- flevel %>%
  select(companyid, year, quarter, dateper, sic_group, rating,
         emp, age, revenue, ebitda, ebit2a, assets, lassets, pledge, levrg,
         CFshare, tot_debt = value)

# joining firm level
dlevel <- left_join(dlevel, flevel, by = c("companyid", "year", "quarter"))
dlevel <- dlevel %>% filter(!is.na(sic_group))

# contract-level rates: average over the contract's life, spreads filled only where missing
dlevel <- dlevel %>%
  group_by(id) %>%
  arrange(year, quarter, .by_group = TRUE) %>%
  fill(intrate, spread, .direction = "downup") %>%
  ungroup()


rm(list = setdiff(ls(), c("flevel", "dlevel", "s_tab")))

haven::write_dta(dlevel, "dlevel_full.dta", version = 13)
saveRDS(dlevel, "dlevel_full.RDS")


