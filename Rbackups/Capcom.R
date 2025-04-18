# COMPUSTAT THAT CAPIQ PANEL! DATA - 2010- 2022
# I RUN INTO MEMORY ISSUES, SO A LOT OF COMPUTATION IS SOLVED WITH FOR LOOPS

library(tidyverse)
library(tictoc)
library(RPostgreSQL)
library(RPostgres)
library(tictoc)

rm(list = ls())

setwd("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V3")
source("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V1//s_tab.R")

# CapitalIQ query -----------------------------------------------------------------

# OPENING THE WRDS CONNECTION
wrds <- dbConnect(Postgres(),
                  host='wrds-pgdata.wharton.upenn.edu',
                  port=9737,
                  dbname='wrds',
                  sslmode='require',
                  user='szbarnabas',
                  password = 'Fsztznenezze410!!')


res <- dbSendQuery(wrds, "select distinct table_schema
                   from information_schema.tables
                   where table_type ='VIEW'
                   or table_type ='FOREIGN TABLE'
                   order by table_schema")
data <- dbFetch(res, n=-1)
dbClearResult(res)
data

tic()

start_year <- 2000
end_year <- 2022
year_increment <- 2


res <- dbSendQuery(wrds, "SELECT currencyid,
                          EXTRACT(YEAR FROM pricedate) AS year,
                          EXTRACT(QUARTER FROM pricedate) AS quarter,
                          AVG(priceclose) AS avg_priceclose
                          FROM ciq.ciqexchangerate
                          GROUP BY currencyid, year, quarter
                          HAVING EXTRACT(YEAR FROM pricedate) > 2000")
exrate <- dbFetch(res, n = -1)
dbClearResult(res)

write_csv(exrate, "exrate.csv")

capiq_from_to <- function(start_date, end_date) {
  query <- paste("SELECT  a.companyname, a.companyid, a.componentid, a.periodenddate, a.filingdate, a.capitalstructuredescription,  
                                  a.descriptiontext, a.dataitemvalue, a.unittypeid, a.issuedcurrencyid, a.maturityhigh, a.interestratehighvalue, 
                                  a.interestratebenchmarktypeid, a.securedtypeid, a.gvkey, a.capitalstructuresubtypeid,
                                  
                                  b.countryid, b.yearfounded, b.zipcode,
                                  
                                  c.country,  c.isocountry3, c.region
                                  
                   FROM ciq.wrds_debt a LEFT JOIN ciq.ciqcompany b
                   ON a.companyid = b.companyid
                   
                   LEFT JOIN ciq.ciqcountrygeo c
                   ON b.countryid = c.countryid
                   
                   WHERE periodenddate BETWEEN '", start_date, "' AND '", end_date, "'", sep = "")
  
  res <- dbSendQuery(wrds, query)
  capiq <- dbFetch(res, n = -1)
  dbClearResult(res)
  
file_name <- paste0("C:/Users/szjud/OneDrive/Asztali gép/EBCs/Data/R/V3/raw_dat/", "capiq_", substr(start_date, 1, 4), "_", substr(end_date, 1, 4), ".csv")
write.csv(capiq, file = file_name, row.names = FALSE, fileEncoding = "UTF-8")

}


tic()
for (i in seq(start_year, end_year, year_increment)) {
  start_date <- paste0(i, "-01-01")
  end_date <- paste0(i + year_increment - 1, "-12-31")
  capiq_from_to(start_date, end_date)  # Call your function here
  print(i)
}
toc()


# Compustat Query --------------------------------------------------------------
comps_from_to <- function(start_date, end_date) {
  
  query <- paste("SELECT a.conm, a.datadate, a.gvkey,  a.curcdq, a.ppegtq, a.sppey, a.dpq, a.dlttq, a.dlcq, a.dvpq, a.cdvcy, a.pdvcy, a.prstkcy, 
                 a.cheq, a.saley, a.prccq, a.cshoq, a.actq, a.atq, a.lctq, a.ltq, a.oibdpq, a.xintq, a.invtq, a.rectq, 
                 a.ppentq, a.ibq, a.cogsq, a.xsgaq, a.dvy, a.dpactq, a.chq, a.capxy, a.revtq, a.indfmt,
                 
                 b.sic, b.spcsrc, b.loc
                 
                 FROM comp.fundq a LEFT JOIN comp.company b
                 ON a.gvkey = b.gvkey
                 
                 WHERE datadate BETWEEN '", start_date, "' AND '", end_date, "'", sep = "")
# quarterly data
res <- dbSendQuery(wrds, query)
compsq <- dbFetch(res, n = -1)
dbClearResult(res)

query2 <- paste("SELECT gvkey, datadate, emp, dm
                                 
                 FROM comp.co_afnd1 
                          
                 WHERE datadate BETWEEN '", start_date, "' AND '", end_date, "'", sep = "")

res <- dbSendQuery(wrds, query2)
compsy <- dbFetch(res, n = -1)
dbClearResult(res)

compsy <- compsy %>%  filter(!is.na(emp) & !is.na(dm)) %>%  mutate(year = year(datadate)) %>% select(-c("datadate"))
comps <-  compsq %>%  mutate(year = year(datadate), qtr = quarter(datadate)) %>% left_join(compsy,  by = c("gvkey", "year"))

file_name <- paste0("C:/Users/szjud/OneDrive/Asztali gép//EBCs/Data/R/V3//raw_dat/", "comps_", substr(start_date, 1, 4), "_", substr(end_date, 1, 4), ".csv")
write.csv(comps, file = file_name, row.names = FALSE, fileEncoding = "UTF-8")

}

tic()
for (i in seq(start_year, end_year, year_increment)) {
  start_date <- paste0(i, "-01-01")
  end_date <- paste0(i + year_increment - 1, "-12-31")
  comps_from_to(start_date, end_date)  # Call your function here
  print(i)
}
toc()

# Capital IQ cleaning and duplication control ---------------------------------------

exrate <- read_csv("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V3//raw_dat//exrate.csv")
matchf <- read_csv("connected_firms.csv") # you can only read this in if it is the second time running the whole code

start_year <- 2000
end_year <- 2022
year_increment <- 2

 for (i in seq(start_year, end_year, year_increment)) {
     start_date <- i
     end_date <- i + year_increment - 1
  
     comps <- read_csv( paste0("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V3//raw_dat//", "comps_", start_date, "_", end_date, ".csv"))
     capiq <- read_csv( paste0("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V3//raw_dat//", "capiq_", start_date, "_", end_date, ".csv"))
   
  # Only keeping firms that can be matched to compustat later - this is to avoid memory issues   
  capiq <- filter(capiq, companyid %in% unique(matchf$value))  # you can only read this in if it is the second time running the whole code
     
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
  
  # ADD ASSETS TO CONNECT
    capiq <- left_join(capiq, comps %>% mutate(year = year(datadate), quarter = quarter(datadate)) %>% select(c("atq", "gvkey", "year", "quarter")), by = c("gvkey", "year", "quarter"))
  
  # Managing duplications: 2nd round -- gvkeys sometimes change for the same firm
  # this is due to subsidiaries: drop duplicates with the lowest assets keep highest consolidation level
  
  capiq <- capiq %>% 
    group_by(companyid, year, quarter) %>% 
    mutate(aux = length(unique(gvkey))) %>% 
    arrange(companyid, id, atq) %>%  
    group_by(id, companyid, year, quarter) %>%
    mutate(dup2 = ifelse(aux == 2, row_number(), 0)) %>%
    ungroup() %>% filter(dup2 != 2)
  
  # dropping remainging dulicates these are usually just the same debt filed twice in a quarter
  capiq <- capiq %>% group_by(id, year, quarter) %>%
    mutate(dup = ifelse(row_number() > 1, row_number(), 0)) %>% filter(dup < 1)
  
  write_csv(capiq, paste0("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V3//raw_dat//", "capiq2_", start_date, "_", end_date, ".csv"))

}

# merging the raw datasets for capiq into one dataset - same for compustat
capiq0 <- read_csv("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V3//raw_dat//capiq2_2010_2011.csv")  %>% filter(FALSE)
comps0 <- read_csv("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V3//raw_dat//comps_2010_2011.csv")  %>% filter(FALSE)
for (i in seq(start_year, end_year, year_increment)) {
  start_date <- i
  end_date <- i + year_increment - 1
  
  capiq <- read_csv(paste0("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V3//raw_dat//", "capiq2_", start_date, "_", end_date, ".csv"))
  comps <- read_csv( paste0("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V3//raw_dat//", "comps_", start_date, "_", end_date, ".csv"))
  
  
  comps0 <- rbind(comps0, comps)
  capiq0 <- rbind(capiq0, capiq)
}
capiq <- capiq0
comps <- comps0  %>% rename(quarter = qtr)
dlevel <- capiq


# COMPUSTAT: some firms are reported in canadian dollars --- bringing all compustat values to USD
comps <- left_join(comps, filter(exrate, currencyid == 27), by = c("year", "quarter"))
monetary_vars <- setdiff(comps0 %>% select_if(is.numeric) %>% names(), c("cshoq", "emp", "year", "qtr", "currencyid"))  
comps0 <- comps %>% filter(curcdq == "CAD") 
comps0[monetary_vars] <- comps0[monetary_vars] / comps0$avg_priceclose 
comps <- rbind(comps %>% filter(curcdq == "USD"), comps0) 

write_csv(comps, "comps.csv")
write_csv(capiq, "capiq.csv")

rm(list = setdiff(ls(), c("comps", "dlevel", "s_tab")))

# Classifying ---------------------------------------
dlevel <- read_csv("capiq.csv")

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


# Asset-based loans
dlevel <- dlevel %>%
  mutate(AB = ifelse(mortgage_flag == 1 | caplease_flag == 1 | ABkey == 1, 1, 0)) %>%
  mutate(AB = ifelse(AB != 1, 0, AB))

# Cash-flow based loans
dlevel <- dlevel %>%
  mutate(CF = ifelse((bond_flag == 1 | debent_flag == 1 | unsec_flag == 1 | CFkey == 1) & mortgage_flag != 1, 1, 0)) %>%
  mutate(CF = ifelse(CF != 1, 0, CF))

# Second round of classification
# Commercial Papers
dlevel <- dlevel %>%
  mutate(CF = ifelse(subtypeid == 1, 1, CF),                      # Commercial Papers
         AB = ifelse(subtypeid == 2 & secured != 3, 1, AB),       # Revolving credit
         AB = ifelse(subtypeid == 3 & secured != 3, 1, AB),       # Term loans
         CF = ifelse(subtypeid == 4 & mortgage_flag != 1, 1, CF), # Bonds and Notes
         AB = ifelse(subtypeid == 7 & secured != 3, 1, AB))       # Other borrowings

# Look at overlaps
table(dlevel$CF, dlevel$AB)

# Fix overlaps - If both 
dlevel <- dlevel %>%
  mutate(AB = ifelse(CF == 1 & AB == 1, 0, AB))

# dropping debttypes I am not interested in (these are less than 0.1% off all contracts)
debt_types <- c("Trust Preferred Securities", "Bank overdraft", "Bills payable", "FEDERAL RESERVE BANK BORROWINGS", "FHLB borrowings",
                "Federal Funds Purchased", "Letter of Credit Outstanding", "Securities loaned", "Securities sold under agreement to repurchase")
dlevel <- dlevel %>% filter(!debttext %in% debt_types )


# dropping classification flags and useless vars
dlevel <- dlevel %>%
  select(-c("mortgage_flag", "caplease_flag", "debent_flag", "unsec_flag", "ABkey", "bond_flag", "CFkey", "aux", "dup", "dup2", "filingdate", "subtypeid", "zipcode", "region"))

write_csv(dlevel, "dlevel.csv")

# CapitalIQ Outlier control ---------------------------------------------------------
dlevel <- read_csv("dlevel.csv")
comps <- read_csv("comps.csv")

# The share of CF based credit over the entire sample
sum_dlevel <- function(dat){
  
  CFsize <- dlevel %>% mutate(ones = 1) %>% group_by(CF) %>%
    summarise(value = sum(value, na.rm = T),
              number = sum(ones)) 
  print(paste0("The share of credit identitfied as CF - by value: ",  CFsize[2,2]/(CFsize[1,2]+CFsize[2,2]))) %>% return()
  print(paste0("The share of credit identitfied as CF - by number: ",  CFsize[2,3]/(CFsize[1,3]+CFsize[2,3]))) %>% return()
  
  CFsize <- dlevel %>% mutate(AB_val = value * AB,
                              CF_val = value * CF) %>% group_by(year, quarter) %>% 
    summarise(AB_val = sum(AB_val, na.rm = T),
              CF_val = sum(CF_val, na.rm = T),
              nobs = max(row_number())) %>% 
    mutate(CFshare =  CF_val / (CF_val + AB_val)) %>% print(n = 100) %>% return()
}
sum_dlevel(dlevel)

# Longitudinal checks ---
# Debt sizes that change a lot are prly incorrect
dlevel <- dlevel %>% group_by(id) %>%
  arrange(companyid, id, year, quarter) %>%
  filter(!(value == 0 | is.na(value))) %>%
  mutate(aux_value = (value - lag(value)) / lag(value)) %>%
  filter(aux_value < 2 | is.na(aux_value)) %>% # drops around 30000 obs
  select(-aux_value)
# not making much difference.. not even sure if it is needed

# for checking new inflows and outflows - this is a bit flawed, cover the end periods
dlevel <- dlevel %>% group_by(id) %>% arrange(id, year, quarter) %>% 
  mutate(inflow = ifelse(row_number() == 1, 1 ,0),
         outflow = ifelse(row_number() == max(row_number()), 1,0))

# 1) Outlier control - GET BACK TO THIS ONCE YOU HAVE FULL COMPUSTAT DATA
dlevel <- dlevel %>%  filter(country == "United States" | country == "Canada")  %>% 
  mutate(d2a = value / atq / 1000000,
             value = if_else(d2a > 3, NA, value))

# Linear interpolation  -------------------------------------------------
dlevel <- dlevel %>% group_by(id) %>% arrange(id, year, quarter) %>% 
  mutate(dateper = (year - 2010)*4 + quarter,    # need a consistent date variable
         dategap = dateper - lag(dateper), dategap  = ifelse(is.na(dategap), 0, dategap),
         aux = lead(dategap), aux  = ifelse(is.na(aux), 0, aux))

# linear interpolation - need to fill up missing rows before I could do that!
tic()
new_rows <- dlevel %>% filter(FALSE)
 for (i in 1:nrow(dlevel)) {
   if (dlevel[[i, "aux"]] > 1 && dlevel[[i, "aux"]] < 9) {
     for (j in 1:(dlevel[[i, "aux"]] - 1)) {

       new_rows <- rbind(new_rows, dlevel[i, ] )
       last <- nrow(new_rows)

       new_rows[last, "dateper"] <- new_rows[last, "dateper"] + j

       new_rows[last, "value"] <- dlevel[[i, "value"]]  + ((dlevel[[i + 1, "value"]] - dlevel[[i, "value"]]) / dlevel[[i, "aux"]])*j
       new_rows[last, "intrate"] <- dlevel[[i, "intrate"]]  + ((dlevel[[i + 1, "intrate"]] - dlevel[[i, "intrate"]]) / dlevel[[i, "aux"]])*j
       new_rows[last, "intbenchmark"] <- dlevel[[i, "intbenchmark"]]  + ((dlevel[[i + 1, "intbenchmark"]] - dlevel[[i, "intbenchmark"]]) / dlevel[[i, "aux"]])*j

     }
   }
   
   if (i %% 10000 == 0) {
     cat("Iteration:", i/nrow(dlevel)*100, "\n")
   }
}
toc()

dlevel <- dlevel %>% mutate(imputed = 0) 
new_rows <- new_rows %>%  mutate(year = 2010 + floor((dateper - 1) / 4),
                                 quarter = (dateper - 1) %% 4 + 1,
                                 imputed = 1) 
dlevel <- rbind(dlevel, new_rows) %>% arrange(id, dateper)  %>% select(-c("aux", "dategap"))


# ONLY RUN IF you want to check non-legacy debt ---------------------------

# 
# dlevel <- read_csv("dlevel_imputed.csv") # CHECK DUPLICATIONS HERE
# dlevel <- dlevel %>% group_by(id) %>% mutate(obs_per = row_number(),
#                                              dateper_firstobs = if_else(obs_per == 1, dateper, NA_real_),
#                                              dateper_firstobs = max(dateper_firstobs, na.rm = T),
#                                              loanage = dateper - dateper_firstobs + 1) %>%  ungroup()
# dlevel <- dlevel %>%  filter(year != 2010, loanage <= 4)


# write_csv(dlevel, "dlevel_imputed.csv") 

# Creating firm level data ------------------------------------------------
dlevel <- read_csv("dlevel_imputed.csv") %>% filter(year >= 2010)

# make a "new outflows of ABL and CFL chart!!!!"
flevel <- dlevel %>% group_by(companyid, companyname, gvkey, yearfounded, country, CF, year, quarter) %>%  # strings, same within firms
  summarise(value = sum(value, na.rm = T)) %>% spread(CF, value) %>% 
  rename("AB_val" = "0",
         "CF_val" = "1") %>%  
  mutate(AB_val = if_else(AB_val == 0, NaN , AB_val),
         CF_val = if_else(CF_val == 0, NaN , CF_val),
         CFshare = CF_val/(AB_val+CF_val),
         CFshare = if_else((is.na(AB_val) & !is.na(CF_val)), 1 , CFshare),
         CFshare = if_else((is.na(CF_val) & !is.na(AB_val)), 0 , CFshare),
         CF_val = replace_na(CF_val, 0),
         AB_val = replace_na(AB_val, 0)) 

# # If you ever want to pick up maturity and intrate - they seem to be cursed somehow
# flevel_aux <- dlevel %>% group_by(companyid, companyname, gvkey, yearfounded, country, year, quarter) %>%
#   summarise(intrate = mean(intrate, na.rm = T),
#             maturity = mean(maturity, na.rm = T))
# flevel <- flevel %>% left_join(select(flevel_aux, c("companyid", "year", "quarter", "intrate", "maturity")), by = c("companyid", "year", "quarter"))
# rm(list = flevel_aux)

write_csv(flevel, "flevel.csv") 

# # the dynamics of new debt - CF or AB
# dlevel <- read_csv("dlevel_imputed.csv") %>% # filter(year != 2010) %>% 
#   mutate(CFinflow = CF*inflow*value,
#          ABinflow = AB*inflow*value) %>% group_by(year) %>% 
#   summarise(AB_val = sum(ABinflow, na.rm = T),
#             CF_val = sum(CFinflow, na.rm = T))
# 
# ggplot(dlevel, aes(x = year)) + geom_line(aes(y = CF_val), color = "blue") + geom_line(aes(y = AB_val), color = "red")
# # not very interesting overall



# Liq vs Reo --------------------------------------------------------------
read.csv()


