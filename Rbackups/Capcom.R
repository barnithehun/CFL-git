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
                  password = 'Fsztznenezze410')


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
  
  file_name <- paste0("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V3//raw_dat//", "capiq_", substr(start_date, 1, 4), "_", substr(end_date, 1, 4), ".csv")
  write.csv(capiq, file = file_name, row.names = FALSE)
  
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
tic()
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

file_name <- paste0("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V3//raw_dat//", "comps_", substr(start_date, 1, 4), "_", substr(end_date, 1, 4), ".csv")
write.csv(comps, file = file_name, row.names = FALSE)

}


for (i in seq(start_year, end_year, year_increment)) {
  start_date <- paste0(i, "-01-01")
  end_date <- paste0(i + year_increment - 1, "-12-31")
  comps_from_to(start_date, end_date)  # Call your function here
  print(i)
}


# Capital IQ cleaning and duplication control ---------------------------------------

exrate <- read_csv("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V3//raw_dat//exrate.csv")
matchf <- read_csv("connected_firms.csv") # you can only read this in if it is the second time running the whole code

start_year <- 2010
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

# merging capiq and comsutat data into one-one dataset
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
monetary_vars <- setdiff(comps0 %>% select_if(is.numeric) %>% names(), c("cshoq", "emp", "year", "quarter", "qtr", "currencyid"))  
comps0 <- comps %>% filter(curcdq == "CAD") 
comps0[monetary_vars] <- comps0[monetary_vars] / comps0$avg_priceclose 
comps <- rbind(comps %>% filter(curcdq == "USD"), comps0) 

rm(list = setdiff(ls(), c("comps", "dlevel", "s_tab")))

write_csv(comps, "comps.csv")

# Classifying ---------------------------------------

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

# Fix overlaps
dlevel <- dlevel %>%
  mutate(AB = ifelse(CF == 1 & AB == 1, 0, AB))

# dropping classification flags and useless vars
dlevel <- dlevel %>%
  select(-c("mortgage_flag", "caplease_flag", "debent_flag", "unsec_flag", "ABkey", "bond_flag", "CFkey", "aux", "dup", "dup2", "secured", "filingdate", "subtypeid", "zipcode", "region"))

write_csv(dlevel, "dlevel.csv")

# CapitalIQ Outlier control ---------------------------------------------------------

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

dlevel <- read_csv("dlevel.csv")
comps <- read_csv("comps.csv")

# Longitudinal checks ---
# Debt sizes that change a lot are prly incorrect
dlevel <- dlevel %>% group_by(id) %>%
  arrange(companyid, id, year, quarter) %>%
  filter(!(value == 0 | is.na(value))) %>%
  mutate(aux_value = (value - lag(value)) / lag(value)) %>%
  filter(aux_value < 2 | is.na(aux_value)) %>% # drops around 30000 obs
  select(-aux_value)
# not making much difference.. not even sure if it is needed

# for checking new CF inflow and outflows
dlevel <- dlevel %>% group_by(id) %>% arrange(id, year, quarter) %>% 
  mutate(inflow = ifelse(row_number() == 1, 1 ,0),
         outflow = ifelse(row_number() == max(row_number()), 1,0))

# 1) Outlier control - GET BACK TO THIS ONCE YOU HAVE FULL COMPUSTAT DATA
dlevel <- dlevel %>%  filter(country == "United States" | country == "Canada")  %>% 
  mutate(d2a = value / atq / 1000000,
             value = if_else(d2a > 3, NA, value))

# Linear interpolation  -------------------------------------------------
# must be super inefficient - takes around 20 mins to run - probably there is a better way of doing this 
dlevel <- dlevel %>% group_by(id) %>% arrange(id, year, quarter) %>% 
  mutate(dateper = (year - 2010)*4 + quarter,                          # need a consistent date variable
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


dlevelimputed <- read_csv("dlevel_imputed.csv")

# dlevel  %>% select(c("companyname", "companyid", "id", "year", "quarter", "value", "intrate", "dateper", "imputed", "CF")) %>% arrange(desc(value)) %>%  print(n = 200)


# Creating firm level data ------------------------------------------------

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
# write_csv(dlevel, "dlevel_imputed.csv") 

# # the dynamics of new debt - CF or AB
# dlevel <- read_csv("dlevel_imputed.csv") %>% # filter(year != 2010) %>% 
#   mutate(CFinflow = CF*inflow*value,
#          ABinflow = AB*inflow*value) %>% group_by(year) %>% 
#   summarise(AB_val = sum(ABinflow, na.rm = T),
#             CF_val = sum(CFinflow, na.rm = T))
# 
# ggplot(dlevel, aes(x = year)) + geom_line(aes(y = CF_val), color = "blue") + geom_line(aes(y = AB_val), color = "red")
# # not very interesting overall

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
         collat = ppentq + invtq + rectq,
         pledge = (ppentq + invtq + rectq)/atq,
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
         inv_rate = (capxq - sppeq)/ lag(atq)) 
# the correct definition would be using 'ppegtq' instead of assets but that has to many missing values 

comps <- comps %>% rename(ebitda = oibdpq, assets = atq, sec_debt = dm, rating = spcsrc, revenue = revtq) %>%  
  select(c("datadate", "revenue","sic","rating","year","quarter","emp","sic_group", "ebitda", "sec_debt",
           "tot_debt","net_debt","levrg","collat","pledge","intcov","eq","liq", "assets", "inv", "inv_rate", "prod"))

# joining flevel and comps
comciq <- flevel %>% left_join(comps, by = c("gvkey", "year", "quarter"))
rm(list = c("comps", "flevel"))

# this is for earlier - so that i would not have to do all observatoins in firm_level
# write_csv(unique(comciq$companyid) %>% as_tibble(), "connected_firms.csv")

# dropping all observations that are not in CAPIQ
comciq <-comciq %>% filter(!is.na(datadate))

# "new variables that could be generated "nd round of new vars after connecting CAPIQ
comciq <- comciq %>%   rename(compsdebt = tot_debt) %>% mutate(
  emp = if_else(emp == 0, NaN , emp),
  sec_share = sec_debt/compsdebt,
  age = year - yearfounded + 1,
  value = (AB_val+CF_val)/1000000,
  d2a = value / assets,
  lassets = log10(assets*1000000),
  consist_comciq = value/compsdebt)

# Only accepting observations where compustat and capiq reports similar data
sum(0.7 < comciq$consist_comciq & 1.3 > comciq$consist_comciq, na.rm = T) / sum(!is.na(comciq$consist_comciq))
# hist(filter(comciq, consist_comciq < 3)$consist_comciq, nclass = 40)
# comciq <- comciq %>% filter(0.8 < consist_comciq & 1.2 > consist_comciq)
 comciq <- comciq %>%  select(-c("compsdebt", "yearfounded", "datadate", "gvkey", "sic"))

write_csv(comciq,"comciq.csv") 

# CLeaning and trimming unlikely observations ----------------------------------------------------------------
truncate_p <- function(x, p) {

  lower_limit <- quantile(x, probs = p, na.rm = T)
  upper_limit <- quantile(x, probs = 1-p , na.rm = T)
  
  x[x < lower_limit] <- NA
  x[x > upper_limit] <- NA
  
  return(x)
}
comciq <- read_csv("comciq.csv") 

# filter if the match quality - compustat and capiq - is poor
comciq <- comciq %>% filter(0.5 < consist_comciq & 2 > consist_comciq,
                            companyname!= "Berkshire Hathaway Inc.",
                            sic_group != ("Nonclassifiable"))

# filtering unlikely observations 
comciq <- comciq %>% 
  mutate(dateper = (year - 2010)*4 + quarter,
         levrg = ifelse(levrg < 0 | levrg > 3, NA, levrg),
         pledge = ifelse(pledge < 0 | pledge > 1, NA, pledge),
         intcov = ifelse(intcov == -Inf | intcov > 500, NA, intcov),
         liq = ifelse(liq < 0, NA, liq),
         inv_rate = ifelse(inv_rate > 0.5 | inv_rate < -0.5, NA, inv_rate),
         inv = ifelse(inv_rate > 0.5 | inv_rate < -0.5, NA, inv),
         age = ifelse(age < 0, NA, age),
         dateper = (year - 2010)*4 + quarter) 

# Further variables for the multivariate analysis
comciq <- comciq %>% group_by(companyid) %>% 
  mutate(gr_assets = (assets - lag(assets)) / lag(assets),
         gr_emp = (emp - lag(emp)) / lag(emp),
         mean_gr_assets = mean(gr_assets, na.rm = T),
         mean_gr_emp = mean(gr_emp, na.rm = T),
         profitability = ebitda / assets,
         profit_vari = sd(ebitda, na.rm = T)^2,
         profit_z = (ebitda - mean(ebitda, na.rm = T))/sd(ebitda, na.rm = T),
         d2c = (collat / assets))

comciq <- comciq %>% ungroup() %>%
   mutate(gr_assets = truncate_p(gr_assets, 0.02),
          gr_emp = truncate_p(gr_emp, 0.02),
          mean_gr_assets = truncate_p(mean_gr_assets, 0.02),
          mean_gr_emp = truncate_p(mean_gr_emp, 0.02),
          profitability = truncate_p(profitability, 0.02),
          profit_vari = truncate_p(profit_vari, 0.02),
          profit_z = truncate_p(profit_z, 0.02),
          d2c = truncate_p(d2c, 0.01),
          prod = truncate_p(prod, 0.02))
  

# I leave the multivariate anlysis for stata for STATA
# write_csv(comciq, "comciq_stata.csv")
# haven::write_dta(comciq, "comciq_stata.dta", version = 13)

# Tables -----------------------------------------------------------
comciq <- read_csv("comciq_stata.csv")

# share of CF held by the top x% of companies 
binnum = 100
comciq <- comciq %>%  mutate(bin_assets = statar::xtile(assets, binnum))
aux <- filter(comciq,!is.na(bin_assets)) %>%  group_by(bin_assets) %>%  summarize(sumCFL = sum(CF_val)) %>% arrange(bin_assets)
print(paste0("Share of the top ", 100/binnum , " percent: ", aux$sumCFL[binnum] / sum(aux$sumCFL) * 100 ))

# table 1 - summary statistics 
summary_stats <- comciq %>% 
  select(c("assets", "revenue", "emp", "eq", "inv_rate",  "prod", "age", "gr_assets", "profit_z", "value", "d2c", "levrg", "pledge", "liq" )) %>%  
  rename(invrate = inv_rate, grassets = gr_assets, profitz = profit_z) %>% 
  mutate(grassets = grassets*100, invrate  = invrate*100,   pledge  = pledge*100) %>% 
  summarize(across(everything(), list(
    mean = ~mean(., na.rm = T),
    p10 = ~quantile(., 0.1, na.rm = T),
    p25 = ~quantile(., 0.25, na.rm = T),
    median = ~median(., na.rm = T),
    p75 = ~quantile(., 0.75, na.rm = T),
    p90 = ~quantile(., 0.9, na.rm = T)
  ))) %>%
  gather(statistic, value)  %>%  separate(statistic, into = c("variable", "statistic"), sep = "_") %>%
  spread(statistic, value) %>%  
  xtable::xtable( caption = "Summary Statistics") %>%
  print(type = "latex", include.rownames = TRUE, booktabs = TRUE,
        sanitize.text.function = function(x) x)

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

dlevel <- read_csv("dlevel.csv")
# number of firms that hold only a single debt contract
s_tab(Ncontr, dlevel %>% group_by(companyid, year, quarter) %>% summarise(Ncontr = n()))


# Plots --------------------------------------------------------------

my_theme <-  theme(plot.title = element_text(face = "bold", hjust = 0.5, size =12),
                   axis.ticks.length=unit(0, "cm"),
                   axis.text.y = element_text(size = 11, angle = 0, vjust = 0.5, hjust=0.5),
                   axis.text.x = element_text(size = 11, angle = 0, vjust = 0.5, hjust=0.5),
                   axis.title.x=element_text(size = 11), axis.title.y=element_text(size = 11),
                   panel.background = element_rect(fill = "white"),
                   panel.grid.major = element_line(size = 0.2, linetype = 'solid', colour = "grey80"),
                   panel.grid.minor = element_line(size = 0.15, linetype = 'solid', colour = "grey90"),
                   panel.border = element_rect(color = "grey50", fill = NA, size = 0.5))

comciq <- comciq %>%  group_by(companyid) %>% 
  mutate(CFcat = ifelse(CFshare == 0, 0, 
                        ifelse(CFshare > 0 & CFshare < 1, 1, 2)))

aux <- s_tab(CFcat, comciq)[["percent"]]

p1 <- ggplot() + 
  stat_ecdf(dat = comciq, aes(CFshare),  color = "blue3",  geom = "step",  size = 1, alpha = 1) +
  geom_hline(yintercept = aux[1], color = "firebrick", linetype = "dashed",  size = 1) +
  geom_hline(yintercept = aux[2]+aux[1], color = "firebrick", linetype = "dashed",  size = 1) + 
    ggtitle("CFL reliance - Cummulative distribution")  +
  ylab("CFL Reliance") + xlab("")   + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(.01, .01)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(.01, .01)) +
  my_theme +  annotate("text", x = 0.05, y = 0.95, label = "Cash flow based borrowers: 12%", color = "grey30", vjust = 1, hjust = 0, fontface = "bold", size = 3 )  +
  my_theme +  annotate("text", x = 0.05, y = 0.60, label = "`Hybrid' borrowers: 52%", color = "grey30", vjust = 1, hjust = 0, fontface = "bold", size = 3)  +
  my_theme +  annotate("text", x = 0.05, y = 0.2, label = "Asset based borrowers: 36%", color = "grey30", vjust = 1, hjust = 0, fontface = "bold", size = 3) 

p2 <- ggplot(filter(comciq, year == 2018, sic_group != "Nonclassifiable" & sic_group != "Construction" &  lassets >=  5)) + 
  geom_smooth(aes(x = lassets, y = CFshare), fill = "lightblue", color = "blue3") +
  ggtitle("CFL reliance over firm size")  + 
  ylab("CFL Reliance") + xlab("Log of assets")   + 
  scale_x_continuous(expand=c(0,0)) +  scale_y_continuous(expand=c(0.01,0.01), limits = c(0,1)) +
  my_theme

cowplot::plot_grid(p1, p2, ncol = 2)
ggsave("smoothcfd.png", dpi = 400, height = 4, width = 7)


comciq <- comciq %>%  group_by(companyid) %>% 
  mutate(CFcat = ifelse(CFshare == 0, "ABL only", 
                        ifelse(CFshare > 0 & CFshare < 1, "Hybrid", "CFL only")))


p1 <- ggplot(filter(comciq, CFcat == "ABL only"), aes(x = lassets)) + 
  geom_histogram(color = "blue3", fill = "blue3", 
                 alpha = 0.5, position = "identity", bins = 50) + 
  scale_y_continuous(expand = c(0, 100), breaks = seq(0, 5000, 1500)) + 
  scale_x_continuous(limits = c(4,12)) + 
  my_theme

p2 <- ggplot(filter(comciq, CFcat == "Hybrid"), aes(x = lassets)) + 
  geom_histogram(color = "firebrick", fill = "firebrick", 
                 alpha = 0.5, position = "identity", bins = 50) + 
  scale_y_continuous(expand = c(0, 100)) + 
  scale_x_continuous(limits = c(4,12)) + 
  my_theme

p3 <- ggplot(filter(comciq, CFcat == "CFL only"), aes(x = lassets)) + 
  geom_histogram(color = "lightblue", fill = "lightblue", 
                 alpha = 0.8, position = "identity", bins = 50) + 
  scale_y_continuous(expand = c(0, 30)) + 
  scale_x_continuous(limits = c(4,12)) + 
  my_theme


cowplot::plot_grid(p1, p2, p3, ncol = 1)
ggsave("histog.png", dpi = 400, height = 6, width = 6)



ggplot(comciq, aes(lassets, fct_rev(CFcat), fill = CFcat)) + 
  ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = 2, alpha = 0.8, size = 0.6, scale = 1,  bandwidth = 0.1) +
  ggtitle("")  + 
  scale_fill_brewer(palette = "Blues", name = "CFL category") +  # Use the "Set2" palette
  ylab("") + xlab("Log of assets")   + 
  scale_y_discrete(expand = c(0, .1)) +
  my_theme
ggsave("ridges.png", dpi = 400, height = 6, width = 8)



  # t-test ------------------------------------------------------------------
# function
t_test <- function(data, var, gr_var) {
  varname <- deparse(substitute(var))
  
  vec1 <- data %>% filter({{ gr_var }} == 1) %>% pull({{ var }}) %>% na.omit()
  vec2 <- data %>% filter({{ gr_var }} == 0) %>% pull({{ var }}) %>% na.omit()
  
  ttest <- t.test(vec1, vec2)
  
  mean_1 <- mean(vec1)[]
  mean_2 <- mean(vec2)[]
  t_statistic <- ttest$statistic[]
  p_value <- ttest$p.value[]
  
  # Return the results
  res <- data.frame(varname, mean_1, mean_2, t_statistic, p_value  )
  
  return(res)
}

vars <- c("assets", "revenue", "emp", "age", "value", "levrg", "prod")

results_list <- lapply(vars, function(var) {
  t_test(comciq, var, CFcat)   })
results <- do.call(rbind, results_list) 
results$varname <- vars
results %>%  xtable::xtable( caption = "Summary Statistics") %>%
  print(type = "latex", include.rownames = TRUE, booktabs = TRUE,
        sanitize.text.function = function(x) x)



# Create a dataset where variables are average through time - last stata regression  --------------------------------------------
comciq <- read_csv("comciq_stata.csv")

# change in ratings are very rare so I will just assume they do not exist
comciq <- comciq %>%  mutate(aux = ifelse(rating != lag(rating) & !is.na(rating) & !is.na(lag(rating)), 1,0))
s_tab(aux, comciq)

vars <- setdiff(comciq %>% select_if(is.numeric) %>% names(),
                c("companyid", "year", "quarter", "dateper", "d2a", "consist_comciq", "mean_gr_assets", "mean_gr_emp"))

comciq_flevel <- comciq %>% group_by(companyid) %>% 
 summarise_at(vars, ~mean(., na.rm = T))


custom_mode <- function(x) {
  tbl <- table(x)
  tbl[tbl == max(tbl)]
}

comciq_flevel2 <- comciq %>% 
  group_by(companyid) %>%
  summarize(rating = paste(names(custom_mode(rating)), collapse = ""),
            sic_group = paste(names(custom_mode(sic_group)), collapse = ""),
            country = paste(names(custom_mode(country)), collapse = ""))


comciq_flevel <- left_join(comciq_flevel, comciq_flevel2, by = "companyid")
rm(list = c("comciq_flevel2"))

haven::write_dta(comciq_flevel, "comciq_flevel.dta", version = 13)

