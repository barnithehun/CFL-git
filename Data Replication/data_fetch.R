# COMPUSTAT THAT CAPIQ PANEL! DATA - 2010- 2024
# I RUN INTO MEMORY ISSUES, SO A LOT OF COMPUTATION IS SOLVED WITH FOR LOOPS

library(tidyverse)
library(tictoc)
library(RPostgreSQL)
library(RPostgres)
library(tictoc)

rm(list = ls())

setwd("C:/Users/szjud/OneDrive/Asztali gép/EBCs/Data/R/V4")
source("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V1//s_tab.R")

# 1a: CapitalIQ query -----------------------------------------------------------------

# OPENING THE WRDS CONNECTION
wrds <- dbConnect(Postgres(),
                  host='wrds-pgdata.wharton.upenn.edu',
                  port=9737,
                  dbname='wrds',
                  sslmode='require',
                  user='szbarnabas',
                  password = 'Fsdztth410!..')


res <- dbSendQuery(wrds, "select distinct table_schema
                   from information_schema.tables
                   where table_type ='VIEW'
                   or table_type ='FOREIGN TABLE'
                   order by table_schema")
data <- dbFetch(res, n=-1)
dbClearResult(res)
data

tic()

# Need to read in smaller batches otherwise I run into memory issues
start_year <- 2000
end_year <- 2024
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
  
  file_name <- paste0("C:/Users/szjud/OneDrive/Asztali gép/EBCs/Data/R/V4/raw_dat/", "capiq_", substr(start_date, 1, 4), "_", substr(end_date, 1, 4), ".csv")
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


# 1b: Compustat Query --------------------------------------------------------------
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
  
  file_name <- paste0("C:/Users/szjud/OneDrive/Asztali gép/EBCs/Data/R/V4/raw_dat/", "comps_", substr(start_date, 1, 4), "_", substr(end_date, 1, 4), ".csv")
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

