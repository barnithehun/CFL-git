# Linking FJC data to wrds does not work, there are very few observations that could be matched, that are no way representative
library(tidyverse)
library(tictoc)
library(RPostgreSQL)
library(RPostgres)
library(tictoc)
library(lubridate)


rm(list = ls())

setwd("C:/Users/szjud/OneDrive/Asztali gép/EBCs/Data/R/V3")
source("C:/Users/szjud/OneDrive/Asztali gép/EBCs/Data/R/V1/s_tab.R")


# reading in FJC data -----------------------------------------------------

# OPENING THE WRDS CONNECTION
wrds <- dbConnect(Postgres(),
                  host='wrds-pgdata.wharton.upenn.edu',
                  port=9737,
                  dbname='wrds',
                  sslmode='require',
                  user='szbarnabas',
                  password = 'Fsztznenezze410!!')

res <- dbSendQuery(wrds, "SELECT casekey, gvkey, filedate, cik
                                 
                          FROM fjc.wrds_bankruptcy_link
                          
                          WHERE filedate BETWEEN '1990-01-01' AND '2023-12-31'")
linking_fjc <- dbFetch(res, n = -1)
dbClearResult(res)

linking_fjc <- linking_fjc %>%
  mutate(cik = sub("^00", "", cik),
         year = year(filedate),
         quarter = quarter(filedate))


# with compustat ----------------------------------------------------------
master_fjc <- read_csv("C:/Users/szjud/OneDrive/Asztali gép/EBCs/Data/R/IDB dataset/master2.csv") %>% select(c("CASEKEY", "chapter")) %>% rename("casekey" = "CASEKEY")
linking_fjc <- left_join(linking_fjc, master_fjc, by =c("casekey"))


res <- dbSendQuery(wrds, "SELECT gvkey, cik
                                 
                          FROM comp.company")
linking_cik <- dbFetch(res, n = -1)
dbClearResult(res)

# adding CIK to compustat data
comps <- read_csv("comps_fjc.csv") %>% left_join(linking_cik, by = c("gvkey")) 
comps_fjc <- comps %>% left_join(linking_fjc, by = c("cik", "year", "quarter")) 
s_tab(chapter, comps_fjc)

# With capiq ----------------------------------------------------------------
cikiq <- read_csv("ciq_capiq.csv")

masterIDB <- read_csv("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//IDB dataset//master2.csv") %>%
  rename("casekey" = "CASEKEY") %>% mutate(year = year(ORGFLDT), quarter = quarter(ORGFLDT)) %>% select(c("year", "quarter", "casekey", "chapter", "EASST", "EASSTnum"))
# alternative to SNAPSHOT use ORGFLDT for time variable


# linking FJC linking table to IDB dataset
ciq_fjc <- linking_fjc %>%  mutate(year = year(filedate), quarter = quarter(filedate)) %>% left_join(masterIDB, by = c("casekey", "year", "quarter"))

# linking the result to capitalIQ and CIK codes
ciq_fjc <- ciq_fjc %>% left_join(cikiq, by = c("cik")) %>% select(-c("gvkey"))

# linking the result compustat data
flevel <- read_csv("C:/Users/szjud/OneDrive/Asztali gép/EBCs/Data/R/V3/flevel.csv")
flevel_intrate <- read_csv("C:/Users/szjud/OneDrive/Asztali gép/EBCs/Data/R/V3/flevel_intrate.csv") %>%
  mutate(year = 2010 + (dateper - 1) %/% 4, quarter = (dateper - 1) %% 4 + 1)

flevel <- flevel %>% left_join(flevel_intrate, by = c("companyid", "year", "quarter")) 

ciq_fjc <- ciq_fjc %>% left_join(flevel, by = c("companyid", "year", "quarter")) 

s_tab(chapter, ciq_fjc)
is.na(ciq_fjc$gvkey) %>% sum()

