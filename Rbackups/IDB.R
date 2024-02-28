# US BANKRUPTCY DATA FROM JUIDICIAL RECORDS
# source: https://www.fjc.gov/research/idb

library(tidyverse)
library(tictoc)
library(RPostgreSQL)
library(RPostgres)
library(tictoc)
library(haven)
library(dataMaid)
library(DataExplorer)

rm(list = ls())

setwd("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//IDB dataset")
source("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V1//s_tab.R")

my_theme <- theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10.5),
                  panel.spacing = unit(2, "lines"), strip.background=element_rect(fill="white"),
                  axis.ticks.length=unit(0, "cm"),
                  axis.text.y = element_text(size = 9, angle = 90, vjust = 0.5, hjust=0.5),
                  axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust=0.5),
                  axis.title.x=element_text(size = 10), axis.title.y=element_text(size = 10),
                  panel.background = element_rect(fill = "white"),
                  panel.grid.major = element_line(size = 0.2, linetype = 'solid', colour = "grey80"),
                  panel.grid.minor = element_line(size = 0.15, linetype = 'solid', colour = "grey90"),
                  panel.border = element_rect(color = "grey50", fill = NA, size = 0.5))

# assembling the full dataset --------------------------------------------------------------------

# things were originally downloaded in "sas7bdat" format but it took too much disk space
# I formatted them to RDS format.. ser original files at: 
#     https://www.fjc.gov/research/idb/bankruptcy-cases-filed-terminated-and-pending-fy-2008-present

#read in like this due to memory constraints
dat1 <- bind_rows(
  readRDS("DB10.RDS"),
  readRDS("DB11.RDS"),
  readRDS("DB12.RDS"),
  readRDS("DB13.RDS")) %>% filter(DBTRTYP == "u" & NTRDBT == "b") # I only keep the business debts by corporations!
dat2 <- readRDS("DB14to23.RDS") %>% filter(DBTRTYP == "u" & NTRDBT == "b")

# wide set of variables
master <- bind_rows(dat1, dat2) %>% select(c("SNAPSHOT", "SNAPFILE", "SNAPPEND", "SNAPCLOS", "CASEKEY", "ORGFLDT", "CRNTCHP", 
                    "EASST", "ELBLTS", "ASSTCASE", "SMLLBUS", "TOTASSTS", "SECURED", "UNSECPR", "UNSECNPR", "CLOSEDT", "CLCHPT"))
# automated EDA
create_report(master)
makeDataReport(master, output = "html", replace = TRUE)

saveRDS(master, "master.RDS")

# cleaning master ---------------------------------------------------------
master <- readRDS("master.RDS") 

## Firm Size
# TOTASSTS and EASST are not consistent! get rid of TOTASSTS an co.
# ggplot(master, aes(x = EASST, y = TOTASSTS)) +  geom_point() +  labs(x = "Category", y = "Value") 
master <- master %>% select(-c("TOTASSTS", "TOTDBT", "SECURED", "UNSECPR", "UNSECNPR")) 

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

# For each case, check the last time it was observed - I only need these
master <- master %>% arrange(SNAPSHOT) %>% group_by(CASEKEY) %>% 
  mutate(last_obs = ifelse(row_number() == max(row_number()), 1, 0)) %>% ungroup() %>% 
  filter(last_obs == 1)

# how long the case was open
master <- master %>% mutate(case_length = as.numeric(CLOSEDT - ORGFLDT))

# Only CH7 and CH11
master <- master %>% filter(!((CRNTCHP == 7 & CLCHPT == 11) | (CRNTCHP == 11 & CLCHPT == 7))) # weird observations in general

# If case is closed then CLCHPT if not closed then CRNTCHP
master <- master %>% mutate(chapter = ifelse(SNAPCLOS == 1, CLCHPT, CRNTCHP)) %>% 
  filter((chapter == 7 | chapter == 11) & EASST != "M")

# Mapping for values - needed for interpolation later
master <- master %>%
  mutate(EASSTnum = case_when(
    EASST == "A" ~ 50000,
    EASST == "B" ~ 100000,
    EASST == "C" ~ 500000,
    EASST == "D" ~ 1000000,
    EASST == "E" ~ 10000000,
    EASST == "F" ~ 50000000,
    EASST == "G" ~ 100000000,
    EASST == "H" ~ 500000000,
    EASST == "I" ~ 1000000000,
    EASST == "J" ~ 5000000000,
    EASST == "K" ~ 10000000000,
    EASST == "L" ~ 50000000000,
    TRUE ~ NA_real_ ))


# Mapping of labels
value_mapping <- c(
  "A" = "$0 to $50k",
  "B" = "$50k to $100k",
  "C" = "$100k to $500k",
  "D" = "$500k to $1m",
  "E" = "$1 to $10m",
  "F" = "$10 to $50m",
  "G" = "$50 to $100m",
  "H" = "$100 to $500m",
  "I" = "$500 to $1000m",
  "J" = "$1b + (pre 2015) ",
  "K" = "$1-$10b",
  "L" = "$10- $50b")

# Use the mutate function to update the EASST variable
master <- master %>% mutate(EASST = factor(EASST, levels = names(value_mapping), labels = value_mapping))

# summarized dataset for liquidatin probability
sum_master <- master %>% group_by(EASST, chapter) %>% 
  summarise(count = n(), totass = mean(EASSTnum)) %>% spread(chapter, count) %>% mutate(pr_liq = `7` / (`7`+`11`))


p1 <- ggplot(master, aes(x = EASST, fill = chapter)) +
  geom_bar(position = "dodge", stat = "count", alpha = 0.8) +
  labs(x = "Estimated total assets", y = "Number of cases (thousands)") +
  ggtitle("Chapter 7 and Chapter 11 Instances") +
  scale_y_continuous(expand = c(0.02, 0), labels = scales::comma_format(scale = 1e-3)) +
  scale_fill_manual(values = c("firebrick", "blue3")) +
  my_theme +
  theme(legend.position = "left", 
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5, "lines"))


p2 <- ggplot(sum_master, aes(x = EASST, y = pr_liq)) +
  geom_col(fill = "blue3", alpha = 0.8) +
  labs(x = "Estimated total assets", y = "Share of liquidations") +
  ggtitle("Unconditional Liquidation Probability") + 
  scale_y_continuous(expand = c(0.01,0), limits = c(0,0.8)) + my_theme


cowplot::plot_grid(p1, p2, ncol = 2)
ggsave("liqprob.png", dpi = 400, height = 4, width = 8)
# it would be nice too see these things over time as well

# logarithmic interpolation - taking the logs of x-axis before linear interpolation
intp <- approx(log10(sum_master$totass), sum_master$pr_liq, method = "linear",
               xout = seq(4, 11, length.out = 200), rule = 2)  # liq. prob stays the same out of bounds
intp <- cbind(intp$x, intp$y) %>%  as.tibble() %>% mutate(x = 10^V1)

ggplot(intp, aes(x = V1, y = V2)) + geom_line()


