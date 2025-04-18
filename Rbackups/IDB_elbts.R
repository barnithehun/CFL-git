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

my_theme <- theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
                  panel.spacing = unit(2, "lines"), strip.background=element_rect(fill="white"),
                  axis.ticks.length=unit(0, "cm"),
                  axis.text.y = element_text(size = 12, angle = 0, vjust = 0.5, hjust=0.5),
                  axis.text.x = element_text(size = 12, angle = 75, vjust = 0.5, hjust=0.5),
                  axis.title.x=element_text(size = 14), axis.title.y=element_text(size = 14),
                  panel.background = element_rect(fill = "white"),
                  panel.grid.major = element_line(size = 0.2, linetype = 'solid', colour = "grey90"),
                  panel.grid.minor = element_line(size = 0.15, linetype = 'solid', colour = "grey95"),
                  panel.border = element_rect(color = "grey50", fill = NA, size = 0.5))

# assembling the full dataset --------------------------------------------------------------------

# things were originally downloaded in "sas7bdat" format but it took too much disk space
# I formatted them to RDS format.. ser original files at: 
#     https://www.fjc.gov/research/idb/bankruptcy-cases-filed-terminated-and-pending-fy-2008-present

#read in like this due to memory constraints
dat1 <- bind_rows(
  readRDS("DB8.RDS"),
  readRDS("DB9.RDS"),
  readRDS("DB10.RDS"),
  readRDS("DB11.RDS"),
  readRDS("DB12.RDS"),
  readRDS("DB13.RDS")) %>% filter(DBTRTYP == "u" & NTRDBT == "b") # I only keep the business debts by corporations!
dat2 <- readRDS("DB14to23.RDS") %>% filter(DBTRTYP == "u" & NTRDBT == "b")

# wide set of variables
master <- bind_rows(dat1, dat2) %>% select(c("SNAPSHOT", "SNAPFILE", "SNAPPEND", "SNAPCLOS", "CASEKEY", "ORGFLDT", "ORGFLCHP", "CRNTCHP", 
                    "EASST", "ELBLTS", "ASSTCASE", "SMLLBUS", "TOTASSTS", "SECURED", "UNSECPR", "UNSECNPR", "CLOSEDT", "CLCHPT"))
# # automated EDA
# create_report(master)
# makeDataReport(master, output = "html", replace = TRUE)

saveRDS(master, "master.RDS")

# cleaning master ---------------------------------------------------------
master <- readRDS("master.RDS") 

## Firm Size
# TOTASSTS and EASST are not consistent! get rid of TOTASSTS an co.
# ggplot(master, aes(x = EASST, y = TOTASSTS)) +  geom_point() +  labs(x = "Category", y = "Value") 
master <- master %>% select(-c("TOTASSTS", "SECURED", "UNSECPR", "UNSECNPR")) 

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

# How often do firms convert from liquidation to reorg
master_closed <- filter(master, SNAPCLOS == 1)
table(master_closed$chapter, master_closed$ORGFLCHP)
table(master_closed$chapter)

# Mapping for values - needed for interpolation later
master <- master %>%
  mutate(ELBLTSnum = case_when(
    ELBLTS == "A" ~ 50000,
    ELBLTS == "B" ~ 100000,
    ELBLTS == "C" ~ 500000,
    ELBLTS == "D" ~ 1000000,
    ELBLTS == "E" ~ 10000000,
    ELBLTS == "F" ~ 50000000,
    ELBLTS == "G" ~ 100000000,
    ELBLTS == "H" ~ 500000000,
    ELBLTS == "I" ~ 1000000000,
    ELBLTS == "J" ~ 5000000000,
    ELBLTS == "K" ~ 10000000000,
    ELBLTS == "L" ~ 50000000000,
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
  "I" =  "$500 to $1000m",
  "J" = "$1b + (pre 2015) ",
  "K" = "$1-$10b",
  "L" = "$10- $50b")

# Use the mutate function to update the EASST variable
master <- master %>% mutate(EASST = factor(EASST, levels = names(value_mapping), labels = value_mapping))

write.csv(master, "master2.csv")

# Interpreting results ----------------------------------------------------
master <- read_csv("master2.csv")

submaster <- master %>% filter(!is.na(CLCHPT) & SNAPSHOT >= "2010-01-01")
unique(submaster$CASEKEY) %>%  length()


master <- master %>% filter(ELBLTSnum > 100000) %>%  mutate(Less10m = ifelse(ELBLTSnum < 10^7, 1, 0))
# master <- master %>% filter(EASST != "$0 to $50k" & EASST != "$50k to $100k") %>%  mutate(Less10m = ifelse(ELBLTSnum < 10^7, 1, 0))
s_tab(Less10m, master) 
master %>% group_by(Less10m, chapter) %>% summarise(count = n()) %>% spread(chapter, count) %>% mutate(pr_liq = `7` / (`7`+`11`))


# summarized dataset for liquidation probability
sum_master <- master %>% group_by(ELBLTS, chapter) %>% 
  summarise(count = n(), totass = mean(ELBLTSnum)) %>% spread(chapter, count) %>% mutate(pr_liq = `7` / (`7`+`11`)) %>% arrange(totass)


write_csv(sum_master, "sum_master.csv")


sum_master$EASST <- factor(sum_master$ELBLTS, levels = c("level1", "level2", "level3", "level4", "level5", "level6", 
                                                        "level7", "level8", "level9", "level10", "level11", "level12"))

p1 <- ggplot(master %>% filter(year(SNAPSHOT) >= 2010), aes(x = as.factor(ELBLTSnum), fill = as.factor(chapter))) +
  geom_bar(position = "dodge", stat = "count", alpha = 0.8) +
  labs(x = "Estimated total assets", y = "Number of cases (thousands)", fill = "Chapter:") +
  ggtitle("Chapter 7 and Chapter 11 Filings") +
  scale_y_continuous(expand = c(0.02, 0), labels = scales::comma_format(scale = 1e-3)) +
  scale_fill_manual(values = c("#FF5555", "#377EB8")) +
  scale_x_discrete(labels = c("$0 to $50k", "$50k to $100k", "$100k to $500k", "$500k to $1m", "$1 to $10m", "$10 to $50m", "$50 to $100m",
                              "$100 to $500m", "$500 to $1000m", "$1b + (pre 2015)", "$1-$10b", "$10- $50b")) + 
  my_theme +
  theme(legend.position = c(0.8, 0.9),  # Adjust the x and y coordinates as needed
        legend.text = element_text(size = 10),
        legend.key.size = unit(1.1, "lines"))


p2 <- ggplot(sum_master, aes(x = as.factor(totass), y = pr_liq, group = 1)) +
  geom_line(color = "#FF5555", linewidth = 2, alpha = 1) +
  labs(x = "Estimated Total Assets", y = "Chapter 7 to Chapter 11 filings") +
  ggtitle("Share of Liquidations") + 
  scale_y_continuous(expand = c(0.01,0), limits = c(0,0.8)) + 
  scale_x_discrete(labels = c( "$0 to $50k", "$50k to $100k","$100k to $500k", "$500k to $1m","$1 to $10m", "$10 to $50m","$50 to $100m",
                               "$100 to $500m","$500 to $1000m", "$1b + (pre 2015)","$1-$10b", "$10- $50b")) + 
  my_theme


cowplot::plot_grid(p1, p2, ncol = 2)
ggsave("C://Users//szjud//OneDrive//Asztali gép//EBCs//CFL-git//Latex codes//Plots//liqprob2.png", dpi = 400, height = 6, width = 10)
# it would be nice too see these things over time as well

# logarithmic interpolation - taking the logs of x-axis before linear interpolation
intp <- approx(log10(sum_master$totass), sum_master$pr_liq, method = "linear",
               xout = seq(4, 11, length.out = 200), rule = 2)  # liq. prob stays the same out of bounds
intp <- cbind(intp$x, intp$y) %>%  as.tibble() %>% mutate(x = 10^V1)

ggplot(intp, aes(x = V1, y = V2)) + geom_line()

write.csv(sum_master, "liquidation_probs.csv")

# liqprob for SMEs --------------------------------------------------------

master <- master %>% mutate(year = year(SNAPSHOT))
master <- master %>% mutate(sme = ifelse(EASSTnum <= 5*10^7, 1, 0))

sum_master <- master %>% group_by(EASST, year, chapter) %>% 
  summarise(count = n()) %>% spread(chapter, count) %>% mutate(pr_liq = `7` / (`7`+`11`))

ggplot(sum_master) +
  geom_line(aes(x = year, y = pr_liq, color = as.factor(EASST))) +
  labs(x = "Estimated total assets", y = "Share of liquidations", linewidth = 4) +
  ggtitle("Unconditional Liquidation Probability") + 
  scale_y_continuous(expand = c(0.01,0), limits = c(0,0.8)) + my_theme


s_tab(EASSTnum, master)
names(master)




