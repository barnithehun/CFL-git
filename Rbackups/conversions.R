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


# New data ----------------------------------------------------------------
master <- readRDS("master.RDS") %>% select(-c("TOTASSTS", "SECURED", "UNSECPR", "UNSECNPR")) 

master <- master %>%
  filter(SNAPCLOS == 1) %>%
  mutate(case_length = as.numeric(CLOSEDT - ORGFLDT),
         converted = if_else(ORGFLCHP == 11 & CRNTCHP == 7, 1, 0))
s_tab(converted, filter(master,ORGFLCHP == 11))

# Only CH7 and CH11
master <- master %>% filter(CRNTCHP == 7 | CRNTCHP == 11)
master <- master %>% filter(ORGFLCHP == 7 | ORGFLCHP == 11)

# How often do firms convert from liquidation to reorg
table(master$CRNTCHP, master$ORGFLCHP)


# Filter and plot
ggplot(filter(master, ORGFLCHP == 11 & ORGFLDT >= as.Date("2000-01-01")), 
       aes(x = ORGFLDT, y = converted)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = TRUE) +
  labs(title = "Conversions Over Time",
       x = "Filing Date",
       y = "Conversion Rate") +
  theme_minimal()


# Filter and plot
ggplot(filter(master, ORGFLCHP == 11 & ORGFLDT >= as.Date("2010-01-01")), 
       aes(x = ORGFLDT, y = case_length)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = TRUE) +
  labs(title = "Conversions Over Time",
       x = "Filing Date",
       y = "Conversion Rate") +
  theme_minimal()
