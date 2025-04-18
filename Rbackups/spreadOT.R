# Here I check how credits spreads change over time ceteris paribus 
# Does not really work, some specifications show the wanted results 
# but these aren't robust to chaging age or log(age) and the age limit
# similarly, regressions on CF based debt do not have the correct sign either

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
                  axis.title.x=element_text(size = 12), axis.title.y=element_text(size = 12),
                  panel.background = element_rect(fill = "white"),
                  panel.grid.major = element_line(linewidth = 0.2, linetype = 'solid', colour = "grey80"),
                  panel.grid.minor = element_line(linewidth = 0.15, linetype = 'solid', colour = "grey90"),
                  panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.5))

# FLEVEL ---------------------------------------

flevel_intrate <- read_csv("flevel_intrate.csv") %>% select(-c("gvkey.x", "gvkey.y"))
comciq <- read_csv("comciq_stata.csv") %>% select(-c("gvkey"))

comciq <- comciq %>% left_join(flevel_intrate, by=c("companyid", "dateper"))

# haven::write_dta(comciq, "flevel_full.dta", version = 13)
# write_csv(comciq, "flevel_full.csv")


# duplications
comciq <- comciq %>% group_by(companyid, dateper) %>%
  mutate(dup = ifelse(row_number() > 1, row_number(), 0)) 
s_tab(dup, comciq)


# Smooths on spreads
ggplot(data = filter(comciq, age < 50 )) + 
  geom_point(aes(x = age, y = intrateCF ), color = "pink", alpha = 0.03) +
  geom_point(aes(x = age, y = intrateAB ), color = "lightblue", alpha = 0.03) +
  geom_smooth(aes(x = age, y = intrateCF ),  fill = "pink", color = "firebrick", se = TRUE) +
  geom_smooth(aes(x = age, y = intrateAB ),  fill = "lightblue", color = "blue4", se = TRUE) +
  ggtitle("") + 
  ylab("Credit Spread") + xlab("Log of Assets") + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + my_theme


# Log-transformation of variables
comciq <- comciq %>%
  mutate(
    lasset = log(assets),
    lemp = log(emp),
    lage = log(age),
    lvalue = log(value),
    lrevenue = log(revenue),
    # Log-modulus transformation
    lebitda = sign(ebitda) * log(abs(ebitda) + 1))

# Filter and clean Ebitda values
comciq <- comciq %>%
  mutate(
    ebit2a = ifelse(ebit2a < -0.5 | ebit2a > 0.2, NA, ebit2a),
    ebitda = ifelse(ebitda < -146.06 | ebitda > 11002, NA, ebitda)
  )


comciq_clean <- comciq %>%
  filter_all(all_vars(!is.infinite(.))) %>%
  drop_na(spreadCF, ebitda, lasset, pledge, levrg, lage, lemp, CFshare, sic_group, rating, dateper)


agelim = 25

modelCF <- lm(
  spreadCF ~ lebitda + lasset + pledge + levrg + age + lemp + CFshare + 
    factor(sic_group) + factor(rating) + factor(dateper) + factor(companyid),
  data = comciq_clean %>% filter(age < agelim)) %>% summary()

modelAB <- lm(
  spreadAB ~ lebitda + lasset + pledge + levrg + age + lemp + CFshare + 
    factor(sic_group) + factor(rating) + factor(dateper) + factor(companyid),
  data = comciq_clean %>% filter(age < agelim)) %>% summary()


modelAB$coefficients[6,]
modelCF$coefficients[6,]

# Regressions to play around with ---------------------------------------

# duplications
comciq <- comciq %>% group_by(companyid, dateper) %>%
  mutate(dup = ifelse(row_number() > 1, row_number(), 0)) 
s_tab(dup, comciq)


comciq <- comciq %>%  mutate(SME = ifelse(assets <= 50, 1,0))
s_tab(SME, comciq)

# Smooths on spreads
ggplot(data = filter(comciq, assets >  10 )) + 
  # geom_point(aes(x = log10(assets), y = d2a ), color = "lightblue", alpha = 0.03) +
  geom_smooth(aes(x = log10(assets), y = d2a ),  fill = "pink", color = "firebrick", se = TRUE) +
  ggtitle("") + 
  ylab("Credit Spread") + xlab("Log of Assets") + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + my_theme


# Log-transformation of variables
comciq <- comciq %>%
  mutate(
    lasset = log(assets),
    lemp = log(emp),
    lage = log(age),
    lvalue = log(value),
    lrevenue = log(revenue),
    # Log-modulus transformation
    lebitda = sign(ebitda) * log(abs(ebitda) + 1))

# Filter and clean Ebitda values
comciq <- comciq %>%
  mutate(
    ebit2a = ifelse(ebit2a < -0.5 | ebit2a > 0.2, NA, ebit2a),
    ebitda = ifelse(ebitda < -146.06 | ebitda > 11002, NA, ebitda)
  )


comciq_clean <- comciq %>%
  filter(age < 20) %>%
  filter_all(all_vars(!is.infinite(.))) %>%
  drop_na(spreadCF, ebitda, lasset, pledge, levrg, lage, lemp, CFshare, sic_group, rating, dateper)


agelim = 1000
lm(
  CFshare ~ lebitda + lasset + pledge + levrg + lage + lemp + CFshare + 
    factor(sic_group) + factor(rating) + factor(dateper) + factor(companyid),
  data = comciq_clean %>% filter(age < agelim)) %>% summary()

lm(
  spreadAB ~ lebitda + lasset + pledge + levrg + lage + lemp + CFshare + 
    factor(sic_group) + factor(rating) + factor(dateper),
  data = comciq_clean %>% filter(age < agelim)) %>% summary()
  