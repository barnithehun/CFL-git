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
                  axis.text.y = element_text(size = 14, angle = 0, vjust = 0.5, hjust=0.5),
                  axis.text.x = element_text(size = 14, angle = 0, vjust = 0.5, hjust=0.5),
                  axis.title.x=element_text(size = 16), axis.title.y=element_text(size = 14),
                  panel.background = element_rect(fill = "white"),
                  panel.grid.major = element_line(linewidth = 0.2, linetype = 'solid', colour = "grey80"),
                  panel.grid.minor = element_line(linewidth = 0.15, linetype = 'solid', colour = "grey90"),
                  panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.5))


# Import data -------------------------------------------------------------

comciq <- read_csv("comciq_stata.csv") %>% select(-c("gvkey"))
dlevel <- read_csv("dlevel_full.csv") %>% arrange(companyid, dateper) 

summary(dlevel$prod)

dlevel <- dlevel %>% mutate(CF = ifelse(CF == 1, "CF debt", "AB debt"))

ggplot(data = dlevel %>% filter(inflow == 1)) + 
  geom_point(aes(x = age, y = spread, color = as.factor(CF)), alpha = 0.05) +
  geom_smooth(aes(x = age, y = spread, color = as.factor(CF)), method = "lm", se = TRUE) +
  ggtitle("Credit Spreads Over Firm Size") + 
  ylab("Credit Spread") + 
  xlab("Log of Assets") + 
  scale_color_manual(
    name = "Loan-Type",
    values = c("#FF5555", "#377EB8")) +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  my_theme +
  theme(
    legend.position = c(0.95, 0.95), # Moves legend to upper right inside plot
    legend.justification = c(1, 1), 
    legend.text = element_text(size = 15), 
    legend.title = element_text(size = 16)
  )

ggsave("C://Users//szjud//OneDrive//Asztali gép//EBCs//CFL-git//Latex codes//Plots//spreads.png", dpi = 400, height = 6, width = 9)

hist(dlevel$ag)

## FOR HIGH PRODUCTIVITY FIRMS
# ggplot(data = dlevel %>% filter(inflow == 1 & ebit2a  > 0.025 )) + 
#   geom_point(aes(x = lassets, y = spread, color = as.factor(CF)), alpha = 0.05) +
#   geom_smooth(aes(x = lassets, y = spread, color = as.factor(CF)), method = "lm", se = TRUE) +
#   ggtitle("") + 
#   ylab("Credit Spread") + 
#   xlab("Log of Assets") + 
#   scale_x_continuous(limits = c(5, 11), expand = c(0, 0)) + 
#   scale_y_continuous(expand = c(0, 0)) + 
#   my_theme


# CF share only for inflows 
aux_CFshare <- dlevel %>% filter(inflow == 1) %>% group_by(companyid, dateper) %>%  
  mutate(ABval = sum(value * AB), 
         CFval = sum(value * CF), 
         CFshare_inflow = CFval/(ABval+CFval)) %>%   
  summarise(CFshare_inflow = mean(CFshare_inflow, na.rm = T))

comciq <- comciq %>% left_join(aux_CFshare, by=c("companyid", "dateper"))


# make intrate date separate for AB and CF debt
flevel_intrate <- dlevel  %>% filter(inflow == 1) %>% 
  group_by(companyid, dateper, CF) %>% 
  mutate(dweight = value / sum(value)) %>%
  summarise(matur = sum(dweight*matumax, na.rm = T),
            intrate = sum(dweight*intrate, na.rm = T),
            spread = sum(dweight*spread, na.rm = T)) %>% 
  mutate(CF = ifelse(CF == 1, "CF", "AB"))


# Use pivot_wider to spread multiple columns
flevel_intrate <- flevel_intrate %>%
  pivot_wider(names_from = CF, 
              values_from = c(spread, intrate, matur))

comciq <- comciq %>% left_join(flevel_intrate, by=c("companyid", "dateper"))


# DiD ---------------------------------------------------------------------
comciq <- comciq %>% mutate(
  treat = ifelse(value < 3, 1, NA),  # Treatment group: value < 2 million
  treat = ifelse(value > 10, 0, treat), # Control group: value > 10 million
  post = ifelse(year >= 2021, 1, NA),        # Post-treatment: year >= 2021
  post = ifelse(year < 2020, 0, post),
  DID = treat * post)       


did_model <- lm(CFshare_inflow ~ DID + post + treat + as.factor(year), data = comciq)
summary(did_model)


sumdat <- comciq %>% group_by(year, treat) %>% filter(!is.na(treat)) %>% 
  summarise(CFshare_inflow = mean(CFshare_inflow, na.rm = T ))


ggplot(data = sumdat) + 
  geom_line(aes(x = year, y = CFshare_inflow, color = as.factor(treat))) +
  ggtitle("") + 
  ylab("Credit Spread") + xlab("Log of Assets") + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  my_theme



