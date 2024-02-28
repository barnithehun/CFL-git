library(tidyverse)
library(tictoc)
library(cowplot)

rm(list = ls())
# rm(list = "capiiq")

setwd("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V3")
source("C://Users//szjud//OneDrive//Asztali gép//EBCs//Data//R//V1//s_tab.R")

# Plot 1 - Smoothed averages --------------------------------------------------------
 
comciq <- read_csv("comciq.csv") # %>% filter(year == 2015)
comciq_stata <- read_csv("comciq_stata.csv") # %>% filter(year == 2015)

my_theme <- theme(plot.title = element_text(face = "bold", hjust = 0.5, size =18),
      panel.spacing = unit(2, "lines"), strip.background=element_rect(fill="white"),
      axis.ticks.length=unit(0, "cm"),
      axis.text.y = element_text(size = 10, angle = 0, vjust = 0.5, hjust=0.5),
      axis.text.x = element_text(size = 10, angle = 0, vjust = 0.5, hjust=0.5),
      axis.title.x=element_text(size = 10), axis.title.y=element_text(size = 10),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(size = 0.2, linetype = 'solid', colour = "grey80"),
      panel.grid.minor = element_line(size = 0.15, linetype = 'solid', colour = "grey90"),
      panel.border = element_rect(color = "grey50", fill = NA, size = 0.5))


comciq <- comciq %>% mutate(lcollat = log10(collat),
                            lintang = log10(assets-collat),
                            dateper = (year - 2010)*4 + quarter)

# assets
p1 <- ggplot(filter(comciq,  sic_group != "Nonclassifiable" & sic_group != "Construction" &  lassets >=  5)) + 
  geom_smooth(aes(x = lassets, y = CFshare),  color = "blue4", fill = "lightblue") +
  ggtitle("")  + 
  ylab("") + xlab("Log of assets")   + 
  scale_x_continuous(expand=c(0,0)) +  scale_y_continuous(expand=c(0.01,0.01), limits = c(0,1)) +
  my_theme

# revenues
p2 <- ggplot(filter(comciq,  sic_group != "Nonclassifiable" & sic_group != "Construction" &  lassets >=  5)) + 
  geom_smooth(aes(x = log10(revenue), y = CFshare), color = "blue4", fill = "lightblue") +
  ggtitle("")  + 
  ylab("") + xlab("Log of revenues")   + 
  scale_x_continuous(expand=c(0,0)) +  scale_y_continuous(expand=c(0.01,0.01), limits = c(0,1)) +
  my_theme
  
# employment
p3 <- ggplot(filter(comciq,  sic_group != "Nonclassifiable" & sic_group != "Construction" &  log10(emp) < 2.5)) + 
  geom_smooth(aes(x = log10(emp), y = CFshare), color = "blue4", fill = "lightblue") +
  ggtitle("")  + 
  ylab("") + xlab("Log of employees")   + 
  scale_x_continuous(expand=c(0,0)) +  scale_y_continuous(expand=c(0.01,0.01), limits = c(0,1)) +
  my_theme
  
# Pledgeability
p4 <- ggplot(filter(comciq,  sic_group != "Nonclassifiable" & sic_group != "Construction" & pledge > 0 & pledge < 1)) + 
  geom_smooth(aes(x = pledge, y = CFshare), color = "blue4", fill = "lightblue") +
  # facet_wrap(~ sic_group, scales = "free_x") + 
  ggtitle("")  + 
  ylab("") + xlab("Pledgeability")   + 
  scale_x_continuous(expand=c(0,0)) +  scale_y_continuous(expand=c(0.01,0.01), limits = c(0,1)) +
  my_theme

# Leverage
p5 <- ggplot(filter(comciq,  sic_group != "Nonclassifiable" & sic_group != "Construction" & levrg > -0.5 & levrg < 1.5)) + 
  geom_smooth(aes(x = levrg, y = CFshare), color = "blue4", fill = "lightblue") +
  # facet_wrap(~ sic_group, scales = "free_x") + 
  ggtitle("")  + 
  ylab("") + xlab("Leverage")   + 
  scale_x_continuous(expand=c(0,0)) +  scale_y_continuous(expand=c(0.01,0.01), limits = c(0,1)) +
  my_theme


# interest coverage
p6 <- ggplot(filter(comciq,  sic_group != "Nonclassifiable" & sic_group != "Construction" &  -10 <  intcov  & intcov  < 20)) + 
  geom_smooth(aes(x = intcov, y = CFshare), color = "blue4", fill = "lightblue") +
  # facet_wrap(~ sic_group, scales = "free_x") + 
  ggtitle("")  + 
  ylab("") + xlab("Interest coverage")   + 
  scale_x_continuous(expand=c(0,0)) +  scale_y_continuous(expand=c(0.01,0.01), limits = c(0,1)) +
  my_theme


# liquidity
p7 <- ggplot(filter(comciq_stata,  sic_group != "Nonclassifiable" & sic_group != "Construction" & profitability > - 0.5)) + 
  geom_smooth(aes(x = profitability, y = CFshare), color = "blue4", fill = "lightblue") +
  # facet_wrap(~ sic_group, scales = "free_x") + 
  ggtitle("")  + 
  ylab("") + xlab("Profitability")   + 
  scale_x_continuous(expand=c(0,0)) +  scale_y_continuous(expand=c(0.01,0.01), limits = c(0,1)) +
  my_theme


# revenue
p8 <-  ggplot(filter(comciq_stata,  sic_group != "Nonclassifiable" & sic_group != "Construction" & liq < 0.4)) + 
  geom_smooth(aes(x = liq, y = CFshare), color = "blue4", fill = "lightblue") +
  # facet_wrap(~ sic_group, scales = "free_x") + 
  ggtitle("")  + 
  ylab("") + xlab("Liquidity")   + 
  scale_x_continuous(expand=c(0,0)) +  scale_y_continuous(expand=c(0.01,0.01), limits = c(0,1)) +
  my_theme


# age
p10 <- ggplot(filter(comciq,  sic_group != "Nonclassifiable" & sic_group != "Construction" & -1 < age & age < 151)) + 
  geom_smooth(aes(x = age, y = CFshare), color = "blue4", fill = "lightblue") +
  ggtitle("")  + 
  ylab("") + xlab("Firm age")   + 
  scale_x_continuous(expand=c(0,0)) +  scale_y_continuous(expand=c(0.01,0.01), limits = c(0,1)) +
  my_theme

plot_grid(p1, p2, p3, p4, p5, p6, p8, p7, p10, ncol = 3)
ggsave("varCFL.png", dpi = 400, height = 6, width = 8)

# rm(list = c("p1", "p2", "p3", "p4", "p5", "p6", "p8", "p9", "p10"))


# Plot 2 - CDF ------------------------------------------------------------

binnum = 5
comciq <- comciq %>%  mutate(bin_assets = statar::xtile(assets, binnum),
                             bin_assets = as.character(bin_assets))

ggplot(filter(comciq, !is.na(bin_assets))) + 
  stat_ecdf(aes(CFshare, color = bin_assets), geom = "step",  size = 1, alpha = 0.8) +
  stat_ecdf(aes(CFshare),  color = "firebrick",  geom = "step",  size = 1.4, alpha = 1) +
  facet_wrap(~ year, scales = "free_x") + 
  ggtitle("")  +
  ylab("Share of CF based borrowing") + xlab("")   + 
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(.01, .01)) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size =12),
        axis.ticks.length=unit(0, "cm"),
        axis.text.y = element_text(size = 11, angle = 0, vjust = 0.5, hjust=0.5),
        axis.text.x = element_text(size = 11, angle = 0, vjust = 0.5, hjust=0.5),
        axis.title.x=element_text(size = 11), axis.title.y=element_text(size = 11),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid', colour = "grey80"),
        panel.grid.minor = element_line(size = 0.15, linetype = 'solid', colour = "grey90"),
        panel.border = element_rect(color = "grey50", fill = NA, size = 0.5)) +
  scale_color_brewer(palette = "YlGn")

# ggsave("cdfs.png", dpi = 300, height = 6, width = 6)


# Plot 3 - Looking at entries to the CFL market  -----------------------------------
comciq <- comciq %>%  group_by(companyid) %>% 
  mutate(CFcat = ifelse( !is.na(CFshare) & CFshare == 0, 0, 1),
         CFenter = as.numeric((CFcat - lag(CFcat)) == 1),
         CFexit = as.numeric((CFcat - lag(CFcat)) == -1)) # there are some firms that only have CFcat = 1 for one period
# mutate(CFcat = ifelse(CFshare == 0, 0, ifelse(CFshare > 0 & CFshare < 1, 1, 2))

# Asset proportion of hybrid firms
# (comciq %>% group_by(CFcat) %>% summarise(assum = sum(assets, na.rm = T)))[2,"assum"] / sum((comciq %>% group_by(CFcat) %>% summarise(assum = sum(assets, na.rm = T)))["assum"])*100

my_theme <- theme(
  plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
  panel.spacing = unit(2, "lines"),
  strip.background = element_rect(fill = "white"),
  axis.ticks.length = unit(0, "cm"),
  axis.text = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size = 10),
  panel.background = element_rect(fill = "white"),
  panel.grid.major = element_line(size = 0.2, linetype = 'solid', colour = "grey80"),
  panel.grid.minor = element_line(size = 0.15, linetype = 'solid', colour = "grey90"),
  panel.border = element_rect(color = "grey50", fill = NA, size = 0.5))

p1 <- ggplot(filter(comciq, sic_group != "Nonclassifiable", age < 101, age >= 0)) + 
  geom_smooth(aes(x = age, y = CFenter), fill = "pink", color = "firebrick") +
  geom_smooth(aes(x = age, y = CFexit), fill = "lightblue", color = "blue4") +
    ggtitle("") + 
  ylab("Share of firms entering") + xlab("Age") + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +  # Set y-axis limits here
  my_theme

p2 <- ggplot(filter(comciq, sic_group != "Nonclassifiable", lassets > 4.99, lassets < 10.5)) + 
  geom_smooth(aes(x = lassets, y = CFenter), fill = "pink", color = "firebrick") +
  geom_smooth(aes(x = lassets, y = CFexit), fill = "lightblue", color = "blue4") +
    ggtitle("") + 
  ylab("Share of firms exiting") + xlab("Log of assets") + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +  # Set y-axis limits here
  my_theme


cowplot::plot_grid(p1, p2, ncol = 2)
ggsave("exit_entry.png", dpi = 400, height = 5, width = 8)

# Plot 4 - CF share evolution  -----------------------------------------------------

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


dlevel <- read_csv("dlevel_imputed.csv")
forplot <- sum_dlevel()
forplot$year_quarter <- zoo::as.yearqtr(paste(forplot$year, forplot$quarter, sep = "Q"))

ggplot(forplot, aes(x= year_quarter, y = CFshare)) + geom_line(color = "blue4", size = 1) +
  ylab("Share of CF debt") + xlab("Year-Quarter") + 
  scale_y_continuous(limits = c(0.5,1)) + 
  my_theme

ggsave("CFshare.png", dpi = 400, height = 5, width = 8)

  
