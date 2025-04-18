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


# Exit, Default or Continue -------------------------------------------------------------------

# Create some sample data
x <- -3.5:4.5
y <- -3.5:4.5

ggplot() +
  geom_line(aes(x,  y), color = "white", size = 1) +  # Scatter plot with sample data
  geom_segment(aes(x = 0, y = 0, xend =min(y) , yend = 0),  # Vertical line from (0,0) to (0,max(y))
               color = "blue4",linetype = "dashed", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = min(x)),  # Horizontal line from (0,0) to (max(x),0)
               color = "blue4", linetype = "dashed", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = max(x), yend = max(x)),  # 45-degree line from (0,0) to (max(x),max(y))
               color = "blue4", size = 1) +
  annotate("text", x = -1.73, y = -1.65,  label = "Default: \n V ≤ 0 and x ≤ 0", size = 4, color = "firebrick",  fontface = "bold") +  
  annotate("text", x = -1.73, y = 2.3,  label = "Continue: \n V > 0 and V > x", size = 4, color = "firebrick",  fontface = "bold") +  
  annotate("text", x = 2.3, y = -1.65,  label = "Exit: \n x > V and x > 0", size = 4, color = "firebrick",  fontface = "bold") +  
  annotate("text", x = 3, y = 3.5,  label = "45°", size = 3, color = "blue4",  fontface = "bold") +  
  theme_minimal() +  # Set theme to minimal
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10.5),
        panel.background = element_rect(fill = "white"),  # Set background color to white
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black", size = 0.2),  # Highlight x and y axes
        axis.title = element_text(size = 10, face = "bold"),  # Set axis title size
        axis.text = element_text(size = 10),    # Set axis text size
        axis.ticks.length = unit(0.2, "cm"),     # Set length of tick marks
        axis.ticks = element_line(color = "white"),  # Set color of tick marks
        axis.ticks.margin = unit(0.5, "cm")  # Set margin for tick marks
  ) +
  labs(title = "Exit, Default or Continue",
       x = "Cash on Hand",
       y = "Continuation Value") +
  scale_y_continuous(expand = c(0, 0), breaks = 0) +  # Set breaks for y-axis
  scale_x_continuous(expand = c(0, 0), breaks = 0)   # Set breaks for x-axis

ggsave("dec.png", dpi = 400, height = 6, width = 8)  



# Plots from Capcom -------------------------------------------------------

comciq <- comciq %>%  group_by(companyid) %>% 
  mutate(CFcat = ifelse(CFshare == 0, 0, 
                        ifelse(CFshare > 0 & CFshare < 1, 1, 2)))

aux <- s_tab(CFcat, comciq)[["percent"]]

# does not run under the new R version
p1 <- ggplot() + 
  stat_ecdf(data = comciq, aes(CFshare),  color = "blue3",  geom = "step",  linewidth = 1, alpha = 1) +
  geom_hline(yintercept = aux[1], color = "firebrick", linetype = "dashed",  linewidth = 1) +
  geom_hline(yintercept = aux[2]+aux[1], color = "firebrick", linetype = "dashed",  linewidth = 1) + 
  ggtitle("CFL reliance - Cummulative distribution")  +
  ylab("CFL Reliance") + xlab("")   + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(.01, .01)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(.01, .01)) +
  my_theme +  annotate("text", x = 0.05, y = 0.95, label = "Cash flow based borrowers: 12%", color = "grey30", vjust = 1, hjust = 0, fontface = "bold", size = 3 )  +
  annotate("text", x = 0.05, y = 0.60, label = "`Hybrid' borrowers: 52%", color = "grey30", vjust = 1, hjust = 0, fontface = "bold", size = 3)  +
  annotate("text", x = 0.05, y = 0.2, label = "Asset based borrowers: 36%", color = "grey30", vjust = 1, hjust = 0, fontface = "bold", size = 3) 

p2 <- ggplot(filter(comciq, year == 2018, sic_group != "Nonclassifiable" & sic_group != "Construction" &  lassets >=  5)) + 
  geom_smooth(aes(x = lassets, y = CFshare), fill = "lightblue", color = "blue3") +
  ggtitle("CFL reliance over firm size")  + 
  ylab("CFL Reliance") + xlab("Log of assets")   + 
  scale_x_continuous(expand=c(0,0)) +  scale_y_continuous(expand=c(0.01,0.01), limits = c(0,1)) +
  my_theme

cowplot::plot_grid(p2, p1, ncol = 2)
#ggsave("smoothcfd.png", dpi = 400, height = 3.5, width = 7)

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





