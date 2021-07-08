##
setwd("speleothem_IADV/")

library(ggplot2)
library(dplyr)

dat <- read.csv("modelled_IADV/region_pre_sd.csv")

# anomalise
anom_df <- data.frame()
for (i in unique(dat$region)){
  for(j in unique(dat$Model)){
    subdat <- dat %>% filter(region == i & Model == j)
    base_val <- subdat %>% filter(Age >=0 & Age <= 1000) %>% summarise(base_mean = mean(precip_sd, na.rm = T),
                                                                        base_sd = sd(precip_sd, na.rm = T))
    subdat$anom <- (subdat$precip_sd - base_val$base_mean)#/base_val$base_sd
    anom_df <- rbind(anom_df, subdat)
  }
}


## linear regressions
mod_reg_lm <- anom_df %>% 
  filter(Age <= 8000 & Age >= 0) %>%
  group_by(region, Model) %>%
  summarise(coeff = summary(lm(anom ~ Age))$coefficients[2,1],
            st_err = summary(lm(anom ~ Age))$coefficients[2,2],
            P_val = summary(lm(anom ~ Age))$coefficients[2,4])
