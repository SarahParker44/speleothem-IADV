## load model sd trends:

library(dplyr)
library(rstatix)

setwd("speleothem-IADV/")

model_sd <- read.csv("modelled_IADV/region_pre_sd.csv")
brkpts <- read.csv("spel_breakpoints.csv")
brkpts[which(brkpts$region == "SW-SAM"),"region"] <- "SAM"

subdat_fit <- data.frame()
subdat_lm_summ <- data.frame() # linear regression results (coeff, p-values, slope, intercept) of groups in each region-model combination
ttest_out <- data.frame() # t-test results amonst segments in each region-model combination
mean_out <- data.frame() # mean IADV values for groups in each region-model combination
for (i in unique(model_sd$region)){
  
  for (j in unique(model_sd$Model)){
    subdat <- model_sd %>% filter(region == i & Model == j)
    
    sub_brkpts <- brkpts %>% filter(region == i & bp < max(subdat$Age))
    
    brk_vals <-  c((min(subdat$Age)-1),sub_brkpts$bp, max(subdat$Age))
    subdat$grp <- cut(x = subdat$Age, breaks = brk_vals, labels = 1:(length(brk_vals)-1))
    
    sub_ttest <- subdat %>%
      pairwise_t_test(precip_sd ~ grp)
    
    sub_ttest$region <- i; sub_ttest$Model <- j
    
    ttest_out <- rbind(ttest_out, sub_ttest)
    
    sub_mean <- subdat %>% group_by(region, Model, grp) %>%
      summarise(mean_grp = mean(precip_sd))
    
    mean_out <- rbind(mean_out, 
                      sub_mean)
    
    for (k in unique(subdat$grp)){
      grp_subdat <- subdat %>% filter(grp == k)
      grp_lm <- lm(precip_sd ~ Age, data = grp_subdat)
      sub_lm <- as.data.frame(summary(grp_lm)$coefficients); sub_lm$grp <- j
      lm_fit <- fitted(grp_lm)
      grp_subdat$fitted_val <- lm_fit
      
      subdat_fit <- rbind(subdat_fit, grp_subdat)
      sub_lm$region <- i; sub_lm$Model <- j; subdat_lm_summ <- rbind(subdat_lm_summ, sub_lm)
    }
  }
}
