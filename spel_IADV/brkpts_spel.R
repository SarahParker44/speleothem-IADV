#### change point analysis

setwd("speleothem-IADV/")

library(dplyr)
library(strucchange)
library(ggplot2)
library(rstatix)

## load spel composite data
spel_comp <- read.csv("spel_IADV/spel_sd.csv")
spel_comp <- spel_comp %>% filter(win_end <= 12000 & win_end >=50 & region %in% c("EAM","ISM","SW-SAM"))

opt_bpts <- function(x) {
  #x = bpts_sum$RSS["BIC",]
  n <- length(x)
  lowest <- vector("logical", length = n-1)
  lowest[1] <- FALSE
  for (i in 2:n) {
    lowest[i] <- x[i] < x[i-1] & x[i] < x[i+1]
  }
  out <- as.integer(names(x)[lowest])
  return(out)
}

comp_fitted <- data.frame()
comp_bp <- data.frame()
comp_lm_summ <- data.frame()
for (i in unique(spel_comp$region)){
  sub_comp <- filter(spel_comp, region == i)
  
  bp <- breakpoints(sub_comp$med ~ sub_comp$age) #breakpoint analysis
  bpts_sum <- summary(bp) 
  opt_brks <- opt_bpts(bpts_sum$RSS["BIC",]) #optimal no. breakpoints
  
  ci_x <- confint(bp, breaks = opt_brks) #get timings of bp's with conf intervals
  ci_ages <- data.frame(bp = sub_comp$age[ci_x$confint[,2]],
                        CI2_5 = sub_comp$age[ci_x$confint[,1]],
                        CI97_5 = sub_comp$age[ci_x$confint[,3]])
  
  ## model:
  #z <- sub_comp$med
  #x <- sub_comp$age
  #best_brk <- ci_ages$bp
  #fm1 <- lm(z ~ x*(x < best_brk[1]) + x*(x > best_brk[1] & x < best_brk[2]) + x*(x > best_brk[2]))
  
  #fac <- breakfactor(bp, breaks = opt_brks)
  #fm1 <- lm(sub_comp$med ~ fac)
  
  
  brk_vals <-  c(max(sub_comp$age),sub_comp$age[ci_x$confint[,2]], (min(sub_comp$age))-1)
  sub_comp$grp <- cut(x = sub_comp$age, breaks = brk_vals, labels = 1:(length(brk_vals)-1))
  
  subcomp_fit <- data.frame()
  subcomp_lm_summ <- data.frame()
  for (j in unique(sub_comp$grp)){
    grp_subcomp <- sub_comp %>% filter(grp == j)
    grp_lm <- lm(med ~ age, data = grp_subcomp)
    sub_lm <- as.data.frame(summary(grp_lm)$coefficients); sub_lm$grp <- j
    lm_fit <- fitted(grp_lm)
    grp_subcomp$fitted_val <- lm_fit
    
    subcomp_fit <- rbind(subcomp_fit, grp_subcomp)
    subcomp_lm_summ <- rbind(subcomp_lm_summ, sub_lm)
  }
  
  # save data
  comp_fitted <- rbind(comp_fitted, subcomp_fit)
  subcomp_lm_summ$region <- i; comp_lm_summ <- rbind(comp_lm_summ, subcomp_lm_summ)
  ci_ages$region <- i; comp_bp <- rbind(comp_bp, ci_ages)
}

write.csv(comp_bp, "spel_breakpoints.csv", row.names = F)

pdf("breakpoints_slope_plots.pdf", height = 15/2.54, width = 9/2.54)
ggplot() +
  geom_line(data = comp_fitted, aes(x = age, y = med), col = "#AAAAAA") +
  geom_line(data = comp_fitted, aes(x = age, y = fitted_val), col = "orange") +
  geom_vline(data = comp_bp, aes(xintercept = bp), lty = 2, col = "#485063") +
  geom_errorbarh(data = comp_bp, aes(xmin = CI2_5, xmax = CI97_5, y = 0.35), col = "red", height = 0.01) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white", colour = "white"), strip.text = element_text(hjust = 0)) +
  ylab("d18O s.d.") + 
  facet_wrap(.~ region, ncol = 1, nrow = 3)
dev.off()


## bp in mean
comp_fitted <- data.frame()
comp_bp <- data.frame()
comp_seg_mean <- data.frame()
comp_ttest <- data.frame()
for (i in unique(spel_comp$region)){
  sub_comp <- filter(spel_comp, region == i)
  
  bp <- breakpoints(sub_comp$med ~ 1) #breakpoint analysis
  bpts_sum <- summary(bp) 
  opt_brks <- opt_bpts(bpts_sum$RSS["BIC",]) #optimal no. breakpoints
  
  if (i == "ISM"){ opt_brks = 2 }
  
  ci_x <- confint(bp, breaks = opt_brks) #get timings of bp's with conf intervals
  ci_ages <- data.frame(bp = sub_comp$age[ci_x$confint[,2]],
                        CI2_5 = sub_comp$age[ci_x$confint[,1]],
                        CI97_5 = sub_comp$age[ci_x$confint[,3]])
  
  brk_vals <-  c(max(sub_comp$age),sub_comp$age[ci_x$confint[,2]], (min(sub_comp$age))-1)
  sub_comp$grp <- cut(x = sub_comp$age, breaks = brk_vals, labels = 1:(length(brk_vals)-1))
  
  sub_tvals <- sub_comp %>% 
    pairwise_t_test(med ~ grp)
  
  sub_meanvals <- sub_comp %>%
    group_by(grp) %>%
    summarise(mean(med))
  
  x_fac <- breakfactor(bp, breaks = opt_brks)
  fm1 <- lm(sub_comp$med ~ x_fac - 1)
  
  fit_sub <- data.frame(region = i,
                        Age = sub_comp$age,
                        bp_fit = fitted(fm1))


  # save data
  comp_fitted <- rbind(comp_fitted, fit_sub)
  sub_meanvals$region <- i; comp_seg_mean <- rbind(comp_seg_mean, sub_meanvals)
  sub_tvals$rgion <- i; comp_ttest <- rbind(comp_ttest, sub_tvals)
  #subcomp_lm_summ$region <- i; comp_lm_summ <- rbind(comp_lm_summ, subcomp_lm_summ)
  ci_ages$region <- i; comp_bp <- rbind(comp_bp, ci_ages)
}

pdf("breakpoints_mean_plots.pdf", height = 15/2.54, width = 9/2.54)
ggplot() +
  geom_line(data = spel_comp, aes(x = age, y = med), col = "#AAAAAA") +
  geom_line(data = comp_fitted, aes(x = Age, y = bp_fit), col = "#30D5C8") +
  geom_vline(data = comp_bp, aes(xintercept = bp), lty = 2, col = "#485063") +
  geom_errorbarh(data = comp_bp, aes(xmin = CI2_5, xmax = CI97_5, y = 0.35), col = "red", height = 0.01) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white", colour = "white"), strip.text = element_text(hjust = 0)) +
  geom_hline(yintercept = Inf, colour = "black") +
  ylab("d18O s.d.") + 
  facet_wrap(.~ region, ncol = 1, nrow = 3, )
dev.off()

