#### 1) Model relationships between spel. d18O stdev and confounding factors
#### 2) Apply correction to spel d18O s.d.
#### 3) Construct and plot uncorrected and corrected spel d18O s.d. (IADV) trends

setwd("speleothem-IADV/")

library(RMySQL)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)

## Load and filter spel data 
# Connect to SISAL SQL database
mydb <- dbConnect(MySQL(), user = "root", password = "", dbname = "sisalv2", 
                  host = "localhost")

# Select tropical/sub tropical sites and Holocene data
Raw_Data <- dbGetQuery(mydb, "SELECT * FROM site JOIN entity USING (site_id) JOIN sample USING (entity_id) JOIN original_chronology USING (sample_id) JOIN d18O USING (sample_id)
                       WHERE (latitude BETWEEN -35 AND 45) AND (interp_age <= 12000);")

Raw_Data <- Raw_Data %>% filter(entity_status != "superseded")

# Filter to monsoon entities
Raw_Data$region <- with(Raw_Data, ifelse(latitude >= 15 & latitude <= 35 & longitude >= 75 & longitude <= 98, "ISM",
                                         ifelse(latitude >= 20 & latitude <= 45 & longitude >= 100 & longitude <= 125, "EAM", 
                                                ifelse(latitude >= -10 & latitude <= 0 & longitude >= -80 & longitude <= -70 | latitude >= -30 & latitude <= -10 & longitude >= -60 & longitude <= -30, "SW-SAM",
                                                       ifelse(latitude >= -30 & latitude <= 5 & longitude >= 80 & longitude <= 170, "IAM",
                                                              ifelse(latitude >= -30 & latitude <= 0 & longitude >= 0 & longitude <= 50, "SAfM", 
                                                                     ifelse(latitude >= 0 & latitude <= 35 & longitude >= -110 & longitude <= -50, "CAM", 
                                                                            ifelse(latitude >= -10 & latitude <= 0 & longitude >= -70 & longitude <= -30, "NE-SAM", "other"))))))))

Raw_Data <- Raw_Data %>% filter(region != "other")


# filter to records at least 100 yrs long and at temporal res of at least 20yrs
source("spel_IADV/coverage_sampled_vs_gap.R") #loads function that gives entity length, excluding gaps and hiatuses
dat_length <- data.frame()
for (i in unique(Raw_Data$entity_id)){
  length <- get_ent_coverage(entity_id = i, age_start = -50, age_end = 12000)[,2]
  sub_df <- data.frame(entity_id = i, length = length)
  dat_length <- rbind(dat_length, sub_df)
}
dat_length <- dat_length %>% filter(length >= 100)
length_df <- data.frame()
source("spel_IADV/entity_sampling_mean_res.R") #loads function that gives average temporal res of entity
for (i in 1:nrow(dat_length)){
  datlen <- get_ent_sampling(entity_id = dat_length$entity_id[i], age_start = -50, age_end = 12000)$sampling_mean
  dat_len <- data.frame(c(dat_length[i,], datlen))
  colnames(dat_len)[3] <- c("mean_res")
  length_df <- rbind(length_df, dat_len)
}

length_df <- length_df %>% filter(mean_res <= 20)

# Calculate d18O stdev for a running window
bin_size <- 100
bin_hw <- bin_size/2
centres <- seq((0-bin_hw),(12000-bin_hw),bin_hw)

dat_ls <- list()
for (j in 1:length(length_df$entity_id)){ #each entity
  ent <- length_df$entity_id[j]
  dat <- Raw_Data %>% filter(entity_id == ent)

  df_out <- data.frame()
  for (i in centres){
    datsub <- dat %>% filter(interp_age >= i-bin_hw & interp_age <= i+bin_hw) #subset to bin
    
    if (nrow(datsub) <= 1){ next } #skip if one or no samples in this bin, can't calculate s.d.
    
    sub_sd <- sd(datsub$d18O_measurement) #calc stdev for bin
    n <- nrow(datsub) # calc nsamples per window
    mean_d18O <- mean(datsub$d18O_measurement) #calc mean d18O per window
    
    #extract measurement uncertainty for bin
    if (length(unique(datsub$d18O_precision)) == 1){
      meas_uncert <- unique(datsub$d18O_precision)
    } else {
      meas_uncert <- mean(datsub$d18O_precision)
    }
    
    if (all(is.na(datsub$depth_sample)) == T){ # if no sample depth data exists
      #print(paste(ent, i))
      #output
      df_sub <- data.frame(entity_id = ent, region = unique(dat$region), win_start = i-bin_hw, win_end = i+bin_hw, gradient = NA, 
                           stdev = sub_sd, nsamples = n, mean_d18O = mean_d18O, meas_uncert = meas_uncert)
    } else { #if sample depth data exists, calculate growth rate
      sub_reg <- lm(datsub$depth_sample ~ datsub$interp_age)
      growth_rate <- sub_reg$coefficients[2]
      if (unique(dat$depth_ref == "from base")){ growth_rate <- growth_rate*-1 } #if depth is from base rather than from top, flip growth rate
      #output
      df_sub <- data.frame(entity_id = ent, region = unique(dat$region), win_start = i-bin_hw, win_end = i+bin_hw, gradient = growth_rate, 
                           stdev = sub_sd, nsamples = n, mean_d18O = mean_d18O, meas_uncert = meas_uncert)
    }
    rownames(df_sub) <- NULL
    df_out <- rbind(df_out, df_sub)
  }
  
  dat_ls[[j]] <- df_out # data output
  
}

dat_df <- bind_rows(dat_ls)

# Filter to entities with >50% windows with stdev above uncert
x <- dat_df %>%
  group_by(entity_id) %>%
  summarise(all_n = n(),
            below_n = length(stdev[meas_uncert > stdev])) %>%
  mutate(prop = below_n/all_n) %>%
  filter(prop < 0.5) # only 1 entity > 0.5

dat_df <- dat_df %>% filter(entity_id %in% x$entity_id)

# remove windows where stdev < measurement uncertainty
dat_df <- dat_df %>% filter(stdev > meas_uncert)

# add site metadata
df_sitedat <- left_join(dat_df, unique(Raw_Data[,c("site_id","site_name","entity_id","latitude","longitude")]), by = "entity_id")

#output sites for map (fig 1)
sites <- df_sitedat %>% group_by(site_id, latitude, longitude) %>% summarise(n())
write.csv(sites[,-4], "spel_IADV/sites_stdev_analysis.csv", row.names = F)

## colinearity of predictor variables?
cor(df_sitedat$gradient, df_sitedat$nsamples, use = "complete.obs")
cor(df_sitedat$gradient, df_sitedat$mean_d18O, use = "complete.obs")
cor(df_sitedat$nsamples, df_sitedat$mean_d18O, use = "complete.obs")



#### Fit multiple linear regression model

## NH mlr
NH_dat <- df_sitedat %>% filter(latitude > 0)
#NH_mlr <- lm(log(stdev) ~ log(nsamples) + mean_d18O, data = NH_dat, na.action = "na.exclude") ##nsamples instead of growth rate - in supplement 
NH_mlr <- lm(log(stdev) ~ log(gradient) + mean_d18O, data = NH_dat, na.action = "na.exclude")
NH_mlr_summ <- summary(NH_mlr)

## SH mlr
SH_dat <- df_sitedat %>% filter(latitude < 0)
#SH_mlr <- lm(log(stdev) ~ log(nsamples) + mean_d18O, data = SH_dat, na.action = "na.exclude") ##nsamples instead of growth rate - in supplement 
SH_mlr <- lm(log(stdev) ~ log(gradient) + mean_d18O, data = SH_dat, na.action = "na.exclude")
SH_mlr_summ <- summary(SH_mlr)

#### MLR residual plots:

## get residuals
# NH
partial.res <- residuals(NH_mlr, "partial") 
colnames(partial.res) <- paste(colnames(partial.res), "_res", sep = "")
NH_site_metdat <- df_sitedat %>% filter(latitude > 0)
resid_dat_NH <- cbind(NH_site_metdat, partial.res)

# SH
partial.res <- residuals(SH_mlr, "partial") 
colnames(partial.res) <- paste(colnames(partial.res), "_res", sep = "")
SH_site_metdat <- df_sitedat %>% filter(latitude < 0)
resid_dat_SH <- cbind(SH_site_metdat, partial.res)

## plots
#growth rate NH
P1 <- ggplot(data = resid_dat_NH, aes(x = log(gradient), y = `log(gradient)_res`)) + geom_point(size = 0.2, col = "#AAAAAA") +
  geom_abline(slope = lm(`log(gradient)_res` ~ log(gradient), data = resid_dat_NH)$coefficients[2], intercept = lm(`log(gradient)_res` ~ log(gradient), data = resid_dat_NH)$coefficients[1], col = "#464646") +
  theme_bw() +
  #theme(text = element_text(size = 12)) +
  xlab("") + ylab("")

# alternative plot for MLR with n-samples - in supplement
#P1 <- ggplot(data = resid_dat_NH, aes(x = log(nsamples), y = `log(nsamples)_res`)) + geom_point(size = 0.2, col = "#AAAAAA") +
#  geom_abline(slope = lm(`log(nsamples)_res` ~ log(nsamples), data = resid_dat_NH)$coefficients[2], intercept = lm(`log(nsamples)_res` ~ log(nsamples), data = resid_dat_NH)$coefficients[1], col = "#464646") +
#  theme_bw() +
#  #theme(text = element_text(size = 12)) +
#  xlab("") + ylab("")

#mean d18O NH  
P2 <- ggplot(data = resid_dat_NH, aes(x = mean_d18O, y = mean_d18O_res)) + geom_point(size = 0.2, col = "#AAAAAA") +
  geom_abline(slope = lm(resid_dat_NH$mean_d18O_res ~ resid_dat_NH$mean_d18O)$coefficients[2], intercept = lm(mean_d18O_res ~ mean_d18O, data = resid_dat_NH)$coefficients[1], col = "#464646") +
  theme_bw() +
  #theme(text = element_text(size = 12)) +
  xlab("") + ylab("")

#growth rate SH
P3 <- ggplot(data = resid_dat_SH, aes(x = log(gradient), y = `log(gradient)_res`)) + geom_point(size = 0.2, col = "#AAAAAA") +
  geom_abline(slope = lm(`log(gradient)_res` ~ log(gradient), data = resid_dat_SH)$coefficients[2], intercept = lm(`log(gradient)_res` ~ log(gradient), data = resid_dat_SH)$coefficients[1], col = "#464646") +
  theme_bw() +
  #theme(text = element_text(size = 12)) #+
  #xlab("log(growth rate) (mm/year)") + ylab("")
  xlab("") + ylab("")

# alternative plot for MLR with n-samples - in supplement
#P3 <- ggplot(data = resid_dat_SH, aes(x = log(nsamples), y = `log(nsamples)_res`)) + geom_point(size = 0.2, col = "#AAAAAA") +
#  geom_abline(slope = lm(`log(nsamples)_res` ~ log(nsamples), data = resid_dat_SH)$coefficients[2], intercept = lm(`log(nsamples)_res` ~ log(nsamples), data = resid_dat_SH)$coefficients[1], col = "#464646") +
#  theme_bw() +
#  #theme(text = element_text(size = 12)) #+
#  #xlab("log(growth rate) (mm/year)") + ylab("")
#  xlab("") + ylab("")

#mean d18O SH
P4 <- ggplot(data = resid_dat_SH, aes(x = mean_d18O, y = mean_d18O_res)) + geom_point(size = 0.2, col = "#AAAAAA") +
  geom_abline(slope = lm(resid_dat_SH$mean_d18O_res ~ resid_dat_SH$mean_d18O)$coefficients[2], intercept = lm(mean_d18O_res ~ mean_d18O, data = resid_dat_SH)$coefficients[1], col = "#464646") +
  theme_bw() +
  #theme(text = element_text(size = 12)) #+
  #xlab(expression(paste("mean ", delta^{18}, "O (\u2030)"))) + ylab("")
  xlab("") + ylab("")

pdf("MLR_fig2.pdf", width = 12/2.54, height = 8/2.54)
grid.arrange(P1,P2,P3,P4, 
             widths = c(1,1),
             layout_matrix = rbind(c(1,2),
                                   c(3,4)),
             left = textGrob(expression(paste("f(", delta^{18}, "O s.d.)")), rot = 90))
dev.off()


#### apply correction
# predict stdev from model:
NH_dat2 <- NH_dat %>% mutate(predicted_stdev = exp(NH_mlr_summ$coefficients["(Intercept)","Estimate"] +
                                                     #(NH_mlr_summ$coefficients["log(nsamples)","Estimate"]*(log(nsamples))) +
                                                     (NH_mlr_summ$coefficients["log(gradient)","Estimate"]*(log(gradient))) +
                                                     (NH_mlr_summ$coefficients["mean_d18O","Estimate"]*mean_d18O)))

SH_dat2 <- SH_dat %>% mutate(predicted_stdev = exp(SH_mlr_summ$coefficients["(Intercept)","Estimate"] +
                                                     #(SH_mlr_summ$coefficients["log(nsamples)","Estimate"]*(log(nsamples))) +
                                                     (SH_mlr_summ$coefficients["log(gradient)","Estimate"]*(log(gradient))) +
                                                     (SH_mlr_summ$coefficients["mean_d18O","Estimate"]*mean_d18O)))

## combine NH and SH data
dat_all <- rbind(NH_dat2, SH_dat2)

## apply correction
dat_all <- dat_all %>% mutate(corrected_stdev = stdev - predicted_stdev)

# regional composite
comp <- dat_all %>% group_by(region, win_end) %>% summarise(med = median(corrected_stdev, na.rm = T), Q1 = quantile(corrected_stdev, na.rm = T)[2], Q3 = quantile(corrected_stdev, na.rm = T)[4])
comp$age <- comp$win_end - bin_hw
write.csv(comp, "spel_IADV/spel_sd.csv", row.names = F) #output

# plot
comp <- comp %>% filter(! region %in% c("NE-SAM","SAfM","CAM","IAM")); 
comp$grp <- "corrected"

# combine corrected and uncorrected stdev
comp_raw <- df_sitedat %>% 
  filter(entity_id %in% unique(na.omit(dat_all)$entity_id)) %>%
  group_by(region, win_end) %>%
  summarise(med = median(stdev, na.rm = T), Q1 = quantile(stdev, na.rm = T)[2], Q3 = quantile(stdev, na.rm = T)[4]) %>% 
  filter(! region %in% c("NE-SAM","SAfM","CAM","IAM"))
comp_raw$grp = "uncorrected"

raw_corr_comp <- rbind(comp, comp_raw)
raw_corr_comp[which(raw_corr_comp$region == "SW-SAM"),"region"] <- "SAM"

#variable order
raw_corr_comp$grp <- factor(raw_corr_comp$grp, c("uncorrected","corrected"))

raw_corr_comp$win_end <- raw_corr_comp$win_end/1000

png("raw_v_corrected_fig3.png", width = 13, height = 8.5, units = "cm", res = 96)
ggplot() + 
  geom_ribbon(data = raw_corr_comp, aes(x = win_end, ymin = Q1, ymax = Q3),col = NA, fill = "#AAAAAA") +
  geom_line(data = raw_corr_comp, aes(x = win_end, y = med), col = "#464646") +
  scale_x_continuous(breaks = seq(0,12,3)) + 
  scale_y_continuous(breaks = seq(0,1.25,0.25)) +
  ylab(expression(paste(delta^{18}, "O s.d. (\u2030)"))) + xlab("Kyrs BP") +
  facet_grid(grp ~ region, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = NA))
dev.off()


## linear regression through composites (table 2)
raw_corr_comp %>%
  filter(win_end <= 6 & win_end >= 0) %>%
  group_by(region) %>%
  summarise(gradient = summary(lm(med ~ win_end))$coefficients[2,1],
            st_err = summary(lm(med ~ win_end))$coefficients[2,2],
            P_val = summary(lm(med ~ win_end))$coefficients[2,4])
