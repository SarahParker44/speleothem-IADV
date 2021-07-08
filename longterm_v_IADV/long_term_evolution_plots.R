### plot long term evolution trends

library(dplyr)
library(ggplot2)
library(ggpubr)
library(palinsol)

setwd("speleothem-IADV/")

## load regional spel trends
spel <- read.csv("longterm_v_IADV/d18O_comp_zscore.csv")
spel2 <- read.csv("longterm_v_IADV/d18O_bin_zscore.csv")

## load precip trends
mod_pre <- read.csv("modelled_IADV/region_pre_sd.csv")

# anomalise
anom_df <- data.frame()
for (i in unique(mod_pre$region)){
  for(j in unique(mod_pre$Model)){
    subdat <- mod_pre %>% filter(region == i & Model == j)
    base_val <- subdat %>% filter(Age >=0 & Age <= 1000) %>% summarise(base_sd = mean(precip_sd, na.rm = T),
                                                                       base_mean = mean(precip_mean, na.rm = T))
    subdat$sd_anom <- subdat$precip_sd - base_val$base_sd
    subdat$mean_anom <- subdat$precip_mean - base_val$base_mean
    anom_df <- rbind(anom_df, subdat)
  }
}

anom_df <- anom_df %>% filter(Age <= 5950 & Age >= 50)

# summer insolation data

# Function that calculates mean insolation for a specified lat and sequence of months
get_insol <- function(Lat, Months){
  mid_month_days_0ka <- seq(15.5, 345.5, by=30)
  tt_present = 0.0
  orbit_present <- astro(tt_present, ber78, degree = FALSE) # present-day atronomical parameters (using Berger, 1978)
  mid_month_tsl_0ka <- day2l(orbit_present, mid_month_days_0ka) #convert days to true solar longitude
  time_BP <- seq(-6000,0,by=50) # Holocene time intervals
  orbital_params <- data.frame(time_BP, t(sapply(time_BP, function(tt) astro(tt, ber78, degree=FALSE)))) # Holocene astro params
  insol_month <- matrix(0, nrow=length(time_BP), ncol=12)
  for (month in seq(1:12)) {
    tsl <- mid_month_tsl_0ka[month] # month solar longitude
    insol_month[,month] <- as.numeric(Insol(orbital_params, long=tsl, lat=Lat*pi/180, S0=1365)) # calc insolation (convert lat to radians) 
  }
  Months_insol <- data.frame(insol_month[,Months]) # filter to specified months
  Months_insol2 <- data.frame(AGE = time_BP*-1, insol = apply(Months_insol, 1, mean)) # mean insolation across months
}

# Calculate Holocene summer insolation for NH and SH
NH_insol <- get_insol(Lat = 30, Month = c(5,6,7,8,9)) # mean MJJAS insolation (NH summer)
SH_insol <- get_insol(Lat = -20, Month = c(11,12,1,2,3)) # mean NDJFM insolation (SH summer)



## plot

# colours
mod_palette <- c("#648FFF","#DC267F","#FE6100","#FFB000")

# NH isol
NH_insol <- NH_insol %>% filter(AGE >= 50 & AGE <= 5950)
p1 <- ggplot(data = NH_insol, aes(x = AGE, y = insol)) +
  geom_line(col = "grey") +
  
  scale_y_continuous(position = "right") +
  scale_x_continuous(breaks = seq(0,6000,1000), expand = c(0.01,0.01)) +
  
  ylab("30°N insolation") +
  
  theme_bw() +
  theme(panel.border = element_rect(colour = NA),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_line(colour = "black"))


# EAM
p2 <- ggplot(data = filter(spel, region == "EAM"), aes(x = age)) +
  geom_ribbon(aes(ymin = conf5, ymax = conf95), fill = "#AAAAAA") +
  geom_line(aes(y = locfit)) +
  geom_line(data = filter(spel2, region == "EAM"), aes(x = age, y = mean_zscore), col = "#8B0000") +
  
  scale_y_reverse() +
  scale_x_continuous(breaks = seq(0,6000,1000), expand = c(0.01,0.01)) +
  
  ylab("spel z-score") +
  
  theme_bw() +
  theme(panel.border = element_rect(colour = NA),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_line(colour = "black"))

p3 <- ggplot(data = filter(anom_df, region =="EAM"), aes(x = Age, y = mean_anom, col = Model)) +
  geom_line() +
  
  scale_color_manual(values = mod_palette) +
  scale_x_continuous(breaks = seq(0,6000,1000), expand = c(0.01,0.01)) +
  scale_y_continuous(position = "right") +
  ylab("mean precip (mm/d)") +
  
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_rect(colour = NA),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.y = element_line(colour = "black"))

# ISM

p4 <- ggplot(data = filter(spel, region == "ISM"), aes(x = age)) +
  geom_ribbon(aes(ymin = conf5, ymax = conf95), fill = "#AAAAAA") +
  geom_line(aes(y = locfit)) +
  geom_line(data = filter(spel2, region == "ISM"), aes(x = age, y = mean_zscore), col = "#8B0000") +
  
  scale_y_reverse() +
  scale_x_continuous(breaks = seq(0,6000,1000), expand = c(0.01,0.01)) +
  
  ylab("spel z-score") +
  
  theme_bw() +
  theme(panel.border = element_rect(colour = NA),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_line(colour = "black"))


p5 <- ggplot(data = filter(anom_df, region =="ISM"), aes(x = Age, y = mean_anom, col = Model)) +
  geom_line() +
  
  scale_color_manual(values = mod_palette) +
  scale_x_continuous(breaks = seq(0,6000,1000), expand = c(0.01,0.01)) +
  scale_y_continuous(position = "right") +
  ylab("mean precip (mm/d)") +
  
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_rect(colour = NA),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.y = element_line(colour = "black"))

# SH isol
SH_insol <- SH_insol %>% filter(AGE >= 50 & AGE <= 5950)
p6 <- ggplot(data = SH_insol, aes(x = AGE, y = insol)) +
  geom_line(col = "grey") +
  
  scale_y_continuous(position = "right") +
  scale_x_continuous(breaks = seq(0,6000,1000), expand = c(0.01,0.01)) +
  
  ylab("30°N insolation") +
  
  theme_bw() +
  theme(panel.border = element_rect(colour = NA),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_line(colour = "black"))
# SAM

p7 <- ggplot(data = filter(spel, region == "SAM"), aes(x = age)) +
  geom_ribbon(aes(ymin = conf5, ymax = conf95), fill = "#AAAAAA") +
  geom_line(aes(y = locfit)) +
  geom_line(data = filter(spel2, region == "SAM"), aes(x = age, y = mean_zscore), col = "#8B0000") +
  
  scale_y_reverse() +
  scale_x_continuous(breaks = seq(0,6000,1000), expand = c(0.01,0.01)) +
  
  ylab("spel z-score") +
  
  theme_bw() +
  theme(panel.border = element_rect(colour = NA),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_line(colour = "black"))


p8 <- ggplot(data = filter(anom_df, region =="SAM"), aes(x = Age, y = mean_anom, col = Model)) +
  geom_line() +
  
  scale_color_manual(values = mod_palette) +
  scale_x_continuous(breaks = seq(0,6000,1000), expand = c(0.01,0.01)) +
  scale_y_continuous(position = "right", breaks = c(0,-0.4,-0.8)) +
  ylab("mean precip (mm/d)") +
  
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_rect(colour = NA),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"))


# multiplot
pdf("lt_evolution.pdf", height = 15/2.54, width = 15/2.54)
ggarrange(p1,p6,p2,p7,p3,p8,p4,p5, ncol = 2, nrow = 4, align = "hv", common.legend = T, legend = "bottom")
dev.off()
