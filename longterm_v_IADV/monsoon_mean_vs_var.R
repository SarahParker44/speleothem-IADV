### Mean monsoon versus s.d. monsoon

library(dplyr)
library(ggplot2)
library(colorspace)
library(cowplot)

setwd("speleothem-IADV/")

# modelled precip
mod_pre <- read.csv("model_IADV/region_pre_sd.csv")
mod_pre$grp <- paste(mod_pre$region, mod_pre$Model)

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

anom_df <- anom_df %>% filter(Age <= 8000)

## spel data
spel_sd <- read.csv("spel_IADV/spel_sd.csv") 
spel_mean <- read.csv("spel_IADV/d18O_bin_zscore.csv")
spel_sd[which(spel_sd$region == "SW-SAM"),"region"] <- "SAM"

spel_all <- left_join(spel_mean, spel_sd)
spel_all <- spel_all %>% filter(age <= 8000)


## colours
mod_palette <- c("#648FFF","#DC267F","#FE6100","#FFB000")


## 1) EAM spel
p1 <- ggplot(data = filter(spel_all, region == "EAM"), aes(x = mean_zscore, y = med)) + 
  geom_point(col = "dark grey", size = 0.1) +
  scale_x_reverse() +
  coord_cartesian(ylim = c(-0.2,0.4)) +
  
  theme_bw() +
  theme(legend.position = "none") +
  xlab("") + ylab("") +
  
  geom_smooth(method = "lm", se = F, size = 0.5, col = "#464646")


p2 <- ggplot(data = filter(anom_df, region == "EAM"), aes(x = mean_anom, y = sd_anom, col = Model)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = mod_palette) +
  
  theme_bw() +
  theme(legend.position = "none") +
  xlab("") + ylab("") +
  scale_y_continuous(position = "right") +
  expand_limits(y = c(-0.17,0.11)) +
  
  geom_smooth(method = "lm", se = F, size = 0.5)
  


## 2) ISM spel

p3 <- ggplot(data = filter(spel_all, region == "ISM"), aes(x = mean_zscore, y = med)) +
  geom_point(col = "dark grey", size = 0.1) +
  
  scale_x_reverse() +
  coord_cartesian(ylim = c(-0.2,0.4)) +
  
  theme_bw() +
  theme(legend.position = "none") +
  xlab("") + ylab("spel s.d.") +
  
  geom_smooth(method = "lm", se = F, size = 0.5, col = "#464646")
 

 p4 <- ggplot(data = filter(anom_df, region == "ISM"), aes(x = mean_anom, y = sd_anom, col = Model)) +
   geom_point(size = 0.1) +
   scale_color_manual(values = mod_palette) +
   
   theme_bw() +
   theme(legend.position = "none") +
   xlab("") + ylab("precip s.d. (mm/d)") +
   scale_y_continuous(position = "right") +
   expand_limits(y = c(-0.17,0.11)) +
   
   geom_smooth(method = "lm", se = F, size = 0.5)
 


## 3) SAM spel
p5 <- ggplot(data = filter(spel_all, region == "SAM"), aes(x = mean_zscore, y = med)) +
  geom_point(col = "dark grey", size = 0.1) +
  scale_x_reverse() +
  #coord_cartesian(ylim = c(-0.2,0.4)) +
  
  theme_bw() +
  theme(legend.position = "none") +
  xlab("spel z-score") + ylab("") +
  
  geom_smooth(method = "lm", se = F, size = 0.5, col = "#464646")

p6 <- ggplot(data = filter(anom_df, region == "SAM"), aes(x = mean_anom, y = sd_anom, col = Model)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = mod_palette) +
  
  theme_bw() +
  theme(legend.position = "none") +
  xlab("mean precip (mm/d)") + ylab("") +
  scale_y_continuous(position = "right") +
  expand_limits(y = c(-0.17,0.11)) +
  
  geom_smooth(method = "lm", se = F, size = 0.5)

 

## multiplot
pdf("mean_vs_var_multiplot2.pdf", height = 5, width = 5)
ggarrange(p1,p2,p3,p4,p5,p6, ncol = 2, nrow = 3, common.legend = T, legend = "bottom")
dev.off()


## lm analysis
spel_lm <- data.frame()
for (i in unique(spel_all$region)){
  subdat <- spel_all %>% filter(region == i)
  sub_lm <- lm(med ~ mean_zscore, data = subdat)
  lm_summ <- summary(sub_lm)
  
  sub_out <- data.frame(region = i,
                        est = lm_summ$coefficients[2,1],
                        P_val = lm_summ$coefficients[2,4]
                        )
  
  spel_lm <- rbind(spel_lm, sub_out)
}

mod_lm <- data.frame()
for (i in unique(anom_df$region)){
  for(j in unique(anom_df$Model)){
    subdat <- anom_df %>% filter(region == i & Model == j)
    
    sub_lm <- lm(sd_anom ~ mean_anom, data = subdat)
    lm_summ <- summary(sub_lm)
    
    sub_out <- data.frame(region = i,Model = j,
                          est = lm_summ$coefficients[2,1],
                          P_val = lm_summ$coefficients[2,4])
    
    mod_lm <- rbind(mod_lm, sub_out)
  }
}
