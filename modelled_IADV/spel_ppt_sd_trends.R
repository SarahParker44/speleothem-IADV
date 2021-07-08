#### plot of precip versus spel d18O s.d. trends

setwd("speleothem-IADV/")

library(ggplot2)
library(dplyr)
library(ggthemes)
library(cowplot)

## load spel d18O s.d. trends
spel_sd <- read.csv("spel_IADV/spel_sd.csv")
spel_sd <- spel_sd %>% filter(win_end <= 8000 & win_end >=50 & region %in% c("EAM","ISM","SW-SAM"))

## load precip s.d. trends
model_sd <- read.csv("modelled_IADV/region_pre_sd.csv")
# anomalise
anom_df <- data.frame()
for (i in unique(model_sd$region)){
  for(j in unique(model_sd$Model)){
    subdat <- model_sd %>% filter(region == i & Model == j)
    base_val <- subdat %>% filter(Age >=0 & Age <= 1000) %>% summarise(base_mean = mean(precip_sd, na.rm = T))
    subdat$anom <- (subdat$precip_sd - base_val$base_mean)
    anom_df <- rbind(anom_df, subdat)
  }
}

anom_df %>%
  filter(Age <= 6000) %>%
  group_by(region, Model) %>%
  summarise(slope = summary(lm(anom ~ Age))$coefficients[2,1],
            summary(lm(anom ~ Age))$coefficients[2,4])

## apply smoothing curve
loess_window = 500

#spel
spel_smooth <- data.frame()
for (i in unique(spel_sd$region)){
  subdat <- spel_sd %>% filter(region == i)
  span <- (loess_window/50)/nrow(subdat)
  
  loess_mod <- loess(med ~ age, data = subdat, span = span)
  loess_out <- predict(loess_mod, se = T)
  
  sub_df <- data.frame(region = i,
                       Age = subdat$age,
                       fit = loess_out$fit,
                       SE1 = loess_out$fit -loess_out$se.fit,
                       SE2 = loess_out$fit +loess_out$se.fit)
  spel_smooth <- rbind(spel_smooth, sub_df)
}

#precip
mod_smooth <- data.frame()
for (i in unique(anom_df$Model)){
  for (j in unique(anom_df$region)){
    subdat <- anom_df %>% filter(Model == i & region == j)
    span <- (loess_window/50)/nrow(subdat)
    
    loess_mod <- loess(anom ~ Age, data = subdat, span = span)
    loess_out <- predict(loess_mod, se = T)
    
    sub_df <- data.frame(Model = i,
                         region = j,
                         Age = subdat$Age,
                         fit = loess_out$fit,
                         SE1 = loess_out$fit -loess_out$se.fit,
                         SE2 = loess_out$fit +loess_out$se.fit)
    mod_smooth <- rbind(mod_smooth, sub_df)
  }
}




# 1) EAM precip
p1 <- ggplot() +
  
  geom_ribbon(data = filter(mod_smooth, region == "EAM"), aes(x = Age, ymin = SE1, ymax = SE2, fill = Model), alpha = 0.5) +
  geom_line(data = filter(mod_smooth, region == "EAM"), aes(x = Age, y = fit, col = Model)) +
  
  expand_limits(x = c(0,8000), y= c(-0.26,0.13)) +
  scale_x_continuous(breaks = seq(0,8000,2000), expand = c(0.01,0.01)) +
  ylab("EAM ppt s.d.") +
  theme_bw() +
  theme(legend.position = "none", 
        panel.border = element_rect(colour = NA),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.y = element_line(colour = "black"))

# 2) EAM spel d18O
p2 <- ggplot() +
  
  geom_ribbon(data = filter(spel_smooth, region == "EAM"), mapping = aes(x = Age, ymin = SE1, ymax = SE2), fill = "grey", alpha = 0.5) +
  geom_line(data = filter(spel_smooth, region == "EAM"), mapping = aes(x = Age, y = fit), col = "dark grey") +
  
  scale_y_continuous(position = "right", breaks = seq(-0.1,0.1,0.1)) +
  expand_limits(x = c(0,8000)#, y = c(-0.26,0.4)
                ) +
  scale_x_continuous(breaks = seq(0,8000,2000), expand = c(0.01,0.01)) +
  ylab(expression(paste("EAM ", delta^18,"O s.d.", sep = ""))) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.border = element_rect(colour = NA),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_line(colour = "black"))

# 3) ISM precip
p3 <- ggplot() +
  
  geom_ribbon(data = filter(mod_smooth, region == "ISM"), aes(x = Age, ymin = SE1, ymax = SE2, fill = Model), alpha = 0.5) +
  geom_line(data = filter(mod_smooth, region == "ISM"), aes(x = Age, y = fit, col = Model)) +
  
  expand_limits(x = c(0,8000), y = c(-0.26,0.13)
                ) +
  scale_x_continuous(breaks = seq(0,8000,2000), expand = c(0.01,0.01)) +
  ylab("ISM ppt s.d.") +
  theme_bw() +
  theme(legend.position = "none", 
        panel.border = element_rect(colour = NA),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_line(colour = "black"))

# 4) ISM d18O
p4 <- ggplot()+
  
  geom_ribbon(data = filter(spel_smooth, region == "ISM"), mapping = aes(x = Age, ymin = SE1, ymax = SE2), fill = "grey", alpha = 0.5) +
  geom_line(data = filter(spel_smooth, region == "ISM"), mapping = aes(x = Age, y = fit), col = "dark grey") +
  
  
  expand_limits(x = c(0,8000)#, y = c(-0.26,0.4)
                ) +
  scale_y_continuous(position = "right") +
  scale_x_continuous(breaks = seq(0,8000,2000), expand = c(0.01,0.01)) +
  ylab(expression(paste("ISM ", delta^18,"O s.d.", sep = ""))) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.border = element_rect(colour = NA),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.y = element_line(colour = "black"))


# 5) SAM ppt
p5 <- ggplot() +
  
  geom_ribbon(data = filter(mod_smooth, region == "SAM"), aes(x = Age, ymin = SE1, ymax = SE2, fill = Model), alpha = 0.5) +
  geom_line(data = filter(mod_smooth, region == "SAM"), aes(x = Age, y = fit, col = Model)) +
  
  expand_limits(x = c(0,8000), y = c(-0.26,0.13)) +
  scale_x_continuous(breaks = seq(0,8000,2000), expand = c(0.01,0.01)) +
  ylab("SAM ppt s.d.") +
  theme_bw() +
  theme(legend.position = "none", 
        panel.border = element_rect(colour = NA),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_line(colour = "black"))

# 6) SAM d18O
p6 <- ggplot() +
  
  geom_ribbon(data = filter(spel_smooth, region == "SW-SAM"), mapping = aes(x = Age, ymin = SE1, ymax = SE2), fill = "grey", alpha = 0.5) +
  geom_line(data = filter(spel_smooth, region == "SW-SAM"), mapping = aes(x = Age, y = fit), col = "dark grey") +
  
  expand_limits(x = c(0,8000)) +
  scale_x_continuous(breaks = seq(0,8000,2000), expand = c(0.01,0.01)) +
  scale_y_continuous(position = "right", breaks = seq(0,0.4,0.1)) +
  ylab(expression(paste(delta^18,"O s.d.", sep = ""))) +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        panel.border = element_rect(colour = NA),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"))

#legend
px <- ggplot() +
  geom_ribbon(data = filter(mod_smooth, region == "EAM"), aes(x = Age, ymin = SE1, ymax = SE2, fill = Model), alpha = 0.5) +
  geom_line(data = filter(mod_smooth, region == "EAM"), aes(x = Age, y = fit, col = Model)) +
  theme_bw() +
  theme(legend.direction = "horizontal", legend.title = element_blank())

grobs <- ggplotGrob(px)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

pdf("legend_spel_pre_sd.pdf", width = 6.5/2.54, height = 1/2.54)
plot_grid(legend)
dev.off()

#plot and save

rel_height <- spel_smooth %>% group_by(region) %>% summarise(min = min(SE1, na.rm = T), max = max(SE2, na.rm = T), range = max(SE2, na.rm = T)-min(SE1, na.rm = T))

pdf("spel_vs_model_sd.pdf", width = 10/2.54, height = 20/2.54)
plot_grid(p1,p2,p3,p4,p5,p6, ncol = 1, align = "v", 
          rel_heights = c(0.8,
                          as.numeric(rel_height[which(rel_height$region == "EAM"),"range"]/rel_height[which(rel_height$region == "SW-SAM"),"range"]),
                          0.8,
                          as.numeric(rel_height[which(rel_height$region == "ISM"),"range"]/rel_height[which(rel_height$region == "SW-SAM"),"range"]),
                          0.8,
                          1))
dev.off()
