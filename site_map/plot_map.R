#### Plot speleothem sites with summer precip - NINO3.4 correlation coefficients ##

setwd("speleothem-IADV/")

library(ncdf4)
library(rgdal)
library(ggplot2)
library(colorspace)


## Load NINO3.4 - CRU precip corr vals

MJJAS <- read.csv("site_map/MJJAS_precipNINO_corr_coeff.csv") 
NDJFM <- read.csv("site_map/NDJFM_precipNINO_corr_coeff.csv") 

corr_mat <- cbind(NDJFM, MJJAS)

## to df with vars = lon, lat, corr_coeff
# lon and lat vals
ncin <- nc_open("cru_ts4.04.1901.2019.pre.dat.nc")
lon <- ncvar_get(ncin, "lon"); lat <- ncvar_get(ncin, "lat")
lat_lims <- sapply(c(-40,45), function(x) which.min(abs(lat-x)))
lat <- lat[lat_lims[1]:lat_lims[2]]

corr_df <- data.frame(expand.grid(lon = lon, lat = lat),
                      corr_coeff = as.numeric(as.matrix(corr_mat)))

corr_df <- na.exclude(corr_df)

## load speleothem sites
spel_sites <- read.csv("speleothem_monsoon_interannual_var_repo/spel_stdev/sites_stdev_analysis.csv")


## load coast lines
wmap <- readOGR(dsn = "ne_110m_land", layer = "ne_110m_land")
wmap@data$id <- rownames(wmap@data)
worldMap <- fortify(wmap)
wmap_DF <- merge(worldMap, wmap@data, by = "id")

## plot map

p <- ggplot() + 
  geom_tile(data = corr_df, aes(x = lon, y = lat, fill = corr_coeff)) + #correlation coeffs
  geom_polygon(data = wmap, aes(x = long, y = lat, group = id), fill = NA, col = "#616161", size = 0.1) + #world map
  geom_point(data = spel_sites, aes(x = longitude, y = latitude), shape = 21) + #spel sites
  geom_rect(mapping = aes(xmin = c(75,100), xmax = c(95,127), ymin = c(10,20), ymax = c(35,42)), fill = NA, col = "#464646") + #ISM,EAM reg limits
  geom_polygon(mapping = aes(x = c(-80,-64,-64,-40,-40,-68,-68,-80,-80), y = c(0,0,-10,-10,-24,-24,-10,-10,0)), fill = NA, col = "#464646") + #SAM reg lims
  scale_fill_binned_diverging("Blue-Red", breaks = c(-0.5,-0.2,0.2,0.5), guide = guide_colorbar(raster = F, frame.colour = "grey")) + #edit colours + colourbar
  scale_x_continuous(breaks = seq(-120,140,20), expand = c(0,0)) + scale_y_continuous(breaks = seq(-40,40,20), expand = c(0,0)) + #edit axes
  coord_fixed(ylim = c(-40,45), xlim = c(-140,160)) + #crop map
  theme(axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_line(colour = "#616161"),
        panel.border = element_rect(colour = "#616161", fill = NA),
        panel.background = element_rect(fill = "white"),
        legend.title = element_blank())

## save map

pdf("C:/Users/sarah/OneDrive/Documents/PhD/IAV_paper/fig1_map.pdf", width = 15/2.54, height = 5/2.54)
print(p)
dev.off()
