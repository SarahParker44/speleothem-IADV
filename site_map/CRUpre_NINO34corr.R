## processing CRU precip data
setwd("speleothem-IADV/")

library(ncdf4)
library(tidyr)
library(lubridate)

ncin <- nc_open("cru_ts4.05.1901.2020.pre.dat.nc")

lon <- ncvar_get(ncin, "lon"); nlon <- length(lon)
lat <- ncvar_get(ncin, "lat"); nlat <- length(lat)
time <- ncvar_get(ncin, "time")
time <- as.Date(time, origin = "1900-01-01")
x <- which(time > "1950-01-01"); tdim <- length(x)

pre <- ncvar_get(ncin, "pre", start = c(1,1,x[1]), count = c(nlon, nlat, tdim))

# subset to (extra)tropics
lat_lims <- sapply(c(-40,45), function(x) which.min(abs(lat-x)))
pre <- pre[,lat_lims[1]:lat_lims[2],]
lat <- lat[lat_lims[1]:lat_lims[2]]


## load NINO 3.4 SSTs
NINO <- read.table("nina34.data.txt") # downloaded 13/04/21

# subset to 1950 -> 2020
NINO <- NINO[which(NINO[,1] >= 1950 & NINO[,1] <= 2020),]

# reshape NINO 
colnames(NINO) <- c("Year",1:12)
NINO <- NINO %>% gather(key = "Month", value = "Nino34", 2:13)
NINO$date <- make_date(year = NINO$Year, month = NINO$Month)
NINO <- NINO %>% dplyr::arrange(date)


### correlate for NH summer
# subset to NH
pre_NH <- pre[,which(lat >= 0)[1]:length(lat),]

cor_mat <- matrix(data = NA, nrow = dim(pre_NH)[1], ncol = dim(pre_NH)[2]) #lon x lat matrix for corr vals
sig_mat <- matrix(data = NA, nrow = dim(pre_NH)[1], ncol = dim(pre_NH)[2]) #lon x lat matrix for P vals

for (i in 1:dim(pre_NH)[1]){
  for (j in 1:dim(pre_NH)[2]){ #for every NH grid
    # subset to grid
    datsub <- pre_NH[i,j,]
    
    if (all(is.na(datsub) == T)){ next } #ocean grids
    if (any(is.na(datsub)) == T){ dry_grids <- rbind(dry_grids, data.frame(lon_pos = i, lat_pos = j))}
    # select summer months for each year, calc mean precip and NINO3.4
    summ_mean <- numeric()
    summ_mean_NINO <- numeric()
    val <- 1
    for (k in 1:(tdim/12)){
      summ_sub <- datsub[(val+4):(val+8)]
      summ_sub_NINO <- NINO[(val+4):(val+8),3]
      pre_mean <- mean(summ_sub)
      NINO_mean <- mean(summ_sub_NINO)
        
      summ_mean[k] <- pre_mean
      summ_mean_NINO[k] <- NINO_mean
        
      val <- val+12
    }
    # correlation
    cor_mat[i,j] <- cor.test(summ_mean, summ_mean_NINO)$estimate
    sig_mat[i,j] <- cor.test(summ_mean, summ_mean_NINO)$p.val
    
  }
}

write.csv(cor_mat, "site_map/MJJAS_precipNINO_corr_coeff.csv", row.names = F)
write.csv(sig_mat, "site_map/MJJAS_precipNINO_corr_pval.csv", row.names = F)


## same for NDJFM
# subset to NH
pre_SH <- pre[,1:max(which(lat <= 0)),]

cor_mat <- matrix(data = NA, nrow = dim(pre_SH)[1], ncol = dim(pre_SH)[2])
sig_mat <- matrix(data = NA, nrow = dim(pre_SH)[1], ncol = dim(pre_SH)[2])
for (i in 1:dim(pre_SH)[1]){
  for (j in 1:dim(pre_SH)[2]){
    # subset to grid
    datsub <- pre_SH[i,j,]
    
    if (all(is.na(datsub) == T)){ next }
    # select summer months for each year, calc mean precip and NINO3.4
    summ_mean <- numeric()
    summ_mean_NINO <- numeric()
    val <- 1
    for (k in 1:((tdim/12)-2)){
      summ_sub <- datsub[(val+10):(val+14)]
      summ_sub_NINO <- NINO[(val+10):(val+14),3]
      pre_mean <- mean(summ_sub)
      NINO_mean <- mean(summ_sub_NINO)
      
      summ_mean[k] <- pre_mean
      summ_mean_NINO[k] <- NINO_mean
      
      val <- val+12
    }
    # correlation
    cor_mat[i,j] <- cor.test(summ_mean, summ_mean_NINO)$estimate
    sig_mat[i,j] <- cor.test(summ_mean, summ_mean_NINO)$p.val
    
  }
}

write.csv(cor_mat, "site_map/NDJFM_precipNINO_corr_coeff.csv", row.names = F)
write.csv(sig_mat, "site_map/NDJFM_precipNINO_corr_pval.csv", row.names = F)
