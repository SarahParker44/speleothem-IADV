#### Extract regional s.d. trends from Pacmedy runs
## requires "pr_IPSL_Sr.nc", "pr_IPSL_Vlr.nc", "slm_IPSL_Sr.nc", "slm_IPSL_Vlr.nc" files
## requires "pr_slm_AWI.RData" and "pr_slm_MPI.RData" generated in "modelled_IADV/reshape_MPI_AWI.R"

library(ncdf4)
library(ggplot2)

setwd("speleothem-IADV/")


## 1) Load Pacmedy data and 2) apply land mask

ptm <- proc.time()
for (i in c("IPSL_Sr","IPSL_Vlr","MPI","AWI")){
  if (i %in% c("IPSL_Sr","IPSL_Vlr")){
    filename <- paste("R_programming/Data/pr_",i,".nc", sep = "")
    ncin <- nc_open(filename)
    lon <- ncvar_get(ncin, "lon"); lat <- ncvar_get(ncin, "lat")
    time <- ncvar_get(ncin, "time_counter"); tdim = length(time)
    
    latvals <- sapply(c(45,-40), function(x) which.min(abs(lat - x))) #load tropics and extratropics only
    lonvals <- sapply(c(-85,130), function(x) which.min(abs(lon - x)))
    
    lat <- lat[latvals[1]:latvals[2]]; nlat <- length(lat)
    lon <- lon[lonvals[1]:lonvals[2]]; nlon <- length(lon)
    
    precip <- ncvar_get(ncin, "precip", start = c(lonvals[1],latvals[1],1), count = c(nlon,nlat,tdim))
    
    nc_close(ncin); rm(ncin)
    
    
    filename <- paste("R_programming/Data/slm_",i,".nc", sep = "") #load landmask
    ncin <- nc_open(filename)
    
    landfrac <- ncvar_get(ncin, "fract_ter", start = c(lonvals[1],latvals[1],1), count = c(nlon,nlat,1))
    landmask <- ifelse(landfrac > 0.5, 1, NA)
    
    nc_close(ncin); rm(ncin)
  } else {
    filename <- paste("Data/pr_slm_", i, ".RData", sep = "") #load preprocessed MPI and AWI pre and landmask data
    load(filename)
    landmask <- ifelse(landfrac ==1, 1, NA)
  }
  
  precip_land <- array(data=NA, dim = dim(precip)) #ocean grids = NA
  for (j in 1:(dim(precip)[3])){
    precip_land[,,j] <- precip[,,j]*landmask
  }
  
  
  ## output
  if (i == "IPSL_Sr"){
    precip_Sr <- precip_land
    lon_Sr <- lon; lat_Sr <- lat
  } else if (i == "IPSL_Vlr") {
    precip_Vlr <- precip_land
    lon_Vlr <- lon; lat_Vlr <- lat
  } else if (i == "MPI"){
    precip_MPI <- precip_land
    lon_MPI <- lon; lat_MPI <- lat
  } else {
    precip_AWI <- precip_land
    lon_AWI <- lon; lat_AWI <- lat
  }
}
proc.time() - ptm #8 minutes

rm(i, j, filename, lon, lat, lonvals, latvals, precip, landfrac, landmask, precip_land, nlon, nlat, tdim, time, varname, time_varname, pre_varname)

## 3) kg/(s*m2) -> mm/d
precip_Sr <- precip_Sr *86400
precip_Vlr <- precip_Vlr *86400
precip_MPI <- precip_MPI *86400
precip_AWI <- precip_AWI *86400


## 4) Regional area-average
regions <- rbind(data.frame(region = "ISM", min_lon = 75, max_lon = 95, min_lat = 10, max_lat = 35),
                 data.frame(region = "EAM", min_lon = 100, max_lon = 127, min_lat = 20, max_lat = 43),
                 data.frame(region = "SAM1", min_lon = -80, max_lon = -64, min_lat = -10, max_lat = 0),
                 data.frame(region = "SAM2", min_lon = -68, max_lon = -40, min_lat = -24, max_lat = -10))


df_all <- data.frame()
for (i in 1:2){ # EAM and ISM
  for (j in c("IPSL_Sr","IPSL_Vlr","MPI","AWI")){ #each model sim
    #select model lon and lat vals, and precip
    if (j == "IPSL_Sr"){
      lon <- lon_Sr; lat <- lat_Sr
      precip <- precip_Sr
    } else if (j == "IPSL_Vlr"){
      lon <- lon_Vlr; lat <- lat_Vlr
      precip <- precip_Vlr
    } else if (j == "MPI"){
      lon <- lon_MPI; lat <- lat_MPI
      precip <- precip_MPI
    } else {
      lon <- lon_AWI; lat <- lat_AWI
      precip <- precip_AWI
    }
    tdim_yrs <- dim(precip)[3]/12 #length of model sim in yrs
    
    # subset to region
    lonlims <- sapply(regions[i,2:3], function(x) which.min(abs(lon - x)))
    latlims <- sapply(regions[i,4:5], function(x) which.min(abs(lat - x)))
    
    regsub <- precip[lonlims[1]:lonlims[2],latlims[2]:latlims[1],]
    
    # area average
    regsub <- apply(regsub, 3, mean, na.rm = T)
    
    # summ vals
    subdat_summ <- regsub[which(rep(1:12,tdim_yrs) %in% c(5:9))] #subset to MJJAS months
      
    summ <- c()
    val <- 1
    for (k in 1:tdim_yrs){
      summ[k] <- mean(subdat_summ[val:(val+4)]) #MJJAs mean
      val <- val+5
    }
      
    # Calculate sliding window s.d. and mean 
    sd_summ <- c()
    mean_summ <- c()
    val <- 1
    for (k in 1:((tdim_yrs/50)-1)){ #each bin
      sd_summ[k] <- sd(summ[val:(val+99)])
      mean_summ[k] <- mean(summ[val:(val+99)])
      val <- val+50
    }

    if (j == "MPI"){ time_yrs <- rev(seq(150,7900,50)) } else { 
      time_yrs <- rev(seq(50,5950,50)) }
    sub_df <- data.frame(region = regions$region[i],
                         Age = time_yrs,
                         Model = j,
                         precip_sd = sd_summ,
                         precip_mean = mean_summ)
    
    df_all <- rbind(df_all, sub_df)
  }
}

## SAM
df_all_SAM <- data.frame()
for (j in c("IPSL_Sr","IPSL_Vlr","MPI","AWI")){
  #select model lon and lat vals, and precip
  if (j == "IPSL_Sr"){
    lon <- lon_Sr; lat <- lat_Sr
    precip <- precip_Sr
  } else if (j == "IPSL_Vlr"){
    lon <- lon_Vlr; lat <- lat_Vlr
    precip <- precip_Vlr
  } else if (j == "MPI"){
    lon <- lon_MPI; lat <- lat_MPI
    precip <- precip_MPI
  } else {
    lon <- lon_AWI; lat <- lat_AWI
    precip <- precip_AWI
  }
  tdim_yrs <- dim(precip)[3]/12
  
  # subset to region
  lonlims1 <- sapply(regions[3,2:3], function(x) which.min(abs(lon - x)))
  latlims1 <- sapply(regions[3,4:5], function(x) which.min(abs(lat - x)))
  
  lonlims2 <- sapply(regions[4,2:3], function(x) which.min(abs(lon - x)))
  latlims2 <- sapply(regions[4,4:5], function(x) which.min(abs(lat - x)))
  
  regsub1 <- precip[lonlims1[1]:lonlims1[2],latlims1[2]:latlims1[1],]
  regsub2 <- precip[lonlims2[1]:lonlims2[2],latlims2[2]:latlims2[1],]
  
  #combine sub regions into 1 region
  regsub1 <- apply(regsub1, 3, function(x) as.numeric(x))
  regsub2 <- apply(regsub2, 3, function(x) as.numeric(x))
  regsub <- rbind(regsub1, regsub2)
  
  # area average
  regsub <- apply(regsub, 2, mean, na.rm = T)
  
  month_vec <- rep(1:12,tdim_yrs)
  #incomplete_summ_months <- c(1:3,(length(month_vec)-1):length(month_vec))
  subdat_summ <- regsub[which(month_vec %in% c(1:3,11:12))]#[-incomplete_summ_months]
  #remove incomplete months
  subdat_summ <- subdat_summ[-c(1:2)] # Nov and Dec of the final year
  subdat_summ <- subdat_summ[-c((length(subdat_summ)-2):length(subdat_summ))] # Jan->Mar of the first year
  
  summ <- c()
  val <- 1
  for (k in 1:(tdim_yrs - 1)){
    summ[k] <- mean(subdat_summ[val:(val+4)])
    val <- val+5
  }
  
  # Calculate sliding window s.d.
  sd_summ <- c()
  mean_summ <- c()
  val <- 1
  for (k in 1:ceiling((length(summ)/50)-1)){
    if (k == 1){
      sd_summ[k] <- sd(summ[val:(val+98)])
      mean_summ[k] <- mean(summ[val:(val+98)])
    } else {
      sd_summ[k] <- sd(summ[(val-1):(val+98)])
      mean_summ[k] <- mean(summ[(val-1):(val+98)])
    }
    val <- val+50
  }
  
  if (j == "MPI"){ time_yrs <- rev(seq(150,7900,50)) } else { 
    time_yrs <- rev(seq(50,5950,50)) }
  sub_df <- data.frame(region = "SAM",
                       Age = time_yrs,
                       Model = j,
                       precip_sd = sd_summ,
                       precip_mean = mean_summ)
  
  df_all_SAM <- rbind(df_all_SAM, sub_df)
}

df <- rbind(df_all, df_all_SAM)


write.csv(df, "modelled_IADV/region_pre_sd.csv", row.names = F)
