### Reshape AWI and MPI so lat is -180->180, not 0->360
## files required: "pr_AWI.nc" and "pr_MPI.nc" 

library(ncdf4)

setwd("speleothem_IADV")


for (i in c("AWI","MPI")){
  # load precip
  filename <- paste("pr_",i,".nc", sep = "")
  ncin <- nc_open(filename)
  lon <- ncvar_get(ncin, "lon"); lat <- ncvar_get(ncin, "lat")
  
  pre_varname <- ifelse(i == "AWI", "var260", "precip") #precip var has different names
  precip <- ncvar_get(ncin, pre_varname) 
  
  nc_close(ncin); rm(ncin)
  
  
  # load land mask
  filename <- paste("slm_",i,".nc", sep = "")
  ncin <- nc_open(filename)
  
  varname <- ifelse(i == "MPI", "var172", "slm") #slm var has different names
  
  landfrac <- ncvar_get(ncin, varname, start = c(1,1,1), count = c(length(lon),length(lat),1))
  
  nc_close(ncin); rm(ncin)
  
  
  # reshape precip
  lon_matconv <- array(data = NA, dim = dim(precip))
  for (k in 1:dim(precip)[3]){
    a <- precip[which(lon > 180),,k]; b <- precip[which(lon <= 180),,k]
    lon_matconv[,,k] <- rbind(a,b)
  }
  precip <- lon_matconv
  
  # reshape mask
  a <- landfrac[which(lon > 180),]; b <- landfrac[which(lon <= 180),]
  landfrac <- rbind(a,b)
  
  # reshape lon
  a <- lon[which(lon > 180)]-360; b <- lon[which(lon <= 180)]
  lon <- c(a,b)
  
  
  ## subset to tropics and extratropics
  latvals <- sapply(c(46,-40), function(x) which.min(abs(lat - x)))
  lonvals <- sapply(c(-85,130), function(x) which.min(abs(lon - x)))
  
  precip <- precip[lonvals[1]:lonvals[2],latvals[1]:latvals[2],]
  landfrac <- landfrac[lonvals[1]:lonvals[2],latvals[1]:latvals[2]]
  
  lat <- lat[latvals[1]:latvals[2]]
  lon <- lon[lonvals[1]:lonvals[2]]
  
  
  # output
  pr_filename <- paste("pr_slm_",i,".RData", sep ="")
  save(lon, lat, precip, landfrac, file = pr_filename)
}

