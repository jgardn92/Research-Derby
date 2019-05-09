# Import GOES hourly data from netCDF files
# SATCORPS CERES GEO Level 2 Edition-4 GOES-15 Northern Hemisphere Version 1.0, Version 1.2

# Functions to import cloud data and find the ...

require(ncdf4)
require(lubridate)
require(ggplot2)

# source("setup.R")
# 
# if(.Platform$OS.type == "unix") {
#   ext.dir <- "~/"
# } else {
#   if(Sys.info()[4] == "DESKTOP-FPVLIPB") {
#     ext.dir <- "D:"
#   } else {
#     ext.dir <- "F:"
#   }
# }


site.loc <- c(-122.254822, 47.6292991) # latitude an longitude of location, in degrees E and degrees N
lat.tol.8k <- 0.1 # Latitude buffer in degrees
lon.tol.8k <- 0.1 # Longitude buffer in degrees

ncs.list <- paste0(ext.dir, "/Research Derby/CERES/", 
                    list.files(paste0(ext.dir, "/Research Derby/CERES/"), 
                               pattern="nc$", recursive = T))
test <- nc_open(ncs.list[1])
ncvar_get(test, varid = "base_time")
nc_close(test)             
             
get_cloud_dat <- function(filepath, # Filepath
                           lat.range, # Set extent for satellite data
                           lon.range) {
  nc.dat <- nc_open(filepath)
  
  # Import cloud visible OD, latitude, longitude, scan initiation time. Raw data use seconds since 1970-01-01 00:00:00 GMT, so I used a conversion of hourly scan time.
  out <- data.frame(vis.od = as.vector(ncvar_get(nc.dat, varid = "cloud_visible_optical_depth")),
                    lat = as.vector(ncvar_get(nc.dat, varid = "latitude")),
                    lon = as.vector(ncvar_get(nc.dat, varid = "longitude")))
  out <- out[complete.cases(out),]
  out$scan.time <- as.POSIXct("1970-01-01 00:00:00", tz = "GMT", format = "%Y-%m-%d %H:%M:%S") + 
    ncvar_get(nc.dat, varid = "base_time") + 
    ncvar_get(nc.dat, varid = "time_offset")
  out <- out[(out$lat > min(lat.range) & out$lat < max(lat.range)),]
  out <- out[(out$lon > min(lon.range) & out$lon < max(lon.range)),]
  nc_close(nc.dat)
  return(out)
  
}

# Cloud, latitude, longitude

scan.time <- as.POSIXct("1970-01-01 00:00:00", tz= "GMT", format = "%Y-%m-%d %H:%M:%S") + ncvar_get(test, varid = "base_time") + ncvar_get(test, varid = "time_offset")


for(i in 1:length(ncs.list)) {
  if(file.exists(ncs.list[i])) {
  vis.iter <- get_cloud_dat(filepath = ncs.list[i],  lat.range = c(site.loc[2]+lat.tol.8k, site.loc[2]-lat.tol.8k), lon.range = c(site.loc[1]+lon.tol.8k, site.loc[1]-lon.tol.8k))
  if(i == 1) {
    vis.od <- vis.iter
  } else {
    vis.od <- rbind(vis.od, vis.iter)
  }
  }
  #if((i %% 10) == 0) {
    print(i)
  #}
}


# scan.time <- with_tz(scan.time, tzone = "America/Los_Angeles")
vis.od$doy <- yday(vis.od$scan.time) + hour(vis.od$scan.time)/24

vis.od.mean <- aggregate(vis.od ~ scan.time + doy, data = vis.od, FUN = mean)
vis.od.mean <- vis.od.mean[order(vis.od.mean$doy),]



