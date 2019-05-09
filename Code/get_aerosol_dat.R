# Read-in and process HDF5 data from OMSP-NPP NMTO3 L2 Total Column Aerosol data product

require(h5)
require(ncdf4)
require(binaryLogic)
require(ggplot2)
require(lubridate)
require(chron)

if(.Platform$OS.type == "unix") {
  root.dir <- "~/Documents/Projects/OneDrive/"
} else {
  if(Sys.info()[4] == "DESKTOP-FPVLIPB") {
    root.dir <- "C:/Users/seanr/OneDrive/"
  } else {
    root.dir <- "D:/Projects/OneDrive/"
  }
}

site.loc <- c(-122.254822, 47.6292991) # latitude an longitude for data location, in degrees E and degrees N
lat.tol <- 0.5 # Latitude buffer in degrees to find cells
lon.tol <- 0.5 # Longitude buffer in degrees to find cells

# Generate list of HDF5 data in the data directory
hdfs.list <- paste0(root.dir, "/uw/Research Derby/data/OMPS-NPP/", 
                    list.files(paste0(root.dir, "/uw/Research Derby/data/OMPS-NPP/"), 
                               pattern="h5$", recursive = T))

# Function to read-in OMPS-NPP Suomi aerosol data product
get_aerosol_dat <- function(filepath, # Filepath
                            lat.range, # Set extent for satellite data
                            lon.range) { # Set extent
    h5.str <- h5::h5file(filepath, mode = "r")
    out <- data.frame(lat = as.vector(as.matrix(readDataSet(h5.str["/GeolocationData/Latitude"]))), 
               lon = as.vector(as.matrix(readDataSet(h5.str["/GeolocationData/Longitude"]))),
               aerosol = as.vector(as.matrix(readDataSet(h5.str["/ScienceData/UVAerosolIndex"]))),
               quality = as.vector(as.matrix(readDataSet(h5.str["ScienceData/QualityFlags"]))),
               utc = as.vector(as.matrix(readDataSet(h5.str["/GeolocationData/UTC_CCSDS_A"])))) # Time in UTC
    out$utc <- sub("T", " ", substr(out$utc, 1,19)) # Remove unnecessary values
    out$utc <- as.POSIXct(out$utc, tz = "UTC")
    out <- out[(out$lat > min(lat.range) & out$lat < max(lat.range)),]
    out <- out[(out$lon > min(lon.range) & out$lon < max(lon.range)),]
  return(out)
}


for(i in 1:length(hdfs.list)) {
  if(i == 1) {
    omps.npp.aeroind <- get_aerosol_dat(filepath = hdfs.list[i], lat.range = c(site.loc[2]+lat.tol, site.loc[2]-lat.tol), lon.range = c(site.loc[1]+lon.tol, site.loc[1]-lon.tol))
  } else {
    omps.npp.aeroind <- rbind(omps.npp.aeroind,
                              get_aerosol_dat(filepath = hdfs.list[i], lat.range = c(site.loc[2]+lat.tol, site.loc[2]-lat.tol), lon.range = c(site.loc[1]+lon.tol, site.loc[1]-lon.tol)))
  }
  
}

omps.npp.aeroind$year <- year(omps.npp.aeroind$utc)
omps.npp.aeroind$pdt <- lubridate::with_tz(omps.npp.aeroind$utc, tzone = "America/Los_Angeles")

binaryLogic::as.binary(omps.npp.aeroind$quality[1], size = 10)
omps.npp.aeroind <- subset(omps.npp.aeroind, aerosol > -100)

# ggplotly(ggplot() + geom_line(data = subset(omps.npp.aeroind, aerosol > -100), aes(x = yday(pdt), y = aerosol, color = factor(year))) +
#   scale_x_continuous(breaks = c(121, 152, 182, 213, 244, 274, 295)) +
#   facet_wrap(~year) + theme(panel.background = element_blank(),
#                             panel.grid.major.y = element_blank(),
#                             panel.grid.major.x = element_line(color = "grey70"),
#                             panel.border = element_rect(color = "black", fill = NA),
#                             axis.line.x = element_line(color = "black"),
#                             axis.line.y = element_line(color = "black")))

