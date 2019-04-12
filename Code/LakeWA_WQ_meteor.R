########################
# Lake Washington WQ and Meteorological data
########################
# Plan ------------------------

# 1. convert solar radiation from meteorological data file to PAR using conversion factors
# in Holtgrieve et al 2010
#   a. Find if there are other conversion factors for cloudy days or smoke

# 2. QAQC water quality data using website thresholds to make sure there aren't negative values
# for things like chl a (e.g. convert values between 0 and -2 to 0; values less than -2 change to NA)
# see website: https://green2.kingcounty.gov/lake-buoy/Parameters.aspx

# Load libraries ---------------
library(dplyr)
source('Code/RD_functions.R')

library(rLakeAnalyzer) # has function thermo.depth() to find thermocline with some bs check for unique depths
source('Code/mod_thermo_depth_func.R')

# Import data ------------------
WQ <- read.csv("Data/Lake_WA_WQ.csv", header = T) 
meteor <- read.csv("Data/Lake_WA_meteor.csv", header = T)

str(WQ)
str(meteor)

# 1. Convert solar radiation to PAR at surface of lake 
# from watts per meter squared to micro Einsteins per second per meter squared
meteor.QC <- meteor
colnames(meteor.QC) <- c("date", "rel_humidity", "solar_rad_Wsqm", "atm_press", "wind_speed", "wind_dir", "bat_volts",
                         "air_temp", "pamout")
meteor.QC$PAR_surf_E <- meteor.QC$solar_rad_Wsqm*0.45*4.57

# 2. QAQC WQ data
WQ.QC <- WQ

# rename columns to short names
colnames(WQ.QC) <- c("date", "depth", "temp", "sp_cond", "DO_conc", "DO_sat", "pH", "chla", "NTU", "phyco", "prov",
                     "end_check", "parmout", "sonde", "prof_number")
str(WQ.QC)

# QC chla
WQ.QC$chla[WQ.QC$chla <= 0 & WQ.QC$chla >= -2] <- 0
WQ.QC$chla[WQ.QC$chla < (-2)] <- NA
WQ.QC$chla[WQ.QC$chla >= 50]

WQ.QC$pH[WQ.QC$pH < 4] <- 4
WQ.QC$pH[WQ.QC$pH > 10]

WQ.QC$DO_conc[WQ.QC$DO_conc < 0.01] <- NA
WQ.QC$DO_conc[WQ.QC$DO_conc > 22]

WQ.QC$DO_sat[WQ.QC$DO_sat < 0.10] 
WQ.QC$DO_sat[WQ.QC$DO_sat > 190]

WQ.QC$NTU[WQ.QC$NTU < (-5)]
WQ.QC$NTU[WQ.QC$NTU > (50)]

# Find thermocline ---------------------------
WQ.short <- WQ.QC[!is.na(WQ.QC$temp),]
length(unique(WQ.short$prof_number))
WQ.depth.check <- WQ.short %>% group_by(prof_number) %>% summarise(uni = length(length(depth) != length(unique(depth))))
WQ.depth.check$uni[WQ.depth.check$uni>1]
WQ.depth.check <- WQ.short %>% group_by(prof_number) %>% summarise(uni = length(depth))
WQ.depth.check <- WQ.short %>% group_by(prof_number) %>% summarise(uni = length(unique(depth)))
# WQ.depth.check[WQ.depth.check$uni < 57,]
WQ.short$depth[WQ.short$prof_number == 11335]
WQ.short$temp[WQ.short$prof_number == 113352]
thermo.depth(WQ.short$temp[WQ.short$prof_number == 113356], WQ.short$depth[WQ.short$prof_number == 113356])
str(WQ.short)
length(unique(WQ.short$prof_number))

WQ.cline <- WQ.short %>% group_by(prof_number) %>% # profile number represents each cast of depth profile data #[!is.na(WQ.QC$temp),] 
mutate(cline = thermo_depth(temp, depth))
head(WQ.cline$cline)

# Plug and chug ------------------------------
PPmax <- PPmax_func(WQ.QC$chla)
head(PPmax)
