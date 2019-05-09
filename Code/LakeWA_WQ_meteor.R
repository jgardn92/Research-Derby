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

# 3. Calculate productivity using Vadeboncoeur equations from file RD_functions.R

# Load libraries ---------------
library(dplyr)
source('Code/RD_functions.R')

library(rLakeAnalyzer) # has function thermo.depth() to find thermocline with some bs check for unique depths

# packages chron or lubridate function yday
library(lubridate)
?yday
?hour



# Import data ------------------
WQ <- read.csv("Data/Lake_WA_WQ.csv", header = T) 
meteor <- read.csv("Data/Lake_WA_meteor.csv", header = T)
AquaPAR <- read.csv("Data/SeattleAquar_PAR.csv", header = T) # use this because lake surface buoy readings seem strange
  # NOTE: AquarPAR readings are in micro moles while PAR from GlobalIrradiance PAR are micro Einsteins
# but apparently that is the same unit. Also, this is in 15 minute intervals
cone <- read.csv("Data/lw_tenth_mbins.csv", header = T) # Bathymetry data
str(cone)
colnames(cone) <- c("depth", "radius", "surf_area", "volume", "side_area")

AOD <- readRDS("Data/2019-04-05aerosol_dat.rds")
str(AOD)

str(WQ)

str(meteor)


# 1. Convert solar radiation to PAR at surface of lake 
# from watts per meter squared to particle flux (micro Einsteins per second per meter squared)
meteor.QC <- meteor
colnames(meteor.QC) <- c("date", "rel_humidity", "solar_rad_Wsqm", "atm_press", "wind_speed", "wind_dir", "bat_volts",
                         "air_temp", "pamout")
meteor.QC$PAR_surf_E <- meteor.QC$solar_rad_Wsqm*0.45*4.57

#######################################################################################################################
# PAR PAR
#######################################################################################################################
plot(PARcs_df$hhour, PARcs_df$PAR_particleFlux)
str(AquaPAR2)
# Join PAR clear sky to PAR from Lake Buoy ----------------------------------------------------------------------------
# str(PARcs_df) # from code file GlobalIrradiance.R hourly data
# meteor.QC$PAR_surf_E <- meteor.QC$solar_rad_Wsqm*0.45*4.57
# # str(meteor.QC)

# meteor2.QC <- meteor.QC
# meteor2.QC$date <- as.POSIXct(strptime(meteor.QC$date, "%m/%d/%Y %H:%M"))
# head(meteor.QC)
# head(meteor2.QC)
# hour(meteor2.QC$date)
# meteor.dates <- meteor2.QC %>% mutate(day = day(date), month = month(date), year = year(date), hhour = hour(date))
# PARs_df <- meteor.dates %>% mutate(PARs = solar_rad_Wsqm*0.45, PARs_particleFlux = solar_rad_Wsqm*0.45*4.57)
# head(PARs_df$PARs_particleFlux, 10)
# head(PARcs_df$PAR_particleFlux, 10)
# str(PARcs_df)
# str(PARs_df)
# PAR_all <- inner_join(PARs_df, PARcs_df, by = c("year", "month", "day", "hhour"))
# 
# plot(PAR_all$PAR_particleFlux, PAR_all$PARs_particleFlux)

# Join PAR clear sky to PAR from Seattle Aquarium PAR ---------------------------------------------------------------
str(AquaPAR)

AquaPAR2 <- AquaPAR
colnames(AquaPAR2) <- c("Date", "PAR_surf_Flux", "Qual_Surf_PAR")
AquaPAR2$Date <- as.POSIXct(strptime(AquaPAR2$Date, "%m/%d/%Y %H:%M")) #, tz = "UTC")
# AquaPAR2$Date <- with_tz(AquaPAR2$Date, tzone = "America/Los_Angeles")
AquaPAR.dates <- AquaPAR2 %>% mutate(doy = yday(Date), year = year(Date), hour = hour(Date), minute = minute(Date)/60)
str(AquaPAR.dates)
AquaPARs <- AquaPAR.dates %>% mutate(hhour = hour + minute)
str(AquaPARs)

PARcs_df <- PARcs_df %>% mutate(doy = yday(D)) #Need to get more updated file from laptop 

PAR_all <- inner_join(AquaPARs, PARcs_df, by = c("doy", "year", "hhour")) #PARcs_df comes from GlobalIrradiance.R file
plot(PAR_all$PAR_particleFlux, PAR_all$PAR_surf_Flux)

str(PAR_all)
# PAR_all_day <- PAR_all[PAR_all$hhour >= PAR_all$sunrise + 0.5 & PAR_all$hhour <= PAR_all$sunset - 0.5,]
PAR_all_day <- PAR_all[PAR_all$hhour >= 10 & PAR_all$hhour <= 14,]
str(PAR_all_day)

# write.csv(PAR_all_day, file = "Output/PAR_daylight_buoy.csv", row.names = F, col.names = T)
# saveRDS(PAR_all_day, file = "Output/PAR_daylight_buoy.rds")

# 2. QAQC WQ data
WQ.QC <- WQ

# rename columns to short names
colnames(WQ.QC) <- c("date", "depth", "temp", "sp_cond", "DO_conc", "DO_sat", "pH", "chla", "NTU", "phyco", "prov",
                     "end_check", "parmout", "sonde", "prof_number")
str(WQ.QC)

# QC chla
# Replacing negative values with zero because negative values are outside the real range of chlorophyll values
# Changes in values to min and max thresholds comes from the buoy sonde expected range values published on-line
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


# Plug and chug ------------------------------
# Eqn 4 Phytoplankton productivity per sec
WQ.short2 <- WQ.QC[!is.na(WQ.QC$temp),]
WQ.short2 <- WQ.short2 %>% mutate(PPmax = PPmax_func(chla))
head(WQ.short2)

# Eqn 5 Find thermocline ---------------------------

# length(unique(WQ.short$prof_number))
# WQ.depth.check <- WQ.short %>% group_by(prof_number) %>% summarise(uni = length(length(depth) != length(unique(depth))))
# WQ.depth.check$uni[WQ.depth.check$uni>1]
# WQ.depth.check <- WQ.short %>% group_by(prof_number) %>% summarise(uni = length(depth))
# WQ.depth.check <- WQ.short %>% group_by(prof_number) %>% summarise(uni = length(unique(depth)))
# # WQ.depth.check[WQ.depth.check$uni < 57,]
# WQ.short$depth[WQ.short$prof_number == 11335]
# WQ.short$temp[WQ.short$prof_number == 113352]
# thermo.depth(WQ.short$temp[WQ.short$prof_number == 113356], WQ.short$depth[WQ.short$prof_number == 113356])
# str(WQ.short)
# length(unique(WQ.short$prof_number))

# WQ.short <- WQ.short2 %>% group_by(prof_number) %>% # profile number represents each cast of depth profile data #[!is.na(WQ.QC$temp),] 
#   mutate(cline = thermo_depth(temp, depth))
# str(WQ.short)

# Eqn 6 Light attenuation function (Kd) ------------------
WQ.short <- WQ.short2 %>% mutate(Kd = Light_atten_func(chla))

which(names(WQ.short) %in% names(PAR_AOD))

# Eqn 7 Light at depth z time t --------------------------
# !!!! Need to join I0t and WQ data base so the vectors will line up
WQ.short$date <- as.POSIXct(strptime(WQ.short$date, "%m/%d/%Y %H:%M"))
WQ.short <- WQ.short %>% mutate(doy = yday(date), year = year(date), hour = hour(date), minute = minute(date)/60)
WQ.short<- WQ.short %>% mutate(hhour = hour + minute)
table(WQ.short$doy, WQ.short$year) # some years have missing days of observation 
# join AOD to PARcs_df
str(AOD)
dplyr::select(AOD, doy, year, AOD)
head(WQ_PAR)
AOD$hhour <- AOD$hour
# PAR_AOD <- inner_join(AOD, PARcs_df, by = c("hhour", "doy", "year"))
PAR_AOD <- merge(dplyr::select(AOD, doy, year, AOD), PARcs_df)
?select
str(test.PAR)
PAR_AOD$surf_light <- clear_to_surface(AOD = PAR_AOD$AOD, CS_PAR = PAR_AOD$PAR_particleFlux)
str(PAR_AOD)
str(WQ.short)
WQ_PAR <- merge(unique(dplyr::select(PAR_AOD, doy, year, dayLength, surf_light)), unique(dplyr::select(WQ.short, chla, PPmax, doy, year, Kd, depth)))
unique(dplyr::select(PAR_AOD, doy, year, dayLength, surf_light))
unique(dplyr::select(WQ.short, chla, PPmax, doy, year, Kd, depth))
#### dayLength and surf_light is a problem -- should be same values within each day.
#WQ_PAR <- inner_join(WQ.short, PAR_AOD, by = c("doy", "year", "hhour"))
str(WQ_PAR)
# WQ_PAR <- WQ_PAR[!is.na(WQ_PAR$Kd),]
# WQ_PAR <- WQ_PAR[!is.na(WQ_PAR$PAR_surf_Flux),]
# Izt <- Izt_func(I0t = WQ_PAR$PAR_surf_Flux, Kd = WQ_PAR$Kd, z = WQ_PAR$depth)
WQ_PAR <- WQ_PAR %>% mutate(Izt = Izt_func(I0t = surf_light, Kd = Kd, z = depth))

# Eqn 9 daily phytoplankton primary production at depth z
# !!!! Need to join PARcs_df to WQ_PAR by date/time to get daylengths

WQ_PAR_cone <- inner_join(WQ_PAR, cone, by = "depth")
# Bathymetry LW_0.1 m bins for Vzdel and Azdel 
#PP_depth <- PPdepth_func(dayLength = WQ_PAR_cone$dayLength, PPmax = WQ_PAR_cone$PPmax, 
#                         Izt = WQ_PAR_cone$Izt, Vzdel = WQ_PAR_cone$volume)
WQ_PAR_cone <- WQ_PAR_cone %>% mutate(PP_depth = PPdepth_func(dayLength = dayLength, PPmax = PPmax, 
                                      Izt = Izt, Vzdel = volume))

# Eqn 10 daily whole lake phytoplankton production ----------
str(WQ_PAR_cone)
WQ_PAR_cone <- WQ_PAR_cone %>% mutate(TPP = TPP_func(PPdepth = PP_depth))
    # Then summarize TPP by day to get right data



# Eqn 11 daily benthic primary production at depth z
# BP_depth <- BPz_func(dayLength = xx, Izt = yyIzt, Azdel = zzArea)
WQ_PAR_cone <- WQ_PAR_cone %>% mutate(BP_depth = BPz_func(dayLength = dayLength, Izt = Izt, Azdel = side_area))

# Eqn 12 whole lake periphyton production TBP
TBP <- TBP(BPz = BP_depth )
    # THen summarize TBP by day to get right data
WQ_PAR_cone <- WQ_PAR_cone %>% mutate(TBP = TBP_func(BPz = BP_depth))
str(WQ_PAR_cone)

(table(WQ_PAR_cone$year))
