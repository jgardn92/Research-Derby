# Functions
# Convert clear-sky irradiance to surface irradiance as a function of AOD

clear_to_surface <- function(AOD, CS_PAR, B0 = -0.322868, B1 = -0.371085) {
  return(CS_PAR * exp(B0 + AOD * B1))
}

# Literature reference for most equations: Sources/Vadeboncoeur et al. 2008 Ecology

# Eq 4 ------------------------------------------------------------------------------------------------------------
# Phytoplankton productivity (PP; mg Cm-3 s-1)

PPmax_func <- function(chl) {
  (chl*2.2) / 3600 # chl is chlorophyll a input
}

# Eq 5 ------------------------------------------------------------------------------------------------------------
# Thermocline

# Modified thermo.depth() function from R package 'rLakeAnalyzer'
# https://rdrr.io/cran/rLakeAnalyzer/src/R/thermo.depth.R
#
#' @title Calculate depth of the thermocline from a temperature profile.
#' 
#' @description This function calculates the location of the thermocline from a temperature
#' profile.  It uses a special technique to estimate where the thermocline lies
#' even between two temperature measurement depths, giving a potentially
#' finer-scale estimate than usual techniques.
#' 
#' 
#' @param wtr a numeric vector of water temperature in degrees C
#' @param depths a numeric vector corresponding to the depths (in m) of the wtr
#' measurements
#' @param Smin Optional paramter defining minimum density gradient for
#' thermocline
#' @param seasonal a logical value indicating whether the seasonal thermocline
#' should be returned. This is fed to thermo.depth, which is used as the
#' starting point.  The seasonal thermocline is defined as the deepest density
#' gradient found in the profile. If \code{FALSE}, the depth of the maximum
#' density gradient is used as the starting point.
#' @param index Boolean value indicated if index of the thermocline depth,
#' instead of the depth value, should be returned.
#' @param mixed.cutoff A cutoff (deg C) where below this threshold,
#' thermo.depth and meta.depths are not calculated (NaN is returned). Defaults
#' to 1 deg C.
#' @return Depth of thermocline. If no thermocline found, value is NaN.
#' @seealso \code{\link{ts.thermo.depth}}, \code{water.density}
#' @keywords manip
#' @examples
#' 
#' # A vector of water temperatures
#' wtr = c(22.51, 22.42, 22.4, 22.4, 22.4, 22.36, 22.3, 22.21, 22.11, 21.23, 16.42, 
#' 		15.15, 14.24, 13.35, 10.94, 10.43, 10.36, 9.94, 9.45, 9.1, 8.91, 8.58, 8.43)
#' 
#'  #A vector defining the depths
#'  depths = c(0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 
#'      17, 18, 19, 20)
#' 
#'  t.d = thermo.depth(wtr, depths, seasonal=FALSE)
#' 
#'  cat('The thermocline depth is:', t.d)
#' 
#' @export
thermo_depth <- function(wtr, depths, Smin = 0.1, seasonal=TRUE, index=FALSE, mixed.cutoff=1){
  
  if(any(is.na(wtr))){
    return(NaN)
  }
  
  if(diff(range(wtr, na.rm=TRUE)) < mixed.cutoff){
    return(NaN)
  }
  
  #We can't determine anything with less than 3 measurements
  # just return lake bottom
  if(length(wtr) < 3){
    return(NaN)
  }
  
  # if(length(depths) != length(unique(depths))){
  #   stop('Depths all must be unique')
  # }
  
  #We need water density, not temperature to do this
  rhoVar = water.density(wtr)
  
  dRhoPerc = 0.15; #in percentage max for unique thermocline step
  numDepths = length(depths);
  drho_dz = rep(NaN, numDepths-1);
  
  #Calculate the first derivative of density
  for(i in 1:(numDepths-1)){
    drho_dz[i] = ( rhoVar[i+1]-rhoVar[i] )/( depths[i+1] - depths[i] );
  }
  
  #look for two distinct maximum slopes, lower one assumed to be seasonal
  thermoInd = which.max(drho_dz)  #Find max slope
  mDrhoZ = drho_dz[thermoInd]
  thermoD = mean( depths[thermoInd:(thermoInd+1)] )
  
  if(thermoInd > 1 && thermoInd < (numDepths-1)){  #if within range
    Sdn = -(depths[thermoInd+1] - depths[thermoInd])/
      (drho_dz[thermoInd+1] - drho_dz[thermoInd])
    
    Sup = (depths[thermoInd]-depths[thermoInd-1])/
      (drho_dz[thermoInd]-drho_dz[thermoInd-1])
    
    upD  = depths[thermoInd];
    dnD  = depths[thermoInd+1];
    
    if( !is.infinite(Sup) & !is.infinite(Sdn) ){
      thermoD = dnD*(Sdn/(Sdn+Sup))+upD*(Sup/(Sdn+Sup));
    }
  }
  
  dRhoCut = max( c(dRhoPerc*mDrhoZ, Smin) )
  locs = findPeaks(drho_dz, dRhoCut)
  pks = drho_dz[locs]
  
  if(length(pks) == 0){
    SthermoD = thermoD
    SthermoInd = thermoInd
  }else{
    mDrhoZ = pks[length(pks)]
    SthermoInd = locs[length(pks)]
    
    if(SthermoInd > (thermoInd + 1)){
      SthermoD = mean(depths[SthermoInd:(SthermoInd+1)])
      
      if(SthermoInd > 1 && SthermoInd < (numDepths - 1)){
        Sdn = -(depths[SthermoInd+1] - depths[SthermoInd])/
          (drho_dz[SthermoInd+1] - drho_dz[SthermoInd])
        
        Sup = (depths[SthermoInd] - depths[SthermoInd-1])/
          (drho_dz[SthermoInd] - drho_dz[SthermoInd-1])
        
        upD  = depths[SthermoInd]
        dnD  = depths[SthermoInd+1]
        
        
        if( !is.infinite(Sup) & !is.infinite(Sdn) ){
          SthermoD = dnD*(Sdn/(Sdn+Sup))+upD*(Sup/(Sdn+Sup))
        }
      }
    }else{
      SthermoD = thermoD
      SthermoInd = thermoInd
    }
  }
  
  if(SthermoD < thermoD){
    SthermoD = thermoD
    SthermoInd = thermoInd
  }
  
  #Ok, which output was requested. Index or value
  # seasonal or non-seasonal
  if(index){
    if(seasonal){
      return(SthermoInd)
    }else{
      return(thermoInd)
    }
  }else{
    if(seasonal){
      return(SthermoD)
    }else{
      return(thermoD)
    }
  }
  
  #list( thermoD, thermoInd, drho_dz, SthermoD, SthermoInd )
}

# Finds the local peaks in a vector. Checks the optionally supplied threshold 
#  for minimum height.
findPeaks <- function(dataIn, thresh=0){
  
  varL = length(dataIn);
  locs = rep(FALSE, varL);
  peaks= rep(NaN, varL);
  
  for(i in 2:varL-1){
    pkI = which.max(dataIn[(i-1):(i+1)])
    posPeak = max(dataIn[(i-1):(i+1)]);
    
    if(pkI == 2) {
      peaks[i] = posPeak;
      locs[i]  = TRUE;
    }
  }
  
  inds = 1:varL;
  locs = inds[locs];
  peaks= peaks[locs];
  
  # remove all below threshold value
  
  useI = peaks > thresh;
  locs = locs[useI];
  
  return(locs)
}

# Eq 6 - Jenny done -----------------------------------------------------------------------------------------------
# Light attenuation coefficient at depth Kd
# Kd = Kb + 0.015[Chla]

#Kb <- 0.8 # based on parameters in Arhonditsis paper

Light_atten_func <- function(chl, kb = 0.8) {
  Kb + 0.015*chl
}

# Using linear regression for surface light equation from AOD -------------------------------------------------------------------------------------
# AOD Aerosol optical depth
# PAR(surface buoy) ~ PAR(clear sky) + AOD
# cloudy sky conversion PAR ratio = 0.41 (amount attenuated by aerosol layer) Kobayashi et al. 2004 JGR Atmostpheres

# Eq 7 --------------------------------------------------------------------------------------------------------------
# Light at depth z, time t; micro Mol per meter squared per second; I0t is surface light
# Izt = I0t*exp(-Kd*z)

Izt_func <- function(I0t, Kd, z) {
  I0t*exp(-Kd*z) # I0t is surface light, Kd is light attenuation coefficient, d is depth
}

# Eq 9 --------------------------------------------------------------------------------------------------------------
# daily phytoplankton primary production at depth z (mg C)
# Ik: 306.8 umol photons/m2/s phytoplankton based on comm composition for lake washington by Arhoditsis et al. 2003
# mix of diatoms, chlorophytes, cyanobacteria -- weighted by relative community composition from literature values
# Citations:
# Cyano: Schuurmans et al. 2015 PLOSONE, Olaizola and Duerr 1990 J App Phycology, Luimstra et al. 2018 Photosynthesis Res
# Diatoms: Heiden et al. 2016 Frontiers in Mar Sci
# Chlorophytes: Xu et al. 2016 Plant phys and biochem

PPdepth_func <- function(dayLength, PPmax, Izt, Ik = 306.8, Vzdel) {
  # daylength must be in seconds between sunrise and sunset
  # PPmax (eqn 4)
  # make sure Izt (eqn 7) only includes times between sunrise and sunset
  # Ik constant for phytoplankton
  # Vzdel is lake Volume as slice between depth z and previous depths.
  dayLength*sum(PPmax*tanh(Izt/Ik)*(Vzdel))
}

# Eqn 10 ------------------------------------------------------------------------------------------------------------
# daily whole-lake phytoplankton production, TPP (mg C/m squared)
# A0 surface area in meters squared
TPP_func <- function(PPdepth, A0 = 75713400) {
  sum(PPdepth)/A0
}

# Eqn 11 ------------------------------------------------------------------------------------------------------------
# Periphyton: Sources/Austin and Hill 1991 Limno and Ocean
#Ik is 217.5 umol photons/m2/s
# BPmax: 50 mg C per meter per hour (convert to seconds) for mesotrophic lakes with unconsolidated sediments. 
# daily benthic primary production, BP, at depth z (mg C)
# WARNING: LAKE WASHINGTON IS NOT A CONE THIS IS ROUGH APPROXIMATION OF BENTHIC PRIMARY PRODUCTIVITY FOR RD 
# b/c don't have bathymetry
BPz_func <- function(dayLength, BPmax = 0.0139, Izt, Ik = 217.5, Azdel) {
  daylength*sum(BPmax*tanh(Izt/Ik)*Azdel)
}

# Eqn 12 --------------------------------------------------------------------------------------------------------------
# daily whole-lake periphyton production, TBP mg C per meter squared
TBP <- function(BPz, A0 = 75713400) {
  sum(BPz)/A0
}

