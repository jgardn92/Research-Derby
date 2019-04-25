# Convert clear-sky irradiance to surface irradiance as a function of AOD

clear_to_surface <- function(AOD, CS_PAR, B0 = -0.322868, B1 = -0.371085) {
  return(CS_PAR * exp(B0 + AOD * B1))
}
