clear_to_surface <- function(AOD, CS_PAR, B0 = -0.322868, B1 = -0.371085){
  return(CS_PAR*exp(B0 + AOD * B1))
}

one_percent <- function(I0,kd = 0.8,z = 13){
  z_oneper <- data.frame(matrix(NA,nrow=0,ncol=2))
  for(i in seq(from = 0, to = 55, by = z)){
    Iz <- I0 * exp(-(kd * i))
    
    if(round(Iz/I0, 2) == 0.01){
      z_oneper <- i
      rowdata <- c(z_oneper,I0)
      z_oneper <- rbind(z_oneper,rowdata)
      }
  }
  return(z_oneper)
}

#How does this translate to primary productivity
PPdepth_func <- function(dayLength, PPmax, Izt, Ik = 306.8, Vzdel,z) {
  # daylength must be in seconds between sunrise and sunset
  # PPmax (eqn 4)
  # make sure Izt (eqn 7) only includes times between sunrise and sunset
  # Ik constant for phytoplankton
  # Vzdel is lake Volume as slice between depth z and previous depths.
  dayLength*sum(PPmax*tanh(Izt/Ik)*(Vzdel))
}

PPmax <- 2.9*2.2
Iz <- 0.713*exp(-(0.8*13))
Izw <- 0.463*exp(-(0.8*13))
ppiz<-PPmax*tanh(Iz/306.8)
ppizw<-PPmax*tanh(Izw/306.8)

ppiz/ppizw

ppizw/ppiz

BPmax = 0.0139

bpiz<-BPmax*tanh(Iz/217.5)
bpizw <- BPmax*tanh(Izw/217.5)

bpizw/bpiz

plot(x=c(Iz,Izw),y=c(ppiz,ppizw))




Iz_function <- function(I0,kd = 0.8,z = 13){
  z_oneper <- data.frame(matrix(NA,nrow=0,ncol=2))
  for(i in seq(from = 0, to = 55, by = z)){
    Iz <- I0 * exp(-(kd * i))
    
    if(round(Iz/I0, 2) == 0.01){
      z_oneper <- i
      rowdata <- c(z_oneper,I0)
      z_oneper <- rbind(z_oneper,rowdata)
    }
  }
  return(z_oneper)
}
round(5.66666,2)
one_percent(0.713,z = 0.1)

-2.303