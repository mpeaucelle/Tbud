############################################################
# Tbud function
# Author: Marc Peaucelle, 2020/06
# Description:
# This function search for bud temperature such as the bud energy balance is at equilibrium
# Call the Ebalance function
############################################################
# inputs: inputs list defined in 1.0_create_database_Tbud.R
#  inputs<-list(tair, in K
#               patm, in kPa
#               qair, in kg kg-1
#               Sw, in W m-2
#               Lw, in W m-2
#               wind, in m s-1
#               abs_sw) 0-1 unitless
############################################################
# outputs: Bud temperature in K
############################################################
 
Tbud<-function(inputs, abs_sw, bud_d){
  
  # computes air relative humidity
  #inputs$RH<-qair2rh(qair = inputs$qair,temp = inputs$tair-273.15,press = inputs$patm*10)
  
  .f <- function(Tbud, inputs) {
    eb <- Ebalance(Tbud = Tbud,
                   T_air = inputs$tair,
                   T_ground = inputs$tground,
                   P = inputs$patm,
                   RH = inputs$RH, 
                   qair = inputs$qair,
                   Sw =  inputs$Sw,
                   Lw = inputs$Lw,
                   wind = inputs$wind,
                   rain =  inputs$rain,
                   bud_d = inputs$bud_d, 
                   abs_sw = inputs$abs_sw, 
                   r=inputs$r, 
                   evap = inputs$evap, 
                   md = inputs$md, 
                   max_dw = inputs$max_dw,
                   tstep = inputs$tstep,
                   kE = 1, 
                   Ebal_only=TRUE)
    return(eb) 
  }
  
  x1<-inputs$tair - 25
  x2<-inputs$tair + 25
  f1<-.f(inputs$tair - 25, inputs = inputs)
  f2<-.f(inputs$tair + 25, inputs = inputs)
  
  if(sign(f1)==sign(f2)){
    print("There is a problem with Tbud")
    print("No root found, we set Tbud == tair")
    print(paste0("Ebal tair-20°C = ",f1))
    print(paste0("Ebal tair+20°C = ",f2))

    return(inputs$tair) 
  } else{
  
  results<-uniroot(
    f = .f, 
    inputs = inputs,
    lower = x1, 
    upper = x2,
    f.lower = f1,
    f.upper = f2,
    trace=10
  )
  return(results$root)
  }
}

# function to compute air relative humidity (rh) from specific air humidity (qair)
qair2rh <- function(qair, temp, press = 1013.25){
  es <-  6.112 * exp((17.67 * temp)/(temp + 243.5))
  e <- qair * press / (0.378 * qair + 0.622)
  rh <- e / es
  rh[rh > 1] <- 1
  rh[rh < 0] <- 0
  return(rh)
}

# function to compute specific air humidity (qair) from air relative humidity (rh)
rh2qair <- function(rh, temp, press = 1013.25){
  es <-  6.112 * exp((17.67 * temp)/(temp + 243.5))
  A <- es * rh
  qair <- ((0.622*A)/press)/(1-(0.378*A/press))
  return(qair)
}

############# Piece of code to include energy storage and thermal time constant ##############
############# Useless with a time resolution > 60 s
# 
# S = 0
# 
# # If storage is computed, we first compute the equilibrium Tbud == Teq
# # We thus apply the thermal time constant to have the real Tbud based on Tbud at the previsous step dT0 = Tbud0
# if (storage){
#   # compute dT/dt
#   Teq = Tbud
#   Tbud = Teq - (Teq-Tbud0)*exp(-1/tau)
#   
#   # change in T for each second
#   dT = Tbud-Tbud0
#   # Tdiff to equilibrium
#   dT0 = Teq-Tbud0
#   # number of seconds to reach the equilibrium
#   n_step = min(dT0/dT,tstep)
#   
#   # calculate energy storage to equilibrium
#   # S = mass * specific heat * dT/dt
#   S = dT * n_step * c_b * rho_b * 4/3*pi*(bud_d/2)**3
# }
