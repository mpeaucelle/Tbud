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

Tbud<-function(inputs){
  
  # computes air relative humidity
  inputs$RH<-qair2rh(qair = inputs$qair,temp = inputs$tair-273.15,press = inputs$patm*10)
  
  .f <- function(Tbud, inputs) {
    eb <- Ebalance(Tbud,inputs$tair,inputs$patm,inputs$RH, inputs$Sw,inputs$Lw,inputs$wind,bud_d=inputs$bud_d,abs_sw=inputs$abs_sw,r=inputs$r,Ebal_only=TRUE)
    return(eb)
  }
  
  x1<-inputs$tair - 20
  x2<-inputs$tair + 20
  f1<-.f(inputs$tair - 20, inputs = inputs)
  f2<-.f(inputs$tair + 20, inputs = inputs)
  
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
