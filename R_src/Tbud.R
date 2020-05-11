Tbud<-function(inputs){
  
  # computes air relative humidity
  inputs$RH<-qair2rh(qair = inputs$qair,temp = inputs$tair-273.15,press = inputs$patm*10)
  
  .f <- function(Tbud, inputs) {
    eb <- Ebalance(Tbud,inputs$tair,inputs$patm,inputs$RH, inputs$Sw,inputs$Lw,inputs$wind,bud_d=0.005,abs_sw=inputs$abs_sw,Ebal_only=TRUE)
    return(eb)
  }
  
  x1<-inputs$tair - 10
  x2<-inputs$tair + 10
  f1<-.f(inputs$tair - 10, inputs = inputs)
  f2<-.f(inputs$tair + 10, inputs = inputs)
  
  if(sign(f1)==sign(f2)){
    x1<-inputs$tair - 15
    x2<-inputs$tair + 15
    f1<-.f(inputs$tair - 15, inputs = inputs)
    f2<-.f(inputs$tair + 15, inputs = inputs)
  }
  
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

qair2rh <- function(qair, temp, press = 1013.25){
  es <-  6.112 * exp((17.67 * temp)/(temp + 243.5))
  e <- qair * press / (0.378 * qair + 0.622)
  rh <- e / es
  rh[rh > 1] <- 1
  rh[rh < 0] <- 0
  return(rh)
}


