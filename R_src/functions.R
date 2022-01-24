############################################################
# Diverse functions
# Author: Marc Peaucelle, 2020/06
# Description:
# Several functions used to average, convert and extract data
############################################################


ave.period<-function(x,fn,xts.v){
  mtx<-fn(xts.v[paste(start_date[x],end_date[x],sep="/")],na.rm=T)
  return(mtx)
}

ave.period_X<-function(x,fn,xts.v,fac=1){
  med<-median(xts.v,na.rm=T)
  iqr<-IQR(xts.v,na.rm=T)
  mini<-which(xts.v>(med+fac*iqr))
  maxi<-which(xts.v<(med-fac*iqr))
  if(length(c(mini,maxi))>0){
   xts.v<-xts.v[-c(mini,maxi)]
  }
  mtx<-fn(xts.v[paste(start_date[x],end_date[x],sep="/")],na.rm=T)
  return(mtx)
}

fmean<-function(x,period){
  #length(x)
  mtx<-mean(x[period[1]:period[2]],na.rm=T)
  return(mtx)
}

fsum<-function(x,period){
  #length(x)
  mtx<-sum(x[period[1]:period[2]],na.rm=T)
  return(mtx)
}

pmean<-function(clim,x){
  mtg<-mean(clim[paste(start_date[x],end_date[x],sep="/")],na.rm=T)      
  return(mtg)
}
pmin<-function(clim,x){
  mtg<-min(clim[paste(start_date[x],end_date[x],sep="/")],na.rm=T)
  return(mtg)
}
pmax<-function(clim,x){
  mtg<-max(clim[paste(start_date[x],end_date[x],sep="/")],na.rm=T)
  return(mtg)
}

psum<-function(clim,x){
  mtg<-sum(clim[paste(start_date[x],end_date[x],sep="/")],na.rm=T)      
  return(mtg)
}


# function to computes bud temperature
fTbud<-function(x,abs_sw,bud_d=0.005,r=0.2,md,tstep,evap=TRUE,max_wd=0.15){
   inputs<-list(tair=tg_tmp[x]+273.15,
               tground=sol_tmp[x],
               patm=pb_tmp[x],
               RH=rh_tmp[x],
               qair=qair_tmp[x],
               Sw=sw_tmp[x],
               Lw=lw_tmp[x],
               wind=wind_tmp[x],
               rain=rain_tmp[x],
               abs_sw=abs_sw,
               bud_d=bud_d,
               r=r,
               evap=evap,
               md=md,
               max_dw=max_wd,
               tstep=tstep)
  
  tmp<-Tbud(inputs,abs_sw = abs_sw,bud_d = bud_d)
  return(tmp)
}

# function to retrieve the corresponding energy budget components at Tbud


fEbal <- function(x,abs_sw,bud_d=0.005,r=0.2,tstep,evap,max_dw=0.15) {
 # print(x)
  inputs<-list(tair=as.numeric(tg_tmp[x])+273.15,
               tground=as.numeric(sol_tmp[x]),
               patm=as.numeric(pb_tmp[x])/1000,
               RH=as.numeric(rh_tmp[x]),
               qair=as.numeric(qair_tmp[x]),
               Sw=as.numeric(sw_tmp[x]),
               Lw=as.numeric(lw_tmp[x]),
               wind=as.numeric(wind_tmp[x]),
               rain=as.numeric(rain_tmp[x]),
               abs_sw=abs_sw,
               bud_d=bud_d,
               r=r,
               evap=evap,
               md=md,
               max_dw=max_dw,
               tstep=tstep)
  
  Tb<-Tbud(inputs,abs_sw = abs_sw,bud_d = bud_d)
  
  eb <- Ebalance(Tb,inputs$tair,inputs$tground,inputs$patm,inputs$RH, inputs$qair, inputs$Sw,inputs$Lw,inputs$wind, inputs$rain,
                 bud_d=bud_d, abs_sw=abs_sw, r=inputs$r, 
                 evap = inputs$evap, md = inputs$md, max_dw = inputs$max_dw, 
                 tstep = inputs$tstep, Ebal_only=FALSE)
  md<<-eb$md
  return(eb) 
}

# computes chilling
fchil<-function(x,xts.v,alti_CRU=NULL,alti_obs=NULL){
  
  if (is.na(alti_CRU)){
     ttg<-xts.v[paste(start_date[x],end_date[x],sep="/")]

  } else {  
    ttg<-xts.v[paste(start_date[x],end_date[x],sep="/")]
    ttg<-ttg+0.0064*(alti_CRU-as.numeric(alti_obs[x]))
  }

  ttg[ttg>5]<-0
  ttg[ttg<0]<-0
  ttg[ttg!=0]<-1
  return(sum(ttg,na.rm=T))
}

# computes forcing

fforc5<-function(x,xts.v,alti_CRU=NULL,alti_obs=NULL){
 if (is.na(alti_CRU)){
     ttg<-xts.v[paste(start_date[x],end_date[x],sep="/")]

  } else {
    ttg<-xts.v[paste(start_date[x],end_date[x],sep="/")]
    ttg<-ttg+0.0064*(alti_CRU-as.numeric(alti_obs[x]))
  }
  
  ttg[ttg<5]<-0
  return(sum(ttg,na.rm=T))
}

# add bisextile days
bis.fn<-function(xts.v,bisextile,tstep=1){
  nb_it=length(bisextile)/tstep

  tmp_v<-c(xts.v[1:bisextile[tstep]])
  
  for(i in 2:nb_it){
    tmp_v<-c(tmp_v,xts.v[bisextile[(i-2)*tstep+1]:bisextile[i*tstep]])
  }
 
 if(bisextile[length(bisextile)]!=length(xts.v)){
    tmp_v<-c(tmp_v,xts.v[bisextile[(i-2)*tstep+1]:length(xts.v)])
  }
  return(tmp_v)
}
