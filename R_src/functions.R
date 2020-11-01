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
fTbud<-function(x,abs_sw){
  inputs<-list(tair=tg_tmp[x]+273.15,
               patm=pb_tmp[x]/1000,
               qair=qair_tmp[x],
               Sw=sw_tmp[x],
               Lw=lw_tmp[x],
               wind=wind_tmp[x],
               abs_sw=abs_sw,
               bud_d=0.005,
               r=0.2)
  
  tmp<-Tbud(inputs)-273.15
  return(tmp)
}

# function to retrieve the corresponding energy budget components at Tbud
fEbal<-function(x,abs_sw,ind,Tbud){
  
  inputs<-list(tair=tg_tmp[x]+273.15,
               patm=pb_tmp[x]/1000,
               qair=qair_tmp[x],
               Sw=sw_tmp[x],
               Lw=lw_tmp[x],
               wind=wind_tmp[x],
               abs_sw=abs_sw,
               bud_d=0.005,
               r=0.2)
  inputs$RH<-qair2rh(qair = inputs$qair,temp = inputs$tair-273.15,press = inputs$patm*10)
  
  tmp<-Ebalance(Tbud[x]+273.15,inputs$tair,inputs$patm,inputs$RH, inputs$Sw,inputs$Lw,inputs$wind,bud_d=0.005,abs_sw=inputs$abs_sw,Ebal_only=FALSE)
  return(tmp[[ind]])
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
bis.fn<-function(xts.v,bisextile){
  tmp_v<-c(xts.v[1:bisextile[1]])
  for(i in 2:length(bisextile)){
    tmp_v<-c(tmp_v,xts.v[bisextile[i-1]:bisextile[i]])
  }
  tmp_v<-c(tmp_v,xts.v[bisextile[i]:length(xts.v)])
  return(tmp_v)
}
