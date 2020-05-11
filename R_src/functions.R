### fucntions
#require(oce)
#Daylength <- function(t, lon=lon, lat=lat)
#{
#  t <- as.numeric(t)
#  alt <- function(t)
#    sunAngle(t, longitude=lon, latitude=lat)$altitude
#  rise <- uniroot(alt, lower=t-86400/2, upper=t)$root
#  set <- uniroot(alt, lower=t, upper=t+86400/2)$root
#  return(set - rise)
#}

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
psum<-function(clim,x){
  mtg<-sum(clim[paste(start_date[x],end_date[x],sep="/")],na.rm=T)      
  return(mtg)
}

quart.fn<-function(x,fn,fn2){
  vec<-NULL

  for (i in 1:10){
    vec<-c(vec,fn(x[i:c(i+2)]))
  }
  leq<-fn2(vec)
  return(leq)
}

fchil<-function(x,xts.v,alti_CRU=NULL,alti_obs=NULL){
  
  if (is.na(alti_CRU)){
     ttg<-xts.v[paste(start_date[x],end_date[x],sep="/")]

  } else {  
    ttg<-xts.v[paste(start_date[x],end_date[x],sep="/")]
    ttg<-ttg+0.0064*(alti_CRU-as.numeric(alti_obs[x]))
  }
#  ttg<-xts.v[paste(start_date[x],end_date[x],sep="/")]
  ttg[ttg>5]<-0
  ttg[ttg<0]<-0
  ttg[ttg!=0]<-1
  return(sum(ttg,na.rm=T))
}
fforc<-function(x,xts.v,alti_CRU=NULL,alti_obs=NULL){
if (is.na(alti_CRU)){
     ttg<-xts.v[paste(start_date[x],end_date[x],sep="/")]

  } else {
    ttg<-xts.v[paste(start_date[x],end_date[x],sep="/")]
    ttg<-ttg+0.0064*(alti_CRU-as.numeric(alti_obs[x]))
  }
#  ttg<-xts.v[paste(start_date[x],end_date[x],sep="/")]
  #ttg[ttg>0]<-1
  ttg[ttg<0]<-0
  return(sum(ttg,na.rm=T))
}

fforc5<-function(x,xts.v,alti_CRU=NULL,alti_obs=NULL){
 if (is.na(alti_CRU)){
     ttg<-xts.v[paste(start_date[x],end_date[x],sep="/")]

  } else {
    ttg<-xts.v[paste(start_date[x],end_date[x],sep="/")]
    ttg<-ttg+0.0064*(alti_CRU-as.numeric(alti_obs[x]))
  }
# ttg<-xts.v[paste(start_date[x],end_date[x],sep="/")]
  #ttg[ttg>0]<-1
  ttg[ttg<5]<-0
  return(sum(ttg,na.rm=T))
}

bis.fn<-function(xts.v,bisextile){
  tmp_v<-c(xts.v[1:bisextile[1]])
  for(i in 2:length(bisextile)){
    tmp_v<-c(tmp_v,xts.v[bisextile[i-1]:bisextile[i]])
  }
  tmp_v<-c(tmp_v,xts.v[bisextile[i]:length(xts.v)])
  return(tmp_v)
}
