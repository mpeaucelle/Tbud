#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

rm(list=ls())
library(RNetCDF)
library(raster)
library(TTR)
library(xts)
library(geosphere)

source("Tbud.R")
source("Ebalance.R")
source("functions.R")
########### LOAD pheno data ##############
#PEP725
# Test for 1 species only 
workdir<-"/home/orchidee02/mpeau/pheno_PEP/"
i<-paste0(workdir,"db_site/betula_pendula_database.txt")
i<-paste0(workdir,"db_site/",args[1])
spname<-"betula"
first_year<-1970
last_year<-2014

obs_db<-read.table(i,sep=";",header=T)
LU_db<-obs_db[obs_db$BBCH==11,]

LU_db<-LU_db[order(LU_db$PEP_ID,LU_db$YEAR),]
LU_db<-LU_db[LU_db$YEAR>c(first_year-1),]

########### LOAD CRU data #######
# need hourly data
nc_tg<-open.nc(paste0(workdir,"netcdf_forcing/crun_tair_1969_2016.nc"))
nc_sw<-open.nc(paste0(workdir,"netcdf_forcing/crun_swdown_1969_2016.nc"))
nc_lw<-open.nc(paste0(workdir,"netcdf_forcing/crun_lwdown_1969_2016.nc"))

lat<-var.get.nc(nc_tg,"latitude")
lon<-var.get.nc(nc_tg,"longitude")
time<-1:17520 #in days since 1950-01-01 00:00
time<-as.Date(time,origin="1968-12-31")
# CRU data are based on a 365d calendar
bisextile<-which(as.character(time)%in%c("1971-12-31","1975-12-31","1979-12-31","1983-12-31",
                                         "1987-12-31","1991-12-31","1995-12-31","1999-12-31",
                                         "2003-12-31","2007-12-31","2011-12-31","2015-12-31"))
time<-1:17532
time<-as.Date(time,origin="1968-12-31")
themask<-var.get.nc(nc_tg,"mask")

alt<-open.nc(paste0(workdir,"netcdf_forcing/force1994_2002.nc"))
alti<-var.get.nc(alt,"altitude")

rast_mask<-raster(t(themask),xmn=min(lon),xmx=max(lon),ymn=min(lat),ymx=max(lat))
rast_alt<-raster(t(alti),xmn=min(lon),xmx=max(lon),ymn=min(lat),ymx=max(lat))

print("fichiers CRU chargés")

### extract site ref-point

lonlat<-unique(cbind(LU_db$LON,LU_db$LAT))
refp_LU<-NULL
check_ref<-NULL
for (j in 1:length(lonlat[,1])){
  if(is.na(lonlat[j,1])==FALSE){
    print("retrieving refpoints for each site")
    print(paste(j,"/",length(lonlat[,1])))
    refp<-extract(x=rast_mask,y=cbind(as.numeric(lonlat[j,1]),as.numeric(lonlat[j,2])))
    
    if((refp[1]==-1) | is.na(refp[1])){	
      rowu<-rowFromY(rast_mask,y = as.numeric(lonlat[j,2]))
      colu<-colFromX(rast_mask,x = as.numeric(lonlat[j,1]))
      zone<-as.matrix(rast_mask)[(rowu-1):(rowu+1),(colu-1):(colu+1)]
      if(all(zone==-1)==FALSE){
        leq2<-which(is.na(zone)==FALSE,arr.ind=T)
        rowu=(rowu-1)+(leq2[1,1]-1)
        colu=(colu-1)+(leq2[1,2]-1) 
        refp_LU[j]<-as.matrix(rast_mask)[rowu,colu]
        check_ref[j]<-1
      } else {
        refp_LU[j]<--1
        check_ref[j]<-2
      }
    } else{
      refp_LU[j]<-refp
      check_ref[j]<-0
    }	
  } else{
    refp_LU[j]<-NA
    check_ref[j]<-2
  }
}

lonlat2<-sapply(X = 1:dim(LU_db)[1],FUN = function(x){paste(LU_db$LON[x],LU_db$LAT[x])},simplify = T)
lonlat<-sapply(X = 1:dim(lonlat)[1],FUN = function(x){paste(lonlat[x,1],lonlat[x,2])},simplify = T)

obref<-NULL
for (j in 1:length(lonlat)){
  obref[lonlat2%in%lonlat[j]]<-refp_LU[j]
}
LU_db<-cbind(LU_db,obref)
names(LU_db)[10]<-'REF'



########### extract CRU climate data for each site #############

origin<-as.Date(paste(LU_db$YEAR,"-01-01",sep=""))
sos_date<-origin+as.numeric(LU_db$DAY)

# presaison = 90 days before leaf unfolding
pre_date<-sos_date-90 # 3 month, fixed date for preseason

# longterm tg, pp, sw, lw
tg_lt<-sw_lt<-lw_lt<-NULL
#preseason longterm tg, sw, lw
tg_pre<-sw_pre<-lw_pre<-sw_pre2<-lw_pre2<-NULL

tbud_pre<-NULL

#CRU altitude
alti<-NULL
#chilling and forcing
chilling<-forcing<-forcing_bud<-chilling_bud<-NULL

refpoint<-as.character(unique(LU_db$REF))

if(length(which(refpoint=="-1"))>0){
  refpoint<-refpoint[-which(refpoint=="-1")]
}

for (j in 1:length(refpoint)){
  print(paste(j,"/",length(refpoint)))
  leq<-LU_db$REF%in%refpoint[j]
  lequel<-which(leq==TRUE)
  
  cell<-which(t(as.matrix(rast_mask))==as.numeric(refpoint[j]))
  alti[lequel]<-rast_alt[cell]
  # in Kelvin -273.15 --> °C
  # We add +1 to refpoint because refpoint is 0:67419, but var.get.nc start at 1 and not 0.
  tg_tmp<-var.get.nc(nc_tg,"Temperature",start=c(as.numeric(refpoint[j])+1,1),count=c(1,17520))-273.15
  # in J/m²/6h *21600 --> W m-2
  sw_tmp<-var.get.nc(nc_sw,"Incoming_Short_Wave_Radiation",start=c(as.numeric(refpoint[j])+1,1),count=c(1,17520))/21600
  lw_tmp<-var.get.nc(nc_lw,"Incoming_Long_Wave_Radiation",start=c(as.numeric(refpoint[j])+1,1),count=c(1,17520))
  
  
  ## add manually bisextile years unvailbale in CRU
  tg_tmp<-bis.fn(xts.v = tg_tmp,bisextile = bisextile)
  sw_tmp<-bis.fn(xts.v = sw_tmp,bisextile = bisextile)
  lw_tmp<-bis.fn(xts.v = lw_tmp,bisextile = bisextile)
  
  
  ## convert in time-series  
  sw_tmp<-xts(sw_tmp,time)
  lw_tmp<-xts(lw_tmp,time)
  tg_tmp<-xts(tg_tmp,time)
 

  # computes bud temperature
  fTbud<-function(x){
 	 inputs<-list(tair=tg_tmp[x],
	              patm=patm_tmp[x],
                      qair=qair_tmp[x],
	              sw=sw_tmp[x],
                      lw=lw_tmp[x],
                      wind=wind_tmp[x])
  
	 tmp<-Tbud(inputs)
	return(tmp)
  } 
 
  tg_bud<-sapply(X=c(1:length(tg_tmp)),FUN=fTbud,simplify=T)

 
  ####### if obs is not in the sea? XD sometimes at the land border ######
  if( all(is.na(tg_tmp))==FALSE){
    
    ######### longterm climate
    start_date<-paste(first_year,"-01-01",sep="")
    end_date<-paste(last_year,"-01-01",sep="")
    
    tg_lt[lequel]<-sapply(X=c(1:length(start_date)),FUN=ave.period,fn=mean,xts.v=tg_tmp,simplify=T)
    sw_lt[lequel]<-sapply(X=c(1:length(start_date)),FUN=ave.period,fn=mean,xts.v=sw_tmp,simplify=T)
    lw_lt[lequel]<-sapply(X=c(1:length(start_date)),FUN=ave.period,fn=mean,xts.v=lw_tmp,simplify=T)       
    
    
    start_date<-pre_date[lequel]
    end_date<-sos_date[lequel]
    
    # we don't take the mean for the preseaon beacuse the mean will depend of the filtering applied after (homogenous years, etc...)
    tg_pre[lequel]<-sapply(X=c(1:length(start_date)),FUN=pmean,clim=tg_tmp,simplify=T)
    sw_pre[lequel]<-sapply(X=c(1:length(start_date)),FUN=pmean,clim=sw_tmp,simplify=T)
    lw_pre[lequel]<-sapply(X=c(1:length(start_date)),FUN=pmean,clim=lw_tmp,simplify=T)
    
    ######### Chilling and Forcing ###########
    end_date<-sos_date[lequel]
    alti_obs<-LU_db$ALT[lequel]

    start_date<-as.Date(x = paste(c(as.numeric(format(end_date,"%Y"))-1),"-11-01",sep=""),format = "%Y-%m-%d")
    chilling[lequel]<-sapply(X=c(1:length(start_date)),FUN=fchil,xts.v=tg_tmp,alti_CRU=unique(alti[lequel]),alti_obs=alti_obs,simplify=T)
    
    start_date<-pre_date[lequel]
    forcing[lequel]<-sapply(X=c(1:length(start_date)),FUN=fforc5,xts.v=tg_tmp,alti_CRU=unique(alti[lequel]),alti_obs=alti_obs,simplify=T)
    
  }
  
}

################## merge data ######################

new_db<-as.data.frame(cbind(LU_db,tg_lt,sw_lt,lw_lt,
                            tg_pre,tbud_pre,sw_pre,lw_pre,
                            chilling,forcing,chilling_bud,forcing_bud))

names(new_db)<-c(names(LU_db),"TGLT","TG_SP",
                 "SWLT","SW_SP",
                 "LWLT","LW_SP",
                 "TGP","SWP","LWP",
                 "CHIL","FORC","FORC5","DL")

save.image(paste("db_",spname,".RData",sep=""))

