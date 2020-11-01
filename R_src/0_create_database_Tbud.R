############################################################
# Create database script
# Author: Marc Peaucelle, 2020/06
# Description:
# This script prepares the dataset by computing preseason bud temperature from its energy balance and environmental conditions. 
# Bud temperature is estimated for phenological observations from the PEP database: http://www.pep725.eu/
# Meteorological data from CRU-JRA are freely avilable at: https://catalogue.ceda.ac.uk/uuid/13f3635174794bb98cf8ac4b0ee8f4ed
# Raw phenology data and meteorological data are not provided within this script. 
# Both PEP data and meteorological are freely abailable after registration and are protected by a common licence. 
# Still, a template file with the structure of PEP data working with this script is provided at https://github.com/mpeaucelle/Tbud 

# References: 
#     - Landsberg, J. J., Butler, D. R. & Thorpe, M. R. Apple bud and blossom temperatures. Journal of Horticultural Science 49, 227–239; 10.1080/00221589.1974.11514574 (1974).
#     - HAMER, P. The heat balance of apple buds and blossoms. Part I. Heat transfer in the outdoor environment. Agricultural and Forest Meteorology 35, 339–352; 10.1016/0168-1923(85)90094-2 (1985).
#     - Muir, C. D. tealeaves: an R package for modelling leaf temperature using energy budgets. AoB PLANTS 11, plz054; 10.1093/aobpla/plz054 (2019).
#     - Monteith, J. L. & Unsworth, M. H. Principles of environmental physics. Plants, animals, and the atmosphere. 4th ed. (Elsevier/Academic Press, Amsterdam, Boston, 2013).
#     - Jones, H. G. Plants and microclimate. A quantitative approach to environmental plant physiology (Cambridge university press, Cambridge, 2013). --> provides absroptivity values for leaves and needles
#     - University Of East Anglia Climatic Research Unit (CRU) & Harris, I. C. CRU JRA v1.1: A forcings dataset of gridded land surface blend of Climatic Research Unit (CRU) and Japanese reanalysis (JRA) data; Jan.1901 - Dec.2017, 2019.

############################################################
#### if running for several species with a batch script
#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly=TRUE)

rm(list=ls())
library(RNetCDF)
library(raster)
library(xts)

source("Tbud.R")
source("Ebalance.R")
source("functions.R")

########### LOAD pheno data ##############
# PEP725
# Script adapted for B.pendula only 
# Adapt these two lines according to your environment. 
workdir<-"XXX workdir here XXX"
i<-paste0(workdir,"betula_pendula_database.txt") # see PEP_database_template.txt for the structure

#### if running for several species with a batch script
#i<-paste0(workdir,"db_site/",args[1])
#spname<-args[2]
spname<-"betula"
# We restrict the analysis to the 1990-2016 period
first_year<-1970
last_year<-2017

# Load observation database
obs_db<-read.table(i,sep=";",header=T)

# Keep only leaf unfolding observation
LU_db<-obs_db[obs_db$BBCH==11,]

# Kepp the selected period
LU_db<-LU_db[order(LU_db$PEP_ID,LU_db$YEAR),]
LU_db<-LU_db[LU_db$YEAR>c(first_year-1),]

########### LOAD CRU data #######
# Here again a template file is provided at https://github.com/mpeaucelle/Tbud 
# For the sake of vectorization, crujra data are here in 1 dimension. 
nc_tg<-open.nc(paste0(workdir,"netcdf_forcing/crujra_1D/crujra_tair_1969_2018.nc"))
nc_sw<-open.nc(paste0(workdir,"netcdf_forcing/crujra_1D/crujra_swdown_1969_2018.nc"))
nc_lw<-open.nc(paste0(workdir,"netcdf_forcing/crujra_1D/crujra_lwdown_1969_2018.nc"))
nc_qair<-open.nc(paste0(workdir,"netcdf_forcing/crujra_1D/crujra_qair_1969_2018.nc"))
nc_uw<-open.nc(paste0(workdir,"netcdf_forcing/crujra_1D/crujra_uwind_1969_2018.nc"))
nc_vw<-open.nc(paste0(workdir,"netcdf_forcing/crujra_1D/crujra_vwind_1969_2018.nc"))
nc_pb<-open.nc(paste0(workdir,"netcdf_forcing/crujra_1D/crujra_press_1969_2018.nc"))

lat<-var.get.nc(nc_tg,"latitude")
lon<-var.get.nc(nc_tg,"longitude")
time<-1:18250 #in days since 1950-01-01 00:00
time<-as.Date(time,origin="1968-12-31")
# CRU data are based on a 365d calendar, thus we artificially add the missing days. 
bisextile<-which(as.character(time)%in%c("1971-12-31","1975-12-31","1979-12-31","1983-12-31",
                                         "1987-12-31","1991-12-31","1995-12-31","1999-12-31",
                                         "2003-12-31","2007-12-31","2011-12-31","2015-12-31"))
time<-1:18262
time<-as.Date(time,origin="1968-12-31")
# Terrestrial mask to convert 1D data to 2D data
themask<-var.get.nc(nc_tg,"mask")

# mean elevation of CRUJRA pixels for temperature correction
#alt<-open.nc(paste0(workdir,"netcdf_forcing/crujra_alti.nc"))
alt<-open.nc(paste0(workdir,"netcdf_forcing/force1994_2002.nc"))
alti<-var.get.nc(alt,"altitude")

rast_mask<-raster(t(themask),xmn=min(lon),xmx=max(lon),ymn=min(lat),ymx=max(lat))
rast_alt<-raster(t(alti),xmn=min(lon),xmx=max(lon),ymn=min(lat),ymx=max(lat))

print("fichiers CRU chargés")

### extract site ref-point 
# Allow to extract 1D netcdf data 

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

# just for information, computes long-term average values of climate conditions. 
# This part was not used to generate figures. 
# initialize variables
# longterm annual tg, pp, sw, lw
tg_lt<-sw_lt<-lw_lt<-NULL
# preseason  tg, sw, lw
tg_pre<-sw_pre<-lw_pre<-wind_pre<-NULL

# preseason bud temperature and associated energy budget components for an absorptivity of 0.8 and 0.5
tbud08_pre<-tbud05_pre<-NULL
rabs08_pre<-rabs05_pre<-NULL
H08_pre<-H05_pre<-NULL

# preseason min and max temperature
tg_min<-tg_max<-NULL
#preseason min and max temperature and energy budget. 
tbud08_min<-tbud05_min<-tbud08_max<-tbud05_max<-NULL
rabs08_min<-rabs05_min<-rabs08_max<-rabs05_max<-NULL
H08_min<-H05_min<-H08_max<-H05_max<-NULL


#CRU altitude
alti<-NULL
#chilling and forcing
chilling<-forcing<-forcing_bud08<-chilling_bud08<-forcing_bud05<-chilling_bud05<-NULL

refpoint<-as.character(unique(LU_db$REF))

if(length(which(refpoint=="-1"))>0){
  refpoint<-refpoint[-which(refpoint=="-1")]
}


# computes Tbud for each site-year

for (j in 1:length(refpoint)){
  print(paste(j,"/",length(refpoint)))
  leq<-LU_db$REF%in%refpoint[j]
  lequel<-which(leq==TRUE)

  # on which pixel is the plot   
  cell<-which(t(as.matrix(rast_mask))==as.numeric(refpoint[j]))
  alti[lequel]<-rast_alt[cell]
  # extract netcdf data for the correponding site
  # We add +1 to refpoint because refpoint is 0:67419, but var.get.nc start at 1 and not 0.
  # in Kelvin -273.15 --> °C
  tg_tmp<-var.get.nc(nc_tg,"Temperature",start=c(as.numeric(refpoint[j])+1,1),count=c(1,18250))-273.15
  # SW is in J/m²/6h *21600 --> W m-2
  sw_tmp<-var.get.nc(nc_sw,"Incoming_Short_Wave_Radiation",start=c(as.numeric(refpoint[j])+1,1),count=c(1,18250))/21600
  lw_tmp<-var.get.nc(nc_lw,"Incoming_Long_Wave_Radiation",start=c(as.numeric(refpoint[j])+1,1),count=c(1,18250))
  qair_tmp<-var.get.nc(nc_qair,"Air_Specific_Humidity",start=c(as.numeric(refpoint[j])+1,1),count=c(1,18250))
  pb_tmp<-var.get.nc(nc_pb,"Pression",start=c(as.numeric(refpoint[j])+1,1),count=c(1,18250))
  uw_tmp<-var.get.nc(nc_uw,"U_wind_component",start=c(as.numeric(refpoint[j])+1,1),count=c(1,18250))
  vw_tmp<-var.get.nc(nc_vw,"V_wind_component",start=c(as.numeric(refpoint[j])+1,1),count=c(1,18250))
  
  # computes wind moment
  wind_tmp<-sqrt (uw_tmp*uw_tmp + vw_tmp*vw_tmp)
  rm(uw_tmp,vw_tmp)
  
  ## add manually bisextile years unvailbale in CRU
  tg_tmp<-bis.fn(xts.v = tg_tmp,bisextile = bisextile)
  sw_tmp<-bis.fn(xts.v = sw_tmp,bisextile = bisextile)
  lw_tmp<-bis.fn(xts.v = lw_tmp,bisextile = bisextile)
  qair_tmp<-bis.fn(xts.v = qair_tmp,bisextile = bisextile)
  pb_tmp<-bis.fn(xts.v = pb_tmp,bisextile = bisextile)
  wind_tmp<-bis.fn(xts.v = wind_tmp,bisextile = bisextile)
  
  
  ## convert in time-series  
  sw_tmp<-xts(sw_tmp,time)
  lw_tmp<-xts(lw_tmp,time)
  tg_tmp<-xts(tg_tmp,time)
  qair_tmp<-xts(qair_tmp,time)
  pb_tmp<-xts(pb_tmp,time)
  wind_tmp<-xts(wind_tmp,time)
  
  
  # computes preseason mean bud temperature and energy budget
  tg_bud08<-sapply(X=c(1:length(tg_tmp)),FUN=fTbud,abs_sw=0.8,simplify=T)
  tg_bud05<-sapply(X=c(1:length(tg_tmp)),FUN=fTbud,abs_sw=0.5,simplify=T)
  
  rabs08<-sapply(X=c(1:length(tg_tmp)),FUN=fEbal,abs_sw=0.8,ind=1,Tbud=tg_bud08,simplify=T)
  rabs05<-sapply(X=c(1:length(tg_tmp)),FUN=fEbal,abs_sw=0.5,ind=1,Tbud=tg_bud05,simplify=T)
  
  H08<-sapply(X=c(1:length(tg_tmp)),FUN=fEbal,abs_sw=0.8,ind=3,Tbud=tg_bud08,simplify=T)
  H05<-sapply(X=c(1:length(tg_tmp)),FUN=fEbal,abs_sw=0.5,ind=3,Tbud=tg_bud05,simplify=T)
  
  
  tg_bud08<-xts(tg_bud08,time)
  tg_bud05<-xts(tg_bud05,time)
  
  rabs08<-xts(rabs08,time)
  rabs05<-xts(rabs05,time)
  H08<-xts(H08,time)
  H05<-xts(H05,time)
  
  
  ####### if obs is not in the sea? :o) sometimes at the land border in CRU-JRA ######
  if( all(is.na(tg_tmp))==FALSE){
    
    ######### longterm climate
    start_date<-paste(first_year,"-01-01",sep="")
    end_date<-paste(last_year,"-01-01",sep="")
    
    tg_lt[lequel]<-sapply(X=c(1:length(start_date)),FUN=ave.period,fn=mean,xts.v=tg_tmp,simplify=T)
    sw_lt[lequel]<-sapply(X=c(1:length(start_date)),FUN=ave.period,fn=mean,xts.v=sw_tmp,simplify=T)
    lw_lt[lequel]<-sapply(X=c(1:length(start_date)),FUN=ave.period,fn=mean,xts.v=lw_tmp,simplify=T)       
    
    
    start_date<-pre_date[lequel]
    end_date<-sos_date[lequel]
    
    tg_pre[lequel]<-sapply(X=c(1:length(start_date)),FUN=pmean,clim=tg_tmp,simplify=T)
    sw_pre[lequel]<-sapply(X=c(1:length(start_date)),FUN=pmean,clim=sw_tmp,simplify=T)
    lw_pre[lequel]<-sapply(X=c(1:length(start_date)),FUN=pmean,clim=lw_tmp,simplify=T)
    wind_pre[lequel]<-sapply(X=c(1:length(start_date)),FUN=pmean,clim=wind_tmp,simplify=T)
 
    tbud08_pre[lequel]<-sapply(X=c(1:length(start_date)),FUN=pmean,clim=tg_bud08,simplify=T)
    tbud05_pre[lequel]<-sapply(X=c(1:length(start_date)),FUN=pmean,clim=tg_bud05,simplify=T)    
    
    rabs08_pre[lequel]<-sapply(X=c(1:length(start_date)),FUN=pmean,clim=rabs08,simplify=T)
    rabs05_pre[lequel]<-sapply(X=c(1:length(start_date)),FUN=pmean,clim=rabs05,simplify=T)
    
    H08_pre[lequel]<-sapply(X=c(1:length(start_date)),FUN=pmean,clim=H08,simplify=T)
    H05_pre[lequel]<-sapply(X=c(1:length(start_date)),FUN=pmean,clim=H05,simplify=T)
    
    
    #preseason longterm tg, sw, lw
    tg_min[lequel]<-sapply(X=c(1:length(start_date)),FUN=pmin,clim=tg_tmp,simplify=T)
    tg_max[lequel]<-sapply(X=c(1:length(start_date)),FUN=pmax,clim=tg_tmp,simplify=T)
    tbud08_min[lequel]<-sapply(X=c(1:length(start_date)),FUN=pmin,clim=tg_bud08,simplify=T)
    tbud05_min[lequel]<-sapply(X=c(1:length(start_date)),FUN=pmin,clim=tg_bud05,simplify=T)
    tbud08_max[lequel]<-sapply(X=c(1:length(start_date)),FUN=pmax,clim=tg_bud08,simplify=T)
    tbud05_max[lequel]<-sapply(X=c(1:length(start_date)),FUN=pmax,clim=tg_bud05,simplify=T)
    
    rabs08_min[lequel]<-sapply(X=c(1:length(start_date)),FUN=pmin,clim=rabs08,simplify=T)
    rabs05_min[lequel]<-sapply(X=c(1:length(start_date)),FUN=pmin,clim=rabs05,simplify=T)
    rabs08_max[lequel]<-sapply(X=c(1:length(start_date)),FUN=pmax,clim=rabs08,simplify=T)
    rabs05_max[lequel]<-sapply(X=c(1:length(start_date)),FUN=pmax,clim=rabs05,simplify=T)
    
    H08_min[lequel]<-sapply(X=c(1:length(start_date)),FUN=pmin,clim=H08,simplify=T)
    H05_min[lequel]<-sapply(X=c(1:length(start_date)),FUN=pmin,clim=H05,simplify=T)
    H08_max[lequel]<-sapply(X=c(1:length(start_date)),FUN=pmax,clim=H08,simplify=T)
    H05_max[lequel]<-sapply(X=c(1:length(start_date)),FUN=pmax,clim=H05,simplify=T)
    
    ######### Chilling and Forcing ###########
    end_date<-sos_date[lequel]
    alti_obs<-LU_db$ALT[lequel]
    
    start_date<-as.Date(x = paste(c(as.numeric(format(end_date,"%Y"))-1),"-11-01",sep=""),format = "%Y-%m-%d")
    chilling[lequel]<-sapply(X=c(1:length(start_date)),FUN=fchil,xts.v=tg_tmp,alti_CRU=unique(alti[lequel]),alti_obs=alti_obs,simplify=T)
    chilling_bud08[lequel]<-sapply(X=c(1:length(start_date)),FUN=fchil,xts.v=tg_bud08,alti_CRU=unique(alti[lequel]),alti_obs=alti_obs,simplify=T)
    chilling_bud05[lequel]<-sapply(X=c(1:length(start_date)),FUN=fchil,xts.v=tg_bud05,alti_CRU=unique(alti[lequel]),alti_obs=alti_obs,simplify=T)
    
    start_date<-pre_date[lequel]
    forcing[lequel]<-sapply(X=c(1:length(start_date)),FUN=fforc5,xts.v=tg_tmp,alti_CRU=unique(alti[lequel]),alti_obs=alti_obs,simplify=T)
    forcing_bud08[lequel]<-sapply(X=c(1:length(start_date)),FUN=fforc5,xts.v=tg_bud08,alti_CRU=unique(alti[lequel]),alti_obs=alti_obs,simplify=T)
    forcing_bud05[lequel]<-sapply(X=c(1:length(start_date)),FUN=fforc5,xts.v=tg_bud05,alti_CRU=unique(alti[lequel]),alti_obs=alti_obs,simplify=T)
    
  }
  
}

################## merge data ######################

new_db<-as.data.frame(cbind(LU_db,tg_lt,sw_lt,lw_lt,
                            tg_pre,tbud08_pre,tbud05_pre,sw_pre,lw_pre,wind_pre,
                            tg_min,tg_max,
                            tbud08_min,tbud08_max,
                            tbud05_min,tbud05_max,
                            chilling,chilling_bud08,chilling_bud05,
                            forcing,forcing_bud08,forcing_bud05,
                            rabs08_pre,rabs05_pre,H08_pre,H05_pre,
                            rabs08_min,rabs05_min,H08_min,H05_min,
                            rabs08_max,rabs05_max,H08_max,H05_max))

names(new_db)<-c(names(LU_db),"TGLT","SWLT","LWLT",
                 "TGP","TBP08","TBP05",
                 "SWP","LWP","WINDP",
                 "TPMIN","TPMAX",
                 "TBP08_MIN","TBP08_MAX",
                 "TBP05_MIN","TBP05_MAX",
                 "CHIL","CHILB08","CHILB05",
                 "FORC","FORCB08","FORCB05",
                 "RABS08","RABS05","H08","H05",
                 "RABS08_MIN","RABS05_MIN","H08_MIN","H05_MIN",
                 "RABS08_MAX","RABS05_MAX","H08_MAX","H05_MAX")


saveRDS(object = new_db,paste0(spname,"_Tbud.rds"),compress = TRUE)
write.csv(new_db,paste0(spname,"_Tbud.csv"),col.names=T,sep=";")

q(save='n')
