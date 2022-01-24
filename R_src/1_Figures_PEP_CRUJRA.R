############################################################
# 
# Author: Marc Peaucelle, 2020/06
# Description:
# Script to generate Figure3, 4 and supplementary Figures of "Toward better phenology analyses based on bud traits and energy budget" (Peaucelle, Penuelas & Verbeeck, submitted)
# Phenological observations were replaced by -1
# Please register to the PEP website to download phenological observations that are freely available at http://www.pep725.eu/
# Data provided in the databases were generated with the script 0.0_create_database_Tbud.R
# 
############################################################
# function to add error bars with the plot() function
add.error.bars <- function(X,Y,SEX,SEY,w,col=1,lwd=1.5){
  X0 = X; Y0 = (Y-SEY); X1 =X; Y1 = (Y+SEY);
  arrows(X0, Y0, X1, Y1, code=3,angle=90,length=w,col=col,lwd=lwd);
  X0 = (X-SEX); Y0 = Y; X1 =(X+SEX); Y1 = Y;
  arrows(X0, Y0, X1, Y1, code=3,angle=90,length=w,col=col,lwd=lwd);
}

############### First we look at preseason environmental conditions over 1970-2015
# This script runs for Betula pendula. Simply change the dataset to generate figures for other species
db_AG<-readRDS("databases/Alnus_Tbud_1990-2015.rds")
db_AH<-readRDS("databases/Aesculus_Tbud_1990-2015.rds")
db_BP<-readRDS("databases/Betula_Tbud_1990-2015.rds")
db_FS<-readRDS("databases/Fagus_Tbud_1990-2015.rds")
db_FE<-readRDS("databases/Fraxinus_Tbud_1990-2015.rds")
db_QR<-readRDS("databases/Quercus_Tbud_1990-2015.rds")

all_db<-rbind(db_AH,db_AG,db_BP,db_FE,db_FS,db_QR)

#/home/orchidee04/mpeau/Tsol_ERA5

################# Function to compute mean yearly values 
meanVal.fn<-function(db){
  # We remove sites for which climate data (at 0.5°) are NA (i.e. along the coast, sea pixel)
  db<-db[which(!is.na(db$TGP)),]

  # We focus only on the 1990-2015 period. 
  #db<-db[db$YEAR>1989,]
  #db<-db[db$YEAR<2016,]
  # We remove sites with less than 20 years of observations per species. 
 
  ID<-paste0(db$SPECIES,db$PEP_ID)

  leq<-which(table(ID)>20)
  db<-db[ID%in%names(table(ID))[leq],]
 

  db$TBP05<-db$TBP05-273.15
  db$TBP08<-db$TBP08-273.15
 
  # we compute mean value with all sites pooled together
  
  TBP05m<-by(data=db$TBP05,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
  TBP08m<-by(data=db$TBP08,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
  TGPm<-by(data=db$TGP,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
  GDDm<-by(data=db$FORC,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
  GDD08m<-by(data=db$FORCB08,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
  GDD05m<-by(data=db$FORCB05,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
  
  SWPm<-by(data=db$SWP,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
  LWPm<-by(data=db$LWP,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
  TGPm<-by(data=db$TGP,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
  DAYm<-by(data=db$DAY,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
  WINDm<-by(data=db$WINDP,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
  
  dif08<-TBP08m-TGPm
  dif05<-TBP05m-TGPm

  RABS05m<-by(data=db$RABS05,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
  RABS08m<-by(data=db$RABS08,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
  
  H05m<-by(data=db$H05,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
  H08m<-by(data=db$H08,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
  
  Sr05<-db$RABS05-db$H05-db$E05
  Sr08<-db$RABS08-db$H08-db$E08
  
  Sr05m<-by(data=Sr05,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
  Sr08m<-by(data=Sr08,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
 
  E05m<-by(data=db$E05,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
  E08m<-by(data=db$E08,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
 
    
  # We estimate associated sd
  sdBB<-by(data=db$DAY,INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)
  sdSW<-by(data=db$SWP,INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)
  sdLW<-by(data=db$LWP,INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)
  sdTG<-by(data=db$TGP,INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)
  sdWIND<-by(data=db$WINDP,INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T) 
  sd05<-by(data=c(db$TBP05-db$TGP),INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)
  sd08<-by(data=c(db$TBP08-db$TGP),INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)
  sdGDD<-by(data=db$FORC,INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)
  sdGDD08<-by(data=db$FORCB08,INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)
  sdGDD05<-by(data=db$FORCB05,INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)
  
  sdR05<-by(data=c(db$RABS05),INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)
  sdR08<-by(data=c(db$RABS08),INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)
  
  sdH05<-by(data=c(db$H05),INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)
  sdH08<-by(data=c(db$H08),INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)
  
  sdE05<-by(data=c(db$E05),INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)
  sdE08<-by(data=c(db$E08),INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)

  sdS05<-by(data=Sr05,INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)
  sdS08<-by(data=Sr08,INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)
  
  
  return(list(TBP05m=TBP05m,
              TBP08m=TBP08m,
              TGPm=TGPm,
              SWPm=SWPm,
              LWPm=LWPm,
              TGPm=TGPm,
              WINDm=WINDm,
              DAYm=DAYm,
              GDDm=GDDm,
              GDD08m=GDD08m,
              GDD05m=GDD05m,
              dif05=dif05,
              dif08=dif08,
              RABS05m=RABS05m,
              RABS08m=RABS08m,
              H05m=H05m,
              H08m=H08m,
              E05m=E05m,
	      E08m=E08m,
              Sr05m=Sr05m,
              Sr08m=Sr08m,
              sd05=sd05,
              sd08=sd08,
              sdBB=sdBB,
              sdGDD=sdGDD,
              sdGDD08=sdGDD08,
              sdGDD05=sdGDD05,
              sdLW=sdLW,
              sdSW=sdSW,
              sdTG=sdTG,
	      sdWIND=sdWIND,
              sdR05=sdR05,
              sdR08=sdR08,
              sdH05=sdH05,
              sdH08=sdH08,
              sdE05=sdE05,
              sdE08=sdE08,
              sdS05=sdS05,
              sdS08=sdS08))
}


# Analysis performed with all data, can be applied on species data
ALL<-meanVal.fn(all_db)


#################################### Function to plot Figure 3a and supplementary figures
plot.fn<-function(db,axes=F){
  if(axes==FALSE){
    plot(c(db$TBP08m-db$TGPm)~as.numeric(names(db$TGPm)),ylim=c(min(db$dif05-db$sd05/2),max(db$dif08+db$sd08/2)),pch=16,axes=FALSE)
    par(new=T)
    plot(c(db$TBP05m-db$TGPm)~as.numeric(names(db$TGPm)),ylim=c(min(db$dif05-db$sd05/2),max(db$dif08+db$sd08/2)),pch=16,axes=FALSE,col="grey50")
  } else {
    par(mar=c(6,6,2,2))
    plot(c(db$TBP08m-db$TGPm)~as.numeric(names(db$TGPm)),ylim=c(0,max(db$dif08+db$sd08/2)),pch=16,xlab="YEAR",ylab=expression(paste("T"[bud], "-T"[air], " (°C)")),cex=2,cex.axis=3,cex.lab=3)
    par(new=T)
    plot(c(db$TBP05m-db$TGPm)~as.numeric(names(db$TGPm)),ylim=c(0,max(db$dif08+db$sd08/2)),col="grey50",pch=16,cex=2,xlab=NULL,ylab=NULL,xaxt='n',yaxt='n',ann=FALSE)
  }
  add.error.bars(X=as.numeric(names(db$TGPm)),Y=db$dif08,SEX=0,SEY=db$sd08/2,w=0,col="black")
  add.error.bars(X=as.numeric(names(db$TGPm)),Y=db$dif05,SEX=0,SEY=db$sd05/2,w=0,col="grey50")
}



pdf("Figure_3.pdf",width=10,height=10)
plot.fn(ALL,axes = T)
dev.off()


###### We look at site level trends in deltaT
db<-all_db
#db<-db[db$YEAR>1989,]
#db<-db[db$YEAR<2016,]
# We remove sites with less than 20 years of observations. 
db<-db[which(!is.na(db$TGP)),]
ID<-paste0(db$SPECIES,db$PEP_ID)

leq<-which(table(ID)>20)
db<-db[ID%in%names(table(ID))[leq],]

ID<-paste0(db$SPECIES,db$PEP_ID)
# length(unique(db$PEP_ID))
# [1] 2960 sites with >20 years of observations between 1990 and 2015

# we compute yearly average
dbm<-meanVal.fn(all_db)
pdf("Supp_Fig2.pdf",width=10,height=10)
par(mar=c(6,6,2,2))
dev.off()
PEPid<-unique(ID)
dif08<-db$TBP08-db$TGP
lm_list<-list()
leq<-which(ID==PEPid[1])
lm_list[[1]]<-lm(dif08[leq]~db$YEAR[leq])  

# we store the slope, latitude and elevation to see if it explains the observed behavior
sl<-lm_list[[1]]$coefficients[2]
lat<-db$LAT[leq][1]
elv<-db$ALT[leq][1]
for (i in 2:length(PEPid)){
  leq<-which(ID==PEPid[i])
  lm_list[[i]]<-lm(dif08[leq]~db$YEAR[leq]) 
  lat<-c(lat,db$LAT[leq][1])
  elv<-c(elv,db$ALT[leq][1])
  if (summary(lm_list[[i]])$coefficients[8]<0.1){ # we plot only sites with a significant trend, pval<0.1
    sl<-c(sl,lm_list[[i]]$coefficients[2])
  } else {
    sl<-c(sl,NA)
  }
}


# We extract sites with a positive and negative trend in deltaT
posID<-PEPid[which(sl>0)]
negID<-PEPid[which(sl<0)]

posDB<-db[ID%in%posID,]
negDB<-db[ID%in%negID,]

posDBm<-meanVal.fn(posDB)
negDBm<-meanVal.fn(negDB)

pdf("Supp_Fig3.pdf",width=10,height=10)
#ylim=c(0.4,1.8)
ylim=c(0,0.6)

plot(c(dbm$TBP08m-dbm$TGPm)~as.numeric(names(dbm$TGPm)),ylim=ylim,pch=16,xlab=NULL,ylab=NULL,xaxt='n',yaxt='n',ann=FALSE,cex=2)
par(new=T)
plot(c(posDBm$TBP08m-posDBm$TGPm)~as.numeric(names(posDBm$TGPm)),ylim=ylim,pch=16,xlab=NULL,ylab=NULL,xaxt='n',yaxt='n',ann=FALSE,cex=2,col="darkorange")
par(new=T)
plot(c(negDBm$TBP08m-negDBm$TGPm)~as.numeric(names(negDBm$TGPm)),ylim=ylim,pch=16,xlab="YEAR",ylab=expression(paste("T"[bud], "-T"[air], " (°C)")),cex=2,cex.axis=3,cex.lab=3,col="dodgerblue")
add.error.bars(X=c(1990:2015),Y=c(dbm$TBP08m-dbm$TGPm),SEX=0,SEY=dbm$sd08/2,w=0,col="black")
add.error.bars(X=c(1990:2015),Y=c(posDBm$TBP08m-posDBm$TGPm),SEX=0,SEY=posDBm$sd08/2,w=0,col="darkorange")
add.error.bars(X=c(1990:2015),Y=c(negDBm$TBP08m-negDBm$TGPm),SEX=0,SEY=negDBm$sd08/2,w=0,col="dodgerblue")

lm1<-lm(c(dbm$TBP08m-dbm$TGPm)~c(1990:2015))
summary(lm1)
segments(x0=1990,y0=(coef(lm1)[1]+1990*coef(lm1)[2]),x1=2015,y1=(coef(lm1)[1]+2015*coef(lm1)[2]),col="black",lwd=3)
lm1<-lm(c(posDBm$TBP08m-posDBm$TGPm)~c(1990:2015))
summary(lm1)
segments(x0=1990,y0=(coef(lm1)[1]+1990*coef(lm1)[2]),x1=2015,y1=(coef(lm1)[1]+2015*coef(lm1)[2]),col="darkorange",lwd=3)
lm1<-lm(c(negDBm$TBP08m-negDBm$TGPm)~c(1990:2015))
summary(lm1)
segments(x0=1990,y0=(coef(lm1)[1]+1990*coef(lm1)[2]),x1=2015,y1=(coef(lm1)[1]+2015*coef(lm1)[2]),col="dodgerblue",lwd=3)

dev.off()
# number of site*species with >20y of data
#> length(sl)
#[1] 5050
# number of site*species with positive deltaT trend
#> length(which(sl>0))
#[1] 902  # 18% of sites
# number of site*species with negative deltaT trend
#> length(which(sl<0))
#[1] 356 # 7% of sites 


# check effect of latitude and elevation
plot(sl~lat)
plot(sl~elv)
# latitude and elevation do not explain the slope. 

# plot sites Supp Figure 2

library(ggplot2)
library("rnaturalearth")
library("rnaturalearthdata")

coord<-unique(cbind(db$LON,db$LAT))
world <- ne_countries(scale = "medium", returnclass = "sf")


sites <- data.frame(longitude = coord[,1], latitude = coord[,2])

pdf("Supplementary_Figure2.pdf",width=12,height=10)
ggplot(data = world) + theme_bw() +
    geom_sf() +
    geom_point(data = sites, aes(x = longitude, y = latitude), size = 2,
        shape = 23, fill = "darkred") +
    coord_sf(xlim = c(-10, 30), ylim = c(40, 60), expand = FALSE)
dev.off()


########### Finally we look at the components of the energy budget in explaining the variability in deltaT

db<-ALL
db<-negDBm

dif08<-c((db$TBP08m-db$TGPm))
df<-as.data.frame(cbind(DIF=dif08,RABS=db$RABS08m,H=db$H08m,E=db$E08m,SR=db$Sr08m))
###### multiple linear regression 
#lm_m<-glm(DIF~RABS+H+SR+E,data = df)
lm_m<-lm(DIF~RABS+H+SR+E,data = df)

summary(lm_m)
pdf("Supplementary_Figure_4.pdf",width=10,height=10)
par(mar=c(6,6,2,2))
par(mfrow=c(2,2))
ylim=c(0.6,1.4)
#plot(db$RABS08m/500~as.numeric(names(db$TGPm)),ylim=ylim,pch=16,xlab="YEAR",ylab="Value",cex=2,cex.axis=2,cex.lab=2,col="darkred")

plot(scale(db$RABS08m,center=F)~as.numeric(names(db$TGPm)),ylim=ylim,pch=16,xlab="YEAR",ylab="Value",cex=2,cex.axis=2,cex.lab=2,col="darkred")

par(new=T)
#plot(db$H08m/50~as.numeric(names(db$TGPm)),ylim=ylim,col="darkorange",pch=16,cex=2,xlab=NULL,ylab=NULL,xaxt='n',yaxt='n',ann=FALSE)
plot(scale(db$H08m,center=F)~as.numeric(names(db$TGPm)),ylim=ylim,col="darkorange",pch=16,cex=2,xlab=NULL,ylab=NULL,xaxt='n',yaxt='n',ann=FALSE)

par(new=T)
#plot(db$Sr08m/500~as.numeric(names(db$TGPm)),ylim=ylim,col="dodgerblue3",pch=16,cex=2,xlab=NULL,ylab=NULL,xaxt='n',yaxt='n',ann=FALSE)
plot(scale(db$Sr08m,center=F)~as.numeric(names(db$TGPm)),ylim=ylim,col="dodgerblue3",pch=16,cex=2,xlab=NULL,ylab=NULL,xaxt='n',yaxt='n',ann=FALSE)

#par(new=T)
#plot(db$E08m/20~as.numeric(names(db$TGPm)),ylim=ylim,col="darkgreen",pch=16,cex=2,xlab=NULL,ylab=NULL,xaxt='n',yaxt='n',ann=FALSE)
#plot(scale(db$E08m,center=F)~as.numeric(names(db$TGPm)),ylim=ylim,col="darkgreen",pch=16,cex=2,xlab=NULL,ylab=NULL,xaxt='n',yaxt='n',ann=FALSE)

par(new=T)
#plot(db$dif08~as.numeric(names(db$TGPm)),ylim=ylim,col="black",pch=16,cex=2,xlab=NULL,ylab=NULL,xaxt='n',yaxt='n',ann=FALSE)
plot(scale(db$dif08,center=F)~as.numeric(names(db$TGPm)),ylim=ylim,col="black",pch=16,cex=2,xlab=NULL,ylab=NULL,xaxt='n',yaxt='n',ann=FALSE)

add.error.bars(X=c(1990:2015),Y=db$RABS08m/500,SEX=0,SEY=db$sdR08/1000,w=0,col="darkred")
add.error.bars(X=c(1990:2015),Y=db$H08m/50,SEX=0,SEY=db$sdH08/100,w=0,col="darkorange")
add.error.bars(X=c(1990:2015),Y=db$Sr08m/500,SEX=0,SEY=db$sdS08/1000,w=0,col="dodgerblue3")
#add.error.bars(X=c(1990:2015),Y=db$E08m/20,SEX=0,SEY=db$sdE08/40,w=0,col="darkgreen")
add.error.bars(X=c(1990:2015),Y=c((db$TBP08m-db$TGPm)),SEX=0,SEY=db$sd08/2,w=0,col="black")


termplot(lm_m,partial.resid = T,col.res = "grey30", se=T,col.se = "darkred",col.term = "darkred",pch=1,cex=2.5,lwd.se=2.5,lwd.term = 2.5,cex.axis=2,cex.lab=2,xlabs=c("Rabs (W m-2)","H (W m-2)","LWbud (W m-2)","E (W m-2)"),ylabs=rep("Part. res. (°C)"))
dev.off()

