############################################################
# 
# Author: Marc Peaucelle, 2020/06
# Description:
# Script to generate Figure3 of "Plants are warming faster than climate" (Peaucelle, Penuelas & Verbeeck, 2020)
# Phenological observations were replaced in the Betula database were replaced by -1
# Please register to the PEP website to download phenological observations that are freely available at http://www.pep725.eu/
# Data provided in the databases were generated with the script 0.0_create_database_Tbud.R
# 
############################################################
library(ggplot2)
library(gridExtra)
# function to add error bars with the plot() function
add.error.bars <- function(X,Y,SEX,SEY,w,col=1,lwd=1.5){
  X0 = X; Y0 = (Y-SEY); X1 =X; Y1 = (Y+SEY);
  arrows(X0, Y0, X1, Y1, code=3,angle=90,length=w,col=col,lwd=lwd);
  X0 = (X-SEX); Y0 = Y; X1 =(X+SEX); Y1 = Y;
  arrows(X0, Y0, X1, Y1, code=3,angle=90,length=w,col=col,lwd=lwd);
}

############### First we look at preseason environmental conditions over 1970-2015
db<-readRDS("database/Betula_Tbud1970-2016.rds")

table(db$YEAR)
# We see that we only have 30 sites for 2016, we remove it. 
db<-db[db$YEAR<2016,]
# We remove sites with less than 10 years of observations. 
leq<-which(table(db$PEP_ID)>10)
db<-db[db$PEP_ID%in%names(table(db$PEP_ID))[leq],]

# we compute mean value with all sites pooled together

TBP05m<-by(data=db$TBP05,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
TBP08m<-by(data=db$TBP08,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
TGPm<-by(data=db$TGP,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)

SWPm<-by(data=db$SWP,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
LWPm<-by(data=db$LWP,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
TGPm<-by(data=db$TGP,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
DAYm<-by(data=db$DAY,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)

dif08<-TBP08m-TGPm
dif05<-TBP05m-TGPm

# We estimate associated sd
sdBB<-by(data=db$DAY,INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)
sdSW<-by(data=db$SWP,INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)
sdLW<-by(data=db$LWP,INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)
sdTG<-by(data=db$TGP,INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)
sd05<-by(data=c(db$TBP05-db$TGP),INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)
sd08<-by(data=c(db$TBP08-db$TGP),INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)


###### Now we plot environmental conditions and average budburst date with all sites pooled together
pdf("Supp_Fig1.pdf",width = 10,height = 10)
par(mfrow=c(2,2),mar=c(5,5,2,2))
plot(DAYm~as.numeric(names(DAYm)),pch=16,xlab="Year",ylab="Budburst date (d)",cex=2,cex.axis=2,cex.lab=2,col="black",ylim=c(90,130))
add.error.bars(X=as.numeric(names(sdBB)),Y=DAYm,SEX=0,SEY=sdBB/2,w=0,col="black")
plot(SWPm~as.numeric(names(SWPm)),pch=16,xlab="Year",ylab="Preseason SW (W m-2)",cex=2,cex.axis=2,cex.lab=2,col="black",ylim=c(90,150))
add.error.bars(X=as.numeric(names(sdSW)),Y=SWPm,SEX=0,SEY=sdSW/2,w=0,col="black")
plot(LWPm~as.numeric(names(LWPm)),pch=16,xlab="Year",ylab=" Preseason LW (W m-2)",cex=2,cex.axis=2,cex.lab=2,col="black", ylim=c(250,300))
add.error.bars(X=as.numeric(names(sdLW)),Y=LWPm,SEX=0,SEY=sdLW/2,w=0,col="black")
plot(TGPm~as.numeric(names(TGPm)),pch=16,xlab="Year",ylab=" Preseason Temp. (C)",cex=2,cex.axis=2,cex.lab=2,col="black",ylim=c(1,6))
add.error.bars(X=as.numeric(names(sdTG)),Y=TGPm,SEX=0,SEY=sdTG/2,w=0,col="black")
dev.off()

pdf("Supp_Fig2.pdf",width = 10,height = 10)
plot(dif08~as.numeric(names(TGPm)),ylim=c(0,1.3),pch=16,xlab="YEAR",ylab=expression(paste("T"[bud], "-T"[air], " (°C)")),cex=2,cex.axis=3,cex.lab=3)
par(new=T)
plot(dif05~as.numeric(names(TGPm)),ylim=c(0,1.3),col="grey50",pch=16,cex=2,xlab=NULL,ylab=NULL,xaxt='n',yaxt='n',ann=FALSE)
add.error.bars(X=as.numeric(names(dif08)),Y=dif08,SEX=0,SEY=sd08/2,w=0,col="black")
add.error.bars(X=as.numeric(names(dif05)),Y=dif05,SEX=0,SEY=sd05/2,w=0,col="grey50")
dev.off()


###### We see from the above Figure an abrupt change in budburst dynamics between 1970-1990 and 1990-2015, inducing an abrupt change in preseason radiation. 
###### This is the reason why we investigate the 1990-2015 period from now. 
db<-readRDS("database/Betula_Tbud1970-2016.rds")
table(db$YEAR)
# We see that we only have 30 sites for 2016, we remove it. 
db<-db[db$YEAR<2016,]
db<-db[db$YEAR>1989,]
# We remove sites with less than 10 observations and not 15 since we work on a shorter period 
leq<-which(table(db$PEP_ID)>10)
db<-db[db$PEP_ID%in%names(table(db$PEP_ID))[leq],]


# We estimate again yearly mean variables
# we compute mean value with all sites pooled together

TBP05m<-by(data=db$TBP05,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
TBP08m<-by(data=db$TBP08,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
TGPm<-by(data=db$TGP,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
dif08<-TBP08m-TGPm
dif05<-TBP05m-TGPm

RABS05m<-by(data=db$RABS05,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
RABS08m<-by(data=db$RABS08,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)

H05m<-by(data=db$H05,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
H08m<-by(data=db$H08,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)

Sr05<-db$RABS05-db$H05
Sr08<-db$RABS08-db$H08

Sr05m<-by(data=Sr05,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)
Sr08m<-by(data=Sr08,INDICES=as.numeric(db$YEAR),FUN=mean,na.rm=T)

# We estimate associated sd
sd05<-by(data=c(db$TBP05-db$TGP),INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)
sd08<-by(data=c(db$TBP08-db$TGP),INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)

sdR05<-by(data=c(db$RABS05),INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)
sdR08<-by(data=c(db$RABS08),INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)

sdH05<-by(data=c(db$H05),INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)
sdH08<-by(data=c(db$H08),INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)

sdS05<-by(data=Sr05,INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)
sdS08<-by(data=Sr08,INDICES=as.numeric(db$YEAR),FUN=sd,na.rm=T)

##### linear regression
lm1<-lm(dif08~c(1990:2015))
lm2<-lm(dif05~c(1990:2015))

summary(lm1)
summary(lm2)

################# We look at the global relationship with all sites pooled togethers
pdf("Figure_3a.pdf",width=10,height=10)
par(mar=c(6,6,2,2))
plot(c(TBP08m-TGPm)~as.numeric(names(TGPm)),ylim=c(0,1.3),pch=16,xlab="YEAR",ylab=expression(paste("T"[bud], "-T"[air], " (°C)")),cex=2,cex.axis=3,cex.lab=3)
par(new=T)
plot(c(TBP05m-TGPm)~as.numeric(names(TGPm)),ylim=c(0,1.3),col="grey50",pch=16,cex=2,xlab=NULL,ylab=NULL,xaxt='n',yaxt='n',ann=FALSE)
add.error.bars(X=c(1990:2016),Y=dif08,SEX=0,SEY=sd08/2,w=0,col="black")
add.error.bars(X=c(1990:2016),Y=dif05,SEX=0,SEY=sd05/2,w=0,col="grey50")

segments(x0=1990,y0=(coef(lm1)[1]+1990*coef(lm1)[2]),x1=2016,y1=(coef(lm1)[1]+2016*coef(lm1)[2]),col="black",lwd=3)
segments(x0=1990,y0=(coef(lm2)[1]+1990*coef(lm2)[2]),x1=2016,y1=(coef(lm2)[1]+2016*coef(lm2)[2]),col="grey50",lwd=3)
dev.off()

################# We look at the relationship site by site

PEPid<-unique(db$PEP_ID)
dif08<-db$TBP08-db$TGP
lm_list<-list()
leq<-which(db$PEP_ID==PEPid[1])
lm_list[[1]]<-lm(dif08[leq]~db$YEAR[leq])  
pdf("Figure_3b.pdf",width=10,height=10)
par(mar=c(6,6,2,2))
plot(dif08[leq]~db$YEAR[leq], col=NULL,xlab=NULL,ylab=NULL,xaxt='n',yaxt='n',ann=FALSE,ylim=c(-2,5),xlim=c(1970,2015))
abline(lm_list[[1]],col="darkred")
# we store the slope, latitude and elevation to see if it explains the observed behavior
sl<-lm_list[[1]]$coefficients[2]
lat<-db$LAT[leq][1]
elv<-db$ALT[leq][1]
for (i in 2:length(PEPid)){
  leq<-which(db$PEP_ID==PEPid[i])
  lm_list[[i]]<-lm(dif08[leq]~db$YEAR[leq]) 
  lat<-c(lat,db$LAT[leq][1])
  elv<-c(elv,db$ALT[leq][1])
  if (summary(lm_list[[i]])$coefficients[8]<0.1){ # we plot only sites with pval<0.1
    sl<-c(sl,lm_list[[i]]$coefficients[2])
    if(lm_list[[i]]$coefficients[2]<0){
      abline(lm_list[[i]],col="darkblue")  
    }else {
      abline(lm_list[[i]],col="darkred")  
    }
    
 } else {
    sl<-c(sl,NA)
  }
  
}
axis(1,cex.axis=2)
axis(2,cex.axis=2)
title(xlab="Year",
      ylab="deltaT (°C)",cex.lab=2)
##### check slope distribution
hist(sl)
# most site have a positive slope
dev.off()

# check effect of latitude and elevation
plot(sl~lat)
plot(sl~elv)
# latitude and elevation do not explain the slope. 

#### same figure but with all the plots, not only the one with p.val<0.1
PEPid<-unique(db$PEP_ID)
dif08<-db$TBP08-db$TGP
lm_list<-list()
leq<-which(db$PEP_ID==PEPid[1])
lm_list[[1]]<-lm(dif08[leq]~db$YEAR[leq])  
pdf("Supp_Fig3.pdf",width=10,height=10)
par(mar=c(6,6,2,2))
plot(dif08[leq]~db$YEAR[leq], col=NULL,xlab=NULL,ylab=NULL,xaxt='n',yaxt='n',ann=FALSE,ylim=c(-2,5),xlim=c(1970,2015))
abline(lm_list[[1]],col="darkred")
# we store the slope, latitude and elevation to see if it explains the observed behavior
sl<-lm_list[[1]]$coefficients[2]
lat<-db$LAT[leq][1]
elv<-db$ALT[leq][1]
for (i in 2:length(PEPid)){
  leq<-which(db$PEP_ID==PEPid[i])
  lm_list[[i]]<-lm(dif08[leq]~db$YEAR[leq])  
  lat<-c(lat,db$LAT[leq][1])
  elv<-c(elv,db$ALT[leq][1])
  sl<-c(sl,lm_list[[i]]$coefficients[2])
    if(lm_list[[i]]$coefficients[2]<0){
      abline(lm_list[[i]],col="darkblue")  
    }else {
      abline(lm_list[[i]],col="darkred")  
    }
}

axis(1,cex.axis=2)
axis(2,cex.axis=2)
title(xlab="Year",
      ylab="deltaT (°C)",cex.lab=2)

##### check slope distribution
hist(sl)
# Most sites have a positive slope
dev.off()

# check effect of latitude and elevation
plot(sl~lat)
plot(sl~elv)
# latitude and elevation do not explain the slope. 


########### Finally we look at the components of the energy budget in explaining the variability in deltaT

dif08<-c((TBP08m-TGPm))
###### multiple linear regression 
lm_m<-glm(dif08~RABS08m+H08m+Sr08m+TGPm)
lm_m<-glm(TBP08m~RABS08m+H08m+Sr08m+TGPm)

summary(lm_m)
pdf("Figure_4.pdf",width=10,height=10)
par(mar=c(6,6,2,2))

plot(RABS08m/500~as.numeric(names(TGPm)),ylim=c(0,1.3),pch=16,xlab="YEAR",ylab="Value",cex=2,cex.axis=2,cex.lab=2,col="darkred")
par(new=T)
plot(H08m/100~as.numeric(names(TGPm)),ylim=c(0,1.3),col="darkorange",pch=16,cex=2,xlab=NULL,ylab=NULL,xaxt='n',yaxt='n',ann=FALSE)
par(new=T)
plot(Sr08m/500~as.numeric(names(TGPm)),ylim=c(0,1.3),col="dodgerblue3",pch=16,cex=2,xlab=NULL,ylab=NULL,xaxt='n',yaxt='n',ann=FALSE)
par(new=T)
plot(dif08~as.numeric(names(TGPm)),ylim=c(0,1.3),col="black",pch=16,cex=2,xlab=NULL,ylab=NULL,xaxt='n',yaxt='n',ann=FALSE)
add.error.bars(X=c(1990:2015),Y=RABS08m/500,SEX=0,SEY=sdR08/1000,w=0,col="darkred")
add.error.bars(X=c(1990:2015),Y=H08m/100,SEX=0,SEY=sdH08/200,w=0,col="darkorange")
add.error.bars(X=c(1990:2015),Y=Sr08m/500,SEX=0,SEY=sdS08/1000,w=0,col="dodgerblue3")
add.error.bars(X=c(1990:2015),Y=c((TBP08m-TGPm)),SEX=0,SEY=sd08/2,w=0,col="black")

par(mfrow=c(2,2))

termplot(lm_m,partial.resid = T, se=T,col.se = "darkred",col.term = "darkred",pch=16,cex=1.5,lwd.se=2,lwd.term = 2,cex.axis=2,cex.lab=2,xlabs=c("Rabs (W m-2)","H (W m-2)","LWbud (W m-2)"),ylabs=rep("Part. res. (°C)"))
dev.off()


df<-data.frame(x=c(1990:2015),Rabs=residuals(lm_m,type = "partial")[,1],H=residuals(lm_m,type = "partial")[,2],Sr=residuals(lm_m,type = "partial")[,3],TG=residuals(lm_m,type = "partial")[,4])
summary(lm(df$Rabs~df$x))
summary(lm(df$H~df$x))
summary(lm(df$Sr~df$x))
summary(lm(df$TG~df$x))
# no relationships between partial residuals and time. All the temporal variability is captured by TG, Rabs, Sw and H. 

p1<-ggplot(df, aes(x=x, y=Rabs)) +
  geom_point(size=4) +
  geom_smooth(method = "lm")+theme_classic(base_size = 20) +
  xlab("Year")+
  ylab("Part. res. (°C")
p2<-ggplot(df, aes(x=x, y=H)) +
  geom_point(size=4) + geom_smooth(method = "lm")+
  theme_classic(base_size = 20) +
  xlab("Year")+
  ylab("Part. res. (°C")
p3<-ggplot(df, aes(x=x, y=Sr)) +
  geom_point(size=4) +
  geom_smooth(method = "lm")+theme_classic(base_size = 20) +
  xlab("Year")+
  ylab("Part. res. (°C")
p4<-ggplot(df, aes(x=x, y=TG)) +
  geom_point(size=4) +
  geom_smooth(method = "lm")+theme_classic(base_size = 20) +
  xlab("Year")+
  ylab("Part. res. (°C")
pdf("Supp_Fig5.pdf",width=10,height=10)
print(grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2))
dev.off()

