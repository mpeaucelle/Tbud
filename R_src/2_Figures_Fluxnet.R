############################################################
# 
# Author: Marc Peaucelle, 2020/06
# Description:
# Script to generate Figure 1 & 2 of "Plants phenological sensitivity to warming lies in their energy budget" (Peaucelle, Penuelas & Verbeeck, 2020)
# FLUXNET data are not provided within this script but can be downloaded after registration at https://fluxnet.fluxdata.org/data/fluxnet2015-dataset/
# 
############################################################
library(RNetCDF)
library(xts)
library(ggplot2)
source("Ebalance.R")
source("Tbud.R")

# Load the netcdf forcinf file. A template of the file structure is provided in the template folder. 
forcing<-"FR-Hes_1997-2006.nc"
# forcing<-"DK-Sor_1996-2006.nc"
# forcing<-"IT-Col_1996-2006.nc"
# forcing<-"US-WCr_1999-2006.nc"
# forcing<-"US-Ha1_1991-2006.nc"
# forcing<-"FI-Hyy_1996-2006.nc"

nc<-open.nc(forcing)
tair<-var.get.nc(nc,"Tair") #K
qair<-var.get.nc(nc,"Qair") #kg/kg
patm<-var.get.nc(nc,"PSurf")/1000 #Pa --> kPa
u<-var.get.nc(nc,"Wind_N") #m/s
v<-var.get.nc(nc,"Wind_E") #m/s
sw<-var.get.nc(nc,"SWdown") #W/m2
lw<-var.get.nc(nc,"LWdown") #W/m2
gpp<-var.get.nc(nc,"GPP") #W/m2

# computes wind momentum
wind= sqrt (u*u + v*v)
dates<-var.get.nc(nc,"tstep")
dates<-as.POSIXct(dates,origin="1996-12-31 23:00:00")
# change origin for other sites
# dates<-as.POSIXct(dates,origin="1995-12-31 23:00:00")


# function to compute bud temperature
fTbud<-function(x,abs_sw,bud_d,r){
  inputs<-list(tair=tair[x],
               patm=patm[x],
               qair=qair[x],
               Sw=sw[x],
               Lw=lw[x],
               wind=wind[x],
               abs_sw=abs_sw,
               bud_d=bud_d,
               r=r)
  
  tmp<-Tbud(inputs)
  return(tmp)
}

######################################################
# Figure 2a with the illustration of the daily course of bud T for 1 day
# We selected a random day, corresponding to time steps 21000:22000

tair1<-xts(tair[21000:22000],dates[21000:22000])

# solar absorptivity to SW = 0.5
tg_bud05<-sapply(X=c(21000:22000),FUN=fTbud,simplify=T,abs_sw=0.5,bud_d=0.005,r=0.2)
tg_bud05<-xts(tg_bud05,dates[21000:22000])

# solar absorptivity to SW = 0.8
tg_bud08<-sapply(X=c(21000:22000),FUN=fTbud,simplify=T,abs_sw=0.8,bud_d=0.005,r=0.2)
tg_bud08<-xts(tg_bud08,dates[21000:22000])

dT05<-(tg_bud05-tair1)["1998-03-24/1998-03-24"]
dT08<-(tg_bud08-tair1)["1998-03-24/1998-03-24"]


dT<-data.frame(Time=seq(0,23.5,0.5),
               dT05=as.numeric(dT05),
               dT08=as.numeric(dT08),
               Sw=as.numeric(xts(sw[21000:22000],dates[21000:22000])["1998-03-24/1998-03-24"])/100,
               Lw=as.numeric(xts(lw[21000:22000],dates[21000:22000])["1998-03-24/1998-03-24"])/100)


colors <- c("dT05" = "grey50", "dT08"="black", "Sw" = "orange", "Lw" = "darkorange")

p<-ggplot(dT, aes(x = Time)) + theme_classic() +
  geom_line(aes(y = Sw, color="Sw"), linetype="twodash",size=1.1)+
  geom_line(aes(y = Lw, color="Lw"), linetype="dashed",size=1.1) +
  geom_line(aes(y = dT05, color = "dT05"),size=1.5) +
  geom_line(aes(y = dT08, color = "dT08"),size=1.5) +
  labs(x="Time",
       color = "Legend") +
   scale_y_continuous(

    # Features of the first axis
    name = "dT (°C)",

    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*100, name=expression("Sw, Lw (W m"^-2*")"))
  ) +
  scale_color_manual(values = colors) +
  theme(axis.title.x = element_text( size=14, face="bold"),
        axis.title.y = element_text( size=14, face="bold"),
        axis.text.x = element_text( size=14),
        axis.text.y = element_text(size=14),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14))

p<- p + theme(panel.border=element_blank(), axis.line=element_blank(), panel.grid.major.y=element_line(colour="grey80"), axis.line.x.bottom = element_line(colour="black",size=0.5),axis.line.y.left = element_line(colour="black",size=0.5))

pdf("Figure_2a.pdf",width=9,height=7)
print(p)
dev.off()



#######################################################
# We compute bud temperature by calling fTbud
# Figure 2b is generated with abs_sw=0.8
# The code is not automated to run on multiple FLUXNET sites, replace forcing path at the begining of the script to generate figures for others sites. 

# We selected the winter and spring preseaon of 1998, from November to May/June
start_date<-"1997-11-01"
end_date<-"1998-06-01"
period<-paste0(start_date,"/",end_date)
preseason<-which((dates>start_date) & (dates<end_date))
tg_bud<-sapply(X=preseason,FUN=fTbud,simplify=T,abs_sw=0.8,bud_d=0.005,r=0.2)
tair1<-xts(tair[preseason],dates[preseason])
tg_bud1<-xts(tg_bud,dates[preseason])

dT<-(tg_bud1-tair1)

#####
gpp1<-xts(gpp[preseason],dates[preseason])
plot(gpp1)
# We see from GPP that budburst is around the 25 of April for the year 1998 for Hesse
budburst<-"1998-04-25"
# budburst<-"1998-05-01" # DK-Sor
# budburst<-"1998-05-01" # IT-Col
# budburst<-"1998-05-01" # US-Ha1
# budburst<-"1998-04-15" # FI-Hyy


# computes daily min, mean and max
tair_mean<-apply.daily(tair1[period]-273.15,mean)
tair_max<-apply.daily(tair1[period]-273.15,max)
tair_min<-apply.daily(tair1[period]-273.15,min)

tbud_mean<-apply.daily(tg_bud1[period]-273.15,mean)
tbud_max<-apply.daily(tg_bud1[period]-273.15,max)
tbud_min<-apply.daily(tg_bud1[period]-273.15,min)


# computes forcing and chilling at budburst
# Note that we use here the simpliest formulation of GDD and NCD
# Better calibration would be requiered to really study budburst

period<-paste0(as.Date(budburst)-60,"/",budburst)
NCD_air<-GDD_air<-apply.daily(tair1[period]-273.15,mean)
NCD_bud<-GDD_bud<-apply.daily(tg_bud1[period]-273.15,mean)

GDD_air[GDD_air<5]<-0
GDD_bud[GDD_bud<5]<-0

GDD_air<-cumsum(GDD_air)
GDD_bud<-cumsum(GDD_bud)

NCD_air[NCD_air>5]<-0
NCD_air[NCD_air<0]<-0
NCD_air[NCD_air!=0]<-1

NCD_bud[NCD_bud>5]<-0
NCD_bud[NCD_bud<0]<-0
NCD_bud[NCD_bud!=0]<-1

NCD_air<-cumsum(NCD_air)
NCD_bud<-cumsum(NCD_bud)

GDD_air[budburst]
GDD_bud[budburst]
NCD_air[budburst]
NCD_bud[budburst]

pdf("Figure2b.pdf",width=10,height=10)
# Change name for other FLUXNET sites
# pdf("Figure_FI-Hyy.pdf",width=10,height=10)
# plot rolling mean 
events <- xts("budburst",
              as.Date(budburst))

plot(na.approx(cbind(rollmean(c(tbud_mean-tair_mean),10),rollmean(c(tbud_max-tair_max),10),rollmean(c(tbud_min-tair_min),10),dT)),col=c("black","darkred","darkblue","grey65"),lwd=c(2.5,2.5,2.5,0.5),main="")
addEventLines(events,col="red",lwd=1.5,on=0)
dev.off()

######################################################
# Role of bud size and ground albedo on simulated bud temperature
# We assessed the role of bud size and soil albedo on simulated Temperature at FLUXNET site
# Because bud size influcences wind and turbulences, we applied the model over the preseason and not only 1 day as in Figure 2a. 

########## Bud size ############
# We selected the winter and spring preseaon of 1998, from November to May
start_date<-"1997-11-01"
end_date<-"1998-05-01"
period<-paste0(start_date,"/",end_date)
preseason<-which((dates>start_date) & (dates<end_date))
tair1<-xts(tair[preseason],dates[preseason])
# bud diameter of 5mm, soil albedo and abs_sw are constant at 0.2 and 0.8 respectively
tg_bud05<-sapply(X=preseason,FUN=fTbud,simplify=T,abs_sw=0.8,bud_d=0.005,r=0.2)
tg_bud05<-xts(tg_bud05,dates[preseason])

# bud diameter of 7mm
tg_bud07<-sapply(X=preseason,FUN=fTbud,simplify=T,abs_sw=0.8,bud_d=0.007,r=0.2)
tg_bud07<-xts(tg_bud07,dates[preseason])

# bud diameter of 10mm
tg_bud10<-sapply(X=preseason,FUN=fTbud,simplify=T,abs_sw=0.8,bud_d=0.01,r=0.2)
tg_bud10<-xts(tg_bud10,dates[preseason])

# bud diameter of 13mm
tg_bud13<-sapply(X=preseason,FUN=fTbud,simplify=T,abs_sw=0.8,bud_d=0.013,r=0.2)
tg_bud13<-xts(tg_bud13,dates[preseason])

dT02<-(tg_bud07-tg_bud05)
dT05<-(tg_bud10-tg_bud05)
dT08<-(tg_bud13-tg_bud05)

# We see that budburst is around the 25 of April for the year 1998
budburst<-"1998-04-25"
# computes daily min, mean and max
tair_mean<-apply.daily(tair1[period]-273.15,mean)
tair_max<-apply.daily(tair1[period]-273.15,max)
tair_min<-apply.daily(tair1[period]-273.15,min)

tbud05_mean<-apply.daily(tg_bud05[period]-273.15,mean)
tbud05_max<-apply.daily(tg_bud05[period]-273.15,max)
tbud05_min<-apply.daily(tg_bud05[period]-273.15,min)

tbud07_mean<-apply.daily(tg_bud07[period]-273.15,mean)
tbud07_max<-apply.daily(tg_bud07[period]-273.15,max)
tbud07_min<-apply.daily(tg_bud07[period]-273.15,min)

tbud10_mean<-apply.daily(tg_bud10[period]-273.15,mean)
tbud10_max<-apply.daily(tg_bud10[period]-273.15,max)
tbud10_min<-apply.daily(tg_bud10[period]-273.15,min)

tbud13_mean<-apply.daily(tg_bud13[period]-273.15,mean)
tbud13_max<-apply.daily(tg_bud13[period]-273.15,max)
tbud13_min<-apply.daily(tg_bud13[period]-273.15,min)


#color palette
colfunc_ave <- colorRampPalette(c("grey50", "black"))
colfunc_min <- colorRampPalette(c("lightblue", "darkblue"))
colfunc_max <- colorRampPalette(c("darkorange", "darkred"))
pdf("Supplementary_Figure_6.pdf",width=10,height=10)
# plot rolling mean for Tave
plot(na.approx(cbind(rollmean(c(tbud07_mean-tbud05_mean),10),rollmean(c(tbud10_mean-tbud05_mean),10),rollmean(c(tbud13_mean-tbud05_mean),10),
rollmean(c(tbud07_min-tbud05_min),10),rollmean(c(tbud10_min-tbud05_min),10),rollmean(c(tbud13_min-tbud05_min),10),
rollmean(c(tbud07_max-tbud05_max),10),rollmean(c(tbud10_max-tbud05_max),10),rollmean(c(tbud13_max-tbud05_max),10))),col=c(colfunc_ave(3),colfunc_min(3),colfunc_max(3)),lwd=rep(2.5,9),main="")

dev.off()


########## Albedo ############
# We selected the winter and spring preseaon of 1998, from November to May
start_date<-"1997-11-01"
end_date<-"1998-05-01"
period<-paste0(start_date,"/",end_date)
preseason<-which((dates>start_date) & (dates<end_date))
tair1<-xts(tair[preseason],dates[preseason])
# bud diameter and abs_sw are constant at 0.05 and 0.8 respectively
# albedo of 0.1 (~ wet bare ground)
tg_bud01<-sapply(X=preseason,FUN=fTbud,simplify=T,abs_sw=0.8,bud_d=0.005,r=0.1)
tg_bud01<-xts(tg_bud01,dates[preseason])

# albedo of 0.5
tg_bud05<-sapply(X=preseason,FUN=fTbud,simplify=T,abs_sw=0.8,bud_d=0.005,r=0.5)
tg_bud05<-xts(tg_bud05,dates[preseason])

# albedo of 0.9 (~ snow)
tg_bud09<-sapply(X=preseason,FUN=fTbud,simplify=T,abs_sw=0.8,bud_d=0.005,r=0.9)
tg_bud09<-xts(tg_bud09,dates[preseason])

# We see that budburst is around the 25 of April for the year 1998
budburst<-"1998-04-25"
# computes daily min, mean and max
tair_mean<-apply.daily(tair1[period]-273.15,mean)
tair_max<-apply.daily(tair1[period]-273.15,max)
tair_min<-apply.daily(tair1[period]-273.15,min)

tbud01_mean<-apply.daily(tg_bud01[period]-273.15,mean)
tbud01_max<-apply.daily(tg_bud01[period]-273.15,max)
tbud01_min<-apply.daily(tg_bud01[period]-273.15,min)

tbud05_mean<-apply.daily(tg_bud05[period]-273.15,mean)
tbud05_max<-apply.daily(tg_bud05[period]-273.15,max)
tbud05_min<-apply.daily(tg_bud05[period]-273.15,min)

tbud09_mean<-apply.daily(tg_bud09[period]-273.15,mean)
tbud09_max<-apply.daily(tg_bud09[period]-273.15,max)
tbud09_min<-apply.daily(tg_bud09[period]-273.15,min)

#color palette
colfunc_ave <- colorRampPalette(c("grey50", "black"))
colfunc_min <- colorRampPalette(c("lightblue", "darkblue"))
colfunc_max <- colorRampPalette(c("darkorange", "darkred"))
pdf("Supplementary_Figure_5.pdf",width=10,height=10)
# plot rolling mean for Tave
plot(na.approx(cbind(rollmean(c(tbud05_mean-tbud01_mean),10),rollmean(c(tbud09_mean-tbud01_mean),10),
rollmean(c(tbud05_min-tbud01_min),10),rollmean(c(tbud09_min-tbud01_min),10),
rollmean(c(tbud05_max-tbud01_max),10),rollmean(c(tbud09_max-tbud01_max),10))),col=c(colfunc_ave(2),colfunc_min(2),colfunc_max(2)),lwd=rep(2.5,6),main="")

dev.off()




