library(RNetCDF)
library(xts)
library(ggplot2)
library(dplyr)
library(tidyr)
source("Ebalance.R")
source("Tbud.R")

forcing<-"/home/orchideeshare/igcmg/IGCM/SRF/METEO/FLUXNET/FR-Hes_1997-2006.nc"
nc<-open.nc(forcing)
tair<-var.get.nc(nc,"Tair") #K
qair<-var.get.nc(nc,"Qair") #kg/kg
patm<-var.get.nc(nc,"PSurf")/1000 #Pa
u<-var.get.nc(nc,"Wind_N") #m/s
v<-var.get.nc(nc,"Wind_E") #m/s
sw<-var.get.nc(nc,"SWdown") #W/m2
lw<-var.get.nc(nc,"LWdown") #W/m2
gpp<-var.get.nc(nc,"GPP") #W/m2

wind= sqrt (u*u + v*v)
dates<-var.get.nc(nc,"tstep")
dates<-as.POSIXct(dates,origin="1996-12-31 23:00:00")
  # computes bud temperature
  fTbud<-function(x){
         inputs<-list(tair=tair[x],
                      patm=patm[x],
                      qair=qair[x],
                      Sw=sw[x],
                      Lw=lw[x],
                      wind=wind[x],
		      abs_sw=0.5)

         tmp<-Tbud(inputs)
        return(tmp)
  }
#######################################################
tg_bud<-sapply(X=c(14400:23500),FUN=fTbud,simplify=T)
gpp1<-xts(gpp[14400:23500],dates[14400:23500])
tair1<-xts(tair[14400:23500],dates[14400:23500])
tg_bud1<-xts(tg_bud,dates[14400:23500])
dT<-(tg_bud1-tair1)

tair_mean<-apply.daily(tair1["1997-11-01/1998-05-04"]-273.15,mean)
tair_max<-apply.daily(tair1["1997-11-01/1998-05-04"]-273.15,max)
tair_min<-apply.daily(tair1["1997-11-01/1998-05-04"]-273.15,min)

tbud_mean<-apply.daily(tg_bud1["1997-11-01/1998-05-04"]-273.15,mean)
tbud_max<-apply.daily(tg_bud1["1997-11-01/1998-05-04"]-273.15,max)
tbud_min<-apply.daily(tg_bud1["1997-11-01/1998-05-04"]-273.15,min)


NCD_air<-GDD_air<-apply.daily(tair1["1997-11-01/1998-05-04"]-273.15,mean)
NCD_bud<-GDD_bud<-apply.daily(tg_bud1["1997-11-01/1998-05-04"]-273.15,mean)

# forcing
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


pdf("dTbud_spring.pdf",width=10,height=10)
plot(na.approx(cbind(rollmean(c(tbud_mean-tair_mean),10),rollmean(c(tbud_max-tair_max),10),rollmean(c(tbud_min-tair_min),10),dT)),col=c("black","darkred","darkblue","grey65"),lwd=c(2.5,2.5,2.5,0.5),main="")
events <- xts("budburst", 
              as.Date(c("1998-04-25")))
addEventLines(events,col="red",lwd=1.5)
dev.off()
######################################################
tg_bud<-sapply(X=c(21000:22000),FUN=fTbud,simplify=T)
tg_bud1<-xts(tg_bud,dates[21000:22000])
tair1<-xts(tair[21000:22000],dates[21000:22000])


dT05<-(tg_bud1-tair1)["1998-03-24/1998-03-24"]
dT08<-(tg_bud1-tair1)["1998-03-24/1998-03-24"]


dT<-data.frame(Time=seq(0,23.5,0.5),
               dT05=as.numeric(dT05),
               dT08=as.numeric(dT08),
               Sw=as.numeric(xts(sw[21000:22000],dates[21000:22000])["1998-03-24/1998-03-24"])/100,
               Lw=as.numeric(xts(lw[21000:22000],dates[21000:22000])["1998-03-24/1998-03-24"])/100)


#df<-dT %>%
#    dplyr::select(Time, dT, Sw, Lw) %>%
#    gather(key = "variable", value = "value", -Time)
#
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

pdf("daily_course_dTbud.pdf",width=9,height=7)
print(p)
dev.off()




