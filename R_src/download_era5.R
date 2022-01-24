Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
devtools::install_github("https://github.com/ErikKusch/KrigR")
library(KrigR)


Extent <- extent(5,20,42,56) # roughly the extent of Saxony
API_User <- "put your API user ID here"
API_Key <- "put your API key here"


datas<-"era5-land"
datas<-"derived-near-surface-meteorological-variables"



vari<-"Near-surface_air_temperature "
vari<-"Near-surface_specific_humidity"
vari<-"Surface_downwelling_longwave_radiation"
vari<-"Near-surface_wind_speed" 
vari<-"Surface_air_pressure" 
vari<-"Surface_downwelling_shortwave_radiation"
vari<-"Rainfall_flux"

State_Raw <- download_ERA(
  Variable = vari,
  DataSet = datas,
  DateStart = "1990-01-01",
  DateStop = "2015-12-30",
  TResolution = "hour",
  TStep = 6,
  Extent = Extent,
  API_User = API_User,
  API_Key = API_Key
)
