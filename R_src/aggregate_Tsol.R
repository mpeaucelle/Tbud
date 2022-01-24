library(raster)
library(RNetCDF)

# load Tsol data
Tsol01<-brick("/home/orchidee04/mpeau/Tsol_ERA5/Tsol_Era5_6h.nc")

# load CRUJra mask
nc_tg<-open.nc("/home/orchidee02/mpeau/pheno_PEP/netcdf_forcing/crujra_1D/crujra_tair_1969_2018.nc")
themask<-var.get.nc(nc_tg,"mask")
lat<-var.get.nc(nc_tg,"latitude")
lon<-var.get.nc(nc_tg,"longitude")
rast_mask<-raster(t(themask),xmn=min(lon),xmx=max(lon),ymn=min(lat),ymx=max(lat))

Tsol025<-resample(Tsol01,rast_mask)

writeRaster(Tsol025,"/home/orchidee04/mpeau/Tsol_ERA5/Tsol025.grd")

q(save="no")


