
library(ncdf4)
gridmet<-nc_open('http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_tmmx_1979_CurrentYear_CONUS.nc')
names(gridmet$dim)
lat_poi <- 42
lon_poi <- -115 
lat_index <- which.min(abs(lat_poi-gridmet$dim$lat$vals))
lon_index <- which.min(abs(lon_poi-gridmet$dim$lon$vals))
tmax_1pixel<-ncvar_get(nc=gridmet,varid="daily_maximum_temperature",start=c(lat_poi,lon_poi,1),count=c(1,1,-1))

lat_nw <- 42
lon_nw <- -124.5 
lat_se <- 32
lon_se <- -115 
lat_index1 <- which.min(abs(lat_nw-gridmet$dim$lat$vals))
lon_index1 <- which.min(abs(lon_nw-gridmet$dim$lon$vals))
lat_index2 <- which.min(abs(lat_se-gridmet$dim$lat$vals))
lon_index2 <- which.min(abs(lon_se-gridmet$dim$lon$vals))

tmax_grid<-ncvar_get(nc=gridmet,varid="daily_maximum_temperature",start=c(lat_index1,lon_index1,1),count=c(lat_index2+1-lat_index1,lon_index2-lon_index1+1,-1))