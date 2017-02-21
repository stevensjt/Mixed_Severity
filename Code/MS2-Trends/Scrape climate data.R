#Scrape climate data?
 
library(tidyverse)
library(rgeos)# For gCentroid
library(rnoaa)
library(sp) #For proj4string
library(stringr) #for str_sub
library(raster)
library(ncdf4)
source(file="./Code/MS2-Trends/Functions.R")

####1. Data management####
#Read in  data
fire.list=read.csv("./Data/Derived/all_fires_ForAnalysis.csv")
fires_usfs=read.csv("./Data/Raw/fires_usfs.csv")
hs_patches=readRDS("./Data/hs_patches.RDS") #CRS EPSG:3310, NAD83 CA Albers
station_data <- ghcnd_stations() #Run once to get list of all the stations (takes a minute)

#Additional data processing
fire.list=merge.data.frame(fires_usfs,fire.list[,c("VB_ID","q")],by="VB_ID",all=T) #This is temporary, add q from fires I already calculated it for. Eventually recalculate everything from the new list.
fire.list=fire.list[-unique(c(grep("0000",fire.list$IGNITION_DATE), 
                              grep("0000",fire.list$CONT_DATE))),] #Remove fires where either the start or end date is unknown (can't get weather) (n=11).
#If the end date is unknown but the start date is known, make it 7 days after the start date.
fire.list[fire.list$CONT_DATE==0,"CONT_DATE"] <- 
  fire.list[fire.list$CONT_DATE==0,"IGNITION_DATE"] + 7
#If the end date is past the 30th of the mnth, add 70 to get to the next month.
fire.list[which(as.numeric(str_sub(fire.list$CONT_DATE,-2))>30),"CONT_DATE"]=
  fire.list[which(as.numeric(str_sub(fire.list$CONT_DATE,-2))>30),"CONT_DATE"] + 70 #
#If the end month is known but the end day is unknown, make the end day the first day of the month:
fire.list[grep("0800",fire.list$CONT_DATE), "CONT_DATE"] = fire.list[grep("0800",fire.list$CONT_DATE), "CONT_DATE"]+1
fire.list[grep("0900",fire.list$CONT_DATE), "CONT_DATE"] = fire.list[grep("0900",fire.list$CONT_DATE), "CONT_DATE"]+1

#fire.list= read.csv("./Data/Derived/fire_list.csv")

####2. Loop through all fires using NOAA data####
#Bad cases involving projection/spTransform/shapefile doesn't exist: 348, 355
for(f in 437:nrow(fire.list)){
#Get centroids of polygons
fires.to.sample=as.character(fire.list$VB_ID)
hs_fire=hs_patches[hs_patches$VB_ID==fires.to.sample[f],]

if(nrow(hs_fire@data)>0){ #CHECK Shapefile: If there is a shapefile in the database, proceed
hs_fire_ll=spTransform(hs_fire,CRS("+proj=longlat +datum=WGS84")) #Reproject to LL to use as NOAA Arg
#proj4string(hs_fire)

#plot(hs_fire_ll,axes=T)
#points(gCentroid(hs_fire_ll))

#Possibilities for accessing noaa data:
#ftp://cran.r-project.org/pub/R/web/packages/rnoaa/rnoaa.pdf
#https://ropensci.org/tutorials/rnoaa_tutorial.html
#https://www.ncdc.noaa.gov/cdo-web/token
#Station acronyms: https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt

#ncdc_datasets(token ="FTEIRBcwLUCArmZqFHwVomZhjQDOCqMv") #GHCND are daily summaries
#ncdc_locs(datasetid='GHCND',token ="FTEIRBcwLUCArmZqFHwVomZhjQDOCqMv")

#Get all station X variable combinations within 100 km of the fire centroid:
stations <- meteo_distance(station_data, lat=gCentroid(hs_fire_ll)@coords[2], 
                           lon=gCentroid(hs_fire_ll)@coords[1], radius = 100)
stations <-stations[-which(is.na(stations$id)),] #Previous call returns a line of NA's; delete that line.
#parms=c("TMAX","TMIN","WSF5","WSFM","EVAP") #Max Temp, Min Temp, Max 5-minute windspeed, Max 1-mile windspeed, Min RH apparently doesn't exist. EVAP= evaporation of water

#Find the nearest station with temperature data and get that data. 
T_Station=stations[stations$element=="TMAX" & 
                     stations$first_year<fire.list[f,"FIRE_YEAR"] & 
                     stations$last_year>=fire.list[f,"FIRE_YEAR"],] 
T_Station=T_Station[which.min(T_Station$distance),]

#Get dates
sd=fire.list[f,"IGNITION_DATE"]
sd=gsub('^(.{4})(.*)$', '\\1-\\2', sd) #Insert dash (crazy regex)
sd=gsub('^(.{7})(.*)$', '\\1-\\2', sd) #Insert dash
ed=paste(fire.list[f,"CONT_DATE"])
ed=gsub('^(.{4})(.*)$', '\\1-\\2', ed) #Insert dash (crazy regex)
ed=gsub('^(.{7})(.*)$', '\\1-\\2', ed) #Insert dash

out <- ncdc(datasetid='GHCND', stationid=paste('GHCND',T_Station$id,sep=':'), 
            datatypeid=c('TMAX','TMIN'), startdate = sd, enddate = ed, 
            limit = 300,
            token ="FTEIRBcwLUCArmZqFHwVomZhjQDOCqMv")


if(nrow(out$data)>0){#CHECK temp data: If temp data exists, proceed
fire.list[f,"TMAX"]=max(out$data[out$data$datatype=="TMAX","value"])/10 #Convert to degrees C
fire.list[f,"MAXTMIN"]=max(out$data[out$data$datatype=="TMIN","value"])/10 #Convert to degrees C
fire.list[f,"T_dist"]=round(T_Station$distance,1) #In KM
}#End CHECK temp data


#Find the nearest station with wind data and get that data:
if(any(stations$element=="WSF5" & 
       stations$first_year<fire.list[f,"FIRE_YEAR"] & 
       stations$last_year>fire.list[f,"FIRE_YEAR"])){ #CHECK wind data: If wind data exists, proceed
  
W_Station=stations[stations$element=="WSF5" & 
                     stations$first_year<fire.list[f,"FIRE_YEAR"] & 
                     stations$last_year>fire.list[f,"FIRE_YEAR"],]
W_Station=W_Station[which.min(W_Station$distance),]

#Get dates
sd=fire.list[f,"IGNITION_DATE"]
sd=gsub('^(.{4})(.*)$', '\\1-\\2', sd) #Insert dash (crazy regex)
sd=gsub('^(.{7})(.*)$', '\\1-\\2', sd) #Insert dash
ed=paste(fire.list[f,"CONT_DATE"])
ed=gsub('^(.{4})(.*)$', '\\1-\\2', ed) #Insert dash (crazy regex)
ed=gsub('^(.{7})(.*)$', '\\1-\\2', ed) #Insert dash

#out <- ncdc(datasetid='GHCND', stationid=paste('GHCND',W_Station$id,sep=':'), 
#            datatypeid=c('WSF5'), startdate = sd, enddate = ed, 
#            limit = 300,
#            token ="FTEIRBcwLUCArmZqFHwVomZhjQDOCqMv")
#fire.list[f,"MAX_WSF5"]=max(out$data[out$data$datatype=="WSF5","value"])/10 #Convert to m/s
#fire.list[f,"W_dist"]=round(W_Station$distance,1) #In KM
}#End CHECK wind data
}#End CHECK Shapefile
gc()
print(f)
}#End for loop

#write_csv(fire.list,'./Data/Derived/fire_list2.csv')

####3. Look into Abatzoglou gridded data####
#fire.list= read.csv("./Data/Derived/fire_list.csv")
#Sample data request page: https://www.reacchpna.org/thredds/ncss/grid/MET/tmmx/tmmx_2016.nc/dataset.html
#Sample url: "https://www.reacchpna.org/thredds/ncss/MET/tmmx/tmmx_1984.nc?var=air_temperature&north=42.15&west=-124.43&east=-114.50&south=32.66&disableProjSubset=on&horizStride=1&time_start=1984-07-18T00%3A00%3A00Z&time_end=1984-07-25T00%3A00%3A00Z&timeStride=1&accept=netcdf"
#To determine full variable, go to data source, select a year of variable of interest, and click "NetcdfSubset" to study example link like what is above.
#Data source: https://www.reacchpna.org/thredds/reacch_climate_MET_catalog.html

#Bad fires: 248 (Gondola; but works if you run independently), 
for(f in 1:507){ #Fire for loop
  #Get centroids of polygons
  fires.to.sample=as.character(fire.list$VB_ID)
  hs_fire=hs_patches[hs_patches$VB_ID==fires.to.sample[f],]
  
  if(nrow(hs_fire@data)>0){ #CHECK Shapefile: If there is a shapefile in the database, proceed.
    hs_fire_ll=spTransform(hs_fire,CRS("+proj=longlat +datum=WGS84")) #Reproject to LL to use as NOAA Arg
    #Get dates
    sd=fire.list[f,"IGNITION_DATE"]
    sd=gsub('^(.{4})(.*)$', '\\1-\\2', sd) #Insert dash (crazy regex)
    sd=gsub('^(.{7})(.*)$', '\\1-\\2', sd) #Insert dash
    ed=paste(fire.list[f,"CONT_DATE"])
    ed=gsub('^(.{4})(.*)$', '\\1-\\2', ed) #Insert dash (crazy regex)
    ed=gsub('^(.{7})(.*)$', '\\1-\\2', ed) #Insert dash
    
    vs=c("tmmx","tmmn","rmax","bi") #Download on "erc" not working properly, so skipping for now.
    vs_full=c("air_temperature","air_temperature","relative_humidity","burning_index_g") #skipping "energy_release_component-g"
    for(index in 1:length(vs)){ #Variable for loop
      v=vs[index]
      v_full=vs_full[index]
      v_link=paste0("https://www.reacchpna.org/thredds/ncss/MET/",v,
                    "/",v,"_",fire.list[f,"FIRE_YEAR"],".nc?var=",v_full,
                    "&north=",gCentroid(hs_fire_ll)$y+0.2,
                    "&west=",gCentroid(hs_fire_ll)$x-0.2,
                    "&east=",gCentroid(hs_fire_ll)$x+0.2,
                    "&south=",gCentroid(hs_fire_ll)$y-0.2,"&disableProjSubset=on&horizStride=1",
                    "&time_start=",sd,"T00%3A00%3A00Z",
                    "&time_end=",ed,"T00%3A00%3A00Z&timeStride=1&accept=netcdf")
      dest <-  paste0("./Data/climate_nc/",v,".nc" )
      #fire.list[f,"lat"]=gCentroid(hs_fire_ll)$y
      #fire.list[f,"long"]=gCentroid(hs_fire_ll)$x
      tmp=try(download.file(url=v_link,destfile=dest), silent=T)
      if(class(tmp)!="try-error"){ #CHECK download error: 
        #if the download produced an error, will have to do it manually so skip the next bit.
      ncin <- nc_open(dest)
      lat <- ncvar_get(ncin,"lat",verbose=F)
      lat_target <- which.min(abs(lat - gCentroid(hs_fire_ll)$y)) #Find closest latitude to fire centroid
      lon <- ncvar_get(ncin,"lon",verbose=F)
      lon_target <- which.min(abs(lon - gCentroid(hs_fire_ll)$x)) #Find closest latitude to fire centroid
      v_array <- ncvar_get(ncin,v_full)
      if(class(v_array)=="matrix"){ #If there was only one burn day so there's a matrix instead of an array
        #Duplicate the arrayso there's no error produced in maximum calculation
        v_array=replicate(2,v_array,simplify="array")
      }
      if(v=="tmmx"){
        #Get maximum high temperature during the burn window
        fire.list[f,paste0("max_",v)]=max(v_array[lat_target,lon_target,])-273.15 #Convert from K to C
      }
      if(v=="tmmn"){
        #Get maximum low temperature during the burn window
        fire.list[f,paste0("max_",v)]=max(v_array[lat_target,lon_target,])-273.15 #Convert from K to C
      }
      #Get minimum high RH during the burn window
      if(v=="rmax"){
        fire.list[f,paste0("min_",v)]=min(v_array[lat_target,lon_target,]) #Not sure of units
      }
      if(v=="bi"){
        #Get max burn index during burn window
        fire.list[f,paste0("max_",v)]=max(v_array[lat_target,lon_target,]) #Not sure of units
      }
      if(v=="erc"){
        #Get max erc index during burn window
        fire.list[f,paste0("max_",v)]=max(v_array[lat_target,lon_target,]) #Not sure of units
      }
      } #END CHECK download error
    } #END variable for loop
    gc()
    
  } #END CHECK shapefile
  print(f)
  gc()
} #END fire for loop

fire.list[which(fire.list$max_bi<0),"max_bi"]=NA #Deer fire is wierd; very small and has negative BI

#write_csv(fire.list,'./Data/Derived/fire_list.csv')

