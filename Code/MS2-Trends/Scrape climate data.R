#Scrape climate data?
 
library(tidyverse)
library(rgeos)# For gCentroid
library(rnoaa)
library(sp) #For proj4string
library(stringr) #for str_sub
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
fire.list[fire.list$CONT_DATE==0,"CONT_DATE"] <- 
  fire.list[fire.list$CONT_DATE==0,"IGNITION_DATE"] + 7
#If the end date is past the 30th of the mnth, add 70 to get to the next month.
fire.list[which(as.numeric(str_sub(fire.list$CONT_DATE,-2))>30),"CONT_DATE"]=
  fire.list[which(as.numeric(str_sub(fire.list$CONT_DATE,-2))>30),"CONT_DATE"] + 70 #

####2. Loop through all fires:
#Bad cases involving projection/spTransform/shapefile doesn't exist: 348, 355, 476, 477, 478, 479
#Checkme: Breaks down starting around f=476. Not sure why.
for(f in 479:nrow(fire.list)){
#Get centroids of polygons
fires.to.sample=as.character(fire.list$VB_ID)
hs_fire=hs_patches[hs_patches$VB_ID==fires.to.sample[f],]
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

#Find the nearest station with temperature data and get that data. Sometimes towards the end of their lifespan the data doesn't exist even though the station does, so we add the >5 yr caveat.
T_Station=stations[stations$element=="TMAX" & 
                     stations$first_year<fire.list[f,"FIRE_YEAR"] & 
                     stations$last_year>fire.list[f,"FIRE_YEAR"]+5,] 
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

if(!is.na(out$data)){
if(nrow(out$data)>0){
fire.list[f,"TMAX"]=max(out$data[out$data$datatype=="TMAX","value"])/10 #Convert to degrees C
fire.list[f,"MAXTMIN"]=max(out$data[out$data$datatype=="TMIN","value"])/10 #Convert to degrees C
fire.list[f,"T_dist"]=round(T_Station$distance,1) #In KM
}
}

#Find the nearest station with wind data and get that data:
if(any(stations$element=="WSF5" & 
       stations$first_year<fire.list[f,"FIRE_YEAR"] & 
       stations$last_year>fire.list[f,"FIRE_YEAR"])){ #If there is wind data to be had
  
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
}
gc()
print(f)
}

#write_csv(fire.list,'./Data/Derived/fire_list.csv')