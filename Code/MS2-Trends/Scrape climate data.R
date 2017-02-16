#Scrape climate data?
 
library(ggplot2)
library(plyr)
library(rgeos)# For gCentroid
library(rnoaa)
source(file="./Analyses/Functions.R")

####1. Read in data####
#Read in processed data
fire.list=read.csv("./Analyses/Processed Data/all_fires_ForAnalysis.csv")
hs_patches=readRDS("./hs_patches.RDS")

#Get centroids of polygons
f=1
fires.to.sample=as.character(fire.list$VB_ID)
hs_fire=hs_patches[hs_patches$VB_ID==fires.to.sample[f],]
#Remove holes
hs_fire=fill_holes(hs_fire=hs_fire)
#proj4string(hs_fire)

plot(hs_fire,axes=T)
points(gCentroid(hs_fire))

#Possibilities for accessing noaa data:
#ftp://cran.r-project.org/pub/R/web/packages/rnoaa/rnoaa.pdf
#https://ropensci.org/tutorials/rnoaa_tutorial.html
#https://www.ncdc.noaa.gov/cdo-web/token

ncdc_datasets(token ="FTEIRBcwLUCArmZqFHwVomZhjQDOCqMv")
ncdc_locs(datasetid='GHCND',token ="FTEIRBcwLUCArmZqFHwVomZhjQDOCqMv")
station_data <- ghcnd_stations()
View(meteo_distance(station_data, 37.8375105,-122.2542293, radius = 10, limit = 30))
#Just need to get LL coords from centroid code, extend radius, and get dates of each fire.