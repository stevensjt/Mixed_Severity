#Scrape climate data?
 
library(ggplot2)
library(plyr)
library(rgeos)# For gCentroid
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

