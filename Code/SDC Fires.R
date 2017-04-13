##This code reads in mapped stand-replacing patches for a r-defined set of fires, and calculates SDC and other spatial statistics on those patches##
##Jens Stevens; stevensjt@gmail.com
####0. load libraries####
library(rgdal)
library(dismo)
library(ggplot2)
library(ggrepel)
require(rgeos)
library(gridExtra)
library(plyr)
library(dplyr)
source("./Code/Functions.R")

####1a. Load burn severity data from scratch (optional)####
#fires <- readOGR("../Large Files/GIS/BurnSev/Current", layer="VegBurnSeverity84-17") 
#Full layer. Very large, takes a prohibitively long time.

#hs_patches <- fires[fires$BURNSEV==4&fires$BEST_ASSES=="YES",] 
#Extract only the high-severity patches, which shrinks down the files size. Takes a while in R, so doing it in QGIS, then loading below in #1b.

#saveRDS(hs_patches,"../Large Files/GIS/hs_patches.RDS")  #takes a while

#rm(fires) #Clear space

####1b. Load high severity patches (if already created)####
hs_patches=readOGR("../Large Files/GIS/BurnSev/Current/", layer="hs_patches")
gc()#Clear space

####1c. Define fires of interest####
#fire.list=read.csv("/Users/Jens/Documents/Davis/Post-Doc/Side Projects/Mixed Severity/Jay's exploratory analysis/Patch_analysis_Results.csv") #Used for first draft of concepts paper
#fire.list=fire.list[(as.character(fire.list$AGENCY)!="USF"&as.character(fire.list$VB_ID)!="2014KING")&as.character(fire.list$Veg)=="Forest",] #Used for first draft; refine to NPS/CDF (optional; used for first draft of concepts paper)
#fire.list=fire.list[(as.character(fire.list$AGENCY)=="USF"|as.character(fire.list$VB_ID)=="2014KING")&as.character(fire.list$Veg)=="Forest",] #Used for first draft; refine to USFS (optional)

#***WHERE IS WHITES FIRE? It's in HS_patches but not in fire.list when read from Data/Derived; Might need to cross-check the current fire.list with fires_usfs and pull in any that are missing to get the SDC values#
#fire.list= read.csv("./Data/Derived/fire_list.csv")

#Create character vector of fires to sample
fires.to.sample <- as.character(fire.list$VB_ID) #Sample all fires
#fires.to.sample <- c("1987EAST","2008CARIBOU") #Sample specific fires

####2. Run Internal Buffering, create "Long" dataset####
Sys.time() #Takes ~12 hours depending on sample size. Should only need to do once. Could be made more efficient in parallel.
for(f in c(1:length(fires.to.sample))){
  cancel=!fires.to.sample[f]%in%hs_patches$VB_ID#If the fire name in question does not have a corresponding shapefile, set cancel to T and bypass the analysis for that fire.
  if(!cancel){#Proceed if you have a valid fire to work with (don't cancel)
    hs_fire=hs_patches[hs_patches$VB_ID==fires.to.sample[f],]
    #Remove holes
    hs_fire=fill_holes(hs_fire=hs_fire)
    #plot(hs_fire,col="darkred",border="transparent")
    Sys.time()
    decay.table=decay(hs_fire=hs_fire,buf_max=1000,buf_inc=10,name=fires.to.sample[f])
    Sys.time()
    ifelse(f==1,fires_long=decay.table,fires_long=rbind(fires_long,decay.table))
  }
  print(paste(fires.to.sample[f],Sys.time()))
  gc()
}

#Save as "Long" dataset, which has the buffering done for each fire (rows=buffer widths).
#Code is split into USFS and NPS/CDF so it takes less time.
#write.csv(fires_long,file="./Data/Derived/Long_Form/NPS_CDF_Long.csv") #Specify whether NPS_CDF or USFS.


####3. Calculate SDC####
#Read in processed full geospatial data (optional; if you saved fires_long earlier):
fires_long=rbind(read.csv("./Data/Derived/Long_Form/USFS_King_Long.csv"),
                      read.csv("./Data/Derived/Long_Form/NPS_CDF_Long.csv")) #All fires; First draft

#Subset to specific fires (optional):
#fires_long=fires_long[as.character(fires_long$name)%in%c("1987EAST","2008CARIBOU"),]
#fires_long$name=as.character(fires_long$name)

#Analyze the fires_long table to calculate sdc
fires_with_stats <- calculate.sdc(fires_long) #Fast!

#Save as "Full" dataset, which has the statistic calculated for each row (rows=buffer widths).
#write.csv(fires_with_stats,file="./Data/Derived/all_fires_Full.csv")

#Summarize to a single row per fire 
summary_fires=ddply(fires_with_stats,.(name),summarize,sdc=mean(sdc))
names(summary_fires)[1] <- "VB_ID"
#Add the sdc parameter to the fire.list file, with a single row per fire ("ForAnalysis" dataset)
fire.list <- merge(fire.list[,c(1:17)],summary_fires[,c("VB_ID","sdc")])
#write.csv(fire.list,file="./Data/Derived/all_fires_ForAnalysis.csv")

####4. Calculate Fragstats####
##4a. Get mean shape index (MSI), area-weighted mean shape index (AWMSI), mean patch fractal dimension (MPFD) and area-weighted mean patch fractal dimension (AWMPFD) for all fires, or specific fires.

fires.to.sample <- as.character(fire.list$VB_ID) #Sample all fires
#fires.to.sample <- c("1987EAST","2008CARIBOU") #Sample specific fires

Sys.time() 
for(f in c(1:length(fires.to.sample))){
  cancel=!fires.to.sample[f]%in%hs_patches$VB_ID#If the fire name in question does not have a corresponding shapefile, set cancel to T and bypass the analysis for that fire.
  if(!cancel){#Proceed if you have a valid fire to work with (don't cancel)
    hs_fire=hs_patches[hs_patches$VB_ID==fires.to.sample[f],]
    #Remove holes
    hs_fire=fill_holes(hs_fire=hs_fire)
    Sys.time()
    SI=get_fragstats(hs_fire=hs_fire)
    fire.list[fire.list$VB_ID==fires.to.sample[f],"MSI"] <- SI[1]
    fire.list[fire.list$VB_ID==fires.to.sample[f],"AWMSI"] <- SI[2]
    fire.list[fire.list$VB_ID==fires.to.sample[f],"MPFD"] <- SI[3]
    fire.list[fire.list$VB_ID==fires.to.sample[f],"AWMPFD"] <- SI[4]
    Sys.time()
  }
  print(paste(fires.to.sample[f],Sys.time()))
  gc()
}
Sys.time()
#write.csv(fire.list,file="./Data/Derived/all_fires_ForAnalysis.csv")

##4b. 