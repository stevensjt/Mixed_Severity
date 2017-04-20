##This code reads in mapped stand-replacing patches for a specified set of fires, and calculates SDC and other spatial statistics on those patches##
##Jens Stevens; stevensjt@gmail.com
####0. load libraries####
library(rgdal) #readOGR
library(tidyverse)

library(dismo)
library(ggrepel)
require(rgeos)
library(gridExtra)
source("./Code/Functions.R")

####1b. Load high severity patches####
#The process of extracting the high severity patches from the raw USFS severity shapefiles in R is very time consuming. I'm doing this in ArcGIS (on PC) and creating shapefile called "hs_patches". The high-severity patches are categorized by setting the RdNBR value associated with 90% basal area mortality as the minimum threshold value.
hs_patches=readOGR("../Large Files/GIS/BurnSev/Current/", layer="hs_patches")
gc()#Clear space

####1c. Define fires of interest####
#fire.list=read.csv("/Users/Jens/Documents/Davis/Post-Doc/Side Projects/Mixed Severity/Jay's exploratory analysis/Patch_analysis_Results.csv") #Used for first draft of concepts paper
#fire.list=fire.list[(as.character(fire.list$AGENCY)!="USF"&as.character(fire.list$VB_ID)!="2014KING")&as.character(fire.list$Veg)=="Forest",] #Used for first draft; refine to NPS/CDF (optional; used for first draft of concepts paper)
#fire.list=fire.list[(as.character(fire.list$AGENCY)=="USF"|as.character(fire.list$VB_ID)=="2014KING")&as.character(fire.list$Veg)=="Forest",] #Used for first draft; refine to USFS (optional)

#Current Concepts MS code below
#Read in raw fire list file from Jay Miller
fire.list <- 
  read.csv("./Data/Raw/fires_usfs.csv") %>%
  filter(Veg == "Forest",(PCT_FS > 0.5 | (PCT_FS < 0.5 & AGENCY == "NPS") ) )
fire.list <- fire.list[fire.list$VB_ID %in% hs_patches$VB_ID,]
#This dataset now has fires >80 ha from 1984 through 2016 that burned predominantly through conifer forest vegetation, and were either at least 50% on Forest Service lands or were at least 50% on Park Service lands (managed by NPS), and have a corresponding shapefile of high-severity patches. Missing 10 fires from 2016 that haven't been mapped yet.

#Create character vector of fires to sample
fires.to.sample <- as.character(fire.list$VB_ID) #Sample all fires
#fires.to.sample <- c("1987EAST","2008CARIBOU") #Sample specific fires

####2. Run Internal Buffering ['decay()'], create "Long" dataset####
Sys.time() #Takes ~11 hours in parallel (Rim takes 2 hours)

for(f in c(122:length(fires.to.sample))){
  hs_fire=hs_patches[hs_patches$VB_ID==fires.to.sample[f],]
  #Remove holes
  hs_fire=fill_holes(hs_fire=hs_fire)
  #CHECKME this step below is necessary to buffer from internal "holes" within the patch. It was missing from the first draft of the paper, could account for some differences.
  hs_fire2 <- createSPComment(hs_fire)
  #plot(hs_fire,col="darkred",border="transparent")
  #Sys.time()
  decay.table=decay(hs_fire=hs_fire2,buf_max=1000,buf_inc=10,name=fires.to.sample[f])
  #Sys.time()
  ifelse(f==1, fires_long <- decay.table, fires_long <- rbind(fires_long,decay.table) )
  print(paste(fires.to.sample[f],Sys.time()))
  gc()
}

#Save as "Long" dataset, which has the buffering done for each fire (rows=buffer widths).
write_csv(fires_long,"./Data/Derived/Long_Form/all_fires_Long.csv")

####3. Calculate SDC####
#Read in processed full geospatial data (optional; if you saved fires_long earlier):
#fires_long=rbind(read.csv("./Data/Derived/Long_Form/USFS_King_Long.csv"), read.csv("./Data/Derived/Long_Form/NPS_CDF_Long.csv")) #All fires; First draft

#Read in long dataset (optional; if you saved fires_long earlier)
#fires_long <- read.csv("./Data/Derived/Long_Form/all_fires_Long.csv")

#Subset to specific fires (optional):
#fires_long=fires_long[as.character(fires_long$name)%in%c("1987EAST","2008CARIBOU"),]
#fires_long$name=as.character(fires_long$name)

#Calculate SDC for each fire (fast), and summarize to a single row per fire 
fires_long <- calculate.sdc(fires_long)
summary_fires <- 
  fires_long %>%
  group_by(name) %>%
  summarise(sdc = mean(sdc))
names(summary_fires)[1] <- "VB_ID"
#Add the sdc parameter to the fire.list file, with a single row per fire ("ForAnalysis" dataset)
fire.list <- merge(fire.list,summary_fires[,c("VB_ID","sdc")])

#Save as "Long" dataset, including the sdc values (need for plotting curves).
write_csv(fires_long,"./Data/Derived/Long_Form/all_fires_Long.csv")

####4. Calculate Fragstats####
##4a. Get mean shape index (MSI), area-weighted mean shape index (AWMSI), mean patch fractal dimension (MPFD) and area-weighted mean patch fractal dimension (AWMPFD) for all fires, or specific fires.

fires.to.sample <- as.character(fire.list$VB_ID) #Sample all fires
#fires.to.sample <- c("1987EAST","2008CARIBOU") #Sample specific fires

for(f in c(1:length(fires.to.sample))){ #Takes ~3 minutes
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
write_csv(fire.list,path="./Data/Derived/all_fires_ForAnalysis.csv")
