library(tidyverse)
library(plyr)
library(zoo) #For moving window
library(rgdal)
library(dismo)
library(ggrepel)
require(rgeos)
library(gridExtra)
source("./Code/Functions.R")

####1. Process data and calculate SDC####

#1a. Load and Filter Data####
#fires=readOGR("/Users/Jens/Documents/Davis/Post-Doc/Side Projects/Mixed Severity/Large Files/GIS/BurnSev/Current/", layer="VegBurnSeverity84-17") #Full layer. CRS EPSG:3310, NAD83 CA Albers
#fires=readOGR("/Users/Jens/Documents/Davis/Post-Doc/Side Projects/Mixed Severity/GIS/BurnSev/", layer="VegBurnSeverity_Sierra_85-15")#Shapefile with all Sierra fires since 2000 (smaller data file for EDA)



#hs_patches=fires[fires$BURNSEV==4&fires$BEST_ASSES=="YES",] #Extract only the high-severity patches, which shrinks down the files size
#saveRDS(hs_patches,"./hs_patches.RDS")
#hs_patches=ReadRDS("./hs_patches.RDS")
hs_patches=readOGR("../Large Files/GIS/BurnSev/Current/", layer="hs_patches")
#rm(fires)#Clear space
gc()#Clear space

#Get fires of interest
#***WHERE IS WHITES FIRE? It's in HS_patches but not in fire.list when read from Data/Derived; Might need to cross-check the current fire.list with fires_usfs and pull in any that are missing to get the SDC values#

fire.list= read.csv("./Data/Derived/fire_list.csv")
#fire.list=fire.list[(as.character(fire.list$AGENCY)!="USF"&as.character(fire.list$VB_ID)!="2014KING")&as.character(fire.list$Veg)=="Forest",] 
fires.to.sample=as.character(fire.list$VB_ID)

#1b. Run Geospatial Analysis - Internal Buffering####
Sys.time() #Takes ~12 hours depending on sample size. Should only need to do once.
for(f in c(1,100,200,300,400,which(is.na(fire.list$q)))){ #
  cancel=!fires.to.sample[f]%in%hs_patches$VB_ID#If the fire name in question does not have a corresponding shapefile, set cancel to T and bypass the analysis for that fire.
  if(!cancel){#Proceed if you have a valid fire to work with (don't cancel)
    hs_fire=hs_patches[hs_patches$VB_ID==fires.to.sample[f],]
    #Remove holes
    hs_fire=fill_holes(hs_fire=hs_fire)
    #plot(hs_fire,col="darkred",border="transparent")
    Sys.time()
    decay.table=decay(hs_fire=hs_fire,buf_max=1000,buf_inc=10,name=fires.to.sample[f])
    Sys.time()
    if(f==1){all_fires=decay.table}else{
      all_fires=rbind(all_fires,decay.table)
    }
  }
  print(paste(fires.to.sample[f],Sys.time()))
  gc()
}
#Save as "Long" dataset, which has the buffering done for each fire (rows=buffer widths).
#write.csv(all_fires,file="./Data/Derived/Long_Form/2015_16_Long.csv")

#1c. Run Statistical Analysis - calculate metric####
fires_for_stats=all_fires #If running on data you just did geospatial analyses

#Analyze the fires_for_stats table to calculate the parameters of interest
fires_with_stats=calculate.sdc(fires_for_stats) #Fast!

#Save as "Full" dataset, which has the statistic calculated for each row (rows=buffer widths).
#write.csv(fires_with_stats,file="./Analyses/Processed Data/all_fires_Full.csv")

#Summarize to a single row per fire 
summary_fires=ddply(fires_with_stats,.(name),summarize,sdc=mean(sdc))
names(summary_fires)[1]="VB_ID"
tmp=merge.data.frame(fire.list,summary_fires,by="VB_ID",all=T)
tmp[which(is.na(tmp$sdc)),"sdc"]=tmp[which(is.na(tmp$sdc)),"q"]

#Add the relevant parameters to the fire.list file, and save (optional)
fire.list$sdc=tmp$sdc
#write_csv(fire.list,path='./Data/Dervied/fire_list2.csv')
#hs_patches2=hs_patches[which(hs_patches$VB_ID%in%fire.list$VB_ID),]
#saveRDS(hs_patches2,"./Data/hs_patches.RDS")

####3. Exploratory analyses####
#3aa. Read in data
#Read in processed data
fire.list= read.csv("./Data/Derived/fire_list.csv")

hs_patches2=readRDS("./Data/hs_patches.RDS")

#3a: Check for normality
hist(fire.list$q)
fire.list$logq=log(fire.list$q)
hist(fire.list$logq)

#3b: Moving-window estimates (averaged per year)
fire.list.annual=ddply(fire.list,.(FIRE_YEAR),summarize,q=mean(q),BA90_PCT=mean(BA90_PCT))

zoo.df=zoo(fire.list.annual$q)
mvw=rollapply(zoo.df, width = 5, by = 1, FUN = mean, align = "center")
fire.list.annual[rownames(as.data.frame(mvw)),"q.mvw5"]=mvw

zoo.df=zoo(fire.list.annual$BA90_PCT)
mvw=rollapply(zoo.df, width = 5, by = 1, FUN = mean, align = "center")
fire.list.annual[rownames(as.data.frame(mvw)),"BA90_PCT.mvw5"]=mvw

####4. Trends over time####
#4a: is q decreasing over time?
m4a=lm(logq~FIRE_YEAR,data=fire.list)
summary(m4a)
ggplot(fire.list,aes(x=FIRE_YEAR,y=logq))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="All Fires")+
  theme_bw()
#conclusion: yes, slightly and significantly

#4b: is percent high severity increasing over time?
m4b=lm(BA90_PCT~FIRE_YEAR,data=fire.list)
summary(m4b)
ggplot(fire.list,aes(x=FIRE_YEAR,y=BA90_PCT))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="All Fires")+
  theme_bw()
#conclusion: yes, slightly but not significantly

#4c: is q decreasing over time in moving-window analysis?
m4c=lm(log(q.mvw5)~FIRE_YEAR,data=fire.list.annual)
summary(m4c)
ggplot(fire.list.annual,aes(x=FIRE_YEAR,y=log(q.mvw5)))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="All Fires")+
  theme_bw()
#conclusion: no, but looks like some weather-driven patterns.

#4d: is percent high severity increasing over time in moving-window analysis?
m4d=lm(BA90_PCT.mvw5~FIRE_YEAR,data=fire.list.annual)
summary(m4d)
ggplot(fire.list.annual,aes(x=FIRE_YEAR,y=BA90_PCT.mvw5))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="All Fires")+
  theme_bw()
#conclusion: no, but looks like some weather-driven patterns.

#4e: percentile analysis
fire.list$BA90_PCT_BIN=data.frame(fire.list$BA90_PCT, bin=cut(fire.list$BA90_PCT, c(seq(from=0,to=10),15,20,25,30,35,40,45,50,60,70,80), include.lowest=TRUE))$bin
for(l in levels(fire.list$BA90_PCT_BIN)){
  fire.list[grep(l,as.character(fire.list$BA90_PCT_BIN),fixed=T),"q.pct"]=
    ecdf(fire.list[grep(l,as.character(fire.list$BA90_PCT_BIN),fixed=T),"q"]) (fire.list[grep(l,as.character(fire.list$BA90_PCT_BIN),fixed=T),"q"])
}
hist(fire.list$q.pct)
ggplot(fire.list[fire.list$BA90_PCT>20,],aes(x=FIRE_YEAR,y=q.pct))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()

#5: 

ggplot(fire.list,aes(x=BA90_PCT,y=logq,col=AGENCY))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="All Fires")+
  theme_bw()

ggplot(fire.list,aes(x=BA90_PCT,y=logq,col=WFU))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="All Fires")+
  theme_bw()