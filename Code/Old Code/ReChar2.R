library(rgdal)
library(ggplot2)
library(raster)
library(rasterVis)
library(rgeos)
library(dplyr)

####1. Read and process data####
#fires=readOGR("/Users/Jens/Documents/Davis/Post-Doc/Side Projects/Mixed Severity/GIS/BurnSev/", layer="VegBurnSeverity85-15")
fires=readOGR("/Users/Jens/Documents/Davis/Post-Doc/Side Projects/Mixed Severity/GIS/BurnSev/", layer="VegBurnSeverity_Sierra_85-15")#Shapefile with all Sierra fires since 2000 (smaller data file for EDA)
fire.list=read.csv("/Users/Jens/Documents/Davis/Post-Doc/Side Projects/Mixed Severity/Jay's exploratory analysis/fires_84-14.csv") #Jay's file
#View(fires@data)

hs_patches=fires[fires$BURNSEV==4&fires$BEST_ASSES=="YES",] #Extract only the high-severity patches, which shrinks down the files size
rm(fires)#Clear space
gc()#Clear space

####1b. Define functions####
fill_holes=function(hs_fire){
  hs_fire_p=slot(hs_fire, "polygons")
  holes <- lapply(hs_fire_p, function(x) sapply(slot(x, "Polygons"), slot, "hole"))
  areas <- lapply(hs_fire_p, function(x) sapply(slot(x, "Polygons"), slot, "area"))
  res <- lapply(1:length(hs_fire_p), function(i) slot(hs_fire_p[[i]], "Polygons")[!(holes[[1]]&areas[[1]]<8100)])#Select the polygons that are not holes. The "i" here is an artifact of the example code; the fires here only have one polygon ID so i=1 (it's a multipart polygon).
  IDs <- row.names(hs_fire)
  hs_fire_fill <- SpatialPolygons(lapply(1:length(res), function(i)
    Polygons(res[[i]], ID=IDs[i])), proj4string=CRS(proj4string(hs_fire)))
  return(hs_fire_fill)
}

####2. Patch distribution analysis####
#2a: Proof of concept on King fire
fires.to.sample="2014KING"
buf_inc= 10 #Buffer increment in m
buf_max= 1000 #Max buffer width is 1000m; needs to be a multiple of buf_inc
dist.table=data.frame(name=rep(fires.to.sample,(buf_max/buf_inc)+1),buf_width=NA,area_ha=NA)
hs_fire=hs_patches[hs_patches$VB_ID==fires.to.sample[1],]
dist.table[dist.table$name==fires.to.sample[[1]],"buf_width"]=seq(0,buf_max,by=buf_inc)
#hs_fire=gBuffer(hs_fire,width=0) #Initial buffer is the entire hs area

Sys.time()
for(w in seq(0,buf_max,by=buf_inc)){ #w = buffer width in m. Takes 1:15 to run for 1000 with inc=5
buf=gBuffer(hs_fire,width=-w) #Buffer the previous buffer layer by an additional w meters
dist.table[dist.table$buf_width==w,"area_ha"]=area(buf)*0.0001
print(w)
gc() #Super important step to improve efficiency; memory dump after every run. Not working as well now.
}
Sys.time()
dist.table$prop.hs=dist.table$area_ha/dist.table$area_ha[1]
plot(prop.hs~buf_width,data=dist.table,main="King Fire",ylab="Prop. high-severity",xlab="Distance from patch edge (m)")
king.table=dist.table

#2b: Trying other big fires
#fires.to.sample.1=c("2014KING","2012CHIPS","2004MEADOW","2001HOOVER")
#fires.to.sample.2=c("2004POWER","2007ANGORA","2007ANTELOPE_CMPLX","2007MOONLIGHT","2008RICH","2008PIUTE","2012CHIPS","2013RIM","2014KING","2015BUTTE")
fires.to.sample.3=
  as.character(fire.list[fire.list$WFU=="yes" & between(fire.list$FIRESIZE_HA,100,1000), "VB_ID"])
fires.to.sample.4=
  as.character(fire.list[fire.list$WFU=="no" & between(fire.list$FIRESIZE_HA,100,1000), "VB_ID"])
fires.to.sample.5=  
  as.character(fire.list[fire.list$AGENCY=="USF" & between(fire.list$FIRESIZE_HA,1000,4000),"VB_ID"])
fires.to.sample.6=  
  as.character(fire.list[fire.list$AGENCY=="NPS" & between(fire.list$FIRESIZE_HA,1000,4000),"VB_ID"])
#Fires that take a long time: Moonlight (0:23), Piute (0:19), Chips (0:15) Rim (1:15), King (0:23), Butte (0:45).
#Fires that are not clean: "2015ROUGH" (might need to re-load newest data; old version hs streaks)
fires.to.sample=fires.to.sample.6

#Initialize
buf_inc= 10 #Buffer increment in m
buf_max= 1000 #Max buffer width is 1000m; needs to be a multiple of buf_inc
dist.table=data.frame(name=rep(fires.to.sample,(buf_max/buf_inc)+1),buf_width=NA,area_ha=NA)

#Run
for(f in c(1:length(fires.to.sample))){
  cancel=!fires.to.sample[f]%in%hs_patches$VB_ID#If the fire name in question does not have a corresponding shapefile, set cancel to T and bypass the analysis for that fire.
  if(!cancel){#Proceed if you have a valid fire to work with (don't cancel)
  hs_fire=hs_patches[hs_patches$VB_ID==fires.to.sample[f],]
  #Remove holes
  hs_fire=fill_holes(hs_fire=hs_fire)
  #START HERE; the problem with filling holes is it's no longer a sp data frame, and there are some warnings, and possibly buf is large too. Need to investigate more.
  #plot(hs_fire,col="darkred",border="transparent")
  dist.table[dist.table$name==fires.to.sample[f],"buf_width"]=seq(0,buf_max,by=buf_inc)
  print(c(fires.to.sample[f],Sys.time()))
  for(w in seq(0,buf_max,by=buf_inc)){ #w = buffer width in m. Takes 1:15 to run for 1000 with inc=5
    if(cancel==F){
      buf=gBuffer(hs_fire,width=-w) #Buffer the previous buffer layer by an additional w meters
      if(is.null(buf)){cancel = T;print("maxed out")} #When you run out of core area because your internal buffer is too wide, cancel the inner loop.
    }
    if(cancel==F){
      dist.table[which(dist.table$buf_width==w & dist.table$name==fires.to.sample[f]),"area_ha"]=area(buf)*0.0001
      print(w) 
    }
    gc() #Super important step to improve efficiency; memory dump after every run. Not working as well now.
  } #end of w buffer for-loop
  } #end of "if(!cancel)
  dist.table[dist.table$name==fires.to.sample[f],"prop.hs"]=
    dist.table[dist.table$name==fires.to.sample[f],"area_ha"]/
    dist.table[dist.table$name==fires.to.sample[f],"area_ha"][1]
  Sys.time()
} #end of f fire for-loop

#write.csv(dist.table,"./Analyses/dist.table_NPS_1k-5kHA.csv")

####3: Analyzing distributions####
dBig=read.csv("./Analyses/dist.table_BigFires.csv")
dBig=dBig[-which(is.na(dBig$buf_width)),]
dBig$category="Big Fires"
dWFU=read.csv("./Analyses/dist.table_WFU_1k-5kHA.csv")
dWFU=WFU[-which(is.na(dWFU$buf_width)),]
dWFU$category="WFU Fires 1-5kHA"
dnWFU=read.csv("./Analyses/dist.table_nonWFU_1k-5kHA.csv")
dnWFU=dnWFU[-which(is.na(dnWFU$buf_width)),]
dnWFU$category="non-WFU Fires 1-5kHA"

#http://stats.stackexchange.com/questions/30975/how-to-add-non-linear-trend-line-to-a-scatter-plot-in-r
#http://stackoverflow.com/questions/26560849/exponential-regression-with-nls-in-r
m.out=list()
q.out=list()

d=dBig
d=d[-which(is.na(d$buf_width)),]

#3a models to calculate q
for(f in 1:length(unique(d$name))){
  md=d[d$name==fires.to.sample[f],]
  md$buf_width[is.na(md$prop.hs)]=NA #Set missing buf_width values to NA.
  m1=lm(log(prop.hs)~buf_width,data=md)
  m2=nls(prop.hs~exp(q*buf_width),data=md,start=list(q=-0.01))
  #test viz:
  plot(prop.hs~buf_width,data=md)
  lines(na.exclude(md$buf_width), predict(m2), col = "red")
  q.out[[f]]=coef(m2)["q"][[1]]
}


#3b plots
pd=rbind(dBig,dWFU,dnWFU)

ggplot(data=pd,aes(x=buf_width,y=prop.hs))+
  #geom_point(aes(col=name))+
  #geom_smooth(formula=y~exp(q.out[[1]]*x))+
  #geom_smooth(aes(col=category),method="glm",family=binomial)+
  geom_smooth(aes(col=category))+
  labs(y="proportion high-severity",x="minimum distance to patch edge (m)")+
  theme_bw()