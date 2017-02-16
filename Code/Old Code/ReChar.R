library(rgdal)
library(ggplot2)
library(raster)
library(rasterVis)
library(rgeos)
####1. California Data####
fires=readOGR("/Users/Jens/Documents/Davis/Post-Doc/Side Projects/Mixed Severity/GIS/BurnSev/", layer="VegBurnSeverity_Sierra_00-15")

#View(fires@data)

hs_patches=fires[fires$BURNSEV==4&fires$BEST_ASSES=="YES",]
#hs_patches=readRDS("/Users/Jens/Documents/Davis/Post-Doc/Side Projects/Mixed Severity/GIS/BurnSev/hs_patches_NorCal_85-15.RDS")
#View(hs_patches@data)

fires.to.sample=c("2004POWER","2007ANGORA","2007ANTELOPE_CMPLX","2007MOONLIGHT","2008RICH","2008PIUTE","2012CHIPS","2013RIM","2014KING","2015BUTTE","2015ROUGH")
#Options include 2012CHIPS, 2001HOOVER, 2004MEADOW

for(f in c(1:length(fires.to.sample))){
fire=hs_patches[hs_patches$VB_ID==fires.to.sample[f],]
#fire=hs_patches[hs_patches$VB_ID=="2001HOOVER",] #If running individually
#plot(fire)

#Calculate patch areas
d.fire=data.frame(name=fire@data$VB_ID,
                  p.area.m2=round(sapply(slot(fire@polygons[[1]], "Polygons"), slot, "area"),0)
)
d.fire$p.area.ha=d.fire$ p.area.m2*0.0001
d.fire$p.area.ha.bin=cut(d.fire$p.area.ha,
                         breaks=c(seq(0, 50, 5),100, 500,1000,Inf),
                         labels=c(seq(0,45, 5),">50",">100",">500",">1000"), right=F
)
d.fire$prop.contrib=d.fire$p.area.ha/sum(d.fire$p.area.ha)

if(f==1){
  d.fires=d.fire
}else{
  d.fires=rbind(d.fires,d.fire)
}
}

#Put NA's in missing combinations
d.fires2=rbind(d.fires[,c("name","p.area.ha.bin","prop.contrib")], cbind(expand.grid(
  name=unique(d.fires$name), 
  p.area.ha.bin=unique(d.fires$p.area.ha.bin)), 
  prop.contrib=NA))

ggplot(d.fires2)+
  geom_bar(aes(x=p.area.ha.bin,weight=prop.contrib,fill=name),position="dodge")+
  labs(x="Patch size class (ha)",y="Proportion of total high-severity area")+
  theme_bw()

#Can I cut HS patches at pinch points?
str(fire@polygons[[1]]@Polygons[[1]])
which.max(sapply(slot(fire@polygons[[1]], "Polygons"), slot, "area"))
#Subset? http://stackoverflow.com/questions/29978826/r-subsetting-and-plotting-a-spatialpoints-object
test=fire@polygons[[1]]@Polygons[[9]]
str(test)
Srs1=Polygons(list(test),"s1")
SpP=SpatialPolygons(list(Srs1))
SpP=readRDS("/Users/Jens/Documents/Davis/Post-Doc/Side Projects/Mixed Severity/SpP.RDS")
plot(test@coords,cex=0.2)
plot(SpP[1],add=F)
#Reproducible?
Sr1 = Polygon(cbind(c(2,4,4,1,2),c(2,3,5,4,2)))
#Simplify?
library(rgeos)
simp=gSimplify(SpP,tol=60,topologyPreserve = F)
plot(simp,add=T)
lpi <- gIntersection(poly, clines)# intersect your line withthe polygon
blpi <- gBuffer(lpi, width = 0.000001) # create a very thin polygon buffer of the intersected line
dpi <- gDifference(poly, blpi) # split using gDifference

#Minimum distance between edges:
#https://stat.ethz.ch/pipermail/r-sig-geo/2009-March/005189.html
#https://stat.ethz.ch/pipermail/r-sig-geo/2009-March/005216.html

#Internal buffer:
#https://cran.r-project.org/web/packages/rgeos/rgeos.pdf
buf50=gBuffer(SpP,width=-50) #byid=F ?
buf100=gBuffer(SpP,width=-100)
buf500=gBuffer(SpP,width=-500)
plot(SpP)
plot(buf50,col=rgb(1,0,0,alpha=0.3),add=T)
plot(buf100,col=rgb(0,1,0,alpha=0.2),add=T)
#plot(buf500,col=rgb(0,0,1,alpha=0.1),add=T)
gArea(buf100) * 0.0001 #30 ha
gArea(buf50) * 0.0001 #49 ha

#ggplot approach (not working)
#SpP_f = fortify(SpP)
#p=ggplot()+
#  geom_path(SpP_f,aes(long,lat))

####2. Other Data (Parks)####
list.dirs <- function(path=".", pattern=NULL, all.dirs=FALSE,
                      full.names=FALSE, ignore.case=FALSE) {
  # use full.names=TRUE to pass to file.info
  all <- list.files(path, pattern, all.dirs,
                    full.names=TRUE, recursive=FALSE, ignore.case)
  dirs <- all[file.info(all)$isdir]
  # determine whether to return full names or just dir names
  if(isTRUE(full.names))
    return(dirs)
  else
    return(basename(dirs))
}

fire.names.Parks=list.dirs(path="./GIS/BurnSev/RDS-2015-0021/Data/severity",pattern='.20')
fire.names.Parks=fire.names.Parks[-grep('SBFC',fire.names.Parks)]
dNBR.paths=paste("./GIS/BurnSev/RDS-2015-FCW/Data/severity/",fire.names.Parks,"/dNBR/dNBR.tif",sep="")
perim.paths=paste("./GIS/BurnSev/RDS-2015-FCW/Data/severity/",fire.names.Parks,"/Perimeter/fire_perimeter.tif",sep="")

#Work with a single raster
for(f in c(6:10)){
perim=rasterToPolygons(raster(perim.paths[f]),dissolve=T)
r=raster(dNBR.paths[f])
rm=mask(r,perim)
#rc=reclassify(rm,c(-Inf,0,0,0,200,1,200,400,2,400,Inf,3))
#plot(rc,col=c('forest green','darkgoldenrod1','darkorange3','darkred'))
#rcv=rasterToPolygons(rc,dissolve=T)
#plot(rcv,col=c('forest green','darkgoldenrod','orange','darkred'))
rc.hs=reclassify(rm,c(-Inf,400,NA,400,Inf,1))
rcv.hs=rasterToPolygons(rc.hs,dissolve=T) #A vector with separate polygons for HS patches

d.fire=data.frame(name=fire.names.Parks[f],
                  p.area.m2=round(sapply(slot(rcv.hs@polygons[[1]], "Polygons"), slot, "area"),0)
)
d.fire$p.area.ha=d.fire$ p.area.m2*0.0001
d.fire$p.area.ha.bin=cut(d.fire$p.area.ha,
                         breaks=c(seq(0, 50, 5),100, 500,1000,Inf),
                         labels=c(seq(0,45, 5),">50",">100",">500",">1000"), right=F
)
d.fire$prop.contrib=d.fire$p.area.ha/sum(d.fire$p.area.ha)
if(f==1){
  d.fires=d.fire
}else{
  d.fires=rbind(d.fires,d.fire)
}
}

d.fires2=rbind(d.fires[,c("name","p.area.ha.bin","prop.contrib")], cbind(expand.grid(
  name=unique(d.fires$name), 
  p.area.ha.bin=unique(d.fires$p.area.ha.bin)), 
  prop.contrib=NA))

ggplot(d.fires2)+
  geom_bar(aes(x=p.area.ha.bin,weight=prop.contrib,fill=name),position="dodge")+
  labs(x="Patch size class (ha)",y="Proportion of total high-severity area")+
  theme_bw()

####3. NorCal vs Klamath, generate statistics####
#START HERE
fires=readOGR("/Users/Jens/Documents/Davis/Post-Doc/Side Projects/Mixed Severity/GIS/BurnSev/", layer="VegBurnSeverity_NorCal_85-15")

#hs_patches=fires[fires$BURNSEV==4&fires$BEST_ASSES=="YES",]
#View(hs_patches@data)
#saveRDS(hs_patches,"/Users/Jens/Documents/Davis/Post-Doc/Side Projects/Mixed Severity/GIS/BurnSev/hs_patches_NorCal_85-15.RDS")
hs_patches=readRDS("/Users/Jens/Documents/Davis/Post-Doc/Side Projects/Mixed Severity/GIS/BurnSev/hs_patches_NorCal_85-15.RDS")
NW=readOGR("/Users/Jens/Documents/Davis/Post-Doc/Side Projects/Mixed Severity/GIS/Boundaries/NAD83 CA Albers/", layer="NW_Forests")
hs_patches_poly=SpatialPolygons(hs_patches@polygons,
  proj4string = hs_patches@proj4string,
  match.ID='6')
NW_poly=SpatialPolygons(NW@polygons,proj4string = NW@proj4string)
#http://gis.stackexchange.com/questions/37503/rename-a-spatialpolygon-class-object-in-r
require(maptools)

gI=gIntersection(hs_patches_poly,NW_poly,byid=F)
