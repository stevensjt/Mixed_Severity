#This code was used to generate Collins et al. 2017 Landscape Ecology paper
#Contact Jens Stevens stevensjt@gmail.com

####0. Load required libraries####
####0. Load libraries####
library(dismo)
library(tidyverse)
library(rgeos)
library(gridExtra)
library(rgdal)
library(cowplot)
library(ggrepel)



####1. This section contains functions for recharacterizing fire regimes####

##Fill holes that are less than 9 pixels large (9*900m2=8100m2, or 0.81 ha)
fill_holes=function(hs_fire){
  hs_fire_p <- slot(hs_fire, "polygons")
  holes <- lapply(hs_fire_p, function(x) sapply(slot(x, "Polygons"), slot, "hole"))
  areas <- lapply(hs_fire_p, function(x) sapply(slot(x, "Polygons"), slot, "area"))
  res <- #Select the polygons that are not holes. The "i" here is an artifact of the example code; the fires here only have one polygon ID so i=1 (it's a multipart polygon).
    lapply(1:length(hs_fire_p), 
           function(i) slot(hs_fire_p[[i]], "Polygons")[!(holes[[1]]&areas[[1]]<8100)])
  IDs <- row.names(hs_fire)
  hs_fire_fill <- SpatialPolygons(lapply(1:length(res), function(i)
    Polygons(res[[i]], ID=IDs[i])), proj4string=CRS(proj4string(hs_fire)))
  return(hs_fire_fill)
  ##One consequence of this is that it's no longer a sp data frame, and there are some warnings that are ok. 
}

##Decay profile: Implement internal buffering on high severity patches, create "Long" dataset
decay=function(hs_fire,buf_max,buf_inc,name,cancel=F){
  require(parallel) #for "mclapply"
  require(rgeos) #for gBuffer
  require(raster) #for "area"
  
  #Set up long data frame:
  dist.table.sub <-
    data.frame(name=rep(name,(buf_max/buf_inc)+1),width=seq(0,buf_max,by=buf_inc),area_ha=NA)
  
  #Core operation:
  buf.list <- #Apply the buffer at every width in X; returns list of length X.
    mclapply(X=-dist.table.sub$width,
             FUN=gBuffer,
             spgeom=hs_fire, 
             byid = FALSE, id = NULL, quadsegs = 5, capStyle = "ROUND", joinStyle = "ROUND", mitreLimit = 1,
             mc.cores = 4) 
  
  #Post-processing of long data frame:
  buf.list <- #Remove any NULL values (when buffer eliminated all HS areas)
    Filter(length,buf.list) 
  dist.table.sub$area_ha <- #Calculate area remaining at each buffer distance
    (as.vector(sapply(buf.list,area))*0.0001)[1:nrow(dist.table.sub)] #0.0001 converts m2 to ha
  if(any(is.na(dist.table.sub$area_ha))){
    #If there are NULL values for the area because the buffer got too wide, remove them.
    dist.table.sub$area_ha[which(is.na(dist.table.sub$area_ha))[1]] <- 0 #Set the first null value to 0
    dist.table.sub <- #Delete the rest of the table rows with NA in area.
      dist.table.sub[-which(is.na(dist.table.sub$area_ha)),] 
  }
  dist.table.sub$prop.hs <- #Calculate the proportion of the original HS area remaining at each buffer distance
    dist.table.sub$area_ha/dist.table.sub$area_ha[1]
  
  #Return the long data frame for the fire in question
  return(dist.table.sub)
}

##Calculate stand-replacing decay coefficient (sdc)
calculate.sdc=function(decay.table){
  #https://stat.ethz.ch/R-manual/R-devel/library/base/html/by.html
  m.list=with(decay.table,
              by(decay.table,name,
                 function(x)
                   nls(prop.hs~1/(10^(sdc*width)),data=x,start=list(sdc=0.01))
              )
  )
  sdc.table=data.frame(name=unique(as.character(decay.table$name)),sdc=sapply(m.list,coef))
  out.table=merge(decay.table,sdc.table,sort=F)
  out.table$sdc.name=as.character(format(round(out.table$sdc,4),scientific=F))
  return(out.table)
}

##Calculate fragstats: mean shape index (MSI), area-weighted mean shape index (AWMSI), mean patch fractal dimension (MPFD) and area-weighted mean patch fractal dimension (AWMPFD). 

get_fragstats <- function (hs_fire){
  EAR=data.frame(poly=NA,perim=NA,area=NA)
  holes=which(sapply(hs_fire@polygons[[1]]@Polygons , slot , "hole"))
  for(p in 1:length(hs_fire@polygons[[1]]@Polygons)){
    if(!p%in%holes){
      #If target polygon not a hole, get perimeter and area, then check for holes within
      df=hs_fire@polygons[[1]]@Polygons[[p]]@coords
      out <- sapply(1:(nrow(df)-1), function(i) {
        d <- dist(df[i:(i+1),1:2])
        c(d)
      })
      EAR[p,"poly"]=p
      EAR[p,"perim"] <-sum(out)
      EAR[p,"area"] <- hs_fire@polygons[[1]]@Polygons[[p]]@area
      
      if((p+1)%in%holes){
        #If target polygon isnt a hole but next polygon is, target is parent. Store info.
        parent.p=p
        coords.storage=hs_fire@polygons[[1]]@Polygons[[p]]@coords
      }
    } else {
      #If target polygon IS a hole (works if there are multiple holes in a row)
      df=hs_fire@polygons[[1]]@Polygons[[p]]@coords
      out <- sapply(1:(nrow(df)-1), function(i) {
        d <- dist(df[i:(i+1),1:2])
        c(d)
      })
      EAR[parent.p,"perim"] <- EAR[parent.p,"perim"] + sum(out) 
      #Add hole's perimeter to parent's perimeter
      
      EAR[parent.p,"area"] <- 
        EAR[parent.p,"area"] - hs_fire@polygons[[1]]@Polygons[[p]]@area
      #Subtract hole's area from parent's area
    }
  }
  EAR$ear <- EAR$perim/EAR$area #Edge:area ratio (= perimeter:area ration = Shape Index)
  EAR$pfd <- 2*log(EAR$perim)/log(EAR$area) #Patch fractal dimension; equation from McGarigal and Marks 1995 p 36
  MSI <- mean(EAR$ear,na.rm=T) #Mean Shape Index
  AWMSI <- weighted.mean(EAR$ear,w = EAR$area,na.rm=T) #Area-weighted Mean Shape Index
  MPFD <- mean(EAR$pfd,na.rm=T) #Mean Patch Fractal Dimension
  AWMPFD <- weighted.mean(EAR$pfd,w = EAR$area,na.rm=T) #Area-weighted Mean Patch Fractal Dimension
  return(c(MSI, AWMSI, MPFD, AWMPFD))
}



####2. This section implements the calculation and plotting of SDC for hypothetical stand-replacing patches####
####2.1. Generate fake data for circles####

##2.1a: many small circles, ~1 ha each, 1000 ha total
fire.edge.length <- 6000 #Side length of a square containing the HS patches
gridsize <- 32
spacing <- fire.edge.length/gridsize
radius=55.75388
df <- expand.grid(lon=seq(from=spacing/2,by=spacing,length.out=gridsize),
               lat=seq(from=spacing/2,by=spacing,length.out=gridsize))
df$radius <- radius
p <- df[,1:2]
rad <- df[,3]
circles.a <- polygons(circles(p, rad, lonlat=FALSE, dissolve=FALSE))

##2.1b: fewer medium circles, 10 ha each, 1000 ha total
gridsize <- 10
spacing <- fire.edge.length/gridsize
radius <- 178.4124
df <- expand.grid(lon=seq(from=spacing/2,by=spacing,length.out=gridsize),
               lat=seq(from=spacing/2,by=spacing,length.out=gridsize))
df$radius <- radius
p <- df[,1:2]
rad <- df[,3]
circles.b <- polygons(circles(p, rad, lonlat=FALSE, dissolve=FALSE))

##2.1c: few large circles, ~100 ha each, 1000 ha total
gridsize <- 3
spacing <- fire.edge.length/gridsize
radius <- 594.708
df <- expand.grid(lon=seq(from=spacing/2,by=spacing,length.out=gridsize),
               lat=seq(from=spacing/2,by=spacing,length.out=gridsize))
df$radius <- radius
p <- df[,1:2]
rad <- df[,3]
circles.c <- polygons(circles(p, rad, lonlat=FALSE, dissolve=FALSE))

##2.1d: one big circle, 1000 ha in total
radius <- 1784.124
df <- data.frame(lon=fire.edge.length/2, lat=fire.edge.length/2, radius=radius) #lat long adjustment to match other plots, e.g. (162/2)-55.75388
p <- df[,1:2]
rad <- df[,3]
circles.d <- polygons(circles(p, rad, lonlat=FALSE, dissolve=FALSE))

####2.2. Create decay profile for each set of circles####

#Initialize
buf_inc <- 10 #Buffer increment in m
buf_max <- 1000 #Max buffer width is 1000m; needs to be a multiple of buf_inc
#Run the decay function on each "fire"
table.a <- decay(hs_fire=circles.a,buf_max=buf_max,buf_inc=buf_inc,name="1ha patches")
table.b <- decay(hs_fire=circles.b,buf_max=buf_max,buf_inc=buf_inc,name="10ha patches")
table.c <- decay(hs_fire=circles.c,buf_max=buf_max,buf_inc=buf_inc,name="100ha patches")
table.d <- decay(hs_fire=circles.d,buf_max=buf_max,buf_inc=buf_inc,name="1000ha patches")
all_fires <- rbind(table.a,table.b,table.c,table.d)
#Add a slight bit of variation so that the model can fit the SDC parameter
all_fires$prop.hs.jitter <- all_fires$prop.hs
all_fires$prop.hs.jitter[!all_fires$prop.hs%in%c(0,1)] <- 
  abs(all_fires$prop.hs[!all_fires$prop.hs%in%c(0,1)]+
        rnorm(all_fires$prop.hs[!all_fires$prop.hs%in%c(0,1)],sd=0.01))

####2.3. Calculate sdc for circles####

#Statistical solution, solved with nonlinear least squares regression
all_fires.sdc <- calculate.sdc(all_fires)

####2.4. Generate plots for circles####

plot.colors=c("#377eb8","#4daf4a","#984ea3","#e41a1c")
p.a=ggplot()+
  geom_path(data=circles.a,aes(x=long,y=lat,group=group),col=plot.colors[1])+
  #lims(x=c(0,5500),y=c(0,5500))+
  xlim(0,fire.edge.length)+ylim(0,fire.edge.length)+
  labs(x="meters",y="meters")+
  theme_bw()+
  labs(title="a")+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title=element_text(hjust=0))

p.b=ggplot()+
  geom_path(data=circles.b,aes(x=long,y=lat,group=group),col=plot.colors[2])+
  xlim(0,fire.edge.length)+ylim(0,fire.edge.length)+
  labs(x="meters",y="meters")+
  theme_bw()+
  labs(title="b")+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title=element_text(hjust=0))
p.c=ggplot()+
  geom_path(data=circles.c,aes(x=long,y=lat,group=group),col=plot.colors[3])+
  xlim(0,fire.edge.length)+ylim(0,fire.edge.length)+
  labs(x="meters",y="meters")+
  theme_bw()+
  labs(title="c")+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title=element_text(hjust=0))
p.d=ggplot()+
  geom_path(data=circles.d,aes(x=long,y=lat,group=group),col=plot.colors[4])+
  xlim(0,fire.edge.length)+ylim(0,fire.edge.length)+
  labs(x="meters",y="meters")+
  theme_bw()+
  labs(title="d")+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title=element_text(hjust=0))

p.fits=ggplot(data=all_fires.sdc,aes(x=width,y=prop.hs))+
  geom_point(aes(col=name),size=1)+
  scale_color_manual(values=plot.colors)+
  geom_path(aes(x=width,y=1/(10^(sdc*width)),group=name,col=name))+
  annotate("text", x = c(0,200,350,950), y = c(0.1,0.1,0.1,0.1), label = unique(all_fires.sdc$sdc.name),size=3)+
  labs(y="proportion stand-replacing",x="distance to edge (m)",col="patch size")+
  theme_bw()+
  labs(title="e")+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title=element_text(hjust=0))
#p.fits

#Fig2=
grid.arrange(arrangeGrob(p.a,p.b,p.c,p.d,ncol=2),p.fits,ncol=1,heights=c(2,1),clip=F)
#dev.copy2pdf(file=paste0("./Figures/Fig2_",Sys.Date(),".pdf"),width=8,height=12)
png(file = paste0("./Figures/Fig2_",Sys.Date(),".png"),width=8,height=12,units="in",res=200)
grid.arrange(arrangeGrob(p.a,p.b,p.c,p.d,ncol=2),p.fits,ncol=1,heights=c(2,1),clip=F)
dev.off()

####2.5. Generate fake data for non-circles####

##2.5a: one big circle, 1000 ha in total
#1000 ha is 1000/0.0001 = 10,000,000 m2. Radius to get that is sqrt(10000000/pi) = 1784.12
radius=1784.124
df <- data.frame(lon=fire.edge.length/2, lat=fire.edge.length/2, radius=radius) #lat long adjustment to match other plots, e.g. (162/2)-55.75388
p <- df[,1:2]
rad <- df[,3]
circles.d <- polygons(circles(p, rad, lonlat=FALSE, dissolve=FALSE))
target.area=area(circles.d) #9,999,491 sq m, or ~10 million m2, or ~1000 ha

##2.5b. An ellipse with area = 1000 ha
library(ellipse)
e=ellipse(0.8,npoints=200) #;plot(e,type='l') #area is 11.28606 
e.factor=target.area/11.28606 #expansion factor to scale to same area as big circle is 886003.7
e=e*sqrt(e.factor)
e=e+(max(e)*1.2) #Add offset to move to center of panel
e.list=list()
e.list$ellipse1=e
eps=list()
eps$ellipse1=Polygon(e.list)
eps1=Polygons(eps,ID=1)
e.spatialpoly=SpatialPolygons(list(eps1),proj4string = CRS("+init=epsg:3310"))
area(e.spatialpoly);eps1@area #Two ways of checking area; should be the same
eps1@area/10000 #convert to ha, target is 1000 ha

##2.5c Make an irregular ellipse
f.spatialpoly=e.spatialpoly
tmp=as.data.frame(f.spatialpoly@polygons[[1]]@Polygons$ellipse1@coords)
names(tmp)=c("x","y")
tmp$zx=NA; tmp$zy=NA
#Create new x-y points either inward or outward from the edge
d=350 #Initialize distance in/out from edge along perpendicular line
for(i in c(seq(from=2,to=nrow(tmp),by=10),seq(from=3,to=nrow(tmp),by=10))){
  m=(tmp[i+1,"y"]-tmp[i-1,"y"])/(tmp[i+1,"x"]-tmp[i-1,"x"])#Calculate slope of two adjacent points
  m2=-(tmp[i+1,"x"]-tmp[i-1,"x"])/(tmp[i+1,"y"]-tmp[i-1,"y"])#Calculate slope of the perpendicular line
  k=d/sqrt(1+m2^2) #Calculate constant k
  tmp[i,"zx"]=tmp[i,"x"]+k
  tmp[i,"zy"]=tmp[i,"y"]+(k*m2)
}#End addition loop
for(i in c(seq(from=7,to=nrow(tmp),by=10),seq(from=8,to=nrow(tmp),by=10))){
  m=(tmp[i+1,"y"]-tmp[i-1,"y"])/(tmp[i+1,"x"]-tmp[i-1,"x"])#Calculate slope of two adjacent points
  m2=-(tmp[i+1,"x"]-tmp[i-1,"x"])/(tmp[i+1,"y"]-tmp[i-1,"y"])#Calculate slope of the perpendicular line
  k=d/sqrt(1+m2^2) #Calculate constant k
  tmp[i,"zx"]=tmp[i,"x"]-k
  tmp[i,"zy"]=tmp[i,"y"]-(k*m2)
}#End subtraction loop
tmp[!is.na(tmp$zx),"x"]=tmp[!is.na(tmp$zx),"zx"]
tmp[!is.na(tmp$zy),"y"]=tmp[!is.na(tmp$zx),"zy"]
tmp.mat=as.matrix(tmp[,-c(3:4)])
f.spatialpoly@polygons[[1]]@Polygons$ellipse1@coords=tmp.mat
area(f.spatialpoly)/10000 #Should match and be ~1000 ha (999.9495)
plot(f.spatialpoly) #Check and see that it looks good and irregular


####2.6. Create decay profile for each set of non-circles####
#Initialize
buf_inc= 10 #Buffer increment in m
buf_max= 1000 #Max buffer width is 1000m; needs to be a multiple of buf_inc
#Run the decay function on each "fire"
table.d=decay(hs_fire=circles.d,buf_max=1000,buf_inc=10,name="circle")
table.e=decay(hs_fire=e.spatialpoly,buf_max=1000,buf_inc=10,name="ellipse")
table.f=decay(hs_fire=f.spatialpoly,buf_max=1000,buf_inc=10,name="irregular ellipse")
#table.g=decay(hs_fire=g.spatialpoly,buf_max=1000,buf_inc=10,name="irregular circle") #Not doing irregular circle for now.
all_fires=rbind(table.d,table.e,table.f)
#Add a slight bit of variation so that the model can fit the R parameter
all_fires$prop.hs.jitter=all_fires$prop.hs
all_fires$prop.hs.jitter[!all_fires$prop.hs%in%c(0,1)]=
  abs(all_fires$prop.hs[!all_fires$prop.hs%in%c(0,1)]+
        rnorm(all_fires$prop.hs[!all_fires$prop.hs%in%c(0,1)],sd=0.01))

####2.7. Calculate sdc for non-circles####
all_fires.sdc.supp=calculate.sdc(all_fires)

####2.8. Generate plots for non-circles####
plot.colors=c("#e41a1c","#ff7f00","#a65628") #,"#df65b0" #Optional color for the irregular circle.
p.d=ggplot()+
  geom_path(data=circles.d,aes(x=long,y=lat,group=group),col=plot.colors[1])+
  #lims(x=c(0,5500),y=c(0,5500))+
  xlim(0,fire.edge.length)+ylim(0,fire.edge.length)+coord_fixed()+
  labs(x="meters",y="meters")+
  theme_bw()+
  labs(title="a")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),
        plot.title=element_text(hjust=0))
p.e=ggplot()+
  geom_path(data=e.spatialpoly,aes(x=long,y=lat,group=group),col=plot.colors[2])+
  xlim(0,fire.edge.length)+ylim(0,fire.edge.length)+coord_fixed()+
  labs(x="meters",y="meters")+
  theme_bw()+
  labs(title="b")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),
        plot.title=element_text(hjust=0))
p.f=ggplot()+
  geom_path(data=f.spatialpoly,aes(x=long,y=lat,group=group),col=plot.colors[3])+
  xlim(0,fire.edge.length)+ylim(0,fire.edge.length)+coord_fixed()+
  labs(x="meters",y="meters")+
  theme_bw()+
  labs(title="c")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),
        plot.title=element_text(hjust=0))
p.fits=ggplot(data=all_fires.sdc.supp,aes(x=width,y=prop.hs))+
  geom_point(aes(col=name),size=1)+
  scale_color_manual(values=plot.colors)+
  geom_path(aes(x=width,y=1/(10^(sdc*width)),group=name,col=name))+
  annotate("text", x = c(950,700,400), y = c(0.35,0.1,0.1), label = unique(all_fires.sdc.supp$sdc.name),size=3)+
  labs(y="proportion \nstand-replacing",x="distance to edge (m)",col="patch shape")+
  theme_bw()+
  labs(title="d")+
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=18,face="bold"),
        plot.title=element_text(hjust=0),
        legend.position=c(0.75,0.75))
#p.fits
p.empty=ggplot(data=all_fires.sdc.supp,aes(x=width,y=prop.hs))+geom_blank()+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        panel.background = element_blank())
grid.arrange(p.d,p.e,p.f,p.fits,ncol=2,heights=c(1,1),clip=F)
#dev.copy2pdf(file=paste0("./Figures/Fig3_",Sys.Date(),".pdf"),width=8,height=8)

####3. This section implements the calculation of SDC for actual stand-replacing patches. Long.####
##This code reads in mapped stand-replacing patches for a specified set of fires, and calculates SDC and other spatial statistics on those patches. It takes a long time when running for 477 fires; advised not to run section 3.2##

####3.1 Load data####
#The process of extracting the high severity patches from the raw USFS severity shapefiles in R is very time consuming. I'm doing this in ArcGIS (on PC) and creating shapefile called "hs_patches". The high-severity patches are categorized by setting the RdNBR value associated with 90% basal area mortality as the minimum threshold value.
hs_patches=readOGR("../Large Files/GIS/BurnSev/Current/", layer="hs_patches")
gc()#Clear space

#Read in raw fire list file from Jay Miller
fire.list <- 
  read.csv("./Data/Raw/fires_usfs.csv") %>%
  filter(Veg == "Forest",(PCT_FS > 0.5 | (PCT_FS < 0.5 & AGENCY == "NPS") ) )
fire.list <- fire.list[fire.list$VB_ID %in% hs_patches$VB_ID,]
#This dataset now has fires >80 ha from 1984 through 2016 that burned predominantly through conifer forest vegetation (>400 ha in NW CA), and were either at least 50% on Forest Service lands or were at least 50% on Park Service lands (managed by NPS), and have a corresponding shapefile of high-severity patches. Missing 10 fires from 2016 that haven't been mapped yet. We run this for 2016 fires but don't include them in plots.

#Create character vector of fires to sample
fires.to.sample <- as.character(fire.list$VB_ID) #Sample all fires

####3.2. Run Internal Buffering ['decay()'], create "Long" dataset. SKIP unless recalculating SDC####
Sys.time() #Takes ~11 hours in parallel (Rim takes 2 hours)

for(f in c(122:length(fires.to.sample))){
  hs_fire=hs_patches[hs_patches$VB_ID==fires.to.sample[f],] #Get a specific fire to work on
  hs_fire=fill_holes(hs_fire=hs_fire) #Remove holes <0.81 ha
  hs_fire2 <- createSPComment(hs_fire) #This step is necessary to buffer from larger internal "holes" within the patch; otherwise gBuffer won't recognize them.
  decay.table=decay(hs_fire=hs_fire2,buf_max=1000,buf_inc=10,name=fires.to.sample[f])
  ifelse(f==1, fires_long <- decay.table, fires_long <- rbind(fires_long,decay.table) )
  print(paste(fires.to.sample[f],Sys.time()))
  gc()
}

#Save as "Long" dataset, which has the buffering done for each fire (rows=buffer widths).
#write_csv(fires_long,"./Data/Derived/Long_Form/all_fires_Long.csv")

####3.3. Calculate SDC for actual fires####

#Read in long dataset (optional; if you saved fires_long earlier)
fires_long <- read.csv("./Data/Derived/Long_Form/all_fires_Long.csv")

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
#write_csv(fires_long,"./Data/Derived/Long_Form/all_fires_Long.csv")

####3.4. Calculate Fragstats####
##3.4a. Get mean shape index (MSI), area-weighted mean shape index (AWMSI), mean patch fractal dimension (MPFD) and area-weighted mean patch fractal dimension (AWMPFD) for all fires, or specific fires.

fires.to.sample <- as.character(fire.list$VB_ID) #Sample all fires

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
#write_csv(fire.list,path="./Data/Derived/all_fires_ForAnalysis.csv")

####4. This section implements plotting of SDC values for actual stand-replacing patches####
####4.1. Load data for analysis####
fire.list <- read_csv("./Data/Derived/all_fires_ForAnalysis.csv")
fires_long <- read_csv("./Data/Derived/Long_Form/all_fires_Long.csv")
hs_patches=readOGR("../Large Files/GIS/BurnSev/Current/", layer="hs_patches")

####4.2.Plot Specific fires####
#Gives a nice comparison of the East and Caribou fires.
fires.to.plot.names=c("1987EAST","2008CARIBOU")
fires.to.plot=fires_long[as.character(fires_long$name)%in%fires.to.plot.names,]
fires.to.plot$name=factor(fires.to.plot$name,labels=c("East","Caribou"))
p.a=ggplot()+
  geom_polygon(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[1],]),
               aes(x=long-min(long),y=lat-min(lat),group=group),col='darkred',fill='darkred')+
  xlim(0,13110) + ylim(0,11280)+ coord_fixed()+
  labs(title="East Fire (1987)",x=" ",y="meters") +
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=18, hjust=0.5))

p.b=ggplot()+
  geom_polygon(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[2],]),
               aes(x=long-min(long),y=lat-min(lat),group=group),fill='darkblue')+
  xlim(0,13110) + ylim(0,11280) + coord_fixed()+
  labs(title="Caribou Fire (2008)",x=" ", y=element_blank()) +
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=18, hjust=0.5))

p.fits=ggplot(data=fires.to.plot,aes(x=width,y=prop.hs))+ 
  geom_point(aes(col=name))+
  scale_color_manual(values=c("darkred","darkblue"),labels=c("East Fire","Caribou Fire")) +
  geom_path(aes(x=width,y=1/(10^(sdc*width)),group=name,col=name))+
  annotate("text", x = c(125,245), y = c(0.25,0.35), label = rev(unique(fires.to.plot$sdc.name)),size=3) +
  labs(y="proportion \n stand-replacing",x="Patch interior buffer distance (m)",title=" ") +
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.position = c(0.9,0.5),legend.justification = c(1, 0.5), legend.title = element_blank())

ggdraw()+
  draw_plot(p.a, x=0, y=.5, width=.5, height=.5) +
  draw_plot(p.b, x=.5, y=.5, width=.5, height=.5) +
  draw_plot(p.fits, x=0, y=0, width=1, height=.5) +
  draw_plot_label(c("A", "B", "C", "meters"), c(0, 0.5, 0, 0.40), c(1, 1, 0.5, 0.55), 
                  size = c(12, 12, 12, 16), fontface=c("bold","bold","bold","plain")) 

#dev.copy2pdf(file=paste0("./Figures/Fig4_",Sys.Date(),".pdf"),width=8,height=8)

####4.3. Summary Analyses of SDC for all fires####
fire.list=fire.list[-which(fire.list$FIRE_YEAR==2016),] #Removing 2016 to be consistent with second paper. Only 6 fires mapped at time of publication. N=477 now.
#6.1 Histograms
p5a=ggplot(fire.list,aes(log(sdc)))+
  geom_histogram(binwidth=0.1,fill="white",col="black")+
  geom_vline(xintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1.2)+
  #labs(x = expression(ln(c[d])))+
  labs(x="ln(SDC)",y="Number of fires",title=" ")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15) )
p5a
#dev.copy2pdf(file=paste0("./Figures/FigS1_",Sys.Date(),".pdf"),width=8,height=8)


p5b=ggplot(fire.list,aes(x=log(FIRESIZE_HA),y=log(sdc)))+ #Optional: add ,col=BA90_PCT
  geom_hline(yintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1)+
  geom_point()+
  geom_smooth(method="lm",col="gray")+
  labs(x ="fire size (ln[ha])\n", y= "ln(SDC)", title = "")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16) )
p5b
#dev.copy2pdf(file=paste0("./Figures/FigS4_",Sys.Date(),".pdf"),width=8,height=6)

p5c=ggplot(fire.list,aes(x=BA90_PCT,y=log(sdc)))+ 
  geom_hline(yintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1)+
  geom_point()+
  geom_smooth(method="lm",col="gray")+
  labs(x ="percentage mapped \n as stand-replacing", y= "ln(SDC)", title = "")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16) )
p5c
#dev.copy2pdf(file=paste0("./Figures/FigS5_",Sys.Date(),".pdf"),width=8,height=6)

ggdraw()+
  draw_plot(p5a, x=0, y=0.6, width=1, height=.4) +
  draw_plot(p5b, x=0, y=0, width=.5, height=.6) +
  draw_plot(p5c, x=0.5, y=0, width=.5, height=.6) +
  draw_plot_label(c("A", "B", "C"), c(0.05, 0.05, 0.55), c(1, 0.6, 0.6), 
                  size = c(12, 12, 12), fontface=c("bold","bold","bold")) 

#dev.copy2pdf(file=paste0("./Figures/Fig5_",Sys.Date(),".pdf"),width=8,height=6)

####4. Summary Analyses of Fragstats for all fires####

#hist(log(fire.list$AWMSI))
fire.list$VB_ID_Special = ifelse (fire.list$VB_ID %in% c("2015CASTLE", "2008VENTURE"),as.character(fire.list$VB_ID),"")
fire.list$VB_ID_SpecialCASTLE = ifelse (fire.list$VB_ID %in% "2015CASTLE",as.character(fire.list$VB_ID),"")
fire.list$VB_ID_SpecialVENTURE = ifelse (fire.list$VB_ID %in% "2008VENTURE",as.character(fire.list$VB_ID),"")
pS2a=ggplot(fire.list,aes(x=log(AWMSI),y=log(sdc)))+ 
  #geom_hline(yintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1)+
  geom_point(aes(col=VB_ID_Special))+
  scale_color_manual(values=c("black","orange","orange"),guide=FALSE)+
  geom_text(aes(label=VB_ID_SpecialCASTLE),hjust=0, vjust=1,size=3)+
  geom_text(aes(label=VB_ID_SpecialVENTURE),hjust=1, vjust=1,size=3)+
  geom_smooth(method="lm",col="gray")+
  labs(x ="ln(AWMSI) ", y= "ln(SDC)", title = "")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))
pS2a 

pS2b=ggplot(fire.list,aes(x=(AWMPFD),y=log(sdc)))+ 
  #geom_hline(yintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1)+
  geom_point()+
  geom_smooth(method="lm",col="gray")+
  labs(x ="AWMPFD ", y= "ln(SDC)", title = "")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16) )
pS2b #Use this one


ggdraw()+
  draw_plot(pS2a, x=0, y=0, width=.5, height=1) +
  draw_plot(pS2b, x=0.5, y=0, width=.5, height=1) +
  draw_plot_label(c("A", "B"), c(0, 0.5), c(1, 1), 
                  size = c(12, 12)) 

#dev.copy2pdf(file=paste0("./Figures/FigS2_",Sys.Date(),".pdf"),width=8,height=8)