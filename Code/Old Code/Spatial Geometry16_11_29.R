
####0. Load libraries####
library(dismo)
library(ggplot2)
require(rgeos)
library(gridExtra)
source("./Analyses/Functions.R")


####1. Generate fake data####

##1a: many small circles, ~1 ha each, 1000 ha total
#1 ha patches = 10,000 m2. 1000 1ha patches = 1000 ha of high severity, sqrt(1000)=31.62, or a 32x32 grid.
#Patch size and radius: To get exactly 1000 ha, 1024 patches (b/c 32x32 grid) * x =1000, x=1000/1024 = 0.9765625 ha = 9765.625 m2. Radius to get 9765.625 m2 is sqrt(9765.625/pi) = 55.75388 m. 
#Spacing: Going with spacing of (55.75388*2)+50 = 162, just to visualize more easily
#Panel a: 1024 patches, patch size = 0.9766 ha, cumulative area = 1000 ha

fire.edge.length=6000 #Side length of a square containing the HS patches
gridsize=32
spacing=fire.edge.length/gridsize
radius=55.75388
df=expand.grid(lon=seq(from=spacing/2,by=spacing,length.out=gridsize),
               lat=seq(from=spacing/2,by=spacing,length.out=gridsize))
df$radius=radius
p <- df[,1:2]
rad <- df[,3]
circles.a <- polygons(circles(p, rad, lonlat=FALSE, dissolve=FALSE))

##1b: fewer medium circles, 10 ha each, 1000 ha total
#10 ha patches = 100,000 m2. 100 10 ha patches= 1000 ha of high severity, 10 x 10 grid.
#Patch size and radius: Patch size = 10 ha = 100000 m2. Radius to get 100,000 m2 is sqrt(100000/pi) = 178.4124 m
#Spacing of centroids should be (178.4124*2)+50 = 407
#Panel b: 100 patches, patch size = 10 ha, cumulative area = 1000 ha

gridsize=10
spacing=fire.edge.length/gridsize
radius=178.4124
df=expand.grid(lon=seq(from=spacing/2,by=spacing,length.out=gridsize),
               lat=seq(from=spacing/2,by=spacing,length.out=gridsize))
df$radius=radius
p <- df[,1:2]
rad <- df[,3]
circles.b <- polygons(circles(p, rad, lonlat=FALSE, dissolve=FALSE))

##1c: few large circles, ~100 ha each, 1000 ha total
#100 ha patches = 1,000,000 m2. 10 100ha patches = 1000 ha of high severity, sqrt(10)=3.1622, or a 3x3 grid.
#Patch size and radius: To get exactly 1000 ha, 9 patches (b/c 3x3 grid) * x =1000, x=1000/9 = 111.1111 ha = 1,111,111 m2. Radius to get 1,111,111 m2 is sqrt(1111111/pi) = 594.708 m. 
#Spacing: Going with spacing of (594.708*2)+50 = 1239 m
#Panel c: 9 patches, patch size = 111 ha, cumulative area = 1000 ha

gridsize=3
spacing=fire.edge.length/gridsize
radius=594.708
df=expand.grid(lon=seq(from=spacing/2,by=spacing,length.out=gridsize),
               lat=seq(from=spacing/2,by=spacing,length.out=gridsize))
df$radius=radius
p <- df[,1:2]
rad <- df[,3]
circles.c <- polygons(circles(p, rad, lonlat=FALSE, dissolve=FALSE))

##1d: one big circle, 1000 ha in total
#1000 ha is 1000/0.0001 = 10,000,000 m2. Radius to get that is sqrt(10000000/pi) = 1784.12
#Panel d: 1 patch, patch size = 1000 ha, cumulative area = 1000 ha

radius=1784.124
df <- data.frame(lon=fire.edge.length/2, lat=fire.edge.length/2, radius=radius) #lat long adjustment to match other plots, e.g. (162/2)-55.75388
p <- df[,1:2]
rad <- df[,3]
circles.d <- polygons(circles(p, rad, lonlat=FALSE, dissolve=FALSE))

#1e: plot to check
#plot(circles.1big)
#spplot(circles.1big)
ggplot()+
  geom_path(data=circles.a,aes(x=long,y=lat,group=group),col="black")+
  geom_path(data=circles.b,aes(x=long,y=lat,group=group),col="black")+
  geom_path(data=circles.c,aes(x=long,y=lat,group=group),col="black")+
  geom_path(data=circles.d,aes(x=long,y=lat,group=group),col="black")+
  theme_bw()

####2. Create decay profile for each set of circles####
#Initialize
buf_inc= 10 #Buffer increment in m
buf_max= 1000 #Max buffer width is 1000m; needs to be a multiple of buf_inc
#Run the decay function on each "fire"
table.a=decay(hs_fire=circles.a,buf_max=1000,buf_inc=10,name="1ha patches")
table.b=decay(hs_fire=circles.b,buf_max=1000,buf_inc=10,name="10ha patches")
table.c=decay(hs_fire=circles.c,buf_max=1000,buf_inc=10,name="100ha patches")
table.d=decay(hs_fire=circles.d,buf_max=1000,buf_inc=10,name="1000ha patches")
all_fires=rbind(table.a,table.b,table.c,table.d)
#Add a slight bit of variation so that the model can fit the R parameter
all_fires$prop.hs.jitter=all_fires$prop.hs
all_fires$prop.hs.jitter[!all_fires$prop.hs%in%c(0,1)]=
  abs(all_fires$prop.hs[!all_fires$prop.hs%in%c(0,1)]+
  rnorm(all_fires$prop.hs[!all_fires$prop.hs%in%c(0,1)],sd=0.01))

####3. Calculate characteristic shape R####
#Statistical solution
all_fires.q=calculate.q2(all_fires)
#Analytical solution (deprecated)
#all_fires.R[grep("1ha",as.character(all_fires.R$name)),"R.analytical"]=round(sqrt((1*10000)/pi),1)
#all_fires.R[grep("10ha",as.character(all_fires.R$name)),"R.analytical"]=round(sqrt((10*10000)/pi),1)
#all_fires.R[grep("100ha",as.character(all_fires.R$name)),"R.analytical"]=round(sqrt((100*10000)/pi),1)
#all_fires.R[grep("1000ha",as.character(all_fires.R$name)),"R.analytical"]=round(sqrt((1000*10000)/pi),1)

####4. Generate plot####
plot.colors=c("#377eb8","#4daf4a","#984ea3","#e41a1c")
p.a=ggplot()+
  geom_path(data=circles.a,aes(x=long,y=lat,group=group),col=plot.colors[1])+
  #lims(x=c(0,5500),y=c(0,5500))+
  xlim(0,fire.edge.length)+ylim(0,fire.edge.length)+
  labs(x="meters",y="meters")+
  theme_bw()
p.b=ggplot()+
  geom_path(data=circles.b,aes(x=long,y=lat,group=group),col=plot.colors[2])+
  xlim(0,fire.edge.length)+ylim(0,fire.edge.length)+
  labs(x="meters",y="meters")+
  theme_bw()
p.c=ggplot()+
  geom_path(data=circles.c,aes(x=long,y=lat,group=group),col=plot.colors[3])+
  xlim(0,fire.edge.length)+ylim(0,fire.edge.length)+
  labs(x="meters",y="meters")+
  theme_bw()
p.d=ggplot()+
  geom_path(data=circles.d,aes(x=long,y=lat,group=group),col=plot.colors[4])+
  xlim(0,fire.edge.length)+ylim(0,fire.edge.length)+
  labs(x="meters",y="meters")+
  theme_bw()

p.fits=ggplot(data=all_fires.q,aes(x=width,y=prop.hs))+
  geom_point(aes(col=name),size=1)+
  scale_color_manual(values=plot.colors)+
  geom_path(aes(x=width,y=1/(10^(q*width)),group=name,col=name))+
  annotate("text", x = c(0,170,350,950), y = c(0.1,0.1,0.1,0.1), label = unique(all_fires.q$q.name),size=3)+
  labs(y="proportion high-severity",x="distance to edge (m)",col="patch size")+
  theme_bw()
#p.fits

Fig1=
  grid.arrange(arrangeGrob(p.a,p.b,p.c,p.d,ncol=2),p.fits,ncol=1,heights=c(2,1),clip=F)
#dev.copy2pdf(file=paste0("./Figures/Fig1_",Sys.Date(),".pdf"),width=8,height=12)

####S1. Supplementary Shape Figure Create Shapes####
##1d: one big circle, 1000 ha in total
#1000 ha is 1000/0.0001 = 10,000,000 m2. Radius to get that is sqrt(10000000/pi) = 1784.12
#Panel d: 1 patch, patch size = 1000 ha, cumulative area = 1000 ha

radius=1784.124
df <- data.frame(lon=fire.edge.length/2, lat=fire.edge.length/2, radius=radius) #lat long adjustment to match other plots, e.g. (162/2)-55.75388
p <- df[,1:2]
rad <- df[,3]
circles.d <- polygons(circles(p, rad, lonlat=FALSE, dissolve=FALSE))
target.area=area(circles.d) #9,999,491 sq m, or ~10 million m2, or ~1000 ha

####S1b. An ellipse with area = 1000 ha####
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

#e.spatialpoly is now set for comparison

####S1.c Make an irregular ellipse####
f.spatialpoly=e.spatialpoly
tmp=as.data.frame(f.spatialpoly@polygons[[1]]@Polygons$ellipse1@coords)
names(tmp)=c("x","y")
tmp$zx=NA; tmp$zy=NA

#Create new x-y points either inward or outward from the edge
#http://math.stackexchange.com/questions/9365/endpoint-of-a-line-knowing-slope-start-and-distance
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

####S1d. Make an irregular circle####
g.spatialpoly=circles.d
tmp=as.data.frame(g.spatialpoly@polygons[[1]]@Polygons[[1]]@coords)
names(tmp)=c("x","y")
tmp$zx=NA; tmp$zy=NA
tmp$x=jitter(tmp$x,0.01); tmp$y=jitter(tmp$y,0.01) #Need to add a little noise to the circle data so that you don't end up with issues calculating the slope of vertical lines.

#Create new x-y points either inward or outward from the edge
d=350 #Initialize distance in/out from edge along perpendicular line
for(i in c(seq(from=2,to=nrow(tmp),by=10),seq(from=3,to=nrow(tmp),by=10))){
  m=(tmp[i+1,"y"]-tmp[i-1,"y"])/(tmp[i+1,"x"]-tmp[i-1,"x"])#Calculate slope of two adjacent points
  m2=-(tmp[i+1,"x"]-tmp[i-1,"x"])/(tmp[i+1,"y"]-tmp[i-1,"y"])#Calculate slope of the perpendicular line
  k=rnorm(1,mean=d,sd=100)/sqrt(1+m2^2) #Calculate constant k
  tmp[i,"zx"]=tmp[i,"x"]+k
  tmp[i,"zy"]=tmp[i,"y"]+(k*m2)
}#End addition loop
for(i in c(seq(from=7,to=nrow(tmp),by=10),seq(from=8,to=nrow(tmp),by=10))){
  m=(tmp[i+1,"y"]-tmp[i-1,"y"])/(tmp[i+1,"x"]-tmp[i-1,"x"])#Calculate slope of two adjacent points
  m2=-(tmp[i+1,"x"]-tmp[i-1,"x"])/(tmp[i+1,"y"]-tmp[i-1,"y"])#Calculate slope of the perpendicular line
  k=rnorm(1,mean=d,sd=100)/sqrt(1+m2^2) #Calculate constant k
  tmp[i,"zx"]=tmp[i,"x"]-k
  tmp[i,"zy"]=tmp[i,"y"]-(k*m2)
}#End subtraction loop
tmp[nrow(tmp),]=tmp[1,] #Set last line of data equal to first in order to "close the loop"; need to do this because of the jittering done earlier which created slight inequalities in the coordinates of the origin point.
tmp[!is.na(tmp$zx),"x"]=tmp[!is.na(tmp$zx),"zx"]
tmp[!is.na(tmp$zy),"y"]=tmp[!is.na(tmp$zx),"zy"]
tmp.mat=as.matrix(tmp[,-c(3:4)])

g.spatialpoly@polygons[[1]]@Polygons[[1]]@coords=tmp.mat
area(g.spatialpoly)/10000 #Should match and be ~1000 ha (999.9495)
plot(g.spatialpoly) #Check and see that it looks good and irregular

####S2. Create decay profile for each set of circles####
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

####3. Calculate characteristic shape R####
#Statistical solution
all_fires.q.supp=calculate.q2(all_fires)

####4. Generate plot####
plot.colors=c("#e41a1c","#ff7f00","#a65628") #,"#df65b0" #Optional color for the irregular circle.
p.d=ggplot()+
  geom_path(data=circles.d,aes(x=long,y=lat,group=group),col=plot.colors[1])+
  #lims(x=c(0,5500),y=c(0,5500))+
  xlim(0,fire.edge.length)+ylim(0,fire.edge.length)+
  labs(x="meters",y="meters")+
  theme_bw()
p.e=ggplot()+
  geom_path(data=e.spatialpoly,aes(x=long,y=lat,group=group),col=plot.colors[2])+
  xlim(0,fire.edge.length)+ylim(0,fire.edge.length)+
  labs(x="meters",y="meters")+
  theme_bw()
p.f=ggplot()+
  geom_path(data=f.spatialpoly,aes(x=long,y=lat,group=group),col=plot.colors[3])+
  xlim(0,fire.edge.length)+ylim(0,fire.edge.length)+
  labs(x="meters",y="meters")+
  theme_bw()
#p.g=ggplot()+
#  geom_path(data=g.spatialpoly,aes(x=long,y=lat,group=group),col=plot.colors[4])+
#  xlim(0,fire.edge.length)+ylim(0,fire.edge.length)+
#  labs(x="meters",y="meters")+
#  theme_bw()

p.fits=ggplot(data=all_fires.q.supp,aes(x=width,y=prop.hs))+
  geom_point(aes(col=name),size=1)+
  scale_color_manual(values=plot.colors)+
  geom_path(aes(x=width,y=1/(10^(q*width)),group=name,col=name))+
  annotate("text", x = c(950,700,400), y = c(0.35,0.1,0.1), label = unique(all_fires.q.supp$q.name),size=3)+
  labs(y="proportion high-severity",x="distance to edge (m)",col="patch shape")+
  theme_bw()
p.fits
p.empty=ggplot(data=all_fires.q.supp,aes(x=width,y=prop.hs))+geom_blank()+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        panel.background = element_blank())

FigS6=
  grid.arrange(arrangeGrob(p.d,p.e,p.empty,p.f,ncol=2),p.fits,ncol=1,heights=c(2,1),clip=F)

dev.copy2pdf(file=paste0("./Figures/FigS6_",Sys.Date(),".pdf"),width=8,height=12)

