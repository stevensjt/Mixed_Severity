
####0. Load libraries####
library(dismo)
library(ggplot2)
require(rgeos)
library(gridExtra)
source("./Analyses/Functions.R")

####1. Rasterizing circles example (skip)####
#http://stackoverflow.com/questions/37829031/creating-raster-with-circles-from-dataframe-in-r

df <- data.frame(lon=c(-0.3,1,1.5,2.7,2.1),
                 lat=c(40.4,42.4,42.4,42.4,42.3), 
                 radius=c(4.4,8.4,11.4,5.4,10.3))
p <- df[,1:2]
rad <- df[,3]*2500
cc <- circles(p, rad, lonlat=TRUE, dissolve=FALSE)
plot(cc)
points(p)
pls <- polygons(cc)
r <- raster(extent(pls), res=.01)
#r=setValues(r,values=1)
pr <- predict(cc, r)
pr2=rasterize(pls,r,field=2,background=1)
#pr3=overlay(pr2,r,fun=max)
plot(pr2)

####2. Fake data####

##2a: many small circles, ~1 ha each, 1000 ha total
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

##2b: fewer medium circles, 10 ha each, 1000 ha total
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

##2c: few large circles, ~100 ha each, 1000 ha total
#100 ha patches = 1,000,000 m2. 10 100ha patches = 1000 ha of high severity, sqrt(10)=3.1622, or a 3x3 grid.
#Patch size and radius: To get exactly 1000 ha, 9 patches (b/c 3x3 grid) * x =1000, x=1000/9 = 111.1111 ha = 1,111,111 m2. Radius to get 1,111,111 m2 is sqrt(1111111/pi) = 594.708 m. 
#Spacing: Going with spacing of (594.708*2)+50 = 1239 m
#Panel c: 1024 patches, patch size = 0.9766 ha, cumulative area = 1000 ha

gridsize=3
spacing=fire.edge.length/gridsize
radius=594.708
df=expand.grid(lon=seq(from=spacing/2,by=spacing,length.out=gridsize),
               lat=seq(from=spacing/2,by=spacing,length.out=gridsize))
df$radius=radius
p <- df[,1:2]
rad <- df[,3]
circles.c <- polygons(circles(p, rad, lonlat=FALSE, dissolve=FALSE))

##2d: one big circle, 1000 ha in total
#1000 ha is 1000/0.0001 = 10,000,000 m2. Radius to get that is sqrt(10000000/pi) = 1784.12
#Panel d: 1 patch, patch size = 1000 ha, cumulative area = 1000 ha

radius=1784.124
df <- data.frame(lon=fire.edge.length/2, lat=fire.edge.length/2, radius=radius) #lat long adjustment to match other plots, e.g. (162/2)-55.75388
p <- df[,1:2]
rad <- df[,3]
circles.d <- polygons(circles(p, rad, lonlat=FALSE, dissolve=FALSE))

#2e: plot to check
#plot(circles.1big)
#spplot(circles.1big)
ggplot()+
  geom_path(data=circles.a,aes(x=long,y=lat,group=group),col="black")+
  geom_path(data=circles.b,aes(x=long,y=lat,group=group),col="black")+
  geom_path(data=circles.c,aes(x=long,y=lat,group=group),col="black")+
  geom_path(data=circles.d,aes(x=long,y=lat,group=group),col="black")+
  theme_bw()

####3. Testing R model from ground up with fake data (skip)####
#Circle with 100 m radius; R is radius
tmp=data.frame(width=seq(0,100,by=10))
tmp$prop.hs=(1-tmp$width/100)^2
#tmp now describes a 31415.93 m2 (or 3.14159 ha) circle
tmp$prop.hs2=c(1,tmp$prop.hs[2:10]+rnorm(tmp$prop.hs[2:10],sd=0.01),0)
nls(prop.hs2~(1-width/q)^2,data=tmp,start=list(q=1))
#This only works when there is some residual error in the data (i.e. prop.hs2); it will not work on perfect circles for which there is an analytical solution

####4. Create decay profile####
#Initialize
buf_inc= 10 #Buffer increment in m
buf_max= 1000 #Max buffer width is 1000m; needs to be a multiple of buf_inc
table.a=decay(hs_fire=circles.a,buf_max=1000,buf_inc=10,name="1ha patches")
table.b=decay(hs_fire=circles.b,buf_max=1000,buf_inc=10,name="10ha patches")
table.c=decay(hs_fire=circles.c,buf_max=1000,buf_inc=10,name="100ha patches")
table.d=decay(hs_fire=circles.d,buf_max=1000,buf_inc=10,name="1000ha patches")
all_fires=rbind(table.a,table.b,table.c,table.d)
all_fires$prop.hs.jitter=all_fires$prop.hs
all_fires$prop.hs.jitter[!all_fires$prop.hs%in%c(0,1)]=
  abs(all_fires$prop.hs[!all_fires$prop.hs%in%c(0,1)]+
  rnorm(all_fires$prop.hs[!all_fires$prop.hs%in%c(0,1)],sd=0.01))

####5. Calculate q####
all_fires.R=calculate.R(all_fires)

####6. Plot EDA####
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

p.fits=ggplot(data=all_fires.R,aes(x=width,y=prop.hs))+
  geom_point(aes(col=name),size=1)+
  scale_color_manual(values=plot.colors)+
  geom_path(aes(x=width,y=(1-width/R)^2,group=name,col=name))+
  annotate("text", x = c(0,170,350,950), y = c(0.1,0.1,0.1,0.1), label = unique(all_fires.R$R.name),size=3)+
  labs(y="proportion high-severity",x="distance to edge (m)",col="patch size")+
  theme_bw()
p.fits

Fig1=
  grid.arrange(arrangeGrob(p.a,p.b,p.c,p.d,ncol=2),p.fits,ncol=1,heights=c(2,1),clip=F)
dev.copy2pdf(file=paste0("./Figures/Fig1_",Sys.Date(),".pdf"),width=8,height=12)
             