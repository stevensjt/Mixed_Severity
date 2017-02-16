
####1. Rasterizing circles example (don't run)####
library(dismo)
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
##2a. 1 big circle; deprecated and repeated below
#http://stackoverflow.com/questions/37829031/creating-raster-with-circles-from-dataframe-in-r
#1000 ha fire = 10 million square meters 10000000; or 10 square km.
#sqrt(10000000)=3162.278 m per side; round to 3150 m per side
#3150/30=105 pixels per side
#The high-severity patch in that fire is (pi*(3150/2)^2)*0.0001 = 779 ha
#e=extent(0,3150,0,3150) #Deprecated unless I start working with rasters.

#df <- data.frame(lon=c(3150/2),lat=c(3150/2), radius=c(3150/2))
#p <- df[,1:2]
#rad <- df[,3]
#circles.1big <- polygons(circles(p, rad, lonlat=FALSE, dissolve=FALSE))

##2b: many small circles, 1000 ha in total
#900 m2 (1 pixel) is 0.9 ha, 1000/0.9 = 1111 individual circles with radius sqrt(900/pi) = 16.93 m; spacing of centroids should be > 17*2 = 34 m. sqrt(1111) = 33 x 33 grid (n=1089)
#To get exactly 1000 ha, 1089 (b/c 33x33 grid) * x =1000, x=1000/1089 = 0.9182736 ha = 9182.736 m2; radius to get 9182.736 m2 is sqrt(9182.736/pi) = 54.06 m. Going with spacing of (54.06*2)+50 = 158 just to visualize more easily

df=expand.grid(lon=seq(from=1,by=158,length.out=33),
               lat=seq(from=1,by=158,length.out=33))
df$radius=54.06
p <- df[,1:2]
rad <- df[,3]
circles.smallest <- polygons(circles(p, rad, lonlat=FALSE, dissolve=FALSE))

##2c: one big circle, 1000 ha in total
#1000 ha is 1000/0.0001 = 10,000,000 m2. Radius to get that is sqrt(10000000/pi) = 1784.12

df <- data.frame(lon=2529, lat=2529, radius=1784.12)
p <- df[,1:2]
rad <- df[,3]
circles.1big <- polygons(circles(p, rad, lonlat=FALSE, dissolve=FALSE))

#2d: plot to check
#plot(circles.1big)
#spplot(circles.1big)
ggplot()+
  geom_polygon(data=circles.1big,aes(x=long,y=lat),fill="transparent",col="black")+
  geom_polygon(data=circles.smallest,aes(x=long,y=lat),fill="transparent",col="black")+
  theme_bw()

circles.smallest.2=aggregate(circles.smallest)
####4. Create decay profile####
#Initialize
buf_inc= 10 #Buffer increment in m
buf_max= 1000 #Max buffer width is 1000m; needs to be a multiple of buf_inc
hs_fire=big1
big1.table=decay(hs_fire=hs_fire,buf_max=1000,buf_inc=10,name="big1")
