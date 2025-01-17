##Code for Concepts paper
####0. load libraries####
library(rgdal) #readOGR
library(cowplot)
library(tidyverse)
library(rgeos) #for gBuffer in section 5; v 0.3-26
library(raster) #for extent() in section 5; v 2.6-7
#library(profvis) #See MCM talk DRUG 1/13/18, function to rip apart code chunks and see time reqs. 
#Also profile -> profile selected lines
#Also library(microbenchmark) and library(rio) for faster imports (but not for shapefiles)

source("./Code/Functions.R")

####1a. Load data for analysis####
#tmp<-profvis({
fire.list <- read_csv("./Data/Derived/all_fires_ForAnalysis.csv")
fires_long <- read_csv("./Data/Derived/Long_Form/all_fires_Long.csv")
hs_patches=readOGR("../Large Files/GIS/BurnSev/Current/", layer="hs_patches")
#})

####2.Plot Specific fires####
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

dev.copy2pdf(file=paste0("./Figures/Fig4_",Sys.Date(),".pdf"),width=8,height=8)

####3. Summary Analyses of SDC for all fires####
fire.list=fire.list[-which(fire.list$FIRE_YEAR==2016),] #Removing 2016 to be consistent with second paper. Only 6 fires mapped at time of publication. N=477 now.
#6.1 Histograms
h=ggplot(fire.list,aes(log(sdc)))+
  geom_histogram(binwidth=0.1,fill="white",col="black")+
  geom_vline(xintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1.2)+
  #labs(x = expression(ln(c[d])))+
  labs(x="ln(SDC)",y="Number of fires",title=" ")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15) )
h
#dev.copy2pdf(file=paste0("./Figures/FigS1_",Sys.Date(),".pdf"),width=8,height=8)

pS2=ggplot(fire.list[fire.list$FIRESIZE_HA>1000,],aes(x=log(FIRESIZE_HA),y=log(sdc)))+
  geom_point()+
  geom_smooth(method="lm",col="black")+
  geom_hline(yintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1)+
  labs(title="Fires > 1000 ha")+
  theme_bw()
pS2
#dev.copy2pdf(file=paste0("./Figures/FigS2_",Sys.Date(),".pdf"),width=8,height=6)

pS3=ggplot(fire.list[fire.list$FIRESIZE_HA<1000,],aes(x=log(FIRESIZE_HA),y=log(sdc)))+
  geom_point()+
  geom_smooth(method="lm",col="black")+
  geom_hline(yintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1)+
  labs(title="Fires < 1000 ha")+
  theme_bw()
pS3
#dev.copy2pdf(file=paste0("./Figures/FigS3_",Sys.Date(),".pdf"),width=8,height=6)

pS4=ggplot(fire.list,aes(x=log(FIRESIZE_HA),y=log(sdc)))+ #Optional: add ,col=BA90_PCT
  geom_hline(yintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1)+
  geom_point()+
  geom_smooth(method="lm",col="gray")+
  labs(x ="fire size (ln[ha])\n", y= "ln(SDC)", title = "")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16) )
pS4
#dev.copy2pdf(file=paste0("./Figures/FigS4_",Sys.Date(),".pdf"),width=8,height=6)

pS5=ggplot(fire.list,aes(x=BA90_PCT,y=log(sdc)))+ 
  geom_hline(yintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1)+
  geom_point()+
  geom_smooth(method="lm",col="gray")+
  labs(x ="percentage mapped \n as stand-replacing", y= "ln(SDC)", title = "")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16) )
pS5
#dev.copy2pdf(file=paste0("./Figures/FigS5_",Sys.Date(),".pdf"),width=8,height=6)

pS6=ggplot(fire.list,aes(x=log(BA90_HA),y=log(sdc),col=FIRE_YEAR))+ 
  geom_hline(yintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1)+
  geom_point()+
  geom_smooth(method="lm",col="gray")+
  labs(x ="high-severity area (ln[ha])\n", y= "ln(SDC)")+
  theme_bw()
pS6
#dev.copy2pdf(file=paste0("./Figures/FigS6_",Sys.Date(),".pdf"),width=8,height=6)

ggdraw()+
  draw_plot(h, x=0, y=0.6, width=1, height=.4) +
  draw_plot(pS4, x=0, y=0, width=.5, height=.6) +
  draw_plot(pS5, x=0.5, y=0, width=.5, height=.6) +
  draw_plot_label(c("A", "B", "C"), c(0.05, 0.05, 0.55), c(1, 0.6, 0.6), 
                  size = c(12, 12, 12), fontface=c("bold","bold","bold")) 

#dev.copy2pdf(file=paste0("./Figures/Fig5_",Sys.Date(),".pdf"),width=8,height=6)

####4. Summary Analyses of Fragstats for all fires####

#hist(log(fire.list$AWMSI))
fire.list$VB_ID_Special = ifelse (fire.list$VB_ID %in% c("2015CASTLE","2008VENTURE"),as.character(fire.list$VB_ID),"")
pS7=ggplot(fire.list,aes(x=log(AWMSI),y=log(sdc)))+ 
  #geom_hline(yintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1)+
  geom_point()+
  geom_text(aes(label=VB_ID_Special),hjust=1, vjust=1,size=3)+
  geom_smooth(method="lm",col="gray")+
  labs(x ="ln(AWMSI) ", y= "ln(SDC)", title = "")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16) )
pS7 #Use this one

pS8=ggplot(fire.list,aes(x=(MSI),y=log(sdc)))+ 
  #geom_hline(yintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1)+
  geom_point()+
  geom_smooth(method="lm",col="gray")+
  labs(x ="MSI ", y= expression(ln(c[d])))+
  theme_bw()
pS8

#hist(fire.list$AWMPFD)
pS9=ggplot(fire.list,aes(x=(AWMPFD),y=log(sdc)))+ 
  #geom_hline(yintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1)+
  geom_point()+
  geom_smooth(method="lm",col="gray")+
  labs(x ="AWMPFD ", y= "ln(SDC)")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16) )
pS9 #Use this one

pS10=ggplot(fire.list,aes(x=(MPFD),y=log(sdc)))+ 
  #geom_hline(yintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1)+
  geom_point()+
  geom_smooth(method="lm",col="gray")+
  labs(x ="MPFD ", y= expression(ln(c[d])))+
  theme_bw()
pS10

#plot(AWMPFD~log(AWMSI), data=fire.list)
ggdraw()+
  draw_plot(pS7, x=0, y=0, width=.5, height=1) +
  draw_plot(pS9, x=0.5, y=0, width=.5, height=1) +
  draw_plot_label(c("A", "B"), c(0, 0.5), c(1, 1), 
                  size = c(12, 12)) 

#dev.copy2pdf(file=paste0("./Figures/FigS1_",Sys.Date(),".pdf"),width=8,height=8)

####5. Example for post-fire GTR####
#Install ggpolypath to deal with issue of holes:
#https://github.com/tidyverse/ggplot2/issues/752
#https://cran.r-project.org/web/packages/ggpolypath/index.html
library(ggpolypath)

##START HERE; swap out Caribou fire for King fire
fires.to.plot.names=c("2015ROUGH","2014KING")
fires.to.plot=fires_long[as.character(fires_long$name)%in%fires.to.plot.names,]
fires.to.plot$name=factor(fires.to.plot$name,labels=c("King", "Rough"))
full_shapes <- readOGR("../Large Files/GIS/BurnSev/Specific Fires/", layer="KingRoughFullPerimsGTR") #CRS EPSG:3310, NAD83 CA Albers

##PLots
p.a_hs <- fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[1],])
p.a_full <- full_shapes[full_shapes$VB_ID==fires.to.plot.names[1],]
p.a_120 <- gBuffer(
  createSPComment(fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[1],])), 
  #createSPComment was done in the analysis used to create the long dataset (from Collins et al), 
  #but needs to be re-done here for plotting a gBuffered layer
  width = -120)
p.b_hs <- fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[2],])
p.b_full <- full_shapes[full_shapes$VB_ID==fires.to.plot.names[2],]
p.b_120 <- gBuffer(
  createSPComment(fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[2],])), 
  width = -120)

p.a=ggplot()+
  geom_polygon(data=p.a_full,
               aes(x=(long-min(extent(p.a_full)@xmin))/1000,
                   y=(lat-min(extent(p.a_full)@ymin))/1000,group=group),
               fill="gray",col="gray")+
  geom_polypath(data=p.a_hs,
               aes(x=(long-min(extent(p.a_full)@xmin))/1000,
                   y=(lat-min(extent(p.a_full)@ymin))/1000,group=group),
               col='darkred',fill='darkred',size = 0.1)+
  geom_polypath(data=p.a_120,
                aes(x=(long-min(extent(p.a_full)@xmin))/1000,
                        y=(lat-min(extent(p.a_full)@ymin))/1000,group=group),
                    col='aquamarine',fill='aquamarine', size = 0.1)+
  xlim(0,48) + ylim(0,42)+ 
  coord_fixed()+
  labs(title="Rough Fire (2015)",x="km",y="km") +
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=18, hjust=0.5))

p.b=ggplot()+
  geom_polygon(data=p.b_full,
               aes(x=(long-min(extent(p.b_full)@xmin))/1000,
                   y=(lat-min(extent(p.b_full)@ymin))/1000,group=group),
               fill="gray",col="gray")+
  geom_polypath(data=p.b_hs,
                aes(x=(long-min(extent(p.b_full)@xmin))/1000,
                    y=(lat-min(extent(p.b_full)@ymin))/1000,group=group),
                col='darkblue',fill='darkblue', size = 0.1)+
  geom_polypath(data=p.b_120,
                aes(x=(long-min(extent(p.b_full)@xmin))/1000,
                    y=(lat-min(extent(p.b_full)@ymin))/1000,group=group),
                col='aquamarine',fill='aquamarine', size = 0.1)+
  xlim(0,48) + ylim(0,42)+ 
  coord_fixed()+
  labs(title="King Fire (2014)",x="km",y="km") +
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=18, hjust=0.5))

p.fits=ggplot(data=fires.to.plot,aes(x=width,y=prop.hs))+ 
  geom_point(aes(col=name))+
  geom_point(data = fires.to.plot[fires.to.plot$width == 120,], aes(x=width, y = prop.hs),
             col = "aquamarine", size = 5)+
  scale_color_manual(values=c("darkblue","darkred"),labels=c("King Fire","Rough Fire")) +
  geom_path(aes(x=width,y=1/(10^(sdc*width)),group=name,col=name))+
  annotate("text", x = c(500,750), y = c(0.10,0.35), 
           label = paste("sdc = ",rev(unique(fires.to.plot$sdc.name)),
                         "\nln(sdc) = ",round(log(rev(unique(fires.to.plot$sdc.name))),3)),
           size=2.5) +
  labs(y="proportion \n stand-replacing",x="Patch interior \nbuffer distance (m)",title=" ") +
  theme_bw()+
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=16),
        legend.position = c(0.9,0.7),legend.justification = c(1, 0.5), 
        legend.title = element_blank())

h=ggplot(fire.list,aes(log(sdc)))+
  geom_histogram(binwidth=0.1,fill="white",col="black")+
  geom_vline(xintercept=log(rev(unique(fires.to.plot$sdc.name))),col=c("darkred","darkblue"),size=1.2)+
  labs(x="ln(SDC)",y="Number of fires",title=" ")+
  xlim(-7.5,-3)+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15) )

pdf("/Users/Jens/Documents/Berkeley/Post-Doc UCB/Research Projects/Side Projects/Postfire GTR/SDC_fig.pdf")
ggdraw()+
  draw_plot(p.a, x=0, y=.5, width=.5, height=.5) +
  draw_plot(p.b, x=.5, y=.5, width=.5, height=.5) +
  draw_plot(p.fits, x=0, y=0, width=.5, height=.5) +
  draw_plot(h, x=.5, y=0, width=.5, height=.5) +
  draw_plot_label(c("A", "B", "C", "D"), c(0, 0.5, 0, 0.50), c(1, 1, 0.5, 0.5), 
                  size = 12, fontface="bold") 

dev.off()
#sum(area(full_shapes)[c(1,2,3,5,7,9,11)])*0.0001 #King Area
#sum(area(full_shapes)[c(4,6,8,10,12,13,14)])*0.0001 #Rough Area
#library(stats)
quantile(log(fire.list$sdc),0.19)

