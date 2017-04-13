##Code for Concepts paper
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

####1a. Load data for analysis####
fire.list <- read.csv("./Analyses/Derived Data/all_fires_ForAnalysis.csv")
fires_full <- read.csv("./Analyses/Processed Data/all_fires_Full.csv")

####2.Plot Specific fires####
#Gives a nice comparison of the East and Caribou fires.

fires.to.plot.names=c("1987EAST","2008CARIBOU")
fires.to.plot=fires_full[as.character(fires_full$name)%in%fires.to.plot.names,]
fires.to.plot$name=factor(fires.to.plot$name,labels=c("East","Caribou"))
p.a=ggplot()+
  geom_polygon(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[1],]),
            aes(x=long-min(long),y=lat-min(lat),group=group),col='darkred',fill='darkred')+
  xlim(0,13110) + ylim(0,11280)+ coord_fixed()+
  labs(title="a              East Fire (1987)",x="meters",y="meters") +
  theme_bw()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title=element_text(hjust=0))
p.b=ggplot()+
  geom_polygon(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[2],]),
            aes(x=long-min(long),y=lat-min(lat),group=group),fill='darkblue')+
  xlim(0,13110) + ylim(0,11280) + coord_fixed()+
  labs(title="b          Caribou Fire (2008)",x="meters",y="meters") +
  theme_bw()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title=element_text(hjust=0))

p.fits=ggplot(data=fires.to.plot,aes(x=width,y=prop.hs))+ 
  geom_point(aes(col=name))+
  scale_color_manual(values=c("darkred","darkblue")) +
  geom_path(aes(x=width,y=1/(10^(q*width)),group=name,col=name))+
  annotate("text", x = c(50,500), y = c(0.25,0.25), label = rev(unique(fires.to.plot$q.name)),size=3) +
  labs(y="proportion \n stand-replacing",x="distance to edge (m)",col="fire",title="c") +
  theme_bw()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title=element_text(hjust=0))

Fig4=
  grid.arrange(arrangeGrob(p.a,p.b,ncol=2),p.fits,ncol=1,heights=c(2,2),clip=F)
dev.copy2pdf(file=paste0("./Figures/Fig4_",Sys.Date(),".pdf"),width=8,height=8)

####3. Summary Analyses####
#Read in processed data
fire.list2=read.csv("./Data/Derived Data/all_fires_ForAnalysis.csv")
#fire.list2=fire.list[,c(2:18,50)]#Remove irrelevant columns (if re-doing the whole buffering analysis)
fire.list2=fire.list #Deprecate and just use fire.list? Not sure why we need a fire.list2

#6.1 Histograms
h=ggplot(fire.list2,aes(log(sdc)))+
  geom_histogram(binwidth=0.1,fill="white",col="black")+
  geom_vline(xintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1.2)+
  #labs(x = expression(ln(c[d])))+
  labs(x="ln(SDC)")+
  theme_bw()
h
#dev.copy2pdf(file=paste0("./Figures/FigS1_",Sys.Date(),".pdf"),width=8,height=8)

pS2=ggplot(fire.list2[fire.list2$FIRESIZE_HA>1000,],aes(x=log(FIRESIZE_HA),y=log(sdc)))+
  geom_point()+
  geom_smooth(method="lm",col="black")+
  geom_hline(yintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1)+
  labs(title="Fires > 1000 ha")+
  theme_bw()
pS2
#dev.copy2pdf(file=paste0("./Figures/FigS2_",Sys.Date(),".pdf"),width=8,height=6)

pS3=ggplot(fire.list2[fire.list2$FIRESIZE_HA<1000,],aes(x=log(FIRESIZE_HA),y=log(sdc)))+
  geom_point()+
  geom_smooth(method="lm",col="black")+
  geom_hline(yintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1)+
  labs(title="Fires < 1000 ha")+
  theme_bw()
pS3
#dev.copy2pdf(file=paste0("./Figures/FigS3_",Sys.Date(),".pdf"),width=8,height=6)

pS4=ggplot(fire.list2,aes(x=log(FIRESIZE_HA),y=log(sdc)))+ #Optional: add ,col=BA90_PCT
  geom_hline(yintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1)+
  geom_point()+
  geom_smooth(method="lm",col="gray")+
  labs(x ="fire size (ln[ha])\n", y= "ln(SDC)")+
  theme_bw()
pS4
#dev.copy2pdf(file=paste0("./Figures/FigS4_",Sys.Date(),".pdf"),width=8,height=6)

pS5=ggplot(fire.list2,aes(x=BA90_PCT,y=log(sdc)))+ 
  geom_hline(yintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1)+
  geom_point()+
  geom_smooth(method="lm",col="gray")+
  labs(x ="percentage mapped \n as high-severity", y= "ln(SDC)")+
  theme_bw()
pS5
#dev.copy2pdf(file=paste0("./Figures/FigS5_",Sys.Date(),".pdf"),width=8,height=6)

#plot(BA90_PCT~log(FIRESIZE_HA),data=fire.list2)
Fig5=
  grid.arrange(h,arrangeGrob(pS4,pS5,ncol=2),ncol=1,heights=c(2,2),clip=F)
#dev.copy2pdf(file=paste0("./Figures/Fig5_",Sys.Date(),".pdf"),width=8,height=6)

hist(log(fire.list2$AWMSI))
pS6=ggplot(fire.list2,aes(x=log(AWMSI),y=log(sdc)))+ 
  #geom_hline(yintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1)+
  geom_point()+
  geom_smooth(method="lm",col="gray")+
  labs(x ="ln(AWMSI) ", y= expression(ln(c[d])))+
  theme_bw()
pS6 #Use this one

pS7=ggplot(fire.list2,aes(x=(MSI),y=log(sdc)))+ 
  #geom_hline(yintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1)+
  geom_point()+
  geom_smooth(method="lm",col="gray")+
  labs(x ="MSI ", y= expression(ln(c[d])))+
  theme_bw()
pS7

hist(fire.list2$AWMPFD)
pS8=ggplot(fire.list2,aes(x=(AWMPFD),y=log(sdc)))+ 
  #geom_hline(yintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1)+
  geom_point()+
  geom_smooth(method="lm",col="gray")+
  labs(x ="AWMPFD ", y= expression(ln(c[d])))+
  theme_bw()
pS8 #Use this one

pS9=ggplot(fire.list2,aes(x=(MPFD),y=log(sdc)))+ 
  #geom_hline(yintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1)+
  geom_point()+
  geom_smooth(method="lm",col="gray")+
  labs(x ="MPFD ", y= expression(ln(c[d])))+
  theme_bw()
pS9

plot(AWMPFD~log(AWMSI), data=fire.list2)
#Which fires exceed the threshold of the purple line?
library(Kmisc)
write.cb(fire.list2[which(fire.list2$q<0.002),c("FIRE_NAME","FIRE_YEAR")])
