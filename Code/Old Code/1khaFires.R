##Analyzing fires with 1000 ha of high severity

####0. load libraries####
library(rgdal)
library(dismo)
library(ggplot2)
library(ggrepel)
require(rgeos)
library(gridExtra)
library(dplyr)
library(plyr)
source("./Analyses/Functions.R")

####1. Load and Filter Data####
#fires=readOGR("/Users/Jens/Documents/Davis/Post-Doc/Side Projects/Mixed Severity/GIS/BurnSev/", layer="VegBurnSeverity85-15") #Full layer
#fires=readOGR("/Users/Jens/Documents/Davis/Post-Doc/Side Projects/Mixed Severity/GIS/BurnSev/", layer="VegBurnSeverity_Sierra_85-15")#Shapefile with all Sierra fires since 2000 (smaller data file for EDA)

#hs_patches=fires[fires$BURNSEV==4&fires$BEST_ASSES=="YES",] #Extract only the high-severity patches, which shrinks down the files size
rm(fires)#Clear space
gc()#Clear space

#fire.list=read.csv("/Users/Jens/Documents/Davis/Post-Doc/Side Projects/Mixed Severity/Jay's exploratory analysis/fires_84-14.csv") #Jay's file
fire.list=read.csv("/Users/Jens/Documents/Davis/Post-Doc/Side Projects/Mixed Severity/Jay's exploratory analysis/Patch_analysis_Results.csv") #Jay's better file
#View(fires@data)

#fire.list=fire.list[between(fire.list$BA90_HA,900,1100)&as.character(fire.list$Veg)=="Forest",]
fire.list=fire.list[(as.character(fire.list$AGENCY)!="USF"&as.character(fire.list$VB_ID)!="2014KING")&as.character(fire.list$Veg)=="Forest",] #START HERE FOR NPS FIRES.
fires.to.sample=as.character(fire.list$VB_ID)


####3. Run Analyses####
Sys.time()
for(f in c(1:length(fires.to.sample))){
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

all_fires=calculate.q(all_fires)
write.csv(all_fires,file="./Analyses/Processed Data/USFS_King_Full.csv")
#all_fires=read.csv("./Analyses/Processed Data/NPS_Full.csv")
summary_fires=ddply(all_fires,.(name),summarize,q=mean(q))
fire.list$q=summary_fires$q
write.csv(fire.list,file="./Analyses/Processed Data/NPS_ForAnalysis.csv")


####4.Plot Specific fires####
#1,7 gives a nice comparison of the East and Caribou fires.
fires.to.plot.names=as.character(summary_fires$name)[c(1,7)]
fires.to.plot=all_fires[as.character(all_fires$name)%in%fires.to.plot.names,]
p.a=ggplot()+
  geom_path(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[1],]),
            aes(x=long,y=lat,group=group),fill="transparent",col='darkred')+
  labs(title=fires.to.plot.names[1])+
  theme_bw()
p.b=ggplot()+
  geom_path(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[2],]),
            aes(x=long,y=lat,group=group),fill="transparent",col='darkblue')+
  labs(title=fires.to.plot.names[2])+
  theme_bw()

p.fits=ggplot(data=fires.to.plot,aes(x=width,y=prop.hs))+
  geom_point(aes(col=name))+
  scale_color_manual(values=c("darkred","darkblue"))+
  geom_path(aes(x=width,y=10^(q*width),group=name,col=name))+
  annotate("text", x = c(50,500), y = c(0.25,0.25), label = rev(unique(fires.to.plot$q.name)),size=3)+
  labs(y="proportion high-severity",x="distance to edge (m)",col="patch size")+
  theme_bw()

Fig2=grid.arrange(arrangeGrob(p.a,p.b,ncol=2),p.fits,ncol=1,widths=c(1,1),heights=c(2,1),clip=F)
dev.copy2pdf(file=paste0("./Figures/Fig2_",Sys.Date(),".pdf"),width=8,height=8)

####5. Trends Analyses####
fire.list2=rbind(read.csv("./Analyses/Processed Data/USFS_King_ForAnalysis.csv"),
                 read.csv("./Analyses/Processed Data/NPS_ForAnalysis.csv"))
fire.list2=fire.list2[,c(2:18,50)]#Remove irrelevant columns

p1=ggplot(fire.list2,aes(x=FIRE_YEAR,y=q))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="All Fires")+
  theme_bw()

p2=ggplot(fire.list2[fire.list2$FIRESIZE_HA>1000,],aes(x=FIRE_YEAR,y=q))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="Fires > 1000 ha")+
  theme_bw()

p3=ggplot(fire.list2[fire.list2$FIRESIZE_HA>10000,],aes(x=FIRE_YEAR,y=q))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="Fires > 10000 ha")+
  theme_bw()

p4=ggplot(fire.list2[fire.list2$FIRESIZE_HA>1000,],aes(x=WFU,y=q))+
  geom_bar(stat="summary",fun.y="mean")+
  #geom_smooth(method="lm")+
  labs(title="Fires > 1000 ha")+
  theme_bw()

p5a=ggplot(fire.list2[fire.list2$FIRESIZE_HA>1000,],aes(x=log(FIRESIZE_HA),y=q,col=WFU))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="Fires > 1000 ha")+
  theme_bw()

p5b=ggplot(fire.list2[fire.list2$FIRESIZE_HA>1000,],aes(x=log(FIRESIZE_HA),y=q,col=AGENCY))+
  geom_point()+
  scale_color_manual(values=c("darkred","darkgreen","orange"))+
  geom_smooth(method="lm",fill=NA)+
  labs(title="Fires > 1000 ha")+
  theme_bw()

p6a=ggplot(fire.list2[fire.list2$FIRESIZE_HA>1000,],aes(x=BA90_PCT,y=q,col=WFU))+
  geom_point()+
  geom_text(data=subset(fire.list2[fire.list2$FIRESIZE_HA>1000,],q>-0.002),
            aes(BA90_PCT,q,label=VB_ID),size=3,hjust=0,angle=65,col="black")+
  geom_smooth(method="lm")+
  labs(title="Fires > 1000 ha")+
  theme_bw()
p6b=ggplot(fire.list2[fire.list2$FIRESIZE_HA>1000,],aes(x=BA90_PCT,y=q,col=WFU))+
  geom_point()+
  geom_text_repel(data=subset(fire.list2[fire.list2$FIRESIZE_HA>1000,],q>-0.002| (as.character(VB_ID)%in%c("2004MEADOW","2001HOOVER","2014CHIPS","2014KING"))),
            aes(BA90_PCT,q,label=VB_ID),size=3,angle=65,col="black")+
  geom_smooth(method="lm")+
  labs(title="Fires > 1000 ha")+
  theme_bw()
