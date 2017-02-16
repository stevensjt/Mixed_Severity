##Analyzing fires with 1000 ha of high severity
#Formerly Fires1kha.R
####0. load libraries####
library(rgdal)
library(dismo)
library(ggplot2)
library(ggrepel)
require(rgeos)
library(gridExtra)
library(plyr)
library(dplyr)
source("./Analyses/Functions.R")

####1. Load and Filter Data####
#fires=readOGR("/Users/Jens/Documents/Davis/Post-Doc/Side Projects/Mixed Severity/GIS/BurnSev/", layer="VegBurnSeverity85-15") #Full layer
#fires=readOGR("/Users/Jens/Documents/Davis/Post-Doc/Side Projects/Mixed Severity/GIS/BurnSev/", layer="VegBurnSeverity_Sierra_85-15")#Shapefile with all Sierra fires since 2000 (smaller data file for EDA)

#hs_patches=fires[fires$BURNSEV==4&fires$BEST_ASSES=="YES",] #Extract only the high-severity patches, which shrinks down the files size
#saveRDS(hs_patches,"./hs_patches.RDS")
hs_patches=ReadRDS("./hs_patches.RDS")
#rm(fires)#Clear space
gc()#Clear space

#fire.list=read.csv("/Users/Jens/Documents/Davis/Post-Doc/Side Projects/Mixed Severity/Jay's exploratory analysis/fires_84-14.csv") #Jay's file
fire.list=read.csv("/Users/Jens/Documents/Davis/Post-Doc/Side Projects/Mixed Severity/Jay's exploratory analysis/Patch_analysis_Results.csv") #Jay's better file
#View(fires@data)

#fire.list=fire.list[between(fire.list$BA90_HA,900,1100)&as.character(fire.list$Veg)=="Forest",]
fire.list=fire.list[(as.character(fire.list$AGENCY)!="USF"&as.character(fire.list$VB_ID)!="2014KING")&as.character(fire.list$Veg)=="Forest",] 
fires.to.sample=as.character(fire.list$VB_ID)

#To get internal buffers for specific fires
#fires.to.sample=c("1987EAST","2008CARIBOU")

####2. Run Geospatial Analysis - Internal Buffering####
Sys.time() #Takes ~12 hours depending on sample size. Should only need to do once.
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
#Save as "Long" dataset, which has the buffering done for each fire (rows=buffer widths).
#write.csv(all_fires,file="./Analyses/Processed Data/NPS_CDF_Long.csv")


####3. Run Statistical Analysis - calculate metric####
#fires_for_stats=all_fires #If running on data you just did geospatial analyses

#Read in processed full geospatial data (optional):
#fires_for_stats=read.csv("./Analyses/Processed Data/USFS_King_Long.csv")
fires_for_stats=rbind(read.csv("./Analyses/Processed Data/USFS_King_Long.csv"),
                      read.csv("./Analyses/Processed Data/NPS_CDF_Long.csv")) #All fires

#Subset to specific fires (optional):
#fires_for_stats=fires_for_stats[as.character(fires_for_stats$name)%in%c("1987EAST","2008CARIBOU"),]
#fires_for_stats$name=as.character(fires_for_stats$name)

#Analyze the fires_for_stats table to calculate the parameters of interest
fires_with_stats=calculate.q2(fires_for_stats) #Fast!

#Save as "Full" dataset, which has the statistic calculated for each row (rows=buffer widths).
#write.csv(fires_with_stats,file="./Analyses/Processed Data/all_fires_Full.csv")

#Summarize to a single row per fire 
summary_fires=ddply(fires_with_stats,.(name),summarize,q=mean(q))

#Add the relevant parameters to the fire.list file, with a single row per fire ("ForAnalysis" dataset)
fire.list$q=summary_fires$q
#write.csv(fire.list,file="./Analyses/Processed Data/all_fires_ForAnalysis.csv")


####4.Plot Specific fires####
#Gives a nice comparison of the East and Caribou fires.
#Load if you haven't yet:
#fires_with_stats=read.csv("./Analyses/Processed Data/all_fires_Full.csv")
fires.to.plot.names=c("1987EAST","2008CARIBOU")
fires.to.plot=fires_with_stats[as.character(fires_with_stats$name)%in%fires.to.plot.names,]
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

####5. Equations and Thresholds####
#Tom's suggestion
"y = 10^(q*x)"
"0.25 = 10^(q*250)"
"log10(0.25)=q*250"
log10(0.25)/250 #=q; q=-0.0024
"A=10^q"
10^-0.0024 # = A; A=0.994489
"q=log10(A)"
"y=10^(log10(A)*x)"
"y=A^x"
0.994489^250 # = y; y=0.2512
#But A is still not normally distributed
hist(10^fire.list2$q)

#My first attempt
"y = 10^(q*x)"
"0.25 = 10^(q*250)"
"log10(0.25)=q*250"
log10(0.25)/250 #=q; q=-0.0024
"R=log10(-q)" #**Define R, which normalizes q; long tail to left are q values of concern.
hist(log10(-fire.list2$q))
log10(0.0024) #=-2.619 = R
"-q=10^R"
"q = -10^R"
-10^(-2.619) #=q; q=-0.0024. This works.
"y=10^((-10^R)*x)"
10^((-10^-2.619)*250) #=0.2506 = y. This works. But it's a messy equation.

#My take 2:
"y=1/(10^(qx))" #This is it!
"R=log10(q)"
"q=10^R"
"y=1/(10^((10^R)*x))"


####6. Summary Analyses####
#Read in processed data
fire.list2=read.csv("./Analyses/Processed Data/all_fires_ForAnalysis.csv")
#fire.list2=fire.list[,c(2:18,50)]#Remove irrelevant columns (if re-doing the whole buffering analysis)

#6.1 Histograms
h=ggplot(fire.list2,aes(log(q)))+
  geom_histogram(binwidth=0.1,fill="white",col="black")+
  geom_vline(xintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1.2)+
  labs(x = expression(ln(c[d])))+
  theme_bw()
h
#dev.copy2pdf(file=paste0("./Figures/FigS1_",Sys.Date(),".pdf"),width=8,height=8)

pS2=ggplot(fire.list2[fire.list2$FIRESIZE_HA>1000,],aes(x=log(FIRESIZE_HA),y=log(q)))+
  geom_point()+
  geom_smooth(method="lm",col="black")+
  geom_hline(yintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1)+
  labs(title="Fires > 1000 ha")+
  theme_bw()
pS2
#dev.copy2pdf(file=paste0("./Figures/FigS2_",Sys.Date(),".pdf"),width=8,height=6)

pS3=ggplot(fire.list2[fire.list2$FIRESIZE_HA<1000,],aes(x=log(FIRESIZE_HA),y=log(q)))+
  geom_point()+
  geom_smooth(method="lm",col="black")+
  geom_hline(yintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1)+
  labs(title="Fires < 1000 ha")+
  theme_bw()
pS3
#dev.copy2pdf(file=paste0("./Figures/FigS3_",Sys.Date(),".pdf"),width=8,height=6)

pS4=ggplot(fire.list2,aes(x=log(FIRESIZE_HA),y=log(q)))+ #Optional: add ,col=BA90_PCT
  geom_hline(yintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1)+
  geom_point()+
  geom_smooth(method="lm",col="gray")+
  labs(x ="fire size (ln[ha])\n", y= expression(ln(c[d])))+
  theme_bw()
pS4
#dev.copy2pdf(file=paste0("./Figures/FigS4_",Sys.Date(),".pdf"),width=8,height=6)

pS5=ggplot(fire.list2,aes(x=BA90_PCT,y=log(q)))+ 
  geom_hline(yintercept=log(c(0.0219,0.0068,0.0020,0.0006)),col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),size=1)+
  geom_point()+
  geom_smooth(method="lm",col="gray")+
  labs(x ="percentage mapped \n as high-severity", y= expression(ln(c[d])))+
  theme_bw()
pS5
#dev.copy2pdf(file=paste0("./Figures/FigS5_",Sys.Date(),".pdf"),width=8,height=6)

#plot(BA90_PCT~log(FIRESIZE_HA),data=fire.list2)
Fig5=
  grid.arrange(h,arrangeGrob(pS4,pS5,ncol=2),ncol=1,heights=c(2,2),clip=F)
dev.copy2pdf(file=paste0("./Figures/Fig5_",Sys.Date(),".pdf"),width=8,height=6)

#Which fires exceed the threshold of the purple line?
library(Kmisc)
write.cb(fire.list2[which(fire.list2$q<0.002),c("FIRE_NAME","FIRE_YEAR")])
