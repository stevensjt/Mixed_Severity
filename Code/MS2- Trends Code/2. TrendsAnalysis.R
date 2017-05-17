##This code analyzes trends in SDC##
##Jens Stevens; stevensjt@gmail.com

library(tidyverse)
library(zoo) #For moving window
library(rgdal)
library(dismo)
library(ggrepel)
require(rgeos)
library(gridExtra)
library(glmulti)
library(tree)
library(ReporteRs)
source("./Code/Functions.R")

####1. Read in and filter data for analysis####
d <- read.csv("./Data/Derived/all_fires_ForAnalysis_weather.csv")
d <- d[,-which(names(d)%in%c("MSI","MPFD"))] #not using mean shape index or mean patch fractal dimension (spatial stats)
d <- d[-which(d$FIRE_YEAR==2016),] #N=477
names(d) <- c("VB_ID","ID_Num","fire_name","fire_year","class","veg_type","agency","ICS_code","pct_fs", #MGMT variables
              "firesize_ha","BA90_ha","BA90_pct","max_patch_ha","ignition_date","contain_date", #FS stats
              "sdc","awmsi","awmpfd","max_tmmx","max_tmmn","min_rmax","max_bi") #spatial stats & weather


d[d$agency%in%c("BIA","CCO","USF,NPS"), "agency"]= NA #Get rid of a few underrepresented agencies/hybrids
d$agency=factor(d$agency)
d$class <- factor(d$class,levels=c("no","yes"), labels =c("SUP","WFU"))
d[which(d$min_rmax==100),c(19:22)]=NA #Two funky fires that had 100 percent humidity and very low temps; suspect potential bias in extraction of meteorological data so excluding met data.
to_scale <-c("fire_year","max_tmmx","max_tmmn","min_rmax","max_bi")
d[,paste(to_scale,"std",sep="_")]=scale(d[,to_scale])

#CHECKME the 2010 fires had pretty high SDC values but are still included in the year split in the regression tree. Tweak them so this isn't the case
#d[d$fire_year==2010 & d$class=="SUP","fire_year"]=2009

####2. Exploratory analyses####
#2a: Check for normality
hist(d$sdc)
hist(log(d$sdc))




####3a: SDC model selection####
#tmax and tmin are correlated, so just using tmax
#WFU and AGENCY are clearl the most important variables
#Best AIC: FIRE_YEAR, AGENCY, WFU, max_tmmx, max_tmmn (890.12). tmmx and tmmn are correlated, so they offset each other in the model but there's a slightly greater effect size for tmmx
#Equivocal model removes tmmn (and reduce the effect size of tmmx (890.73))
#Equivocal model replaces FIRE_YEAR with burn index (they are correlated) (890.86)
#Big jumps happen from 26:27 (can't really explain; adding FIRE_YEAR to an equation that had all the weather variables except tmax) and 30:31 (can't really explain; adding tmmn and rmax to equation that just had FIRE_YEAR plus WFU and AGENCY). Generally it's a pretty gradual shift.

potential_parms=c("fire_year","class","agency","max_tmmx","max_tmmn","min_rmax","max_bi")
potential_parms_std=c("fire_year_std","class","agency","max_tmmx_std","max_tmmn_std","min_rmax_std","max_bi_std")
#m2c <- glmulti(y=log(sdc)~ 0 + fire_year + class + agency + max_tmmx + max_tmmn + min_rmax + max_bi, 
#               data=d,level=1)
m2c <- glmulti(y="log(sdc)", 
        xr=potential_parms,
        data=d,level=1,method="h") #If running interactions (level 2), try genetic algorithm (method = "g")

#Make table
m_max <- 5
x=matrix(NA,nrow=length(c("AIC",rev(rownames(coef(m2c) ) ) ) ),ncol=m_max+1)
x[,1] <-c("AIC",rev(rownames(coef(m2c))))

for(m in c(1:m_max)){
  c <- round(summary(m2c@objects[[m]])$coefficients[,1],3)
  x[which(x[,1]%in%names(c)),m+1] <- c
  x[1,m+1] <- round(AIC(m2c@objects[[m]]),2)
}
# Set up general table properties and formatting
cell_p = cellProperties(padding.right=3, padding.left=3)
par_p = parProperties(text.align="right")
# Create table
ft = FlexTable(x, header.columns=FALSE, body.cell.props=cell_p, body.par.props=par_p)
ft = addHeaderRow(ft, text.properties=textBold(), c("","Model #"),
                  colspan=c(1,m_max), par.properties=parCenter())
ft = addHeaderRow(ft, text.properties=textBold(), c("Model AIC \n /coefficients",c(1:m_max)),
                  colspan=rep(1,times=m_max+1), par.properties=parCenter())
ft


####3b: Regression tree with best model from above####
library(rpart)
library(rpart.plot)
#START HERE the formula below is a good one; can probably justify a candidate model above without tmmn. I think the weird tmmx<39 result can be explained by the fact that most fires were in the northwest in 1987 when it was probably very hot but also complex topography can give more complex stand-replacing fire dynamics. Need to add a dummy variable for the Northwest vs Sierra.
#http://blog.revolutionanalytics.com/2013/06/plotting-classification-and-regression-trees-with-plotrpart.html
tree.1 <- rpart(formula = formula(m2c@objects[[2]]), data=d)

#png(file = paste0("./Manuscripts/MS2- Trends/Figures/Fig2_",Sys.Date(),".png"),width=4,height=4,units="in",res=200)
prp(tree.1)					# Will plot the tree
#dev.off()




####4. Trends over time####
####4a: Moving-window estimates (averaged per year)####
d.annual <-
  group_by(d,fire_year) %>%
  summarise(sdc=mean(sdc), log_sdc=mean(log(sdc)), BA90_pct = mean(BA90_pct),
            max_bi=mean(max_bi,na.rm=T),max_tmmx=mean(max_tmmx,na.rm=T))
log_sdc_mvw=rollapply(zoo(log(d.annual$sdc)), width = 5, by = 1, FUN = mean, align = "center")
max_bi_mvw=rollapply(zoo(d.annual$max_bi), width = 5, by = 1, FUN = mean, align = "center")
max_tmmx_mvw=rollapply(zoo(d.annual$max_tmmx), width = 5, by = 1, FUN = mean, align = "center")
BA90_pct_mvw=rollapply(zoo(d.annual$BA90_pct), width = 5, by = 1, FUN = mean, align = "center")
d.annual[rownames(as.data.frame(log_sdc_mvw)),"log_sdc_mvw5"]=log_sdc_mvw
d.annual[rownames(as.data.frame(max_bi_mvw)),"max_bi_mvw5"]=max_bi_mvw
d.annual[rownames(as.data.frame(max_tmmx_mvw)),"max_tmmx_mvw5"]=max_tmmx_mvw
d.annual[rownames(as.data.frame(max_tmmx_mvw)),"BA90_pct_mvw5"]=BA90_pct_mvw

sdc <- 
  ggplot(d.annual,aes(x=fire_year,y=log_sdc))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x="fire year", y= "ln(sdc)",title="average SDC")+
  theme_bw()
sdc_mvw <- 
  ggplot(d.annual,aes(x=fire_year,y=log_sdc_mvw5))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x="fire year", y= "ln(sdc)",title="5-year moving window average SDC")+
  theme_bw()
bi <- 
  ggplot(d.annual,aes(x=fire_year,y=max_bi))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x="fire year", y= "burn index",title="average burn index")+
  theme_bw()
bi_mvw <- 
  ggplot(d.annual,aes(x=fire_year,y=max_bi_mvw5))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x="fire year", y= "burn index",title="5-year moving window average burn index")+
  theme_bw()
tmmx <- 
  ggplot(d.annual,aes(x=fire_year,y=max_tmmx))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x="fire year", y= "maximum temperature (C)",title="average maximum temperature (C)")+
  theme_bw()
tmmx_mvw <- 
  ggplot(d.annual,aes(x=fire_year,y=max_tmmx_mvw5))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x="fire year", y= "maximum temperature (C)",title="5-year moving window average maximum temperature (C)")+
  theme_bw()

#png(file = paste0("./Manuscripts/MS2- Trends/Figures/Fig3_",Sys.Date(),".png"),width=6,height=9,units="in",res=200)
grid.arrange(sdc,sdc_mvw,bi,bi_mvw,tmmx,tmmx_mvw,ncol=2)
#dev.off()

#Changes over time in log(sdc)
#Significantly more negative in both the total data and the 5-year average
summary(lm(log(sdc)~fire_year,data=d.annual))
summary(lm(log_sdc_mvw5~fire_year,data=d.annual))

#Changes over time in burn index
#Significantly more positive
summary(lm(max_bi~fire_year,data=d.annual))
summary(lm(max_bi_mvw5~fire_year,data=d.annual))

#Changes over time in max temperature
#No significant trend (significant if you consider the 5-year average)
summary(lm(max_tmmx~fire_year,data=d.annual))
summary(lm(max_tmmx_mvw5~fire_year,data=d.annual))

#Changes over time in percent high-severity
#Significantly more positive, but not significant if you consider the 5-year average.
#A stronger signal in SDC than in percent high severity
summary(lm(BA90_pct~fire_year,data=d.annual))
summary(lm(BA90_pct_mvw5~fire_year,data=d.annual))


#4b: percentile analysis CHECKME this part needs to be updated.
fire.list$BA90_PCT_BIN=data.frame(fire.list$BA90_PCT, bin=cut(fire.list$BA90_PCT, c(seq(from=0,to=10),15,20,25,30,35,40,45,50,60,70,80), include.lowest=TRUE))$bin
for(l in levels(fire.list$BA90_PCT_BIN)){
  fire.list[grep(l,as.character(fire.list$BA90_PCT_BIN),fixed=T),"q.pct"]=
    ecdf(fire.list[grep(l,as.character(fire.list$BA90_PCT_BIN),fixed=T),"q"]) (fire.list[grep(l,as.character(fire.list$BA90_PCT_BIN),fixed=T),"q"])
}
hist(fire.list$q.pct)
ggplot(fire.list[fire.list$BA90_PCT>20,],aes(x=FIRE_YEAR,y=q.pct))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()

#5: 
#Trouble with plotting NA's as gray
#df <- data.frame (V1=factor(c("A","B","A","B",NA)),x=c(1:5),y=c(1:5))
#ggplot(df,aes(x=x,y=y,col=V1))+
#  geom_point()

agency_pct <-
  ggplot(na.omit(d[,c("BA90_pct","sdc","agency")]),aes(x=BA90_pct,y=log(sdc),col=agency))+
  geom_point()+
  scale_color_manual(values=c("darkred","darkgreen","orange"))+
  geom_smooth(method="lm")+
  labs(title="All Fires")+
  theme_bw()
agency_ha <-
  ggplot(na.omit(d[,c("firesize_ha","sdc","agency")]),aes(x=log(firesize_ha),y=log(sdc),col=agency))+
  geom_point()+
  scale_color_manual(values=c("darkred","darkgreen","orange"))+
  geom_smooth(method="lm")+
  labs(title="All Fires")+
  theme_bw()
class_pct <-
  ggplot(na.omit(d[,c("BA90_pct","sdc","class")]),aes(x=BA90_pct,y=log(sdc),col=class))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="All Fires")+
  theme_bw()
class_ha <-
  ggplot(na.omit(d[,c("firesize_ha","sdc","class")]),aes(x=log(firesize_ha),y=log(sdc),col=class))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="All Fires")+
  theme_bw()

summary(lm(log(sdc)~BA90_pct+class,data=d))
summary(lm(log(sdc)~firesize_ha+class,data=d))
summary(lm(log(sdc)~BA90_pct+agency,data=d))
summary(lm(log(sdc)~firesize_ha+agency,data=d))

png(file = paste0("./Manuscripts/MS2- Trends/Figures/Fig4_",Sys.Date(),".png"),width=10,height=10,units="in",res=200)
grid.arrange(class_pct,class_ha,agency_pct,agency_ha,ncol=2)
dev.off()

####5. Identify and plot examples####
#Run the chunk below for BA_range <- c(0,10), c(10,20), c(20,30), c(30-40), and c(45-55)
BA_range <- c(45,55)
BA_name <- paste("BA",BA_range[1],BA_range[2],"candidate",sep="_")
for(r in 1:nrow(d)){ #Find fires with BA90 in a desired range
  if(between(d$BA90_pct[r],BA_range[1],BA_range[2])){
    d[r,BA_name] = TRUE
  } else d[r,BA_name] = FALSE
}
FS_range <- quantile(as.vector(na.exclude(d[d[,BA_name] == TRUE,"firesize_ha"])),c(0.45,0.55))
for(r in 1:nrow(d)){ #Identify those fires of similar size
  if(d[r,BA_name] & between(d$firesize_ha[r],FS_range[1],FS_range[2])){
    d[r,BA_name] = TRUE
  } else {d[r,BA_name] = FALSE}
}
for(r in 1:nrow(d)){ #Identify those fires of similar size
  if(d[r,BA_name] & between(d$firesize_ha[r],FS_range[1],FS_range[2])){
    d[r,BA_name] = TRUE
  } else {d[r,BA_name] = FALSE}
}
d[d$sdc %in% c( max(d[d[,BA_name]==TRUE,"sdc"]), min(d[d[,BA_name]==TRUE,"sdc"]) ),"Plot_candidates"]=TRUE

#Plot examples from 10-55% HS
d[which(d$Plot_candidates),"Plot_candidates"] <- c(5,2,3,1,7,8,6,4)
d[which(is.na(d$Plot_candidates)),"Plot_candidates"] <- ""
hs_patches=readOGR("../Large Files/GIS/BurnSev/Current/", layer="hs_patches") #CRS EPSG:3310, NAD83 CA Albers

####4.2.Plot Specific fires####
#Gives a nice comparison of the East and Caribou fires.
fires.to.plot.names=as.character(d[ order(d[,"Plot_candidates"])[470:477] , "VB_ID"])
#fires.to.plot=fires_long[as.character(fires_long$name)%in%fires.to.plot.names,]
#fires.to.plot$name=factor(fires.to.plot$name,labels=c("East","Caribou"))
p.a=ggplot()+
  geom_polygon(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[1],]),
               aes(x=long-min(long),y=lat-min(lat),group=group),col='darkred',fill='darkred')+
  xlim(0,10000) + ylim(0,8000) + coord_fixed()+
  labs(title="1. Highway Fire (2001)",x=" ",y="meters") +
  annotate("text",x=0,y=6000,label=
             paste("sdc =",round(d[d$VB_ID==fires.to.plot.names[1],"sdc"],4),
                   "\nln(sdc) = ",round(log(d[d$VB_ID==fires.to.plot.names[1],"sdc"]),2) ),
           hjust=0,vjust=0,size=3)+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=18, hjust=0.5))

p.b=ggplot()+
  geom_polygon(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[2],]),
               aes(x=long-min(long),y=lat-min(lat),group=group),col='darkred',fill='darkred')+
  xlim(0,10000) + ylim(0,8000) + coord_fixed()+
  labs(title="2. Recer Fire (1990)",x=" ",y="meters") +
  annotate("text",x=0,y=6000,label=
             paste("sdc =",round(d[d$VB_ID==fires.to.plot.names[2],"sdc"],4),
                   "\nln(sdc) = ",round(log(d[d$VB_ID==fires.to.plot.names[2],"sdc"]),2) ),
           hjust=0,vjust=0,size=3)+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=18, hjust=0.5))

p.c=ggplot()+
  geom_polygon(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[3],]),
               aes(x=long-min(long),y=lat-min(lat),group=group),col='darkred',fill='darkred')+
  xlim(0,10000) + ylim(0,8000) + coord_fixed()+
  labs(title="3. Horton Fire (1999)",x=" ",y="meters") +
  annotate("text",x=0,y=6000,label=
             paste("sdc =",round(d[d$VB_ID==fires.to.plot.names[3],"sdc"],4),
                   "\nln(sdc) = ",round(log(d[d$VB_ID==fires.to.plot.names[3],"sdc"]),2) ),
           hjust=0,vjust=0,size=3)+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=18, hjust=0.5))

p.d=ggplot()+
  geom_polygon(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[4],]),
               aes(x=long-min(long),y=lat-min(lat),group=group),col='darkred',fill='darkred')+
  xlim(0,10000) + ylim(0,8000) + coord_fixed()+
  labs(title="4. Goose Fire (2001)",x=" ",y="meters") +
  annotate("text",x=0,y=6000,label=
             paste("sdc =",round(d[d$VB_ID==fires.to.plot.names[1],"sdc"],4),
                   "\nln(sdc) = ",round(log(d[d$VB_ID==fires.to.plot.names[1],"sdc"]),2) ),
           hjust=0,vjust=0,size=3)+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=18, hjust=0.5))

p.e=ggplot()+
  geom_polygon(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[5],]),
               aes(x=long-min(long),y=lat-min(lat),group=group),col='darkred',fill='darkred')+
  xlim(0,10000) + ylim(0,8000) + coord_fixed()+
  labs(title="5. Rack Fire (1989)",x=" ",y="meters") +
  annotate("text",x=0,y=6000,label=
             paste("sdc =",round(d[d$VB_ID==fires.to.plot.names[5],"sdc"],4),
                   "\nln(sdc) = ",round(log(d[d$VB_ID==fires.to.plot.names[5],"sdc"]),2) ),
           hjust=0,vjust=0,size=3)+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=18, hjust=0.5))

p.f=ggplot()+
  geom_polygon(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[6],]),
               aes(x=long-min(long),y=lat-min(lat),group=group),col='darkred',fill='darkred')+
  xlim(0,10000) + ylim(0,8000) + coord_fixed()+
  labs(title="6. Boulder Complex Fire (2006)",x=" ",y="meters") +
  annotate("text",x=0,y=6000,label=
             paste("sdc =",round(d[d$VB_ID==fires.to.plot.names[6],"sdc"],4),
                   "\nln(sdc) = ",round(log(d[d$VB_ID==fires.to.plot.names[6],"sdc"]),2) ),
           hjust=0,vjust=0,size=3)+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=18, hjust=0.5))

p.g=ggplot()+
  geom_polygon(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[7],]),
               aes(x=long-min(long),y=lat-min(lat),group=group),col='darkred',fill='darkred')+
  xlim(0,10000) + ylim(0,8000) + coord_fixed()+
  labs(title="7. Stream Fire (2001)",x=" ",y="meters") +
  annotate("text",x=0,y=6000,label=
             paste("sdc =",round(d[d$VB_ID==fires.to.plot.names[7],"sdc"],4),
                   "\nln(sdc) = ",round(log(d[d$VB_ID==fires.to.plot.names[7],"sdc"]),2) ),
           hjust=0,vjust=0,size=3)+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=18, hjust=0.5))

p.h=ggplot()+
  geom_polygon(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[8],]),
               aes(x=long-min(long),y=lat-min(lat),group=group),col='darkred',fill='darkred')+
  xlim(0,10000) + ylim(0,8000) + coord_fixed()+
  labs(title="8. Sims Fire (2004)",x=" ",y="meters") +
  annotate("text",x=0,y=6000,label=
             paste("sdc =",round(d[d$VB_ID==fires.to.plot.names[8],"sdc"],4),
                   "\nln(sdc) = ",round(log(d[d$VB_ID==fires.to.plot.names[8],"sdc"]),2) ),
           hjust=0,vjust=0,size=3)+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=18, hjust=0.5))

png(file = paste0("./Manuscripts/MS2- Trends/Figures/Fig5_",Sys.Date(),".png"),width=10,height=10,units="in",res=200)
grid.arrange(p.a,p.b,p.c,p.d,p.e,p.f,p.g,p.h,ncol=2)
dev.off()

#Ultimately better way to plot is below, using "draw_plot" from cowplot.
ggdraw()+
  draw_plot(p.a, x=0, y=.5, width=.5, height=.5) +
  draw_plot(p.b, x=.5, y=.5, width=.5, height=.5) +
  draw_plot(p.fits, x=0, y=0, width=1, height=.5) +
  draw_plot_label(c("A", "B", "C", "meters"), c(0, 0.5, 0, 0.40), c(1, 1, 0.5, 0.55), 
                  size = c(12, 12, 12, 16), fontface=c("bold","bold","bold","plain")) 

#dev.copy2pdf(file=paste0("./Figures/Fig4_",Sys.Date(),".pdf"),width=8,height=8)


#Redo Fig. 5
class_pct <-
  ggplot(na.omit(d[,c("BA90_pct","sdc","class","Plot_candidates")]),aes(x=BA90_pct,y=log(sdc),col=class))+
  geom_point()+
  geom_text(aes(label=Plot_candidates),hjust=1, vjust=1,size=6,col="black",fontface="bold")+
  geom_smooth(method="lm")+
  labs(title="All Fires")+
  theme_bw()
png(file = paste0("./Manuscripts/MS2- Trends/Figures/Fig4_",Sys.Date(),".png"),width=10,height=10,units="in",res=200)
grid.arrange(class_pct,class_ha,agency_pct,agency_ha,ncol=2)
dev.off()
