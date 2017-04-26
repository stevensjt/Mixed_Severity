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
        data=d,level=1)

#Make table
m_max <- 10
x=matrix(NA,nrow=length(c("AIC",rev(rownames(coef(m2c) ) ) ) ),ncol=m_max+1)
x[,1] <-c("AIC",rev(rownames(coef(m2c))))

for(m in c(1:10)){
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
                  colspan=c(1,10), par.properties=parCenter())
ft = addHeaderRow(ft, text.properties=textBold(), c("Model AIC \n /coefficients",c(1:10)),
                  colspan=rep(1,times=11), par.properties=parCenter())
ft


####3b: Regression tree with best model from above####
library(rpart)
library(rpart.plot)
#START HERE the formula below is a good one; can probably justify a candidate model above without tmmn. I think the weird tmmx<39 result can be explained by the fact that most fires were in the northwest in 1987 when it was probably very hot but also complex topography can give more complex stand-replacing fire dynamics.
#http://blog.revolutionanalytics.com/2013/06/plotting-classification-and-regression-trees-with-plotrpart.html
tree.1 <- rpart(formula = formula(m2c@objects[[2]]), data=d)

#plot(m2d); text(m2d)
png(file = paste0("./Manuscripts/MS2- Trends/Figures/Fig2_",Sys.Date(),".png"),width=4,height=4,units="in",res=200)
prp(tree.1)					# Will plot the tree
dev.off()




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
d.annual[rownames(as.data.frame(mvw_log_sdc)),"log_sdc_mvw5"]=log_sdc_mvw
d.annual[rownames(as.data.frame(mvw_max_bi)),"max_bi_mvw5"]=max_bi_mvw
d.annual[rownames(as.data.frame(mvw_max_tmmx)),"max_tmmx_mvw5"]=max_tmmx_mvw
d.annual[rownames(as.data.frame(mvw_max_tmmx)),"BA90_pct_mvw5"]=BA90_pct_mvw

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

png(file = paste0("./Manuscripts/MS2- Trends/Figures/Fig3_",Sys.Date(),".png"),width=6,height=9,units="in",res=200)
grid.arrange(sdc,sdc_mvw,bi,bi_mvw,tmmx,tmmx_mvw,ncol=2)
dev.off()

#Changes over time in log(sdc)
#Significantly more negative in both the total data and the 5-year average
summary(lm(log(sdc)~fire_year,data=d))
summary(lm(log_sdc_mvw5~fire_year,data=d.annual))

#Changes over time in burn index
#Significantly more positive
summary(lm(max_bi~fire_year,data=d))
summary(lm(max_bi_mvw5~fire_year,data=d.annual))

#Changes over time in max temperature
#No significant trend (significant if you consider the 5-year average)
summary(lm(max_tmmx~fire_year,data=d))
summary(lm(max_tmmx_mvw5~fire_year,data=d.annual))

#Changes over time in percent high-severity
#Significantly more positive, but not significant if you consider the 5-year average.
#A stronger signal in SDC than in percent high severity
summary(lm(BA90_pct~fire_year,data=d))
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

png(file = paste0("./Manuscripts/MS2- Trends/Figures/Fig4_",Sys.Date(),".png"),width=10,height=10,units="in",res=200)
grid.arrange(agency_pct,agency_ha,class_pct,class_ha,ncol=2)
dev.off()
