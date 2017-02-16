library(ggplot2)
library(plyr)
library(zoo) #For moving window

####1. Process data and calculate q####

####2. Read in data####
#Read in processed data
fire.list=read.csv("./Analyses/Processed Data/all_fires_ForAnalysis.csv")
hs_patches=readRDS("./hs_patches.RDS")

####3. Exploratory analyses####
#3a: Check for normality
hist(fire.list$q)
fire.list$logq=log(fire.list$q)
hist(fire.list$logq)

#3b: Moving-window estimates (averaged per year)
fire.list.annual=ddply(fire.list,.(FIRE_YEAR),summarize,q=mean(q),BA90_PCT=mean(BA90_PCT))

zoo.df=zoo(fire.list.annual$q)
mvw=rollapply(zoo.df, width = 5, by = 1, FUN = mean, align = "center")
fire.list.annual[rownames(as.data.frame(mvw)),"q.mvw5"]=mvw

zoo.df=zoo(fire.list.annual$BA90_PCT)
mvw=rollapply(zoo.df, width = 5, by = 1, FUN = mean, align = "center")
fire.list.annual[rownames(as.data.frame(mvw)),"BA90_PCT.mvw5"]=mvw

####4. Trends over time####
#4a: is q decreasing over time?
m4a=lm(logq~FIRE_YEAR,data=fire.list)
summary(m4a)
ggplot(fire.list,aes(x=FIRE_YEAR,y=logq))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="All Fires")+
  theme_bw()
#conclusion: yes, slightly and significantly

#4b: is percent high severity increasing over time?
m4b=lm(BA90_PCT~FIRE_YEAR,data=fire.list)
summary(m4b)
ggplot(fire.list,aes(x=FIRE_YEAR,y=BA90_PCT))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="All Fires")+
  theme_bw()
#conclusion: yes, slightly but not significantly

#4c: is q decreasing over time in moving-window analysis?
m4c=lm(log(q.mvw5)~FIRE_YEAR,data=fire.list.annual)
summary(m4c)
ggplot(fire.list.annual,aes(x=FIRE_YEAR,y=log(q.mvw5)))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="All Fires")+
  theme_bw()
#conclusion: no, but looks like some weather-driven patterns.

#4d: is percent high severity increasing over time in moving-window analysis?
m4d=lm(BA90_PCT.mvw5~FIRE_YEAR,data=fire.list.annual)
summary(m4d)
ggplot(fire.list.annual,aes(x=FIRE_YEAR,y=BA90_PCT.mvw5))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="All Fires")+
  theme_bw()
#conclusion: no, but looks like some weather-driven patterns.

#4e: percentile analysis
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

ggplot(fire.list,aes(x=BA90_PCT,y=logq,col=AGENCY))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="All Fires")+
  theme_bw()

ggplot(fire.list,aes(x=BA90_PCT,y=logq,col=WFU))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="All Fires")+
  theme_bw()