#ESA SDC analysis

####EDA PLOTS####
fire.list=fire.list[fire.list$Veg=="Forest",]
#***WHERE IS WHITES FIRE? It's in HS_patches but not in fire.list; Might need to cross-check the current fire.list with fires_usfs and pull in any that are missing to get the SDC values#
for(y in 1984:2016){
  q=quantile(fire.list[fire.list$FIRE_YEAR==y,"sdc"],0.2)
  fire.list[fire.list$FIRE_YEAR==y,"q20"]=ifelse(fire.list[fire.list$FIRE_YEAR==y,"sdc"]<q,1,0)
}
####1. Change over time
#All forest fires:
ggplot(fire.list)+
  geom_point(aes(x=FIRE_YEAR,y=log(sdc),col=BA90_PCT))+
  geom_smooth(aes(x=FIRE_YEAR,y=log(sdc)),method="lm",col="black")+
  labs(col="% High Severity",y="log(SDC)")
summary(lm(log(sdc)~FIRE_YEAR,data=fire.list)) #t=-2.437,df=461,P=0.0152

#Is percent high-severity changing too?
ggplot(fire.list)+
  geom_point(aes(x=FIRE_YEAR,y=log(BA90_PCT),col=BA90_PCT))+
  geom_smooth(aes(x=FIRE_YEAR,y=log(BA90_PCT)),method="lm",col="black")
summary(lm(log(BA90_PCT)~FIRE_YEAR,data=fire.list)) #Yes, P=0.02

#Only bottom 20% of sdc within a given year
d=fire.list[fire.list$q20==1,]
ggplot()+
  geom_point(aes(x=FIRE_YEAR,y=log(sdc),col=BA90_PCT))+
  geom_smooth(aes(x=FIRE_YEAR,y=log(sdc)),method="lm",col="black")
summary(lm(log(sdc)~FIRE_YEAR,data=d)) #NS

#All fires except widespread fire years
d=fire.list[!fire.list$FIRE_YEAR%in%c(1987,2008),]
ggplot(d)+
  geom_point(aes(x=FIRE_YEAR,y=log(sdc),col=BA90_PCT))+
  geom_smooth(aes(x=FIRE_YEAR,y=log(sdc)),method="lm",col="black")
summary(lm(log(sdc)~FIRE_YEAR,data=d)) #t=-3.060,df=345,P=0.002

#What about model residuals?
fire.list[!is.na(fire.list$q),"sdc_BA90_resid"]=resid(lm(log(q)~BA90_PCT,data=fire.list))
fire.list[which(!is.na(fire.list$q)&!is.na(fire.list$max_bi)),"sdc_BA90_maxbi_resid"]=resid(lm(log(q)~BA90_PCT+max_bi,data=fire.list))

ggplot(fire.list)+
  geom_point(aes(x=FIRE_YEAR,y=sdc_BA90_resid,col=BA90_PCT))+
  geom_smooth(aes(x=FIRE_YEAR,y=sdc_BA90_resid),method="lm",col="black")

ggplot(fire.list[-which(fire.list$FIRE_YEAR%in%c(1987, 1999, 2008)),])+
  geom_point(aes(x=FIRE_YEAR,y=sdc_BA90_maxtmmx_resid,col=min_rmax))+
  geom_smooth(aes(x=FIRE_YEAR,y=sdc_BA90_maxtmmx_resid),method="lm",col="black")

lm(sdc_BA90_maxbi_resid~FIRE_YEAR,
   data=fire.list[-which(fire.list$FIRE_YEAR%in%c(1987, 1999, 2008)),]) %>%
  summary()

####2. Managing agency####
d=fire.list[fire.list$AGENCY%in%c("CDF","NPS","USF"),]
d$AGENCY=factor(d$AGENCY,levels=c("NPS","CDF","USF"))
dls=d[d$BA90_PCT<30,]
dls2=d[d$BA90_PCT<max(d[d$AGENCY=="NPS","BA90_PCT"]),]
dhs=d[d$BA90_PCT>=30,]

ggplot(dls2)+
  geom_point(aes(x=BA90_PCT,y=log(sdc),col=AGENCY))+
  geom_smooth(aes(x=BA90_PCT,y=log(sdc),col=AGENCY),method="lm")

summary(lm(log(sdc)~BA90_PCT*AGENCY,data=d)) #Significant interaction
summary(lm(log(sdc)~AGENCY,data=dls2)) #Significantly more negative sdc in USFS and CDF
summary(lm(log(sdc)~AGENCY+BA90_PCT,data=dls2)) #SDC responds independently of %HS

####3. WFU####
ggplot(fire.list)+
  geom_point(aes(x=BA90_PCT,y=log(sdc),col=WFU))+
  geom_smooth(aes(x=BA90_PCT,y=log(sdc),col=WFU),method="lm")

dls2=fire.list[fire.list$BA90_PCT<max(fire.list[fire.list$WFU=="yes","BA90_PCT"]),]
dls2=dls2[dls2$FIRE_YEAR<2010,]
#dls2$WFU=as.character(dls2$WFU)
#dls2[dls2$FIRE_YEAR==2008,"WFU"]="2008"
ggplot(dls2[order(dls2$WFU),])+
  geom_point(aes(x=BA90_PCT,y=log(sdc),col=WFU))+
  geom_smooth(aes(x=BA90_PCT,y=log(sdc),col=WFU),method="lm")

summary(lm(log(sdc)~BA90_PCT*WFU,data=dls2)) #No significant interaction (slopes are same)
summary(lm(log(sdc)~WFU,data=dls2)) #Significantly more negative sdc in WFU than non-WFU; t=5.63,P<0.0001
summary(lm(log(sdc)~BA90_PCT + WFU,data=dls2)) #Significant effect AFTER controlling for mean % high severity

