library(ggplot2)

####2. Trends Analyses####
#Read in processed data
fire.list=read.csv("./Analyses/Processed Data/NPS_ForAnalysis.csv")

fire.list2=rbind(read.csv("./Analyses/Processed Data/USFS_King_ForAnalysis.csv"),
                 read.csv("./Analyses/Processed Data/NPS_ForAnalysis.csv"))
fire.list2=fire.list2[,c(2:18,50)]#Remove irrelevant columns

p1=ggplot(fire.list2,aes(x=FIRE_YEAR,y=log(q)))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="All Fires")+
  theme_bw()

m1=lm(log(q)~FIRE_YEAR,data=fire.list2)
summary(m1)

p2=ggplot(fire.list2[fire.list2$FIRESIZE_HA>1000,],aes(x=FIRE_YEAR,y=log(q)))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="Fires > 1000 ha")+
  theme_bw()

m2a=lm(log(q)~FIRE_YEAR,data=fire.list2[fire.list2$FIRESIZE_HA>1000,])
summary(m2a)

p2b=ggplot(fire.list2[fire.list2$FIRESIZE_HA>1000&!(fire.list2$FIRE_YEAR%in%c(1987,2008)),],aes(x=FIRE_YEAR,y=log(q)))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="Fires > 1000 ha without '87 and '08")+
  theme_bw()

m2b=lm(log(q)~FIRE_YEAR,data=fire.list2[fire.list2$FIRESIZE_HA>1000&!(fire.list2$FIRE_YEAR%in%c(1987,2008)),])
summary(m2b)

p3=ggplot(fire.list2[fire.list2$FIRESIZE_HA>10000,],aes(x=FIRE_YEAR,y=log(q)))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="Fires > 10000 ha")+
  theme_bw()

p4=ggplot(fire.list2[fire.list2$FIRESIZE_HA>1000,],aes(x=WFU,y=log(q)))+
  geom_bar(stat="summary",fun.y="mean")+
  #geom_smooth(method="lm")+
  labs(title="Fires > 1000 ha")+
  theme_bw()

p5a=ggplot(fire.list2[fire.list2$FIRESIZE_HA>1000,],aes(x=log(FIRESIZE_HA),y=log(q),col=WFU))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="Fires > 1000 ha")+
  theme_bw()

p5b=ggplot(fire.list2[fire.list2$FIRESIZE_HA>1000,],aes(x=log(FIRESIZE_HA),y=log(q),col=AGENCY))+
  geom_point()+
  scale_color_manual(values=c("darkred","darkgreen","orange"))+
  geom_smooth(method="lm",fill=NA)+
  labs(title="Fires > 1000 ha")+
  theme_bw()

p6a=ggplot(fire.list2[fire.list2$FIRESIZE_HA>1000,],aes(x=BA90_PCT,y=log(q),col=WFU))+
  geom_point()+
  geom_text(data=subset(fire.list2[fire.list2$FIRESIZE_HA>1000,],log(q)<(-6)),
            aes(BA90_PCT,log(q),label=VB_ID),size=3,hjust=0,angle=65,col="black")+
  geom_smooth(method="lm")+
  labs(title="Fires > 1000 ha")+
  theme_bw()

p6b=ggplot(fire.list2[fire.list2$FIRESIZE_HA>1000,],aes(x=BA90_PCT,y=log(q),col=WFU))+
  geom_point()+
  geom_text_repel(data=subset(fire.list2[fire.list2$FIRESIZE_HA>1000,],q>-0.002| (as.character(VB_ID)%in%c("2004MEADOW","2001HOOVER","2014CHIPS","2014KING"))),
                  aes(BA90_PCT,q,label=VB_ID),size=3,angle=65,col="black")+
  geom_smooth(method="lm")+
  labs(title="Fires > 1000 ha")+
  theme_bw()