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
library(car) #For durbinWatsonTest
library(grid) #For viewport
library(RColorBrewer) #For brewer.pal
source("./Code/Functions.R")

####1. Read in and filter data for analysis####
d <- read.csv("./Data/all_fires_ForAnalysis_weather.csv")
d <- d[,-which(names(d)%in%c("MSI","MPFD"))] #not using mean shape index or mean patch fractal dimension (spatial stats)
d <- d[-which(d$FIRE_YEAR==2016),] #N=477
names(d) <- c("VB_ID","ID_Num","fire_name","fire_year","class","veg_type","agency","ICS_code","pct_fs", #MGMT variables
              "firesize_ha","BA90_ha","BA90_pct","max_patch_ha","ignition_date","contain_date", #FS stats
              "sdc","awmsi","awmpfd","max_tmmx","max_tmmn","min_rmax","max_bi") #spatial stats & weather


d$agency=as.character(d$agency)
d[d$agency%in%c("BIA","CCO"), "agency"]= NA #Get rid of a few underrepresented agencies
d[grep("USF",d$agency), "agency"]= "USFS" #Recode "USF" as "USFS" for consistency in text; this also classifies five fires with agency of "USF/NPS" as "USFS" because they had a National Forest as the ICS code.
d$agency <- factor(d$agency)
d$region <- ifelse(d$ICS_code%in%c("KNF","MEU","MNF","SHF","SHU","SRF"),"NW","SCSN") #Identify fires in NW CA.
d$region <- factor(d$region,levels=c("SCSN","NW"))
d$class <- factor(d$class,levels=c("no","yes"), labels =c("SUP","WFU"))
d[which(d$min_rmax==100),c(19:22)]=NA #Two funky fires that had 100 percent humidity and very low temps; suspect potential bias in extraction of meteorological data so excluding met data for these fires.
#write_csv(d,"./Data/Derived/data_for_SDC_paper.csv")

####2. Exploratory analyses####
#2a: Check for normality
#hist(d$sdc)
#hist(log(d$sdc))

#2b: Summary table (Table 2)
Table2 <- 
  d %>%
  group_by(agency,class) %>%
  summarize(N = length(sdc), min_ha = min(firesize_ha), median_ha = median(firesize_ha), 
            max_ha = max(firesize_ha,na.rm=T), median_fireyear = round(median(fire_year,na.rm=T),0), 
            mean_max_tmmx = mean(max_tmmx,na.rm=T), mean_max_bi = mean(max_bi,na.rm=T), 
            mean_max_tmmn = mean(max_tmmn,na.rm=T), mean_min_rmax = mean(min_rmax,na.rm=T)
  )

Table2[,c(8:11)] <- round(Table2[,c(8:11)],1)
Table2[,"agency"] <- as.character(Table2[,"agency"][[1]])
names(Table2) <- c("agency","class","N","min size \n(ha)","median \nsize (ha)","max size \n(ha)","median \nfire year",
                   "mean maximum \nhigh temperature","mean maximum\nburn index","mean maximum \nlow temperature", 
                   "mean minimum \nhigh humidity")
Table2[is.na(Table2$agency),"agency"] <- "NA"
# Set up general table properties and formatting
cell_p = cellProperties(padding.right=3, padding.left=3)
par_p = parProperties(text.align="right")
# Make Table
ft = FlexTable(Table2, header.columns=FALSE, body.cell.props=cell_p, body.par.props=par_p)
ft = addHeaderRow(ft, text.properties=textBold(), names(Table2),
                   par.properties=parCenter())
ft #Save as HTML Table 2

####3a: SDC model selection####
#tmax and tmin are correlated, so just using tmax
#WFU and AGENCY are clearl the most important variables
#Best AIC: FIRE_YEAR, AGENCY, WFU, max_tmmx, max_tmmn (890.12). tmmx and tmmn are correlated, so they offset each other in the model but there's a slightly greater effect size for tmmx
#Equivocal model removes tmmn (and reduce the effect size of tmmx (890.73))
#Equivocal model replaces FIRE_YEAR with burn index (they are correlated) (890.86)
#Big jumps happen from 26:27 (can't really explain; adding FIRE_YEAR to an equation that had all the weather variables except tmax) and 30:31 (can't really explain; adding tmmn and rmax to equation that just had FIRE_YEAR plus WFU and AGENCY). Generally it's a pretty gradual shift.

potential_parms=c("fire_year","class","agency","region","max_tmmx","max_tmmn","min_rmax","max_bi")
m2c <- glmulti(y="log(sdc)", 
        xr=potential_parms,
        data=d,level=1,method="h") #If running interactions (level 2), try genetic algorithm (method = "g")

summary(m2c@objects[[2]])

#Make table 1: Model summary
m_max <- 5
x=matrix(NA,nrow=length(c("AIC",rev(rownames(coef(m2c) ) ) ) ),ncol=m_max+1)
x[,1] <-c("AIC",sort(rownames(coef(m2c))))

for(m in c(1:m_max)){
  c <- round( summary(m2c@objects[[m]])$coefficients[,1]
              [order(names(summary(m2c@objects[[m]])$coefficients[,1]))], 3)
  x[which(x[,1]%in%names(c)),m+1] <- c
  x[1,m+1] <- round(AIC(m2c@objects[[m]]),3)
}
x=x[c(1,2,3,4,5,6,9,8,11,7,10),] #Order weather variables by their inclusion in the model
# Set up general table properties and formatting
cell_p = cellProperties(padding.right=3, padding.left=3)
par_p = parProperties(text.align="right")
# Create table
ft = FlexTable(x, header.columns=FALSE, body.cell.props=cell_p, body.par.props=par_p)
ft = addHeaderRow(ft, text.properties=textBold(), c("","Model #"),
                  colspan=c(1,m_max), par.properties=parCenter())
ft = addHeaderRow(ft, text.properties=textBold(), c("Model AIC \n /coefficients",c(1:m_max)),
                  colspan=rep(1,times=m_max+1), par.properties=parCenter())
ft #Save as HTML Table 1


####3b: Regression tree with best model from above####
library(rpart)
library(rpart.plot)
#The formula below is a good one; a simple model within one AIC point of the best model. I think  Tried adding a dummy variable for "region" but it wasn't appreciably better than this model.
#http://blog.revolutionanalytics.com/2013/06/plotting-classification-and-regression-trees-with-plotrpart.html
#Recode variables for regression tree
rtree.d <- d
names(rtree.d)[19] <- "max_high_temp"
levels(rtree.d$class) <- c("suppression","wildland fire use")
#formula(m2c@objects[[2]]) #Get formula and modify max_high_temp as below
tree.1 <- rpart(formula = "log(sdc) ~ 1 + class + agency + fire_year + max_high_temp", data=rtree.d)

split.labs <- function(x, labs, digits, varlen, faclen) {
  sapply(labs, function(lab) 
    if (grepl("fire_year", lab)) {
      rhs <- sub(".* ", "", lab);
      lab <- sub(rhs, ceiling(as.numeric(rhs))+1, lab) #+1 here is the modification
    } else lab)
} #Note: The default printing takes the fire_year split at 2010.5 (see print(tree.1)) and incorrectly rounds it to >=2010 instead of >=2011. The "+1" modification above is the only way I could figure out to make this round correctly.

png(file = paste0("./Manuscripts/MS2- Trends/Figures/Fig2_",Sys.Date(),".png"),width=4.5,height=4,units="in",res=500)

h <- #Histogram of all fires to put ln(SDC) in context
  ggplot(d)+
  geom_histogram(aes(log(sdc)), binwidth = 0.333,
                 fill=c("darkred","darkred",brewer.pal(10,"RdYlBu")),
                 col="black")+
  geom_vline(xintercept=c(-3.8,-5.1),col="black", lty=2)+
  labs(x = "ln(sdc)", title = "all fires")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 9),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 7))

prp(tree.1, varlen = 0, faclen = 0, type = 3, extra = 1,
    box.col = brewer.pal(10,"RdYlBu")[c(4,4,4,4,4,4,4,5,5,7,6,6,6,8,7)],
    #box.col is custom to match the histogram; values match tree.1$frame$yval
    #extra = 1, under = TRUE,
    #yesno.yshift=-1, #Doesn't apply if using type = 3
    clip.right.labs = FALSE, #Only applies if using type = 3
    mar = c(2,2,1,2), 
    split.fun = split.labs, main="                            ln(SDC)") # Plot the tree
#varlen: no limit to variable name length
#faclen: no limit to factor name length
#Type: How to place the labels (good types = 0,3)
#yesno.yshift: move yes/no labels further down

print(h, vp=viewport(.8, .25, .4, .45))
dev.off()


#Investigate the wierd max high temp > 39 result. It can be explained by the fact that all (18) of these fires were in the northwest, where complex topography can give more complex stand-replacing fire dynamics, and most (10) were in 1987 which was a warm year but where topography might have regulated stand-replacing effects.
#d_tmp <- d[d$class=="SUP" & d$max_tmmx >= 24 & d$fire_year<2011 & d$max_tmmx >= 39 & !is.na(d$fire_name),]
#d_tmp <- d_tmp[!is.na(d_tmp$fire_name),] 

####4. Trends over time####
####4a: Moving-window estimates (averaged per year); deprecate this sub-section eventually####

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
  labs(x=" ", y= "ln(sdc)",title="a")+
  annotate("text", x=2015, y=-3.4, size = 4, label = paste("R^2 == ", 0.11), parse = TRUE,hjust=1)+
  annotate("text", x=2015, y=-3.7, size = 4, label = paste("P = ", 0.058),hjust=1)+
  theme_bw()+
  theme(plot.title=element_text(size=14))
sdc_mvw <- 
  ggplot(d.annual,aes(x=fire_year,y=log_sdc_mvw5))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x=" ", y= "ln(sdc)",title="b       5-year mean value")+
  xlim(ggplot_build(sdc)$layout$panel_scales_x[[1]]$range$range) + 
  ylim(ggplot_build(sdc)$layout$panel_scales_y[[1]]$range$range) +
  annotate("text", x=2015, y=-3.4, size = 4, label = paste("R^2 == ", 0.14), parse = TRUE,hjust=1)+
  annotate("text", x=2015, y=-3.7, size = 4, label = paste("P = ", 0.047),hjust=1)+
  theme_bw()+
  theme(plot.title=element_text(size=14))
bi <- 
  ggplot(d.annual,aes(x=fire_year,y=max_bi))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x=" ", y= "burn index", title = "b")+
  annotate("text", x=2015, y=55, size = 4, label = paste("R^2 == ", 0.32), parse = TRUE,hjust=1)+
  annotate("text", x=2015, y=52, size = 4, label = paste("P = ", 0.001),hjust=1)+
  theme_bw()
bi_mvw <- 
  ggplot(d.annual,aes(x=fire_year,y=max_bi_mvw5))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x=" ", y= "burn index", title = "d")+
  xlim(ggplot_build(bi)$layout$panel_scales_x[[1]]$range$range) + 
  ylim(ggplot_build(bi)$layout$panel_scales_y[[1]]$range$range) +
  annotate("text", x=2015, y=55, size = 4, label = paste("R^2 == ", 0.69), parse = TRUE,hjust=1)+
  annotate("text", x=2015, y=52, size = 4, label = paste("P < ", 0.001),hjust=1)+
  theme_bw()
tmmx <- 
  ggplot(d.annual,aes(x=fire_year,y=max_tmmx))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x="fire year", y= "maximum temperature (C)", title = "c")+
  annotate("text", x=2015, y=26, size = 4, label = paste("R^2 == ", 0.10), parse = TRUE,hjust=1)+
  annotate("text", x=2015, y=24, size = 4, label = paste("P = ", 0.077),hjust=1)+
  theme_bw()
tmmx_mvw <- 
  ggplot(d.annual,aes(x=fire_year,y=max_tmmx_mvw5))+
  geom_point()+
  geom_smooth(method="lm")+
  #xlim(ggplot_build(tmmx)$layout$panel_ranges[[1]]$x.range) + 
  xlim(ggplot_build(tmmx)$layout$panel_scales_x[[1]]$range$range)+
  ylim(ggplot_build(tmmx)$layout$panel_scales_y[[1]]$range$range) +
  annotate("text", x=2015, y=25, size = 4, label = paste("R^2 == ", 0.29), parse = TRUE,hjust=1)+
  annotate("text", x=2015, y=24, size = 4, label = paste("P = ", 0.003),hjust=1)+
  labs(x="fire year", y= "maximum temperature (C)", title = "f")+
  theme_bw()

png(file = paste0("./Manuscripts/MS2- Trends/Figures/Fig4_",Sys.Date(),".png"),width=3.5,height=9,units="in",res=500)
grid.arrange(sdc,bi,tmmx,ncol=1)
dev.off()

#Changes over time in log(sdc)
#Significantly more negative in both the total data and the 5-year average
m <- lm(log(sdc)~fire_year,data=d.annual)
summary(m)
durbinWatsonTest(m)
m <- lm(log_sdc_mvw5~fire_year,data=d.annual)
summary(m)
durbinWatsonTest(m)


#Changes over time in burn index
#Significantly more positive
m <- lm(max_bi~fire_year,data=d.annual)
summary(m)
durbinWatsonTest(m)
m <- lm(max_bi_mvw5~fire_year,data=d.annual)
summary(m)
durbinWatsonTest(m)

#Changes over time in max temperature
#No significant trend (significant if you consider the 5-year average)
m <- lm(max_tmmx~fire_year,data=d.annual)
summary(m)
durbinWatsonTest(m)
m <- lm(max_tmmx_mvw5~fire_year,data=d.annual)
summary(m)
durbinWatsonTest(m)

####4a.2 Subset by region####
#We're not going to re-create Jay's 2012 analysis at all. The weighting was an interesting idea but not going to pursue it. Just report the regional split for all fires, tell Jay that if we just look at USFS fires we don't see the same trend.
#Maybe apply to suppression fires only (since that's the first split and where most of the action happens)
d.annual.region <-
  #group_by(d[d$agency=="USFS",],fire_year, region) %>% #Jay's original recommendatin; Nothing
  #group_by(d[d$agency!="NPS",],fire_year, region) %>% #Nothing
  #group_by(d[d$class=="SUP"&d$agency=="USFS",],fire_year, region) %>% #Moving window for SCSN marg. sig., positive trend over time
  #group_by(d[d$class=="SUP"&d$agency!="NPS",],fire_year, region) %>% #Nothing
  #group_by(d[d$agency=="USFS"&d$fire_year<2007,],fire_year, region) %>% #Trying to approximate Miller et al 2009; can't do it CHECKME 
  #group_by(d[d$agency=="USFS"&d$fire_year<2011&d$region=="SCSN",],fire_year, region) %>% #Trying to approximate Miller et al 2012; can't do it CHECKME
  group_by(d,fire_year, region) %>% #Go with this one in an appendix, SDC only. START HERE
  summarise(sdc=mean(sdc), log_sdc=mean(log(sdc)),  BA90_pct = mean(BA90_pct), 
            max_bi=mean(max_bi,na.rm=T),max_tmmx=mean(max_tmmx,na.rm=T)) #w = length(BA90_pct), deprecated.
log_sdc_mvw=rollapply(zoo(log(d.annual.region$sdc)), width = 5, by = 1, FUN = mean, align = "center")
max_bi_mvw=rollapply(zoo(d.annual.region$max_bi), width = 5, by = 1, FUN = mean, align = "center")
max_tmmx_mvw=rollapply(zoo(d.annual.region$max_tmmx), width = 5, by = 1, FUN = mean, align = "center")
BA90_pct_mvw=rollapply(zoo(d.annual.region$BA90_pct), width = 5, by = 1, FUN = mean, align = "center")
d.annual.region[rownames(as.data.frame(log_sdc_mvw)),"log_sdc_mvw5"]=log_sdc_mvw
d.annual.region[rownames(as.data.frame(max_bi_mvw)),"max_bi_mvw5"]=max_bi_mvw
d.annual.region[rownames(as.data.frame(max_tmmx_mvw)),"max_tmmx_mvw5"]=max_tmmx_mvw
d.annual.region[rownames(as.data.frame(max_tmmx_mvw)),"BA90_pct_mvw5"]=BA90_pct_mvw
#d.annual.region <- d.annual.region[-which(is.na(d.annual.region$fire_year)),] 

sdc_region <- 
  ggplot(d.annual.region,aes(x=fire_year,y=log_sdc,col=region))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x="fire year", y= "ln(sdc)")+
#  annotate("text", x=2015, y=-3.4, size = 4, label = paste("R^2 == ", 0.14), parse = TRUE,hjust=1)+
#  annotate("text", x=2015, y=-3.7, size = 4, label = paste("P = ", 0.047),hjust=1)+
  theme_bw()+
  theme(plot.title=element_text(size=14),legend.position = "none")
sdc_mvw_region <- 
  ggplot(d.annual.region,aes(x=fire_year,y=log_sdc_mvw5,col=region))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x="fire year", y= "ln(sdc)",title="b    5-year mean value")+
  xlim(ggplot_build(sdc_region)$layout$panel_scales_x[[1]]$range$range) + 
  ylim(ggplot_build(sdc_region)$layout$panel_scales_y[[1]]$range$range) +
#  annotate("text", x=2015, y=-3.4, size = 4, label = paste("R^2 == ", 0.14), parse = TRUE,hjust=1)+
#  annotate("text", x=2015, y=-3.7, size = 4, label = paste("P = ", 0.047),hjust=1)+
  theme_bw()+
  theme(plot.title=element_text(size=14))

BA90_region <- 
  ggplot(d.annual.region,aes(x=fire_year,y=BA90_pct,col=region))+
  geom_point()+
  geom_smooth(method="lm",aes(weight=w))+
  labs(x="fire year", y= "percent high-severity",title="annual mean value")+
  theme_bw()+
  theme(plot.title=element_text(size=14, hjust=0.5))

png(file = paste0("./Manuscripts/MS2- Trends/Figures/FigA3_",Sys.Date(),".png"),width=3,height=3,units="in",res=500)
sdc_region
dev.off()

summary(lm(log_sdc~fire_year,data=d.annual.region[d.annual.region$region=="NW",]))
summary(lm(log_sdc_mvw5~fire_year,data=d.annual.region[d.annual.region$region=="NW",]))
summary(lm(log_sdc~fire_year,data=d.annual.region[d.annual.region$region=="SCSN",]))
summary(lm(log_sdc_mvw5~fire_year,data=d.annual.region[d.annual.region$region=="SCSN",]))

summary(lm(BA90_pct~fire_year,data=d.annual.region[d.annual.region$region=="NW",]))
summary(lm(BA90_pct_mvw5~fire_year,data=d.annual.region[d.annual.region$region=="NW",]))
summary(lm(BA90_pct~fire_year,data=d.annual.region[d.annual.region$region=="SCSN",]))
summary(lm(BA90_pct_mvw5~fire_year,data=d.annual.region[d.annual.region$region=="SCSN",]))


####4b: percentile analysis CHECKME this part is deprecated####
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


####5a. Identify specific example fires####
#Run the chunk below for BA_range <- c(0,10), c(10,20), c(20,30), c(30-40), and c(45-55). Only needed to run once; skip to 5b.
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
#d[which(d$Plot_candidates),"Plot_candidates"] <- c(5,2,3,1,7,8,6,4)
#Above line only run once; candidates now determined so now skip to 5b

####5b.Identify  specific example fires####
#Gives a nice comparison of 8 example fires.
d[as.character(d$VB_ID) %in% 
    c("2001HIGHWAY","1990RECER","1999HORTON2","2009GOOSE","1989RACK","2006BOULDER_CMPLX",
      "2001STREAM","2004SIMS"),"Plot_candidates"] <-  c(5,2,3,1,7,8,6,4)
d[which(is.na(d$Plot_candidates)),"Plot_candidates"] <- ""
####5c.Plot  specific example fires####
#hs_patches <- readOGR("../Large Files/GIS/BurnSev/Current/", layer="hs_patches") #CRS EPSG:3310, NAD83 CA Albers
full_shapes <- readOGR("../Large Files/GIS/BurnSev/Current/", layer="SampleFiresFullPerims") #CRS EPSG:3310, NAD83 CA Albers

fires.to.plot.names=as.character(d[ order(d[,"Plot_candidates"])[470:477] , "VB_ID"])
pd.a <- fortify(full_shapes[full_shapes$VB_ID==fires.to.plot.names[1],])
p.a=
  ggplot()+
  geom_polygon(data=pd.a,
               aes(x=long-min(long),y=lat-min(lat),group=group),fill="gray",col="gray")+
  geom_polygon(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[1],]),
               aes(x=long-min(pd.a$long),y=lat-min(pd.a$lat),group=group),col='darkred',fill='darkred')+
  xlim(0,10000) + ylim(0,10000) + coord_fixed()+
  labs(title="1. Highway Fire (2001)",x=" ",y="meters") +
  annotate("text",x=0,y=8000,label=
             paste("sdc =",round(d[d$VB_ID==fires.to.plot.names[1],"sdc"],4),
                   "\nln(sdc) = ",round(log(d[d$VB_ID==fires.to.plot.names[1],"sdc"]),2) ),
           hjust=0,vjust=0,size=3)+
  theme_bw()+
  theme(axis.text.x=element_text(size=12,angle = 45, hjust = 1),
        axis.text.y=element_text(size=12,angle = 0),
        axis.title=element_text(size=13),
        plot.title=element_text(size=14, hjust=0.5))

pd.b <- fortify(full_shapes[full_shapes$VB_ID==fires.to.plot.names[2],])
p.b=ggplot()+
  geom_polygon(data=pd.b,
               aes(x=long-min(long),y=lat-min(lat),group=group),fill="gray",col="gray")+
  geom_polygon(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[2],]),
               aes(x=long-min(pd.b$long),y=lat-min(pd.b$lat),group=group),col='darkred',fill='darkred')+
  xlim(0,10000) + ylim(0,10000) + coord_fixed()+
  labs(title="2. Recer Fire (1990)",x=" ",y=" ") +
  annotate("text",x=0,y=8000,label=
             paste("sdc =",round(d[d$VB_ID==fires.to.plot.names[2],"sdc"],4),
                   "\nln(sdc) = ",round(log(d[d$VB_ID==fires.to.plot.names[2],"sdc"]),2) ),
           hjust=0,vjust=0,size=3)+
  theme_bw()+
  theme(axis.text.x=element_text(size=12,angle = 45, hjust = 1),
        axis.text.y=element_text(size=12,angle = 0),
        axis.title=element_text(size=13),
        plot.title=element_text(size=14, hjust=0.5))

pd.c <- fortify(full_shapes[full_shapes$VB_ID==fires.to.plot.names[3],])
p.c=ggplot()+
  geom_polygon(data=pd.c,
               aes(x=long-min(long),y=lat-min(lat),group=group),fill="gray",col="gray")+
  geom_polygon(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[3],]),
               aes(x=long-min(pd.c$long),y=lat-min(pd.c$lat),group=group),col='darkred',fill='darkred')+
  xlim(0,10000) + ylim(0,10000) + coord_fixed()+
  labs(title="3. Horton Fire (1999)",x=" ",y="meters") +
  annotate("text",x=0,y=8000,label=
             paste("sdc =",round(d[d$VB_ID==fires.to.plot.names[3],"sdc"],4),
                   "\nln(sdc) = ",round(log(d[d$VB_ID==fires.to.plot.names[3],"sdc"]),2) ),
           hjust=0,vjust=0,size=3)+
  theme_bw()+
  theme(axis.text.x=element_text(size=12,angle = 45, hjust = 1),
        axis.text.y=element_text(size=12,angle = 0),
        axis.title=element_text(size=13),
        plot.title=element_text(size=14, hjust=0.5))

pd.d <- fortify(full_shapes[full_shapes$VB_ID==fires.to.plot.names[4],])
p.d=ggplot()+
  geom_polygon(data=pd.d,
               aes(x=long-min(long),y=lat-min(lat),group=group),fill="gray",col="gray")+
  geom_polygon(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[4],]),
               aes(x=long-min(pd.d$long),y=lat-min(pd.d$lat),group=group),col='darkred',fill='darkred')+
  xlim(0,10000) + ylim(0,10000) + coord_fixed()+
  labs(title="4. Goose Fire (2009)",x=" ",y=" ") +
  annotate("text",x=0,y=8000,label=
             paste("sdc =",round(d[d$VB_ID==fires.to.plot.names[4],"sdc"],4),
                   "\nln(sdc) = ",round(log(d[d$VB_ID==fires.to.plot.names[4],"sdc"]),2) ),
           hjust=0,vjust=0,size=3)+
  theme_bw()+
  theme(axis.text.x=element_text(size=12,angle = 45, hjust = 1),
        axis.text.y=element_text(size=12,angle = 0),
        axis.title=element_text(size=13),
        plot.title=element_text(size=14, hjust=0.5))

pd.e <- fortify(full_shapes[full_shapes$VB_ID==fires.to.plot.names[5],])
p.e=ggplot()+
  geom_polygon(data=pd.e,
               aes(x=long-min(long),y=lat-min(lat),group=group),fill="gray",col="gray")+
  geom_polygon(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[5],]),
               aes(x=long-min(pd.e$long),y=lat-min(pd.e$lat),group=group),col='darkred',fill='darkred')+
  xlim(0,10000) + ylim(0,10000) + coord_fixed()+
  labs(title="5. Rack Fire (1989)",x=" ",y="meters") +
  annotate("text",x=0,y=8000,label=
             paste("sdc =",round(d[d$VB_ID==fires.to.plot.names[5],"sdc"],4),
                   "\nln(sdc) = ",round(log(d[d$VB_ID==fires.to.plot.names[5],"sdc"]),2) ),
           hjust=0,vjust=0,size=3)+
  theme_bw()+
  theme(axis.text.x=element_text(size=12,angle = 45, hjust = 1),
        axis.text.y=element_text(size=12,angle = 0),
        axis.title=element_text(size=13),
        plot.title=element_text(size=14, hjust=0.5))

pd.f <- fortify(full_shapes[full_shapes$VB_ID==fires.to.plot.names[6],])
p.f=ggplot()+
  geom_polygon(data=pd.f,
               aes(x=long-min(long),y=lat-min(lat),group=group),fill="gray",col="gray")+
  geom_polygon(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[6],]),
               aes(x=long-min(pd.f$long),y=lat-min(pd.f$lat),group=group),col='darkred',fill='darkred')+
  xlim(0,10000) + ylim(0,10000) + coord_fixed()+
  labs(title="6. Boulder Complex (2006)",x=" ",y=" ") +
  annotate("text",x=0,y=4000,label=
             paste("sdc =",round(d[d$VB_ID==fires.to.plot.names[6],"sdc"],4),
                   "\nln(sdc) = ",round(log(d[d$VB_ID==fires.to.plot.names[6],"sdc"]),2) ),
           hjust=0,vjust=0,size=3)+
  theme_bw()+
  theme(axis.text.x=element_text(size=12,angle = 45, hjust = 1),
        axis.text.y=element_text(size=12,angle = 0),
        axis.title=element_text(size=13),
        plot.title=element_text(size=14, hjust=0.5))

pd.g <- fortify(full_shapes[full_shapes$VB_ID==fires.to.plot.names[7],])
p.g=ggplot()+
  geom_polygon(data=pd.g,
               aes(x=long-min(long),y=lat-min(lat),group=group),fill="gray",col="gray")+
  geom_polygon(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[7],]),
               aes(x=long-min(pd.g$long),y=lat-min(pd.g$lat),group=group),col='darkred',fill='darkred')+
  xlim(0,10000) + ylim(0,10000) + coord_fixed()+
  labs(title="7. Stream Fire (2001)",x="meters",y="meters") +
  annotate("text",x=0,y=8000,label=
             paste("sdc =",round(d[d$VB_ID==fires.to.plot.names[7],"sdc"],4),
                   "\nln(sdc) = ",round(log(d[d$VB_ID==fires.to.plot.names[7],"sdc"]),2) ),
           hjust=0,vjust=0,size=3)+
  theme_bw()+
  theme(axis.text.x=element_text(size=12,angle = 45, hjust = 1),
        axis.text.y=element_text(size=12,angle = 0),
        axis.title=element_text(size=13),
        plot.title=element_text(size=14, hjust=0.5))

pd.h <- fortify(full_shapes[full_shapes$VB_ID==fires.to.plot.names[8],])
p.h=ggplot()+
  geom_polygon(data=pd.h,
               aes(x=long-min(long),y=lat-min(lat),group=group),fill="gray",col="gray")+
  geom_polygon(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[8],]),
               aes(x=long-min(pd.h$long),y=lat-min(pd.h$lat),group=group),col='darkred',fill='darkred')+
  xlim(0,10000) + ylim(0,10000) + coord_fixed()+
  labs(title="8. Sims Fire (2004)",x="meters",y=" ") +
  annotate("text",x=0,y=8000,label=
             paste("sdc =",round(d[d$VB_ID==fires.to.plot.names[8],"sdc"],4),
                   "\nln(sdc) = ",round(log(d[d$VB_ID==fires.to.plot.names[8],"sdc"]),2) ),
           hjust=0,vjust=0,size=3)+
  theme_bw()+
  theme(axis.text.x=element_text(size=12,angle = 45, hjust = 1),
        axis.text.y=element_text(size=12,angle = 0),
        axis.title=element_text(size=13),
        plot.title=element_text(size=14, hjust=0.5))

png(file = paste0("./Manuscripts/MS2- Trends/Figures/Fig6_",Sys.Date(),".png"),width=6,height=10,units="in",res=500)
grid.arrange(p.a,p.b,p.c,p.d,p.e,p.f,p.g,p.h,ncol=2)
dev.off()

#Export shapefile with just the ones we're using:
#hs_used <- hs_patches[hs_patches$VB_ID%in%d$VB_ID,]
#writeOGR(hs_used, "../Large Files/GIS/BurnSev/Current","hs_patches_used", driver = "ESRI Shapefile") #CRS EPSG:3310, NAD83 CA Albers

####6: Plot SDC vs class, agency, pct_hs, and fire size#### 
#Make sure you've run #5b first.
#Trouble with plotting NA's as gray
#df <- data.frame (V1=factor(c("A","B","A","B",NA)),x=c(1:5),y=c(1:5))
#ggplot(df,aes(x=x,y=y,col=V1))+
#  geom_point()

class_pct <-
  ggplot(na.omit(d[,c("BA90_pct","sdc","class","Plot_candidates")]),aes(x=BA90_pct,y=log(sdc),col=class))+
  geom_point()+
  geom_smooth(method="lm")+
  geom_text(aes(label=Plot_candidates),hjust=1, vjust=1,size=6,col="black",fontface="bold")+
  labs(y= "ln(sdc)",x=" ", title = "            percent high-severity\na")+
  annotate("text", x=60, y=-3, label = paste("R^2 == ", 0.67), parse = TRUE) +
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13),
        legend.position = 'none')
class_ha <-
  ggplot(na.omit(d[,c("firesize_ha","sdc","class")]),aes(x=log(firesize_ha),y=log(sdc),col=class))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(y= " ", x=" ", title = "                       fire area\nb")+
  annotate("text", x=10, y=-3, label = paste("R^2 == ", 0.22), parse = TRUE) +
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13))
agency_pct <-
  ggplot(na.omit(d[,c("BA90_pct","sdc","agency")]),aes(x=BA90_pct,y=log(sdc),col=agency))+
  geom_point()+
  scale_color_manual(values=c("darkred","darkgreen","orange"),guide=FALSE)+
  geom_smooth(method="lm")+
  labs(y= "ln(sdc)", x="% high-severity", title = "c")+
  annotate("text", x=60, y=-3, label = paste("R^2 == ", 0.67), parse = TRUE) +
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13))
agency_ha <-
  ggplot(na.omit(d[,c("firesize_ha","sdc","agency")]),aes(x=log(firesize_ha),y=log(sdc),col=agency))+
  geom_point()+
  scale_color_manual(values=c("darkred","darkgreen","orange"))+
  geom_smooth(method="lm")+
  labs(y= " ", x= "ln(area [ha])", title = "d")+
  annotate("text", x=10, y=-3, label = paste("R^2 == ", 0.21), parse = TRUE) +
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13))

#summary(lm(log(sdc)~BA90_pct+class,data=d))
#summary(lm(log(sdc)~firesize_ha+class,data=d))
#summary(lm(log(sdc)~BA90_pct+agency,data=d))
#summary(lm(log(sdc)~BA90_pct+agency,data=within(d, agency <- relevel(agency, ref = 3) ) ) )
#summary(lm(log(sdc)~firesize_ha+agency,data=d))
#summary(lm(log(sdc)~firesize_ha+agency,data=within(d, agency <- relevel(agency, ref = 3) ) ) )

png(file = paste0("./Manuscripts/MS2- Trends/Figures/Fig3_",Sys.Date(),".png"),width=7.5,height=10,units="in",res=500)
grid.arrange(class_pct,class_ha,agency_pct,agency_ha,ncol=2,widths=c(0.44,0.56))
dev.off()



####7: Forest loss by agency####
#Calculate core patch area 
d$P <- 1/(10^(d$sdc*120))
d$CPA_120 <- d$BA90_ha * d$P

d_CPA <-
  d %>%
  group_by(fire_year,agency) %>%
  summarise(CPA_120 = sum(CPA_120)) %>%
  group_by(agency)%>%
  mutate(CPA_120_cumul = cumsum(CPA_120))
d_CPA <- d_CPA[complete.cases(d_CPA),]
  
p_CPA <-
  ggplot(d_CPA,aes(x=fire_year,y=CPA_120_cumul,col=agency))+
  geom_line()+
  scale_color_manual(values=c("darkred","darkgreen","orange"))+
  labs(x = "fire year", y = "cumulative area >120 m in from patch edge (ha)")+
  scale_x_continuous(minor_breaks = seq(1984, 2015, 1),breaks = seq(1985, 2015, 5))+
  theme_bw()+
  theme(panel.grid.minor = element_line(size=0.5),
        panel.grid.major.x = element_line(color = "black"))

area_burned <- 
  group_by(d,agency) %>%
  summarise(area_burned=sum(firesize_ha))
d_CPA[d_CPA$fire_year==2015,"CPA_120_cumul"]/area_burned[c(1:3),"area_burned"]

png(file = paste0("./Manuscripts/MS2- Trends/Figures/Fig6_",Sys.Date(),".png"),width=7.5,height=4,units="in",res=500)
p_CPA
dev.off()
  