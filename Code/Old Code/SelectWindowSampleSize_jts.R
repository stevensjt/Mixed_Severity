library(maptools)
library(sp)
library(raster)
library(spatstat)
library(plyr)
library(party)
library(caret)

setwd ("C:\\JamieProjects\\Rim_JFSP\\Treatments")

##load files using readShapePoly in package maptools for polygons and raster command in the raster package for rasters

#mask area to select random points within
mask2500 <- readShapePoly("Perim_2295")

#Severity class raster
Severity <- raster("sevcl_adj")

#Treatment class and age rasters
Trts <- raster("trt_reclass")
TrtAge <- raster("trtage_cl")

#BI and ERC rasters
BI <- raster("BI")
ERC <- raster("ERC")

#AET and water deficit rasters
AET <- raster("aet_30")
Def <- raster("def_30")

#Plantation age class
#this is a reclassified raster where plantations <=20 have value 1, >20 have value 2 and not plantation has value 3
plantation <- raster("plantnclass")

#vegetation raster
#reclassified Landfire EVT by EVT_PHYS (see notes in "Steps for preliminary moving window analysis.docx")
veg <- raster("evt_reclass")


####### Some kind of loop for n would be ideal here!!!! #######

#Set number of random points
#I've been running these one at a time
n <-50
n <-60
n <-70
n <-80
n <-90
n <-100
n <-110
n <-120


##JTS Code
sample.sizes=seq(from=50,to=120,by=10) #You can enter manually as c(50,60,70...etc) but this is faster. A sequence of candidate sample sizes ranging from 50 to 120 in increments of 10.
df_n=RF_vars_n=list() #Storage lists can hold an output data frame for each n (used below).
for(n in sample.sizes){ #Run the code chunk between the {brackets} for each unique value of n.
##End JTS Code
  
#generate random points using spsample from sp package
ptsn_2500 <- spsample(mask2500, n, "random")

#OR if want a minimum distance set between points use the rSSI function from the spatstat package:
#ptsn_2500_m <- SpatialPoints(coords(rSSI(2000, n=10, win = mask2500)))

####extract values for windows using package 'raster'######

#Percent in each severity class
Ext_Sev <- extract(Severity, ptsn_2500, buffer=1795)
counts_Sev <- lapply(Ext_Sev,table)
pct_Sev <- lapply(counts_Sev, FUN=function(x){ x / sum(x) } ) 

#convert to data.frame using plyr package
non.null.list <- lapply(pct_Sev, lapply, function(x)ifelse(is.null(x), NA, x))
pct_Sev.df <-rbind.fill(lapply(non.null.list, as.data.frame))
#replace NA values with 0
pct_Sev.df[is.na(pct_Sev.df)] <- 0

Pct_H <- pct_Sev.df[,"X4"]


##Percent in each treatment class
Ext_Trts <- extract(Trts, ptsn_2500, buffer=1795)
counts_Trts <- lapply(Ext_Trts,table)
pct_Trts <- lapply(counts_Trts, FUN=function(x){ x / sum(x) } ) 

#convert list to data frame 
non.null.list <- lapply(pct_Trts, lapply, function(x)ifelse(is.null(x), NA, x))
pct_Trts.df <-rbind.fill(lapply(non.null.list, as.data.frame))
pct_Trts.df[is.na(pct_Trts.df)] <- 0

#calculate total % area treated within each window
Pct_TrtTot <- 1-pct_Trts.df$X9

#calculate % treated excluding High severity
Pct_Trt_noH <- Pct_TrtTot-pct_Trts.df$X4

#create matrix with only proportion values for treatments 
TrtOnly <- subset(pct_Trts.df, select = -X9)

#calculate maximum treatment type for windows with at least 25% or 50% treated area
TrtMaj50 <- ifelse(pct_Trts.df$X9 < 0.5, apply(TrtOnly,1,FUN=function(x){names(which.max(x))}), "X9")
TrtMaj25 <- ifelse(pct_Trts.df$X9 < 0.75, apply(TrtOnly,1,FUN=function(x){names(which.max(x))}), "X9")


#Calculate majority treatment age for windows with at least 25% or 50% treated
Ext_Age <- extract(TrtAge, ptsn_2500, buffer=1795)
counts_Age <- lapply(Ext_Age,table)
pct_Age <- lapply(counts_Age, FUN=function(x){ x / sum(x) } ) 

non.null.list <- lapply(pct_Age, lapply, function(x)ifelse(is.null(x), NA, x))
pct_Age.df <-rbind.fill(lapply(non.null.list, as.data.frame))
pct_Age.df[is.na(pct_Age.df)] <- 0

TrtAgeOnly <- subset(pct_Age.df, select = -X100)

##change column names to numeric
names(TrtAgeOnly) <- gsub("X", "", names(TrtAgeOnly))

AgeMaj50 <- as.numeric(ifelse(pct_Age.df$X100 < 0.5, apply(TrtAgeOnly,1,FUN=function(x){names(which.max(x))}), 100))
AgeMaj25 <- as.numeric(ifelse(pct_Age.df$X100 < 0.75, apply(TrtAgeOnly,1,FUN=function(x){names(which.max(x))}), 100))

##Extract vegetation variables

#Percent plantation <= 20 years, and percent plantation >20 years
Ext_pl <- extract(plantation, ptsn_2500, buffer=1795)
counts_pl <- lapply(Ext_pl,table)
pct_pl <- lapply(counts_pl, FUN=function(x){ x / sum(x) } ) 

non.null.list <- lapply(pct_pl, lapply, function(x)ifelse(is.null(x), NA, x))
pct_pl.df <-rbind.fill(lapply(non.null.list, as.data.frame))
pct_pl.df[is.na(pct_pl.df)] <- 0

Plant_lte20 <- pct_pl.df[,c("X1")]
Plant_gt20 <- pct_pl.df[,c("X2")]


#majority veg type and % in classes of interest

Ext_veg <- extract(veg, ptsn_2500, buffer=1795)
counts_veg <- lapply(Ext_veg,table)
pct_veg <- lapply(counts_veg, FUN=function(x){ x / sum(x) } ) 

non.null.list <- lapply(pct_veg, lapply, function(x)ifelse(is.null(x), NA, x))
pct_veg.df <-rbind.fill(lapply(non.null.list, as.data.frame))
pct_veg.df[is.na(pct_veg.df)] <- 0

#Majority veg type
Maj_veg <- apply(pct_veg.df,1,FUN=function(x){names(which.max(x))})

#Proportion barren, open water and sparsely vegetated
Pct_BOS <- pct_veg.df[,c("X2")]

#Proportion shrubland
Pct_Shrub <- pct_veg.df[,c("X4")]

##Extract mean BI, ERC, AET and water deficit for each window
Mean_BI <- extract(BI, ptsn_2500, buffer=1795, fun=mean)
Mean_ERC <- extract(ERC, ptsn_2500, buffer=1795, fun=mean)
Mean_AET <- extract(AET, ptsn_2500, buffer=1795, fun=mean)
Mean_Def <- extract(Def, ptsn_2500, buffer=1795, fun=mean)


##combine variables in a dataframe and run random forests

##JTS Code
##At this point it seems like you're done doing your extractions. It might make sense to end the loop here and store each resultant data frame as an individual list element. Note that above I defined df_n and RF_vars_n as lists.
df_n[[as.character(n)]] <- data.frame(ptsn_2500, Pct_H, Pct_TrtTot, Pct_Trt_noH, TrtMaj50, TrtMaj25, AgeMaj50, AgeMaj25, 
		 Mean_BI, Mean_ERC, Mean_AET, Mean_Def, Plant_lte20, Plant_gt20, Maj_veg, Pct_BOS,Pct_Shrub)

RF_vars_n[[as.character(n)]] <- data.frame(Pct_Trt_noH, TrtMaj50, TrtMaj25, AgeMaj50, AgeMaj25, 
		 Mean_BI, Mean_ERC, Mean_AET, Mean_Def, Plant_lte20, Plant_gt20, Maj_veg, Pct_BOS,Pct_Shrub)

##Note the "as.character(n)" addition above names each list element with a character string version of the sample size; if you just do df_n[[n]], and the first n is 50, then you'll have 49 empty list elements until you get to the first data frame you want. It's easy to toggle back and forth from characters to numbers if you need to- the converse function is as.numeric(). In any case, once you have the data frames, you can run the RF analysis below manually, or put it in its own for loop following the syntax I've used here and where you call "data=RF_vars_n", just change it to RF_vars_n[[as.character(n)]], for example.
}##Close the for loop
##End JTS Code

#select random start seed
r <- sample(1:n,1)
set.seed(r)

#run RF
forest2500_n = cforest(Pct_H~ ., data=RF_vars_n, controls= cforest_unbiased (mtry = 5, ntree=1000))

#Output Error
error_n <- cforestStats(forest2500_n)

#compile error from each run
e <- data.frame(error_50,error_60,error_70,error_80,error_90,error_100,error_110,error_120)
