library(shiny) #originally 1.0.3
library(rgdal) # 1.2.5
library(parallel) #for "mclapply" version 3.3.2
library(rgeos) #for gBuffer version 0.3.21
library(raster) #for "area" version 2.5.8
library(ggplot2) #version 2.2.1
library(RCurl) #version 1.95.4.8


hs_fire <- 
  readOGR("../Large Files/GIS/BurnSev/Specific Fires/HS Layers/", layer="2013American_CT")
hs_fire <- 
  readOGR("../Large Files/GIS/BurnSev/Specific Fires/HS Layers/", layer="2013AMERICAN")

##Fill holes:
fill_holes=function(hs_fire){
  hs_fire_p=slot(hs_fire, "polygons")
  holes <- lapply(hs_fire_p, function(x) sapply(slot(x, "Polygons"), slot, "hole"))
  areas <- lapply(hs_fire_p, function(x) sapply(slot(x, "Polygons"), slot, "area"))
  res <- lapply(1:length(hs_fire_p), 
                function(i) slot(hs_fire_p[[i]], "Polygons")[!(holes[[1]]&areas[[1]]<8100)]
  )
  #Select the polygons that are not holes. The "i" here is an artifact of the example code; the fires here only have one polygon ID so i=1 (it's a multipart polygon).
  IDs <- row.names(hs_fire)
  hs_fire_fill <- SpatialPolygons(lapply(1:length(res), function(i)
    Polygons(res[[i]], ID=IDs[i])), proj4string=CRS(proj4string(hs_fire)))
  return(hs_fire_fill)
  ##One consequence of this is that it's no longer a sp data frame, and there are some warnings. 
}

##Decay profile: 
##Implement internal buffering on high severity patches, create "Long" dataset
decay=function(hs_fire,buf_max=1000,buf_inc=input$buf_inc,name="your_fire",cancel=F){
  #Set up long data frame:
  dist.table.sub <-
    data.frame(name=rep(name,(buf_max/buf_inc)+1),
               width=seq(0,buf_max,by=buf_inc),
               area_ha=NA)
  #Core operation:
  buf.list <- #Apply the buffer at every width in X; returns list of length X.
    mclapply(X=-dist.table.sub$width[2:length(dist.table.sub$width)],
             FUN=gBuffer,
             spgeom=hs_fire, 
             byid = FALSE, id = NULL, quadsegs = 5, 
             capStyle = "ROUND", joinStyle = "ROUND", mitreLimit = 1,
             mc.cores = 4) 
  
  #Post-processing of long data frame:
  buf.list <- #Remove any NULL values (when buffer eliminated all HS areas)
    Filter(length,buf.list) 
  dist.table.sub$area_ha <- #Calculate area remaining at each buffer distance
    #0.0001 converts m2 to ha
    (as.vector(sapply(buf.list,raster::area))*0.0001)[1:nrow(dist.table.sub)] 
  if(any(is.na(dist.table.sub$area_ha))){
    #If there are NULL values for the area because the buffer got too wide, 
    #remove them. Set the first null value to 0.
    dist.table.sub$area_ha[which(is.na(dist.table.sub$area_ha))[1]] <- 0 
    if(any(is.na(dist.table.sub$area_ha))){
      #If there are STILL NULL values, delete the rest
      dist.table.sub <- #Delete the rest of the table rows with NA in area.
        dist.table.sub[-which(is.na(dist.table.sub$area_ha)),] 
    } 
  }
  dist.table.sub$prop.hs <- 
    #Calculate the proportion of the original HS area remaining,
    #at each buffer distance
    dist.table.sub$area_ha/dist.table.sub$area_ha[1]
  
  #Return the long data frame for the fire in question
  return(dist.table.sub)
}

##Calculate stand-replacing decay coefficient (sdc)
calculate.sdc=function(decay.table){
  #https://stat.ethz.ch/R-manual/R-devel/library/base/html/by.html
  m.list=with(decay.table,
              by(decay.table,name,
                 function(x)
                   nls(prop.hs~1/(10^(sdc*width)),data=x,start=list(sdc=0.01))
              )
  )
  sdc.table=data.frame(name=unique(as.character(decay.table$name)),
                       sdc=sapply(m.list,coef))
  out.table=merge(decay.table,sdc.table,sort=F)
  out.table$sdc.name=as.character(format(round(out.table$sdc,4),scientific=F))
  return(out.table)
}

#Don't copy the below code back to app.R; this just implements the functions.
hs_fire <- fill_holes(hs_fire)

input = data.frame(buf_inc = 100)
Sys.time()
decay.table <- decay(hs_fire)
Sys.time()
sdc.table <- calculate.sdc(decay.table) 
unique(sdc.table$sdc.name)
#print(paste0("your sdc = ",unique(sdc$sdc.name)))

#Time comparisons for 2014WHITES_HS.shp
#1. 10m buffer with no holes, 3:27; sdc = 0.0035
#2. 10m buffer with holes, 4:22
#3. 20m buffer with no holes, 1:42; sdc = 0.0035
#3. 30m buffer with no holes, 1:06; sdc = 0.0035
#4. 40m buffer with no holes, 1:03; sdc = 0.0.0035
#5. 50m buffer with no holes, 0:42; sdc = 0.0035
#6. 100m buffer with no holes, 0:23; sdc = 0.0039