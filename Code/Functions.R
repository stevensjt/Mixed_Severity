##This code contains functions for recharacterizing fire regimes##

##Fill holes that are less than 9 pixels large (9*900m2=8100m2, or 0.81 ha)
fill_holes=function(hs_fire){
  hs_fire_p=slot(hs_fire, "polygons")
  holes <- lapply(hs_fire_p, function(x) sapply(slot(x, "Polygons"), slot, "hole"))
  areas <- lapply(hs_fire_p, function(x) sapply(slot(x, "Polygons"), slot, "area"))
  res <- lapply(1:length(hs_fire_p), function(i) slot(hs_fire_p[[i]], "Polygons")[!(holes[[1]]&areas[[1]]<8100)])#Select the polygons that are not holes. The "i" here is an artifact of the example code; the fires here only have one polygon ID so i=1 (it's a multipart polygon).
  IDs <- row.names(hs_fire)
  hs_fire_fill <- SpatialPolygons(lapply(1:length(res), function(i)
    Polygons(res[[i]], ID=IDs[i])), proj4string=CRS(proj4string(hs_fire)))
  return(hs_fire_fill)
  ##One consequence of this is that it's no longer a sp data frame, and there are some warnings. 
}

##Decay profile: Implement internal buffering on high severity patches, create "Long" dataset
decay=function(hs_fire,buf_max,buf_inc,name,cancel=F){
  require(parallel) #for "mclapply"
  require(rgeos) #for gBuffer
  require(raster) #for "area"
  
  #Set up long data frame:
  dist.table.sub <-
    data.frame(name=rep(name,(buf_max/buf_inc)+1),width=seq(0,buf_max,by=buf_inc),area_ha=NA)
  
  #Core operation:
  buf.list <- #Apply the buffer at every width in X; returns list of length X.
    mclapply(X=-dist.table.sub$width,
             FUN=gBuffer,
             spgeom=hs_fire, 
             byid = FALSE, id = NULL, quadsegs = 5, capStyle = "ROUND", joinStyle = "ROUND", mitreLimit = 1,
             mc.cores = 4) 
  
  #Post-processing of long data frame:
  buf.list <- #Remove any NULL values (when buffer eliminated all HS areas)
    Filter(length,buf.list) 
  dist.table.sub$area_ha <- #Calculate area remaining at each buffer distance
    (as.vector(sapply(buf.list,area))*0.0001)[1:nrow(dist.table.sub)] #0.0001 converts m2 to ha
  if(any(is.na(dist.table.sub$area_ha))){
    #If there are NULL values for the area because the buffer got too wide, remove them.
    dist.table.sub$area_ha[which(is.na(dist.table.sub$area_ha))[1]] <- 0 #Set the first null value to 0
    dist.table.sub <- #Delete the rest of the table rows with NA in area.
      dist.table.sub[-which(is.na(dist.table.sub$area_ha)),] 
  }
  dist.table.sub$prop.hs <- #Calculate the proportion of the original HS area remaining at each buffer distance
    dist.table.sub$area_ha/dist.table.sub$area_ha[1]
  
  #Return the long data frame for the fire in question
  return(dist.table.sub)
}
#


##Calculate stand-replacing decay coefficient (sdc)
calculate.sdc=function(decay.table){
  #https://stat.ethz.ch/R-manual/R-devel/library/base/html/by.html
  m.list=with(decay.table,
               by(decay.table,name,
                  function(x)
                    nls(prop.hs~1/(10^(sdc*width)),data=x,start=list(sdc=0.01))
                  )
               )
  sdc.table=data.frame(name=unique(as.character(decay.table$name)),sdc=sapply(m.list,coef))
  out.table=merge(decay.table,sdc.table,sort=F)
  out.table$sdc.name=as.character(format(round(out.table$sdc,4),scientific=F))
  return(out.table)
}

##Calculate fragstats: mean shape index (MSI), area-weighted mean shape index (AWMSI), mean patch fractal dimension (MPFD) and area-weighted mean patch fractal dimension (AWMPFD). 

get_fragstats <- function (hs_fire){
  EAR=data.frame(poly=NA,perim=NA,area=NA)
  holes=which(sapply(hs_fire@polygons[[1]]@Polygons , slot , "hole"))
  for(p in 1:length(hs_fire@polygons[[1]]@Polygons)){
    if(!p%in%holes){
      #If target polygon not a hole, get perimeter and area, then check for holes within
      df=hs_fire@polygons[[1]]@Polygons[[p]]@coords
      out <- sapply(1:(nrow(df)-1), function(i) {
        d <- dist(df[i:(i+1),1:2])
        c(d)
      })
      EAR[p,"poly"]=p
      EAR[p,"perim"] <-sum(out)
      EAR[p,"area"] <- hs_fire@polygons[[1]]@Polygons[[p]]@area
      
      if((p+1)%in%holes){
        #If target polygon isnt a hole but next polygon is, target is parent. Store info.
        parent.p=p
        coords.storage=hs_fire@polygons[[1]]@Polygons[[p]]@coords
      }
    } else {
      #If target polygon IS a hole (works if there are multiple holes in a row)
      df=hs_fire@polygons[[1]]@Polygons[[p]]@coords
      out <- sapply(1:(nrow(df)-1), function(i) {
        d <- dist(df[i:(i+1),1:2])
        c(d)
      })
      EAR[parent.p,"perim"] <- EAR[parent.p,"perim"] + sum(out) 
      #Add hole's perimeter to parent's perimeter
      
      EAR[parent.p,"area"] <- 
        EAR[parent.p,"area"] - hs_fire@polygons[[1]]@Polygons[[p]]@area
      #Subtract hole's area from parent's area
    }
  }
  EAR$ear <- EAR$perim/EAR$area #Edge:area ratio (= perimeter:area ration = Shape Index)
  EAR$pfd <- 2*log(EAR$perim)/log(EAR$area) #Patch fractal dimension; equation from McGarigal and Marks 1995 p 36
  MSI <- mean(EAR$ear,na.rm=T) #Mean Shape Index
  AWMSI <- weighted.mean(EAR$ear,w = EAR$area,na.rm=T) #Area-weighted Mean Shape Index
  MPFD <- mean(EAR$pfd,na.rm=T) #Mean Patch Fractal Dimension
  AWMPFD <- weighted.mean(EAR$pfd,w = EAR$area,na.rm=T) #Area-weighted Mean Patch Fractal Dimension
  return(c(MSI, AWMSI, MPFD, AWMPFD))
}
