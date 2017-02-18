##Functions for recharacterizing fire regimes

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
  ##One consequence of this is that it's no longer a sp data frame, and there are some warnings, and possibly buf is large too. 
}

##Decay profile
decay=function(hs_fire,buf_max,buf_inc,name,cancel=F){
  dist.table.sub=
    data.frame(name=rep(name,(buf_max/buf_inc)+1),width=seq(0,buf_max,by=buf_inc),area_ha=NA)
  #Core operation:
  buf.list=lapply(X=-dist.table.sub$width,FUN=gBuffer,spgeom=hs_fire, byid = FALSE, id = NULL, quadsegs = 5, capStyle = "ROUND", joinStyle = "ROUND", mitreLimit = 1) #Apply the buffer at every width in X; returns list of length X.
  buf.list=Filter(length,buf.list) #Remove any NULL values (when buffer eliminated all HS areas)
  dist.table.sub$area_ha=(as.vector(sapply(buf.list,area))*0.0001)[1:nrow(dist.table.sub)] #0.0001 converts m2 to ha
  if(any(is.na(dist.table.sub$area_ha))){#If there are NULL values for the area because the buffer got too wide, remove them.
    dist.table.sub$area_ha[which(is.na(dist.table.sub$area_ha))[1]]=0 #Set the first null value to 0
    dist.table.sub=dist.table.sub[-which(is.na(dist.table.sub$area_ha)),] #Delete the rest of the table rows with NA in area.
  }
  dist.table.sub$prop.hs=dist.table.sub$area_ha/dist.table.sub$area_ha[1]
  return(dist.table.sub)
}

##Calculate decay coefficient sdc

#Version 1
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
