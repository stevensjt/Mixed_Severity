#Source for filling holes: https://stat.ethz.ch/pipermail/r-sig-geo/2014-January/020139.html
library(maptools)
NZy <- c(-53,-29)
NZx <- c(160,180)
gshhs.c.b <- system.file("share/gshhs_c.b", package="maptools")
NZ <- Rgshhs(gshhs.c.b, xlim=NZx, ylim=NZy, level=2)$SP
NZp <- slot(NZ, "polygons")
holes <- lapply(NZp, function(x) sapply(slot(x, "Polygons"), slot, "hole"))
res <- lapply(1:length(NZp), function(i) slot(NZp[[i]], "Polygons")[!holes[[i]]])
IDs <- row.names(NZ)
NZfill <- SpatialPolygons(lapply(1:length(res), function(i)
  Polygons(res[[i]], ID=IDs[i])), proj4string=CRS(proj4string(NZ)))

#Using the creation functions should ensure reasonable reconstruction, with 
#two reservations:
  
#  1) this provides no protection against an exterior ring nesting inside 
#another exterior ring, so needing a hole to avoid overlapping, and

#2) this doesnt secure the correct construction of comments, needed by 
#GEOS/rgeos to place interior rings in the correct exterior rings.

#You can handle the latter (and possibly your error with your code), by 
#assigning comments or by removing them:

slot(NZfill, "polygons") <- lapply(slot(NZfill, "polygons"),
   checkPolygonsHoles)

#will assign comments, seen by:

lapply(slot(NZfill, "polygons"), comment)

#or remove them by:

slot(NZfill, "polygons") <- lapply(slot(NZfill, "polygons"),
   "comment<-", NULL)

#In rgeos, a NULL comment forces the internal generation of comments, the 
#correct comments in the Polygons objects avoid using time to do this, but 
#stale comments (like yours) are not corrected, and generate an error.

#It is hard to keep track if you dont learn to like and use *apply.

#Hope this helps,

#Roger

#For our purposes:
#Load an hs_fire
plot(hs_fire,col="darkred",border="transparent")
fill_holes=function(hs_fire){
hs_fire_p=slot(hs_fire, "polygons")
holes <- lapply(hs_fire_p, function(x) sapply(slot(x, "Polygons"), slot, "hole"))
areas <- lapply(hs_fire_p, function(x) sapply(slot(x, "Polygons"), slot, "area"))
res <- lapply(1:length(hs_fire_p), function(i) slot(hs_fire_p[[i]], "Polygons")[!(holes[[1]]&areas[[1]]<8100)])#Select the polygons that are not holes. The "i" here is an artifact of the example code; the fires here only have one polygon ID so i=1 (it's a multipart polygon).
IDs <- row.names(hs_fire)
hs_fire_fill <- SpatialPolygons(lapply(1:length(res), function(i)
  Polygons(res[[i]], ID=IDs[i])), proj4string=CRS(proj4string(hs_fire)))
return(hs_fire_fill)
}
plot(hs_fire_fill,col="darkred",border="transparent")
