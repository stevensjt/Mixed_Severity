#ESA SDC analysis

####EDA PLOTS####
#plot a hypothetical temperature heat map with fire perimeter:
r=raster(v_array[,,8],ymn=min(lat),ymx=max(lat),xmn=min(lon),xmx=max(lon))
plot(r)
plot(hs_fire_ll,axes=T,add=T)
points(gCentroid(hs_fire_ll))

fire.list[!is.na(fire.list$q),"sdc_BA90_resid"]=resid(lm(log(q)~BA90_PCT,data=fire.list))
fire.list[which(!is.na(fire.list$q)&!is.na(fire.list$max_bi)),"sdc_BA90_maxbi_resid"]=resid(lm(log(q)~BA90_PCT+max_bi,data=fire.list))

ggplot(fire.list)+
  geom_point(aes(x=BA90_PCT,y=log(q),col=max_bi))

ggplot(fire.list)+
  geom_point(aes(x=max_bi,y=log(q)))+
  geom_smooth(aes(x=max_bi,y=log(q)),method="lm",col="black")

ggplot(fire.list)+
  geom_point(aes(x=BA90_PCT,y=log(q),col=FIRE_YEAR))+
  geom_smooth(aes(x=BA90_PCT,y=log(q)),method="lm",col="black")

ggplot(fire.list)+
  geom_point(aes(x=FIRE_YEAR,y=log(q),col=BA90_PCT))+
  geom_smooth(aes(x=FIRE_YEAR,y=log(q)),method="lm",col="black")

ggplot(fire.list[!is.na(fire.list$q),])+
  geom_point(aes(x=max_bi,y=BA90_PCT,col=log(q)))+
  geom_smooth(aes(x=max_bi,y=BA90_PCT),method="lm",col="black")

ggplot(fire.list)+
  geom_point(aes(x=FIRE_YEAR,y=sdc_BA90_resid,col=BA90_PCT))+
  geom_smooth(aes(x=FIRE_YEAR,y=sdc_BA90_resid),method="lm",col="black")

ggplot(fire.list[-which(fire.list$FIRE_YEAR%in%c(1987, 1999, 2008)),])+
  geom_point(aes(x=FIRE_YEAR,y=sdc_BA90_maxtmmx_resid,col=min_rmax))+
  geom_smooth(aes(x=FIRE_YEAR,y=sdc_BA90_maxtmmx_resid),method="lm",col="black")

lm(sdc_BA90_maxbi_resid~FIRE_YEAR,
   data=fire.list[-which(fire.list$FIRE_YEAR%in%c(1987, 1999, 2008)),]) %>%
  summary()