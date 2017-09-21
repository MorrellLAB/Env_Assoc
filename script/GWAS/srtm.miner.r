library(raster)
library(sp)
library(rgdal)

## load location data
loc.d = read.delim('s13059-015-0712-3-s1.txt',header=T)
altitudes = c()

for (i in 1:nrow(loc.d)){
  longi = loc.d[i,'Longitude']
  lati = loc.d[i,'Latitude']
  if (abs(lati) > 60){
    altitudes = c(altitudes,NA)
    next
  }
  dat = 0
  while (length(dat)==1) {
    dat = getData('SRTM',lon=longi,lat=lati)
  }
  dat = setMinMax(dat)
  tile.loc = dat@extent
  lonmin = tile.loc@xmin
  lonmax = tile.loc@xmax
  latmin = tile.loc@ymin
  latmax = tile.loc@ymax
  
  ## col is longitude (min to max), row is latitude (max to min)
  
  x.coor = round( (longi - lonmin) / (lonmax - lonmin) * dat@ncols )
  y.coor = round( (latmax - lati) / (latmax - latmin) * dat@nrows )
  
  altitudes = c(altitudes, dat[y.coor,x.coor,1])  
}
loc.d$altitude = altitudes

write.table(loc.d,file='location.altitude.txt',sep='\t',
            quote=F,row.names=F)
#plot(dat)
#points(longi,lati)
