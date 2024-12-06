library(tidyr)
library(dplyr)
library(purrr)
library(raster)
library(automap)
setwd("~/Spatial_interpolation")

###-------- import Enterococcus -----#########
Entero.site<-read.csv("Entero_STX.csv")[,-1]

Entero.data.filt.year=Entero.data.filt %>% 
  select(Site, Year) %>% 
  unique() %>% 
  group_by(Year) %>% 
  summarize(Enteron = n())
Entero.data.filt.month=Entero.data.filt %>% 
  select(Site, Month) %>% 
  unique() %>% 
  group_by(Month) %>% 
  summarize(Enteron = n())


###-------- import Dissolved Oxygen -----#########
DO.site<-read.csv("DO_STX.csv")[,-1]

DO.data.filt.year=DO.data.filt %>% 
  select(Site, Year) %>% 
  unique() %>% 
  group_by(Year) %>% 
  summarize(DOn = n())
DO.data.filt.month=DO.data.filt %>% 
  select(Site, Month) %>% 
  unique() %>% 
  group_by(Month) %>% 
  summarize(DOn = n())


###-------- import E. coli -----#########
Ecoli.site<-read.csv("Ecoli_STX.csv")[,-1]

Ecoli.data.filt.year=Ecoli.data.filt %>% 
  select(Site, Year) %>% 
  unique() %>% 
  group_by(Year) %>% 
  summarize(Ecolin = n())
Ecoli.data.filt.month=Ecoli.data.filt %>% 
  select(Site, Month) %>% 
  unique() %>% 
  group_by(Month) %>% 
  summarize(Ecolin = n())


###-------- import Nitrogen -----#########
Nitrogen.site<-read.csv("Nitrogen_STX.csv")[,-1]

Nitrogen.data.filt.year=Nitrogen.data.filt %>% 
  select(Site, Year) %>% 
  unique() %>% 
  group_by(Year) %>% 
  summarize(Nitrogenn = n())
Nitrogen.data.filt.month=Nitrogen.data.filt %>% 
  select(Site, Month) %>% 
  unique() %>% 
  group_by(Month) %>% 
  summarize(Nitrogenn = n())


###-------- import pH -----#########
pH.site<-read.csv("pH_STX.csv")[,-1]

pH.data.filt.year=pH.data.filt %>% 
  select(Site, Year) %>% 
  unique() %>% 
  group_by(Year) %>% 
  summarize(pHn = n())
pH.data.filt.month=pH.data.filt %>% 
  select(Site, Month) %>% 
  unique() %>% 
  group_by(Month) %>% 
  summarize(pHn = n())


###-------- import Phosphorus -----#########
Phosphorus.site<-read.csv("Phosphorus_STX.csv")[,-1]

Phosphorus.data.filt.year=Phosphorus.data.filt %>% 
  select(Site, Year) %>% 
  unique() %>% 
  group_by(Year) %>% 
  summarize(Phosphorusn = n())
Phosphorus.data.filt.month=Phosphorus.data.filt %>% 
  select(Site, Month) %>% 
  unique() %>% 
  group_by(Month) %>% 
  summarize(Phosphorusn = n())



###-------- import Secchi -----#########
Secchi.site<-read.csv("Secchi_STX.csv")[,-1]

Secchi.data.filt.year=Secchi.data.filt %>% 
  select(Site, Year) %>% 
  unique() %>% 
  group_by(Year) %>% 
  summarize(Secchin = n())
Secchi.data.filt.month=Secchi.data.filt %>% 
  select(Site, Month) %>% 
  unique() %>% 
  group_by(Month) %>% 
  summarize(Secchin = n())


###-------- import Temperature -----#########
Temp.site<-read.csv("Temp_STX.csv")[,-1]

Temp.data.filt.year=Temp.data.filt %>% 
  select(Site, Year) %>% 
  unique() %>% 
  group_by(Year) %>% 
  summarize(Tempn = n())
Temp.data.filt.month=Temp.data.filt %>% 
  select(Site, Month) %>% 
  unique() %>% 
  group_by(Month) %>% 
  summarize(Tempn = n())


# Observations by month
month_n = full_join(Entero.data.filt.month, full_join(DO.data.filt.month, full_join(Ecoli.data.filt.month, full_join(Nitrogen.data.filt.month, full_join(pH.data.filt.month, full_join(Phosphorus.data.filt.month, full_join(Secchi.data.filt.month, Temp.data.filt.month)))))))
write.csv(month_n, "Month_n.csv")
month_n=read.csv("Month_n.csv")

# Observations by year
year_n = full_join(Entero.data.filt.year, full_join(DO.data.filt.year, full_join(Ecoli.data.filt.year, full_join(Nitrogen.data.filt.year, full_join(pH.data.filt.year, full_join(Phosphorus.data.filt.year, full_join(Secchi.data.filt.year, Temp.data.filt.year)))))))
write.csv(year_n, "Year_n.csv")
year_n=read.csv("Year_n.csv")

# Import shoreline map
library(maptools)
setwd("~/Documents/SERC/") 
if (!rgeosStatus()) gpclibPermit()
gshhs.f.b <- "gshhg-bin-2.3.6/gshhs_f.b"
sf1 <- getRgshhsMap(gshhs.f.b, xlim = c(-65,-64.5), ylim = c(17.6,17.9)) %>%
  fortify()

# Example plot
ggplot()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_blank(), plot.title = element_text(hjust = 0.5))+
  geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = "grey70", color='black', lwd = 0.1)+
  geom_point(data=Entero.site, aes(x=Lon, y=Lat, color = Entero_mean))+
  scale_color_viridis_c(na.value = NA)+
  coord_fixed(xlim = c(-65,-64.5), ylim = c(17.6,17.9), expand = 0)



###-------- Spatial interpolation -----#########

# Create convex hull around sampling sites
testpts<-Entero.site[,2:3]
maphull.index<-chull(testpts)
maphull.index <- as.integer(maphull.index)
map.hull <- Entero.site[maphull.index, 3:2]
ggplot()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_blank(), plot.title = element_text(hjust = 0.5))+
  geom_polygon(data = map.hull, aes(x=Lon, y = Lat), fill = NA, color='black', lwd = 0.5)+
  geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = "grey70", color='black', lwd = 0.1)+
  geom_point(data=Entero.site, aes(x=Lon, y=Lat))+
  scale_color_viridis_c(na.value = NA)+
  coord_fixed(xlim = c(-65,-64.5), ylim = c(17.6,17.9), expand = 0)

# Stretch the hull by generating a circle of points that surround the vertices
n <- 6 # number of points you want on the circle
pts.circle <- t(sapply(1:n,function(r)c(cos(2*r*pi/n),sin(2*r*pi/n))))
buffer.data <- list()
for(i in 1:nrow(map.hull)){
  # The higher the denominator, the smaller the circle/buffer
  buffer.data[[i]] <- data.frame(x = map.hull$Lon[i] + pts.circle[,1]/20, 
                                 y = map.hull$Lat[i] + pts.circle[,2]/20)
}
buffer.data.long <- do.call(rbind, buffer.data)
# Generate a new concave hull using the outside of the circles as the reference
map.hull.stretch <- data.frame(concaveman::concaveman(as.matrix(buffer.data.long), concavity = 1, length_threshold = 0.3))
#map.hull.stretch <- read.table('map_hull_stretch.txt', sep = " ", header = FALSE)
map.hull.stretch <- as.data.frame(map.hull.stretch)
ggplot()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_blank(), plot.title = element_text(hjust = 0.5))+
  geom_polygon(data = map.hull.stretch, aes(x=V1, y = V2), fill = NA, color='black', lwd = 0.5)+
  geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = "grey70", color='black', lwd = 0.1)+
  geom_point(data=Entero.site, aes(x=Lon, y=Lat))+
  scale_color_viridis_c(na.value = NA)+
  coord_fixed(xlim = c(-65,-64.5), ylim = c(17.6,17.9), expand = 0)

# Convert hull to Spatial object
library(raster)
hull.poly <- SpatialPolygons(list(Polygons(list(Polygon(cbind(map.hull.stretch$V1, map.hull.stretch$V2))), ID=1)), proj4string = CRS("+proj=longlat +datum=NAD83"))
hull.poly.df <- SpatialPolygonsDataFrame(hull.poly, data=data.frame(ID=1))
hull.grid <- raster(hull.poly, res = 1/250)



###-------- Select coral sampling sites for extracting data -----#########
setwd("~/Documents/STX/SppFiles")
site.locations <- read.csv('sitesSamples.csv', header = T)

# View collection sites on the map
ggplot()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_blank(), plot.title = element_text(hjust = 0.5))+
  geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = "grey70", color='black', lwd = 0.1)+
  geom_point(data = site.locations, aes(x = Longitude, y = Latitude), color = 'red', size = 0.5)+
  scale_fill_viridis_c(na.value = NA)+
  coord_fixed(xlim = c(-65,-64.5), ylim = c(17.6,17.9), expand = 0)

site.coords <- SpatialPoints(site.locations[,c('Longitude', 'Latitude')], proj4string = CRS("+proj=longlat +datum=NAD83"))
site.locations.shifted <- vector()
for (each.site in 1:nrow(site.locations)){
  # geodesic distances
  distances <- fields::rdist.earth(coordinates(site.coords[each.site]), coordinates(krige.mask), miles=FALSE)
  # coordinates of the pixel that is at shorter distance 
  shifted.coordinates<-coordinates(krige.mask)[which.min(distances),]
  site.locations.shifted <- rbind(site.locations.shifted, shifted.coordinates)
}
new.locations <- data.frame(site.locations.shifted) %>%
  cbind(site.locations, .)



###-------- Convert sampling data to Spatial object -----#########

# Choose one variable at a time to send through kriging loop
Temp.only<-Temp.site[,2:8]
wq.sp <- SpatialPoints(Temp.only[,2:1], proj4string = CRS("+proj=longlat +datum=NAD83"))
wq.spdf <- SpatialPointsDataFrame(wq.sp, as.data.frame(Temp.only))

DO.only<-DO.site[,2:8]
wq.sp <- SpatialPoints(DO.only[,2:1], proj4string = CRS("+proj=longlat +datum=NAD83"))
wq.spdf <- SpatialPointsDataFrame(wq.sp, as.data.frame(DO.only))

Ecoli.only<-Ecoli.site[,2:8]
Ecoli.only$Ecoli_min=NULL
wq.sp <- SpatialPoints(Ecoli.only[,2:1], proj4string = CRS("+proj=longlat +datum=NAD83"))
wq.spdf <- SpatialPointsDataFrame(wq.sp, as.data.frame(Ecoli.only))

Entero.only<-Entero.site[,2:9]
wq.sp <- SpatialPoints(Entero.only[,2:1], proj4string = CRS("+proj=longlat +datum=NAD83"))
wq.spdf <- SpatialPointsDataFrame(wq.sp, as.data.frame(Entero.only))

Nitrogen.only<-Nitrogen.site[,2:8]
wq.sp <- SpatialPoints(Nitrogen.only[,2:1], proj4string = CRS("+proj=longlat +datum=NAD83"))
wq.spdf <- SpatialPointsDataFrame(wq.sp, as.data.frame(Nitrogen.only))

pH.only<-pH.site[,2:8]
wq.sp <- SpatialPoints(pH.only[,2:1], proj4string = CRS("+proj=longlat +datum=NAD83"))
wq.spdf <- SpatialPointsDataFrame(wq.sp, as.data.frame(pH.only))

Phosphorus.only<-Phosphorus.site[,2:8]
wq.sp <- SpatialPoints(Phosphorus.only[,2:1], proj4string = CRS("+proj=longlat +datum=NAD83"))
wq.spdf <- SpatialPointsDataFrame(wq.sp, as.data.frame(Phosphorus.only))

Secchi.only<-Secchi.site[,2:8]
wq.sp <- SpatialPoints(Secchi.only[,2:1], proj4string = CRS("+proj=longlat +datum=NAD83"))
wq.spdf <- SpatialPointsDataFrame(wq.sp, as.data.frame(Secchi.only))



###-------- Kriging loop -----#########

setwd("~/Documents/TNC/Data/Waterqualitydata/Interp")
library(automap)

data.out = data.frame(new.locations[,1:3])
for(c in names(wq.spdf)) { 
  interp.var <- c
  test <- na.omit(wq.spdf@data[,c("Lat", "Lon", c)])
  #make test into a new spatial dataframe
  wq.sp.new <- SpatialPoints(test[,2:1], proj4string = CRS("+proj=utm +datum=NAD83"))
  wq.spdf.new <- SpatialPointsDataFrame(wq.sp.new, as.data.frame(test))
  
  # Auto-fit a model variogram
  variogram<-automap::autofitVariogram(get(interp.var)~1, wq.spdf.new)
  
  plot(variogram)
  spatial.points <- SpatialPoints(coords = coordinates(hull.grid), 
                                  proj4string = CRS("+proj=utm +datum=NAD83"))
  spatial.grid <- SpatialPixels(points = spatial.points)
  kriging_result <- automap::autoKrige(get(interp.var)~1, wq.spdf.new, spatial.grid)
  # save variogram
  pdf(file = paste0("variogram_", interp.var, ".pdf"))
  plot(variogram)
  dev.off()
  
  # Create Rasterbrick object and mask to St. Croix polygon
  krige.brick <- brick(kriging_result$krige_output)
  krige.mask <- mask(krige.brick, hull.poly)
  names(krige.mask) <- c('prediction', 'variance', 'stdev')
  # save prediction map
  pdf(file = paste0("prediction_", interp.var, ".pdf"))
  plot(krige.mask)
  dev.off()
  
  krige.predict <- cbind(as.data.frame(krige.mask$prediction), as.data.frame(spatial.grid))
  ggpredict <- ggplot()+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title = element_blank(), plot.title = element_text(hjust = 0.5))+
    geom_tile(data = krige.predict, aes(x = x, y = y, fill = prediction, col = prediction))+
    geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = "grey70", color='black', lwd = 0.1)+
    scale_fill_viridis_c(na.value = NA)+
    scale_color_viridis_c(na.value = NA)+
    coord_fixed(xlim = c(-65,-64.5), ylim = c(17.6,17.9), expand = 0)
  # # save ggplot prediction map
  ggsave(file = paste0("ggpredict_", interp.var, ".pdf"))
  plot(ggpredict)
  dev.off()
  
  # save a raster
  krig.fixed<-data.frame(cbind(krige.predict$x, krige.predict$y, krige.predict$prediction))
  r = rasterFromXYZ(krig.fixed)
  plot(r)
  writeRaster(r, file = paste0(interp.var), format="raster", overwrite=TRUE)
  
  #extract data
  extract.data <- data.frame(raster::extract(krige.mask, new.locations[,4:5])[,1])
  colnames(extract.data) <- interp.var
  data.out <- cbind(data.out, extract.data)
  data.out$Lat=NULL
  data.out$Lon=NULL
  
}

write.csv(data.out, "STXenv_allvars.csv")







