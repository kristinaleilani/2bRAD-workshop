# Install TESS3R
install.packages("devtools",repos="https://cloud.r-project.org")
devtools::install_github("bcm-uga/TESS3_encho_sen")
# You may need to install Rtools if you have not already. If so, follow the instructions in your console.
library(tess3r)
library(maptools)
library(maps)
library(data.table)
library(sp)
library(scatterpie)
library(dplyr)
library(ggplot2)
library(rworldmap)
library(PBSmapping)
library(tidyr)
setwd("~/Documents/STX") 
# Remember to change "STX" to your species abbreviation (AAGA, PAST, PSTR, OFAV, SSID, or MCAV) throughout 


# Plot sites where samples were collected around St. Croix, VI
sites <- read.csv("sitesSamples.csv") # Import site coordinates
plot(sites[,c("Longitude","Latitude")],pch=16,cex=1,col="red")
text(sites[,c("Longitude","Latitude")],labels=sites$Site,cex=0.7,col="red",pos=4)

# Import shoreline map, taken from https://www.ngdc.noaa.gov/mgg/shorelines/
# Put the downloaded file into your STX directory and unzip it
gshhs.f.b <- "gshhg-bin-2.3.7/gshhs_f.b"
sf1 <- importGSHHS(gshhs.f.b, xlim = c(295,295.5), ylim = c(17.6,17.9)) %>%
  fortify()
sf1$X <- (360-sf1$X)*(-1) # Change longitude coordinates form 0-360 to -180-180
names(sf1)[1] <- "group"
names(sf1)[4] <- "long"
names(sf1)[5] <- "lat"

# Map of collection sites around St. Croix (no genetic info included yet)
ggplot()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_blank(), plot.title = element_text(hjust = 0.5))+
  geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = "grey70", color='black', lwd = 0.1)+
  geom_point(data=sites, aes(x=Longitude, y=Latitude))+
  scale_color_viridis_c(na.value = NA)+
  coord_fixed(xlim = c(-65,-64.5), ylim = c(17.6,17.9), expand = 0)




##### PART ONE: Scatterpie ----------------------------------------##########################
# Let's first look at our admixture groups by plotting pie charts on the map

# Import site info for each sample
samples=read.table("bams.nr")
names(samples)[1]<-"Sample"
samples$Sample=paste(sub(".bam","",samples$Sample),sep="")
rownames(samples)<-samples$Sample
# Merge samples with lat/lon coordinates
i2p=read.csv("STX_sitesample.csv")
i2p=i2p[i2p$Sample %like% "STX", ] # Subset your species
mord=left_join(samples, i2p, by="Sample")
mordi=left_join(mord, sites, by="Site")

# Import admixture groups
pies <- read.table('STX2.qopt')
pies1=cbind(pies, mordi)

pies2 <- pies1 %>%
  mutate(Value = 1) %>%
  group_by(Site) %>%
  summarise(
    radius = sum(Value),
  )
pies3=merge(pies2, pies1, by="Site")
pies3$Sample=NULL
pies4 = pies3 %>% 
  group_by(Site) %>% 
  summarize(V1=sum(V1), V2=sum(V2))
pies3$V1=NULL
pies3$V2=NULL
pies5=merge(pies4, pies3, by="Site")
pies6=aggregate(list(numdup=rep(1,nrow(pies5))),pies5, length)
pies6$Longitude=as.numeric(pies6$Longitude)
pies6$Latitude=as.numeric(pies6$Latitude)
pies6$rad=pies6$numdup/300
ggplot() + 
  geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = 'grey70', color='black', lwd = 0.1) +
  geom_scatterpie(data = pies6, aes(x = Longitude, y = Latitude, r=rad), alpha=0.7, cols = c('V1','V2'), sorted_by_radius=T)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c('tomato', 'lightblue', 'wheat'))+
  coord_equal(xlim = c(-65,-64.5), ylim = c(17.6,17.9), expand = 0)
# Save your scatterpie map
ggsave("STX_scatterpie.tiff", units="in", width=9, height=6, dpi=300, compression = 'lzw')



###---- Plotting the same scatterpie map by # individuals instead of % admixture
# Import cluster.admix (from STX_popstructure.R)
load("STX_clusters.RData")
pies12=cbind(pies1, cluster.admix)
pies2 <- pies12 %>%
  mutate(Value = 1) %>%
  group_by(Site) %>%
  dplyr::summarise(
    radius = sum(Value),
  )
pies3=merge(pies2, pies12, by="Site")
pies31 <- pies3 %>%
  pivot_wider(names_from = "cluster.admix",
              names_prefix = "clust",
              values_from = "cluster.admix",
              values_fn = list(cluster.admix = length),
              values_fill=0)
pies4 = pies31 %>% 
  group_by(Site) %>% 
  dplyr::summarize(clust1=sum(clust1), clust2=sum(clust2))
pies31$Sample=NULL
pies31$clust1=NULL
pies31$clust2=NULL
pies31$V1=NULL
pies31$V2=NULL
pies5=merge(pies4, pies31, by="Site")
pies6=aggregate(list(numdup=rep(1,nrow(pies5))),pies5, length)
ggplot() + 
  geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = 'grey70', color='black', lwd = 0.1) +
  geom_scatterpie(data = pies6, aes(x = Longitude, y = Latitude, r=radius/250), alpha=0.7, cols = c('clust1','clust2'), sorted_by_radius=T)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c('tomato', 'lightblue', 'wheat'))+
  coord_equal(xlim = c(-65,-64.5), ylim = c(17.6,17.9), expand = 0)






##### PART TWO: TESS3R ----------------------------------------##########################
# Tess3r can infer spatial boundaries for admixture groups on the map. 
# To do this, we'll first make an outline (or polygon) around our sites to confine our spatial analysis

# Create convex hull around sampling sites
testpts<-sites[,2:3]
maphull.index<-chull(testpts)
maphull.index <- as.integer(maphull.index)
map.hull <- sites[maphull.index, 3:2]
ggplot()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_blank(), plot.title = element_text(hjust = 0.5))+
  geom_polygon(data = map.hull, aes(x=Longitude, y = Latitude), fill = NA, color='black', lwd = 0.5)+
  geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = "grey70", color='black', lwd = 0.1)+
  geom_point(data=sites, aes(x=Longitude, y=Latitude))+
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
map.hull.stretch <- as.data.frame(map.hull.stretch)
ggplot()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_blank(), plot.title = element_text(hjust = 0.5))+
  geom_polygon(data = map.hull.stretch, aes(x=V1, y = V2), fill = NA, color='black', lwd = 0.5)+
  geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = "grey70", color='black', lwd = 0.1)+
  geom_point(data=sites, aes(x=Longitude, y=Latitude))+
  scale_color_viridis_c(na.value = NA)+
  coord_fixed(xlim = c(-65,-64.5), ylim = c(17.6,17.9), expand = 0)

# Save your map outline so you can import it next time and won't have to recreate it
#map.hull.stretch <- write.table("map.hull.stretch")
#map.hull.stretch <- read.table('map_hull_stretch', sep = " ", header = FALSE)



# Plot sites inside the map outline
coordinates = mordi[,c("Longitude","Latitude")]
outline=map.hull.stretch
names(outline)=c("Longitude","Latitude")
margin=0.2
plot(coordinates, 
     pch = 19, 
     cex = .5, 
     xlab = "Longitude (°E)", 
     ylab = "Latitude (°N)",
     asp=1,
     xlim=c(min(mordi$Longitude)-margin,
            max(mordi$Longitude)+margin))
polygon(outline,col="grey90")
points(coordinates, pch = 19, cex = 0.6,col="red")
text(mordi[,c("Longitude","Latitude")],labels=mordi$Site,cex=0.7,col="red",pos=4)

# Set a color palette
my.colors <- c("tomato", "lightblue", "wheat","olivedrab", "cyan3","hotpink","gold","orange")
my.palette <- CreatePalette(my.colors, 8)

# Import admixture groups, and plot as barplot
qm=as.qmatrix(as.matrix(read.table("STX2.qopt")))
bp=barplot(qm,border=NA,space=0,col.palette = my.palette)
axis(1, at = 1:nrow(qm), labels = i2p$Site[bp$order], las = 3, cex.axis = .4) 

Outline=Polygon(outline)
BL=Polygons(list(Outline),"outline")
BLL=SpatialPolygons(list(BL))
out=SpatialPolygonsDataFrame(BLL,data.frame(N=c("one"),row.names=c("outline")))

# Plot admixture boundaries on the map
plot(qm,coordinates,method="map.max",interpol=FieldsKrigModel(10),resolution=c(300,300),cex=1,map.polygon=out,asp=1,col.palette = CreatePalette(my.colors, 8),xlab = "Longitude (°E)", ylab = "Latitude (°N)")
text(mordi[,c("Longitude","Latitude")],labels=mordi$Site,cex=0.7,col="red",pos=4)

# Plot map using ggplot (prettier option)
pl <- ggtess3Q(qm, coordinates, map.polygon = out, interpol=FieldsKrigModel(5),resolution=c(300,300), col.palette = CreatePalette(my.colors, 8))
pl +
  geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = "grey70", color='black', lwd = 0.1)+
  geom_point(data = as.data.frame(coordinates), aes(x = Longitude, y = Latitude), size = 2) + 
  xlab("Longitute") +
  ylab("Latitude") + 
  theme_bw()+
  coord_fixed(xlim = c(-64.95,-64.5), ylim = c(17.65,17.8), expand = 0)
# Save your tess3r plot:
ggsave("STX_tess3r.tiff", units="in", width=9, height=6, dpi=300, compression = 'lzw')

