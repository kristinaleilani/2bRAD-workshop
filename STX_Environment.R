# To install gradient forest, see https://github.com/z0on/RDA-forest
library(ggplot2)
library(vegan)
library(gradientForest)
library(dplyr)
library(maptools)
library(PBSmapping)
setwd("~/Documents/STX") 
source("~/Documents/STX/RDA-forest_functions.R")
# Remember to change "STX" to your species abbreviation (ex. PSTR, AAGA, PAST, MCAV, OFAV, SSID)


###---- Load the shoreline (download from www.ngdc.noaa.gov/mgg/shorelines)
gshhs.f.b <- "gshhg-bin-2.3.7/gshhs_f.b"
sf1 <- importGSHHS(gshhs.f.b, xlim = c(295,295.5), ylim = c(17.6,17.9)) %>%
  fortify()
sf1$X <- (360-sf1$X)*(-1) # Change longitude from 0-360 to -180-180
names(sf1)[1] <- "group"
names(sf1)[4] <- "long"
names(sf1)[5] <- "lat"

###---- Import species files 
IBS=as.matrix(read.table("STX2.ibsMat"))
samples=read.table("bams.nr")
samples$V1=paste(sub(".bam","",samples$V1),sep="")
samples=samples$V1
dimnames(IBS)=list(samples,samples)

inds=as.data.frame(samples)
colnames(inds) = "Sample"
row.names(inds)<-inds$Sample

# If necessary, remove technical replicates. If not, skip.
goods=samples[! samples %in% c('PSTR63a')]
IBS=IBS[goods,goods] 
inds=as.data.frame(inds[goods,] )
colnames(inds) = "Sample"
row.names(inds)<-inds$Sample
samples=goods
names(samples)[1]<-"Sample"



###---- Import site info 
sites=read.csv("STX_sitesample.csv")
coords <- read.csv("sitesSamples.csv") 
sites0 <- merge(inds, sites, by="Sample", sort=F)
sites1 <- merge(sites0, coords, by="Site", sort=F)
rownames(sites1)<-sites1$Sample
latlon=sites1[,c("Latitude","Longitude")]
names(latlon)=c("lat","lon")

###----  Import environmental data
depth=read.csv("STX_Sample_depth.csv")
depth$Site=NULL
env0=read.csv("STX_env.csv")
env1=merge(env0, sites0, by="Site")
env=merge(env1, depth, by="Sample")
rownames(env)=env$Sample
env$Sample=NULL
env$Site=NULL
env$Lat=NULL
env$Lon=NULL






#####--------- Check pop structure of lineages and remove clones if necessary ---------#######
hc=hclust(as.dist(IBS),"ave")
plot(hc,cex=0.5) # clustering of samples by IBS (great to detect clones or closely related individuals)
abline(h=0.15,col="red") # this seems like a reasonable  "low-hang" threshold for calling related groups
# If you don't have clones to remove, continue to spatial predictors

# If you do have clones, remove them here:
cuts=cutree(hc,h=0.15)
goods=c();i=1
for (i in unique(cuts)) {
  goods[i]=names(cuts[cuts==i])[1]
}
length(goods)  # how many samples are left?
# Subsetting all data for only the retained samples
IBS=IBS[goods,goods] 
latlon=latlon[goods,]
env=env[goods,]
inds=inds[goods,]
sites1=sites1[goods,]

# Evaluate PCoA for outliers before continuing
ord.all=capscale(as.dist(IBS)~1)
plot(ord.all,scaling=1)
points(ord.all,scaling=1, pch=16)
summary(ord.all)




#####---------   Explore isolation by distance ---------#######
geo_dist<- dist(latlon)
ibs_dist<- as.dist(IBS)
#Plot IBS vs geo distance
ggplot(NULL, aes(x=geo_dist, y=ibs_dist) ) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis", trans="log") +
  geom_smooth(method="gam", formula = y ~s(x))+
  labs(title = "Isolation by distance") +
  xlab('Geographic distance') +
  ylab('Genetic distance') +
  theme_bw()
ggsave("STX_IBD.tiff", units="in", width=4, height=4, dpi=300, compression = 'lzw')



#####---------  Plot PCoA in ggplot, colored by environmental variables. Look for gradients across PCA space. ---------#######
ord=capscale(IBS~1)
summary(ord) 
ords=scores(ord,display="sites")
axes2plot=c(1:4) # which PCAs to plot
scores=data.frame(ord$CA$u[,axes2plot])
scores=cbind(scores, env)
scores$cluster.admix<-as.factor(scores$cluster.admix)

# Plot a PCoA colored by depth:
ggplot(scores,aes(scores[,1],scores[,2], asp=1, fill=Depth)) + 
  geom_point(aes(size=1)) +
  scale_fill_continuous(trans = 'reverse') +
  theme_bw()+
  coord_fixed()+
  xlab(names(scores)[1])+
  ylab(names(scores)[2])+
  #geom_label(label=rownames(scores))+
  guides(size = "none")
ggsave("STX_PCoA_depth.tiff", units="in", width=9, height=5, dpi=300, compression = 'lzw')

# Plot a PCoA colored by any other environmental variable:
ggplot(scores,aes(scores[,1],scores[,2], asp=1, fill=Temp_mean)) + 
  geom_point(aes(size=1)) +
  theme_bw()+
  coord_fixed()+
  xlab(names(scores)[1])+
  ylab(names(scores)[2])+
  #geom_label(label=rownames(scores))+
  guides(size = "none")
ggsave("STX_PCoA_env.tiff", units="in", width=9, height=5, dpi=300, compression = 'lzw')

# Plot a PCoA colored by any other longitude or latitude:
scores=cbind(scores, latlon)
ggplot(scores,aes(scores[,1],scores[,2], asp=1, fill=lat)) + 
  geom_point(aes(size=1)) +
  theme_bw()+
  coord_fixed()+
  xlab(names(scores)[1])+
  ylab(names(scores)[2])+
  #geom_label(label=rownames(scores))+
  guides(size = "none")
ggsave("STX_PCoA_lat.tiff", units="in", width=9, height=5, dpi=300, compression = 'lzw')

ggplot(scores,aes(scores[,1],scores[,2], asp=1, fill=lon)) + 
  geom_point(aes(size=1)) +
  theme_bw()+
  coord_fixed()+
  xlab(names(scores)[1])+
  ylab(names(scores)[2])+
  #geom_label(label=rownames(scores))+
  guides(size = "none")
ggsave("STX_PCoA_lon.tiff", units="in", width=9, height=5, dpi=300, compression = 'lzw')

       





#############------------------ GRADIENT FOREST: Finding genotype-environment associations  -----------------###############

#####--------- Create spatial predictors  ---------#######
# principal components of space - lat and lon rotated to be uncorrelated, and centered 
xy=scores(capscale(dist(latlon)~1),scaling=1)$sites
colnames(xy)=c("xx","yy")
# Moran eigenvector maps (MEMs) - capturing possible spatial trends
mems=data.frame(pcnm(dist(xy))$vectors)
# adding xy and first 5 MEMs (will be called PCNM1-5) to env
env=cbind(env,xy)
env=cbind(env,mems[,1:5])
colnames(env)
# remembering the names of spatial predictors
space=colnames(env)[grep("PCNM|xx|yy",colnames(env))]



#####--------- Variable selection  ---------#######
# See https://github.com/z0on/RDA-forest for more explanation

mm=mtrySelection(Y=IBS,X=env,nreps=11, prop.positive.cutoff=0.5,top.pcs=25)
# If one of the predictors is constant across sites, remove and rerun mtrySelection
env$Entero_min=NULL
env$Entero_mean_yearly_range=NULL

# boxplot of importance differences at different mtry 
ggplot(mm$delta,aes(var,values))+
  geom_boxplot()+
  coord_flip()+
  geom_hline(yintercept=0,col="red")
ggsave("STX_gf_boxplot.tiff", units="in", width=4, height=4, dpi=300, compression = 'lzw')

# bar chart of proportion of positive change in response to higher mtry, good variables would be the ones above the red line
ggplot(mm$prop.positive,aes(var,prop.positive))+
  geom_bar(stat="identity")+
  coord_flip()+
  geom_hline(yintercept=0.5,col="red")
ggsave("STX_gf_barchart.tiff", units="in", width=4, height=4, dpi=300, compression = 'lzw')


#--------------- Spatial bootstrap of mtry-passing variables
ll=load("rasters_XY.RData") # load table of environmental values for a grid of spatial points where we want to predict how our creatures would adapt.
ll
# "rasters" "XY"
plot(XY,pch=".",asp=1) # view the grid 

# if there are no new points to predict, just skip the newX option in the call to spatialBootstrap; 
# the predictions will be made for original data points then.
sb=spatialBootstrap(Y=IBS,X=env,newX=rasters,nreps=25,top.pcs=25)

# plot importance boxplot including space variables
ggplot(sb$all.importances,aes(variable,importance))+geom_boxplot()+coord_flip()+theme_bw()
sum(sb$median.importance) #proportion of variation explained
ggsave("STX_gf_bootstrap_space.tiff", units="in", width=4, height=4, dpi=300, compression = 'lzw')

# plot importance boxplot without space variables
ggplot(sb$all.importances[!(sb$all.importances$variable %in% space),],aes(variable,importance))+geom_boxplot()+coord_flip()+theme_bw()
sum(sb$median.importance[!(names(sb$median.importance) %in% space)]) #proportion of variation explained without spatial variables
ggsave("STX_gf_bootstrap_nospace.tiff", units="in", width=4, height=4, dpi=300, compression = 'lzw')



#####------------------- Plotting genetic turnover (adaptation) map ---------#######
# see https://gradientforest.r-forge.r-project.org/biodiversity-survey.pdf for more explanation

turnovers=sb$turnovers
write.csv(turnovers, file=paste("STX_turnovers.csv", sep='_'))

# if we want to keep only top 5 predictors (recommended for cases when there are clear "winner" predictors):
raster.vars=names(sb$median.importance[!names(sb$median.importance) %in% space])[1:5]

# principal component analysis of the predicted genetic turnover patterns
pc <- prcomp(turnovers[, raster.vars])
plot(pc$sdev)
# we will try to visualize the first 3 PCs
pcs2show=c(1,2,3)

# color flippage flags - change between -1 and 1 to possibly improve color representation in the final map
flip.pc1=(1)
flip.pc2=(1)
flip.pc3=(1)

flip=""
if(flip.pc1==(-1)) { flip=paste(flip,"1",sep="")}
if(flip.pc2==(-1)) { flip=paste(flip,"2",sep="")}
if(flip.pc3==(-1)) { flip=paste(flip,"3",sep="")}

# magic to color three PCA dimensions on a map
pc1 <- flip.pc1*pc$x[, pcs2show[1]]
pc2 <- flip.pc2*pc$x[, pcs2show[2]]
pc3 <- flip.pc3*pc$x[, pcs2show[3]]
b <- pc1 - pc2
g <- -pc1
r <- pc3 + pc2 - pc1
r <- (r - min(r))/(max(r) - min(r)) * 255
g <- (g - min(g))/(max(g) - min(g)) * 255
b <- (b - min(b))/(max(b) - min(b)) * 255
nvs <- dim(pc$rotation)[pcs2show[1]]
lv <- length(raster.vars)
vind <- rownames(pc$rotation) %in% raster.vars
scal <- 30
xrng <- range(pc$x[, 1], pc$rotation[, pcs2show[1]]/scal) * 1.9
yrng <- range(pc$x[, 2], pc$rotation[, pcs2show[2]]/scal) * 1.9
man.colors=rgb(r, g, b, max = 255)


# -------- plotting the map of predicted adaptive communities

coltxt="coral" 
important=raster.vars[order(sqrt(pc$rotation[raster.vars, 1]^2+pc$rotation[raster.vars, 2]^2),decreasing=T)][1:min(3,length(raster.vars))]
lv=length(important)
par(mfrow=c(1,2))
plot((pc$x[, pcs2show[1:2]]), xlim = xrng, ylim = yrng, pch = ".", cex = 4, col = rgb(r, g, b, max = 255), asp = 1,xaxt="none",yaxt="none",bty="none",xlab="",ylab="")
jit <- rnorm(lv,0.005,0.001)
arrows(rep(0, lv), rep(0, lv), pc$rotation[important, pcs2show[1]]/scal, pc$rotation[important, pcs2show[2]]/scal, length = 0.0625,col=coltxt)
text(pc$rotation[important, 1]/scal + jit * sign(pc$rotation[important, pcs2show[1]]), pc$rotation[important, pcs2show[2]]/scal + jit * sign(pc$rotation[important, pcs2show[2]]), labels = important,cex=0.7,col=coltxt)
plot(XY, pch=15,cex = 0.5, asp = 1, col = man.colors)
map(coasts,add=T,col="grey80",fill=T,border="grey80",lwd=1)

# contrasting colors = habitats requiring differential adaptation, likely driven by factors in the legend.
# you can try plotting any of the environmental variables you think are important
library(viridis)
ggplot(XY,aes(x,y,color=rasters$TEMP_yearly_range))+geom_point()+scale_color_viridis()+coord_equal()
ggplot(XY,aes(x,y,color=rasters$TEMP_mean))+geom_point()+scale_color_viridis()+coord_equal()





