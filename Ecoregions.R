library(dplyr)
library(vegan)
library(ggplot2)
library(data.table)
library(tidyr)
setwd("~/Documents/STX/SppFiles") 
source("RDAforest_functions.R")

########-------- Import lineage assignments for all species ------#######
# Siderastrea siderea
SSIDclust=read.csv("SSID_3clust.csv")[,-1]
SSIDclust$cluster = paste0('SSID_', SSIDclust$cluster)
# Orbicella faveolata
OFAVclust=read.csv("OFAV_2clust.csv")[,-1]
OFAVclust$cluster = paste0('OFAV_', OFAVclust$cluster)
# Montastrea cavernosa
MCAVclust=read.csv("MCAV_6clust.csv")[,-1]
MCAVclust$cluster = paste0('MCAV_', MCAVclust$cluster)
# Pseudodiploria strigosa
PSTRclust=read.csv("PSTR_4clust.csv")[,-1]
PSTRclust$Site=NULL
PSTRclust$Latitude=NULL
PSTRclust$Longitude=NULL
PSTRclust$cluster = paste0('PSTR_', PSTRclust$cluster)
# Porites astreoides
PASTclust=read.csv("PAST_2clust.csv")[,-1]
PASTclust$cluster = paste0('PAST_', PASTclust$cluster)
# Agaricia agaricites
AAGAclust=read.csv("AAGA_2clust.csv")[,-1]
AAGAclust$cluster = paste0('AAGA_', AAGAclust$cluster)
########-------- Combine lineage assignments from all species into a single table of community abundances ------#######
clust=rbind(AAGAclust, MCAVclust, OFAVclust, PASTclust, PSTRclust, SSIDclust)
# Merge samples with lat/lon coordinates
i2p=read.csv("STX_sitesample.csv")
mord=left_join(clust, i2p, by="Sample")
clust0=mord
clust0$Sample=NULL
clust.wide= clust0 %>%
  pivot_wider(names_from = "cluster",
              values_from = "cluster",
              values_fn = list(cluster = length),
              values_fill=0)
clust.wide=as.data.frame(clust.wide)

########-------- Get ecoregions ------#######
ecor1 <- data.frame(Site=c('Isaac Bay', 'Columbus landing', 'Deep end', 'Butler Bay',
                          'Carambola shallow', 'The Palms', 'Fredericksted Pier', 'Cane Bay',
                          'Carambola deep', 'Cane Bay deep', 'WAPA', 'Channel Rock', 
                          'Joes Reef', 'North Star'),
                  ecoregion=c('ecoregionA', 'ecoregionC', 'ecoregionA', 'ecoregionC',
                              'ecoregionC', 'ecoregionB', 'ecoregionD', 'ecoregionC',
                              'ecoregionC', 'ecoregionC', 'ecoregionB', 'ecoregionA', 
                              'ecoregionA', 'ecoregionC'))
ecor1[13,"Site"]="Joe's Reef"
rownames(ecor1)=ecor1$Site
ecor2=ecor1
ecor2$Site=NULL
# Dummify categorical ecoregion variables
ecor=dummify(ecor2)
ecor$Site=ecor1$Site

########-------- Import coordinates for sites ------#######
coords <- read.csv("sitesSamples.csv") 
latlon=coords[,c("Latitude","Longitude")]
names(latlon)=c("lat","lon")
xy=as.data.frame(scores(capscale(dist(latlon)~1),scaling=1)$sites)
colnames(xy)=c("xx","yy")
xy$Site=coords$Site

########-------- Import depth of sites ------#######
depth=read.csv("STX_Sample_depth.csv")
depth$Sample=NULL
depth.s=aggregate(list(numdup=rep(1,nrow(depth))),depth, length)
depth.s2=depth.s
depth.s2$Site=paste(depth.s2$Site, rowid(depth.s$Site), sep = "_")
depth.s3=cbind(depth.s$Site, depth.s2)
names(depth.s3)[2]="Sitenumber"
names(depth.s3)[1]="Site"
depth2=left_join(depth, depth.s3, by=c("Site", "Depth"))
depth2=left_join(depth2, xy, by="Site")
depth=read.csv("STX_Sample_depth.csv")
depth3=cbind(depth$Sample, depth2)
names(depth3)[1]="Sample"
ecod=left_join(depth3, ecor, by="Site")
ecod$numdup=NULL
ecod$Sample=NULL
ecod$Site=NULL
ecod=aggregate(list(numdup=rep(1,nrow(ecod))),ecod, length)
ecod$numdup=NULL
names(ecod)[2]="Site"

########-------- Cleanup X and Y tables for gradient forest ------#######
mord0=left_join(mord, depth3, by=c("Sample", "Site"))
mord0$Site=NULL
names(mord0)[4]="Site"
mord0$Depth=NULL
mord0$numdup=NULL
clust1=mord0
clust1$Sample=NULL
clust.wide= clust1 %>%
  pivot_wider(names_from = "cluster",
              values_from = "cluster",
              values_fn = list(cluster = length),
              values_fill=0)
clust.wide=as.data.frame(clust.wide)

ecod=ecod %>% semi_join(clust.wide, by = "Site") 
rownames(ecod)=ecod$Site
rownames(clust.wide)=clust.wide$Site
ecod$Site=NULL
clust.wide$Site=NULL
# Match order between dataframes
clust.wide=clust.wide[ order(match(rownames(clust.wide), rownames(ecod))), ]
# Turn community adundance table into a Bray-Curtis dissimilarity matrix
clust=as.matrix(vegdist(clust.wide, method="bray"))
clust=data.frame(clust)
env=as.data.frame(ecod)
# Remove nonsense variables and variables that are the same at every site
env$numdup=NULL
env$Entero_mean_yearly_range=NULL
env$Entero_min=NULL
# Regress latitude and longitude out of the ordination of the community matrix
clust0=capscale(clust~1+Condition(as.matrix(env[,c("xx","yy")])))
env0=env[,!names(env) %in% c("xx","yy")]
# Look at PCA of sites
plot(clust0, scaling=1, choices=c(1,2))
plot(clust0$CA$eig) # See how many PCs to keep



#######-------- Gradient forest analysis --------#########

#######-------- Explore associations with ecoregions --------#######
# First cleanup community table
clust.wide= clust1 %>%
  pivot_wider(names_from = "cluster",
              values_from = "cluster",
              values_fn = list(cluster = length),
              values_fill=0)
clust.wide=as.data.frame(clust.wide)
ecod=ecod %>% semi_join(clust.wide, by = "Site") 
rownames(ecod)=ecod$Site
ecod$Site=NULL
# Match order between dataframes
clust.wide=clust.wide[ order(match(rownames(clust.wide), rownames(ecod))), ]
env=as.data.frame(ecod)
# Now cleanup environmental table
# Remove nonsense variables and variables that are the same at every site
env$Entero_mean_yearly_range=NULL
env$Entero_min=NULL
env$Nitrogen_min=NULL
env$xx=NULL
env$yy=NULL
clust.wide=clust.wide[,4:22]

# Implement gradient forest 
maxLevel <- log2(0.368*nrow(env)/2)
gfeco <- gradientForest(cbind(env, clust.wide), predictor.vars=colnames(env), response.vars=colnames(clust.wide), ntree=500, maxLevel=maxLevel, trace=T, corr.threshold=0.50)
plot(gfeco, plot.type = "O")
importance(gfeco) # Depth = 0.02672922  
imp=sum_up_importances(gfeco, ecor2) 
imp # Ecoregions = 0.09793522
by.importance <- names(importance(gfeco))[2]
plot(gfeco, plot.type = "C", imp.vars = by.importance, show.overall = T, legend = T, leg.posn = "topleft", leg.nspecies = 5, cex.lab = 0.7, cex.legend = 0.4, cex.axis = 0.6, line.ylab = 0.9, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0, 0.3, 0, 0)))




#######-------- Explore associations with all environmental variables --------#######
enva=read.csv("STXenv_allvars.csv")
# Note that shifted coordinates are corrected for map projections
env1=merge(enva, ecor, by="Site")
ecod=left_join(depth3, env1, by="Site")
ecod$Sample=NULL
ecod$Site=NULL
ecod=aggregate(list(numdup=rep(1,nrow(ecod))),ecod, length)
ecod$numdup=NULL
names(ecod)[2]="Site"
# Remake input matrices
clust.wide= clust1 %>%
  pivot_wider(names_from = "cluster",
              values_from = "cluster",
              values_fn = list(cluster = length),
              values_fill=0)
clust.wide=as.data.frame(clust.wide)
ecod=ecod %>% semi_join(clust.wide, by = "Site") 
rownames(ecod)=ecod$Site
rownames(clust.wide)=clust.wide$Site
ecod$Site=NULL
clust.wide$Site=NULL
# Match order between dataframes
clust.wide=clust.wide[ order(match(rownames(clust.wide), rownames(ecod))), ]
env=as.data.frame(ecod)
# Cleanup environmental table
# Remove nonsense variables and variables that are the same at every site
env$numdup=NULL
env$X=NULL
env$Entero_mean_yearly_range=NULL
env$Entero_min=NULL
env$Nitrogen_min=NULL
env$Latitude=NULL
env$Longitude=NULL
env=env[,1:40]
# Remove redundant coordinate information
clust.wide$xx=NULL
clust.wide$yy=NULL
env$xx=NULL
env$yy=NULL

# Implement gradient forest 
maxLevel <- log2(0.368*nrow(env)/2)
gfeco <- gradientForest(cbind(env, clust.wide), predictor.vars=colnames(env), response.vars=colnames(clust.wide), ntree=500, maxLevel=maxLevel, trace=T, corr.threshold=0.50)
plot(gfeco, plot.type = "O")
importance(gfeco) # Depth =0.0542375834, pH_mean_yearly_range= 0.0132380840, Entero_max=0.0114337200, Phosphorus_mean= 0.0100938414, pH_max=0.0093445185 
by.importance <- names(importance(gfeco))[1:4]
plot(gfeco, plot.type = "C", imp.vars = by.importance, show.overall = T, legend = T, leg.posn = "topleft", leg.nspecies = 5, cex.lab = 0.7, cex.legend = 0.4, cex.axis = 0.6, line.ylab = 0.9, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0, 0.3, 0, 0)))
