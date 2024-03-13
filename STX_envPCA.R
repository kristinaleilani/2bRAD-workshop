setwd("~/Documents/STX") 
library(ggplot2)
library(vegan)
library(dplyr)
library(caret)
library(corrplot)
###----  Import environmental data
depth=read.csv("STX_Sample_depth.csv")
depth$Sample=NULL
depth=aggregate(list(numdup=rep(1,nrow(depth))),depth, length)
depth$numdup=NULL
env0=read.csv("STX_env.csv")
rownames(env0)=env0$Site
#env1=merge(env0, depth, by="Site")
#rownames(env1) <- make.names(env1$Site, unique = TRUE)
env1=env0
env1$Site=NULL
env1$Lat=NULL
env1$Lon=NULL
env1$Entero_min=NULL
env1$Entero_mean_yearly_range=NULL

# scale variables
env_scale=as.data.frame(scale(env1))
# plot correlation of variables
M0<-cor(env_scale)
corrplot(M0, method="color")

corr=cor(env_scale)
# Hierarchical Clustering with hclust
hc <- hclust(as.dist(1-corr), "ave")
# Plot the result
plot(hc)
abline(h=0.22,col="red") # choose a cutoff threshold 
cuts=cutree(hc,h=0.22)
goods=c();i=1
for (i in unique(cuts)) {
  goods[i]=names(cuts[cuts==i])[1]
}
length(goods)  # how many variables are left?
# subsetting for only the retained samples
envt_scale=as.data.frame(t(env_scale))
envt_scale=envt_scale[goods,] 
# plot correlations of retained variables
env_scale=as.data.frame(t(envt_scale))
M0<-cor(env_scale)
corrplot(M0, method="color")

# Look at environmental correlations among sites
scorr=cor(envt_scale)
hc <- hclust(as.dist(1-scorr), "ave")
# Plot the result as hierarchical clustering tree
plot(hc)
# And plot as PCA
a=capscale(as.dist(1-scorr)~1)
plot(a, scaling=1)
# eigenvalues: how many good eigenvalues (exceeding broken stick expectation) are there:
plot(a$CA$eig,ylab="proportion of variance explained",xlab="eigenvector #")
plot(a$CA$eig/sum(a$CA$eig),ylab="proportion of variance explained",xlab="eigenvector #")
lines(bstick(length(a$CA$eig)),col="red") # "broken stick" model
# Check other MDS axes
ords=scores(a,display="sites")
axes2plot=c(1:4) # which axes to plot
scores=data.frame(a$CA$u[,axes2plot])
scores=cbind(scores, colnames(envt_scale))
colnames(scores)[5]<-"Site"
ggplot(scores,aes(scores[,2],scores[,3], asp=1, label = Site)) + 
  geom_label() +
  theme_bw()+
  xlab(names(scores)[2])+
  ylab(names(scores)[3])+
  guides(size = "none")





