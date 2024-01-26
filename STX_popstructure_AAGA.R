########  AAGA population genetics
library(dplyr)
library(vegan)
library(ggplot2)
library(pheatmap)
library(tess3r)
setwd("~/Documents/STX") 
# Change all species names in this file (Bonnet, Agaricia, or Recruits) to yours (PSTR, MCAV, OFAV, SSID, AAGA, or PAST)

# Import sample names
samples=read.table("bams.nr2")
names(samples)[1]<-"Sample"
samples$Sample=paste(sub(".bam","",samples$Sample),sep="")
rownames(samples)<-samples$Sample

# Import IBS matrix
IBS=as.matrix(read.table("AAGA1.ibsMat"))
dimnames(IBS)=list(rownames(samples),rownames(samples))






#######----- Exploring genetic structure of your population -----#######


#--- Plot a hierarchical tree
hc=hclust(as.dist(IBS),"ave")
png("AAGA_IBS_dendrogram.png", width = 8, height = 8, units = 'in', res = 300)
plot(hc,cex=0.0005) # clustering of samples by IBS (great to detect closely related individuals)
dev.off()



# If you have many clones, we will remove them for creating the PCoA and then project them back in.

# ======= STEP 1: looking for least similar samples, to be removed for good

toomuch=c()
minIBS=c()
IBS1=IBS
diag(IBS1)=1
for(i in 1:round(ncol(IBS)*0.75)){
  message(i)
  mxcor=apply(IBS1,1,min)
  b=sample(which(mxcor==min(mxcor)),1)
  minIBS=append(minIBS,mxcor[b])
  IBS1=IBS1[-b,-b]
  toomuch=append(toomuch,names(b))
}
plot(minIBS,type="l")
# minIBS.0=minIBS[-length(minIBS)]
# minIBS.1=minIBS[-1]
# maxIBSdiff=minIBS.1-minIBS.0
# # plot(scale(maxIBSdiff),type="l")
# abline(h=2,col="red")
#plot(maxIBSdiff,type="l")
# cutname=maxIBSdiff[scale(maxIBSdiff)>2][length(maxIBSdiff[scale(maxIBSdiff)>2])]

# CHOOSE CUTOFF for witholding samples from PCoA contruction (to be projected back later)
cut=25
plot(minIBS,type="l")
abline(v=cut,col="red")
verybads=which(colnames(IBS) %in% toomuch[1:cut])
toodiff=colnames(IBS)[verybads]

IBS.nd=IBS[-verybads,-verybads]
plot(hclust(as.dist((1-IBS.nd)),method="ave"),cex=0.0001,main="dissimilar ones removed")

#======= STEP 2: looking for too-similar samples, to withhold them from PCoA construction (but back-projecting them later)

toolittle=c()
maxIBS=c()
IBS1=IBS.nd
diag(IBS1)=0
for(i in 1:round(ncol(IBS.nd)*0.75)){
  message(i)
  mxcor=apply(IBS1,1,max)
  b=sample(which(mxcor==max(mxcor)),1)
  maxIBS=append(maxIBS,mxcor[b])
  IBS1=IBS1[-b,-b]
  toolittle=append(toolittle,names(b))
}

plot(maxIBS,type="l")
# CHOOSE CUTOFF for witholding samples from PCoA contrustion (to be projected back later)
cut=75
plot(maxIBS,type="l")
abline(v=cut,col="red")
bads=which(colnames(IBS.nd) %in% toolittle[1:cut])

# making PCoA out of remaining samples
IBS.sel=IBS.nd
plot(hclust(as.dist((1-IBS.sel)),method="ave"),cex=0.0001)
caps.sel=capscale((1-IBS.sel)~1)
plot(caps.sel,scaling=1)

# projecting all samples
caps.pred=predict(caps.sel,(1-IBS),type='sp',scaling="sites")
plot(caps.pred,asp=1)


ord=caps.sel
# eigenvalues: how many good eigenvalues (exceeding broken stick expectation) are there:
plot(ord$CA$eig,ylab="proportion of variance explained",xlab="eigenvector #")
plot(ord$CA$eig/sum(ord$CA$eig),ylab="proportion of variance explained",xlab="eigenvector #")
lines(bstick(length(ord$CA$eig)),col="red") # "broken stick" model: expected random distribution of eigenvalues
# looks like we have 2 good eigenvalues


#--- Plot a heatmap 
png("AAGA_IBS_heatmap.png", width = 8, height = 8, units = 'in', res = 300)
pheatmap(IBS)
dev.off()





#######----- Exploring relatedness in your population -----#######

# reading long relatedness table, output of ngsRelate
rel=read.table("AAGA.res",sep="\t",header=T)

# creatig an empty square matrix
relm=matrix(0,nrow=length(unique(rel$a))+1,ncol=length(unique(rel$a))+1)

# filling up the square matrix with entries from "rab" column
for (a in unique(rel$a)) {
  for (b in unique(rel$b)) {
    if (b<=a) {next}
    relm[a,b]=relm[b,a]=rel[rel$a==a & rel$b==b,"rab"]
  }
}
diag(relm)=1
# adding names to columns and rows - assuming ngsRelate was run on angsd result obtained for bams.nr file
bams=scan("bams.qc",what="character") # list of bam files
bams=sub(".bam","",bams)
dimnames(relm)=list(bams,bams)

png("AAGA_relate_dendrogram.png", width = 8, height = 8, units = 'in', res = 300)
plot(hclust(as.dist(1-relm),method="ave"),cex=0.0005) # clustering of samples by IBS (great to detect closely related individuals)
dev.off()



png("AAGA_relate_heatmap.png", width = 8, height = 8, units = 'in', res = 300)
pheatmap(as.dist(relm))
dev.off()







#######----- Exploring admixture in your population -----#######
source("plot_admixture_v5_function.R")
npops <- 2 #choose number of populations
inName <- paste0('AAGA', npops, '.qopt') #inName corresponds to .qopt files. Change to your species name.

i2p=read.csv("STX_sitesample.csv")
i2p=i2p[i2p$Sample %like% "AAGA", ] # Subset your species
pop=left_join(samples, i2p, by="Sample")
names(pop)=c("ind", "pop")

npops=as.numeric(sub("\\D+(\\d+)\\..+","\\1",inName))
q=read.table(paste(inName,sep=""),header=F)
tbl=cbind(q,pop)
rownames(tbl) <- tbl$ind <- sub("(.*?)\\..*$", "\\1", tbl$ind)
head(tbl,20) # this is how the resulting dataset must look
tbl$pop=factor(tbl$pop,levels = c("Fredericksted Pier", "Butler Bay", "Carambola deep", "Carambola shallow",
                                  "North Star", "Cane Bay deep", "Cane Bay", "Columbus landing",
                                  "The Palms", "WAPA", "Deep end", "Joe's Reef", "Channel Rock", "Isaac Bay")) #Set order of sites to match the map.
#png("AAGA_admix_k2.png", width = 10, height = 4, units = 'in', res = 300)
plotAdmixture(data=tbl,npops=npops,grouping.method="distance")
#dev.off()

# You can also take a look at admixture scenarios for 3, 4, 5, and 6 clusters

# If it looks like you have two main admixture groups, let's save these assignments:
npops=ncol(q2)
cluster.admix=apply(q2[,1:npops],1,function(x) {return(paste(which(x>0.5),collapse=".")) })
#save(cluster.admix,file=paste("AAGA_clusters.RData",sep=""))

# Let's plot a PCoA colored by two admixture groups:
summary(ord) # Check how much genetic variation are explained by MDS axes
# Save your summary as a text file
s <- summary(ord)
capture.output(s, file = "AAGA_IBS_capscale.txt")

ords=scores(ord,display="sites")
axes2plot=c(1:4) # which axes to plot
scores=data.frame(ord$CA$u[,axes2plot])
scores=cbind(scores, cluster.admix)
scores$cluster.admix<-as.factor(scores$cluster.admix)

# Plot a PCoA colored by admixture group:
ggplot(scores,aes(scores[,1],scores[,2], asp=1, fill=cluster.admix)) + 
  geom_point(aes(size=1, colour = cluster.admix)) +
  theme_bw()+
  coord_fixed()+
  xlab(names(scores)[1])+
  ylab(names(scores)[2])+
  scale_color_manual(values=c("tomato", "lightblue"))+
  geom_label(label=rownames(scores))+
  guides(size = "none")
# Save your PCoA without labels, colored by admixture group:
ggplot(scores,aes(scores[,1],scores[,2], asp=1, fill=cluster.admix)) + 
  geom_point(aes(size=1, colour = cluster.admix), alpha=0.85) +
  theme_bw()+
  coord_fixed()+
  xlab(names(scores)[1])+
  ylab(names(scores)[2])+
  scale_color_manual(values=c("tomato", "lightblue"))+
  guides(size = "none")
ggsave("AAGA_PCoA_admix.tiff", units="in", width=9, height=5, dpi=300, compression = 'lzw')

# Get separate samples lists of each admixture group
admix1 <- subset(rownames(scores), scores$cluster.admix == 1) #22 samples
admix1 <- paste0(admix1, ".bam")
write.table(admix1, "AAGA_cluster1", sep="\t", col.names = F, row.names = F)
admix2 <- subset(rownames(scores), scores$cluster.admix == 2) #33 samples
admix2 <- paste0(admix2, ".bam")
write.table(admix2, "AAGA_cluster2", sep="\t", col.names = F, row.names = F)
