########  STX population genetics
library(dplyr)
library(vegan)
library(ggplot2)
library(pheatmap)
setwd("~/Documents/STX") 
# If you're using windows and set your working directory by copying the path, then you will need to change the directions of the backslashes in your path here.
# Change all species names in this file (ex. change "STX" to your species: PSTR, MCAV, OFAV, SSID, AAGA, or PAST)

# Import sample names
samples=read.table("bams.nr")
names(samples)[1]<-"Sample"
samples$Sample=paste(sub(".bam","",samples$Sample),sep="")
rownames(samples)<-samples$Sample

# Import IBS matrix
IBS=as.matrix(read.table("STX2.ibsMat"))
dimnames(IBS)=list(rownames(samples),rownames(samples))






#######----- Exploring genetic structure of your population -----#######


#--- Plot a hierarchical tree
hc=hclust(as.dist(IBS),"ave")
plot(hc,cex=0.5) # clustering of samples by IBS (great to detect closely related individuals)
abline(h=0.15,col="red") # this seems like a reasonable  "low-hang" threshold for calling related groups

# Save your dendrogram plot without labels
png("STX_IBS_dendrogram.png", width = 8, height = 8, units = 'in', res = 300)
plot(hc,cex=0.0005) 
dev.off()

# If you tree seems to have clones, remove them here. If not, move onto plotting PCoA.
cuts=cutree(hc,h=0.1)
goods=c();i=1
for (i in unique(cuts)) {
  goods[i]=names(cuts[cuts==i])[1]
}
length(goods)  # how many samples are left?
# subsetting all data for only the retained samples
IBS=IBS[goods,goods] 
samples=as.data.frame(samples[goods,])
names(samples)[1]<-"Sample"


#--- Plot a PCoA
ord=capscale(as.dist(IBS)~1)
plot(ord,scaling=1) # Plot PCoA

# eigenvalues: how many good eigenvalues (exceeding broken stick expectation) are there:
plot(ord$CA$eig,ylab="proportion of variance explained",xlab="eigenvector #")
plot(ord$CA$eig/sum(ord$CA$eig),ylab="proportion of variance explained",xlab="eigenvector #")
lines(bstick(length(ord$CA$eig)),col="red") # "broken stick" model: expected random distribution of eigenvalues
# looks like we have 2 good eigenvalues


#--- Plot a heatmap 
pheatmap(IBS)
# Save your heatmap plot
png("STX_IBS_heatmap.png", width = 8, height = 8, units = 'in', res = 300)
pheatmap(IBS)
dev.off()







#######----- Exploring relatedness in your population -----#######

# reading long relatedness table, output of ngsRelate
rel=read.table("STX.res",sep="\t",header=T)

# creating an empty square matrix
relm=matrix(0,nrow=length(unique(rel$a))+1,ncol=length(unique(rel$a))+1)

# filling up the square matrix with entries from "rab" column
for (a in unique(rel$a)) {
  for (b in unique(rel$b)) {
    if (b<=a) {next}
    relm[a,b]=relm[b,a]=rel[rel$a==a & rel$b==b,"rab"]
  }
}
diag(relm)=1
# adding names to columns and rows - assuming ngsRelate was run on angsd result obtained for bams.qc file
bams=scan("bams.qc",what="character") # list of bam files
bams=sub(".bam","",bams)
dimnames(relm)=list(bams,bams)

# View a dendrogram of relatedness
plot(hclust(as.dist(1-relm),method="ave"))
# Save your dendrogram of relatedness
png("STX_relate_dendrogram.png", width = 8, height = 8, units = 'in', res = 300)
plot(hclust(as.dist(1-relm),method="ave"),cex=0.0005) # clustering of samples by IBS (great to detect closely related individuals)
dev.off()

# View a heatmap of relatedness
pheatmap(as.dist(1-relm))
# Save your heatmap of relatedness
png("STX_relate_heatmap.png", width = 8, height = 8, units = 'in', res = 300)
pheatmap(as.dist(relm))
dev.off()






#######----- Exploring admixture in your population -----#######

q2<-read.table("STX2.qopt") # name of the input file to plot, output of ngsAdmix run
pop<-read.table("pops.txt") # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list.
names(pop)=c("ind","pop")
ord<-order(pop$pop) #order it by population
barplot(t(q2)[,ord],
        col=1:2,
        names=pop$pop[ord],
        las=2,
        space=0,
        border=NA,
        xlab="Site",
        ylab="Admixture proportions for K=2")

# Save your barplot with k=2 clusters
png("STX_admix_k2.png", width = 10, height = 4, units = 'in', res = 300)
barplot(t(q2)[,ord],
        col=1:2,
        names=pop$pop[ord],
        las=2,
        space=0,
        border=NA,
        xlab="Site",
        ylab="Admixture proportions for K=2")
dev.off()

# If it looks like you have two main admixture groups, let's save these assignments:
npops=ncol(q2)
cluster.admix=apply(q2[,1:npops],1,function(x) {return(paste(which(x>0.5),collapse=".")) })
#save(cluster.admix,file=paste("Bonnet_clusters.RData",sep=""))

# Let's plot a PCoA colored by two admixture groups:
ord=capscale(IBS~1)
summary(ord) # Check how much genetic variation are explained by MDS axes
# Save your summary as a text file
s <- summary(ord)
capture.output(s, file = "STX_IBS_capscale.txt")

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
  geom_label(label=rownames(scores))+
  guides(size = "none")
# Save your PCoA without labels, colored by admixture group:
ggplot(scores,aes(scores[,1],scores[,2], asp=1, fill=cluster.admix)) + 
  geom_point(aes(size=1, colour = cluster.admix), alpha=0.85) +
  theme_bw()+
  coord_fixed()+
  xlab(names(scores)[1])+
  ylab(names(scores)[2])+
  guides(size = "none")
ggsave("STX_PCoA_admix.tiff", units="in", width=9, height=5, dpi=300, compression = 'lzw')

# Get separate samples lists of each admixture group
admix1 <- subset(rownames(scores), scores$cluster.admix == 1) #22 samples
admix1 <- paste0(admix1, ".bam")
write.table(admix1, "STX_cluster1", sep="\t", col.names = F, row.names = F)
admix2 <- subset(rownames(scores), scores$cluster.admix == 2) #33 samples
admix2 <- paste0(admix2, ".bam")
write.table(admix2, "STX_cluster2", sep="\t", col.names = F, row.names = F)





#######-----  Next, let's explore a k=3 scenario...
q3<-read.table("STX3.qopt") # name of the input file to plot, output of ngsAdmix run
pop<-read.table("pops.txt") # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list.
names(pop)=c("ind","pop")
ord<-order(pop$pop) #order it by population
barplot(t(q3)[,ord],
        col=1:3,
        names=pop$pop[ord],
        las=2,
        space=0,
        border=NA,
        xlab="Site",
        ylab="Admixture proportions for K=3")
# Save your barplot with k=2 clusters
png("STX_admix_k3.png", width = 10, height = 4, units = 'in', res = 300)
barplot(t(q3)[,ord],
        col=1:2,
        names=pop$pop[ord],
        las=2,
        space=0,
        border=NA,
        xlab="Site",
        ylab="Admixture proportions for K=3")
dev.off()
# If it looks like you have two main admixture groups, let's save these assignments:
npops=ncol(q3)
cluster.admix=apply(q3[,1:npops],1,function(x) {return(paste(which(x>0.33),collapse=".")) })
# You can plot a PCoA colored by 3 admixture groups as well.


# Also take a look at admixture scenarios for 4, 5, and 6 clusters


