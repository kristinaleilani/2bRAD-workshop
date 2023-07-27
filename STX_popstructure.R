########  STX population genetics
library(dplyr)
library(vegan)
library(ggplot2)
library(pheatmap)
setwd("~/Documents/STX") 
# Change all species names in this file (Bonnet, Agaricia, or Recruits) to yours (PSTR, MCAV, OFAV, SSID, AAGA, or PAST)

# Import sample names
samples=read.table("bams.nr")
names(samples)[1]<-"Sample"
samples$Sample=paste(sub(".bam","",samples$Sample),sep="")
rownames(samples)<-samples$Sample

# Import IBS matrix
IBS=as.matrix(read.table("Bonnet.ibsMat"))
dimnames(IBS)=list(rownames(samples),rownames(samples))






#######----- Exploring genetic structure of your population -----#######


#--- Plot a hierarchical tree
hc=hclust(as.dist(IBS),"ave")
plot(hc,cex=0.5) # clustering of samples by IBS (great to detect closely related individuals)
abline(h=0.15,col="red") # this seems like a reasonable  "low-hang" threshold for calling related groups
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






#######----- Exploring relatedness in your population -----#######

# reading long relatedness table, output of ngsRelate
rel=read.table("Recruits.res",sep="\t",header=T)

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

plot(hclust(as.dist(1-relm),method="ave"))
pheatmap(as.dist(1-relm))






#######----- Exploring admixture in your population -----#######


# For K=2 (assuming two genetic clusters)
dir="~/Documents/STX" # path to input files
inName="Agaricia_k2.qopt" # name of the input file to plot, output of ngsAdmix run
pops="pops.txt" # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list.
npops=as.numeric(sub(".+(\\d+)\\..+","\\1",inName))
tbl=read.table(paste(inName,sep=""),header=F)
i2p=read.table(paste(pops,sep=""),fill = T, header=T)
names(i2p)=c("ind","pop")
tbl=cbind(tbl,i2p)
row.names(tbl)=tbl$ind
head(tbl)

qm=as.qmatrix(as.matrix(read.table("Agaricia_k2.qopt")))
bp=barplot(qm,border=NA,space=0)
axis(1, at = 1:nrow(qm), labels = i2p$pop[bp$order], las = 3, cex.axis = .4) 


# For K=3 (assuming three genetic clusters)
dir="~/Documents/STX" # path to input files
inName="Agaricia_k3.qopt" # name of the input file to plot, output of ngsAdmix run
pops="pops.txt" # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list.
npops=as.numeric(sub(".+(\\d+)\\..+","\\1",inName))
tbl=read.table(paste(inName,sep=""),header=F)
i2p=read.table(paste(pops,sep=""),fill = T, header=T)
names(i2p)=c("ind","pop")
tbl=cbind(tbl,i2p)
row.names(tbl)=tbl$ind
head(tbl)


my.colors <- c("blue4","cyan3","hotpink","gold","orange")
my.palette <- CreatePalette(my.colors, 8)
qm=as.qmatrix(as.matrix(read.table("Agaricia_k3.qopt")))
bp=barplot(qm,border=NA,space=0,col.palette = my.palette)
axis(1, at = 1:nrow(qm), labels = i2p$pop[bp$order], las = 3, cex.axis = .4) 


