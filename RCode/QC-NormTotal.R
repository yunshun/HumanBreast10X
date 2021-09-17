smythlab
cd cheny/Project/BPal/10X_Human/ScientificData/NormTotal
module load R/3.6.1
module load python/3.7.0

R

#library(scater)
#library(scran)
library(Seurat)
library(edgeR)
library(vcd)
#library(gplots)
library(ggplot2)
library(pheatmap)

col.pMedium <- scater:::.get_palette("tableau10medium")
col.pDark <- scater:::.get_palette("tableau20")[2*(1:10)-1]
col.pLight <- scater:::.get_palette("tableau20")[2*(1:10)]
col.p <- c(col.pDark, col.pLight)


### Samples to be combined
#SamplesComb <- c("PM0019","PM0233","PM0092","PM0230","MH0169","PM0095Total",
#    "MH0023Total","MH0064Total","PM0342Total","PM0372Total",
#    "MH0021Total","MH275Total","MH288Total")
Samples <- c("N-0019-total","N-0233-total","N-0092-total","N-0230.17-total","N-0169-total","N-0093-total",
    "N-0123-total","N-0064-total","N-0342-total","N-0372-total","N-0021-total","N-0275-total","N-0288-total")
SamplesComb <- gsub("-","_",Samples)


### Readin the data
DGE <- paste0("dge_", SamplesComb)
DD <- paste0("dd_", SamplesComb)
for(i in 1:length(SamplesComb)) {
    eval(call("<-", as.name(DGE[i]), 
        read10X(path=paste0("../Data/",Samples[i]), DGEList=TRUE)))
}

seqDepth <- NGenes <- list()
for(i in 1:length(SamplesComb)) {
    y <- get(DGE[i])
    seqDepth[[i]] <- y$samples$lib.size
    NGenes[[i]] <- colSums(y$counts!=0)
}
names(seqDepth) <- names(NGenes) <- SamplesComb

Labels <- rep(Samples, times=sapply(seqDepth, length))
Labels <- factor(Labels, levels=Samples)
seqDepth_vec <- unlist(seqDepth)
NGenes_vec <- unlist(NGenes)
dat <- data.frame(Sample=Labels, seqDepth=seqDepth_vec, NGenes=NGenes_vec)
NS <- length(Samples)

##### Plots

pdf("QC-NormTotal-LibSize.pdf", height=5, width=7)
par(mar=c(7.5,5,2,1))
boxplot(seqDepth_vec ~ Labels, varwidth=FALSE, xlab="", ylab="Library size\n",log="y", las=2, range=0, 
    col="grey70", main="Normal total samples", frame=FALSE)
dev.off()

pdf("QC-NormTotal-NGenes.pdf", height=5, width=7)
par(mar=c(7.5,5,2,1))
boxplot(NGenes_vec ~ Labels, varwidth=FALSE, xlab="", ylab="Detected genes",log="y", las=2, range=0, 
    col="grey70", main="Normal total samples", frame=FALSE)
dev.off()

#save(seqDepth_vec, NGenes_vec, Labels, file="QC-NormTotal.RData")


