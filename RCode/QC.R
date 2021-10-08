smythlab
cd cheny/Project/BPal/10X_Human/ScientificData/QC
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
Samples <- c("N-0092-total", "N-0019-total", "N-0280-epi", "N-0093-epi", "N-0093-total", "N-1469-epi", "N-0408-epi", 
    "N-1105-epi", "N-0230.17-total", "N-0064-epi", "N-0064-total", "N-0230.16-epi", "N-0233-total", "N-0169-total", 
    "N-0123-epi", "N-0123-total", "N-0342-epi", "N-0342-total", "N-0288-total", "N-0021-total", "N-0275-epi", 
    "N-0275-total", "N-0372-epi", "N-0372-total", "B1-0894", "B1-0033", "B1-0023", "B1-0090", "TN-0126", "TN-0135", 
    "TN-0106", "TN-0114-T2", "TN-B1-4031", "TN-B1-0131", "TN-B1-0554", "TN-B1-0177", "HER2-0308", "HER2-0337", 
    "HER2-0031", "HER2-0069", "HER2-0161", "HER2-0176", "ER-0319", "ER-0001", "ER-0125", "ER-0360", "ER-0114-T3", 
    "ER-0032", "ER-0042", "ER-0025", "ER-0151", "ER-0163", "ER-0029-7C", "ER-0029-9C", "ER-0040-T", "ER-0040-LN", 
    "ER-0043-T", "ER-0043-LN", "ER-0056-T", "ER-0056-LN", "ER-0064-T", "ER-0064-LN", "ER-0167-T", "ER-0167-LN", 
    "ER-0173-T", "ER-0173-LN", "mER-0178", "mER-0068-T", "mER-0068-LN")
SamplesComb <- gsub("-","_",Samples)

Type <- c("Normal", "BRCA1 Preneo", "TN", "TN BRCA1", "HER2+", "ER+")
Group <- rep(Type, times=c(24, 4, 4, 4, 6, 27))
Group <- factor(Group, levels=Type)

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
col <- col.pLight[Group]

pdf("QC.pdf", height=7.5, width=13)
#par(mfrow=c(2,1))
layout(matrix(c(1,2,3), 3, 1, byrow = TRUE),
   heights=c(3,3,1))
par(mar=c(7.5,5,2,1))
boxplot(seqDepth_vec ~ Labels, varwidth=FALSE, xlab="", ylab="",log="y", las=2, range=0, medlwd=1, 
    col=col, main="Library size", frame=FALSE)
boxplot(NGenes_vec ~ Labels, varwidth=FALSE, xlab="", ylab="",log="y", las=2, range=0, medlwd=1, 
    col=col, main="Detected genes", frame=FALSE)
par(mar=c(0,0,0,0))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend = Type, pch=15, pt.cex=3, cex=1.5, bty='n',  col = col.pLight[1:length(Type)], ncol=3)
dev.off()

#save(seqDepth_vec, NGenes_vec, Labels, col, Type, Group, Samples, col.pLight, file="QC.RData")


SampleStats <- read.delim("../Data/SampleStats.txt")
m <- match(Samples, SampleStats$SampleName)
SampleStats <- SampleStats[m,]
nd <- SampleStats$CellNum
nfd <- SampleStats$CellNum - SampleStats$CellNumAfter
names(nd) <- Samples

### CellNumber
pdf("CellNumber.pdf", height=5, width=15)
yupper <- 2.3e4
ylim <- c(0, yupper)
par(mar=c(7.5,5,2,1))
bp <- barplot(rbind(nd-nfd, nfd), space=0.3, beside=FALSE, col=c("lightblue", "salmon"), 
    main="Number of cells", legend.text=c("Cells retained", "Cell filtered"), 
    args.legend=list(x="topright",bty="n"), las=2, border=NA, 
    axes=FALSE, ylab="", ylim=ylim)
axis(side=2, at=5e3*(0:4), las=2)
text(x=bp, y=nd + yupper/20, labels=paste0(round(100*(nfd/nd),digits=1),"%"), xpd=TRUE, cex=0.7, srt=90)
dev.off()



############################################
### Combine altogether

load("../Data/QC.RData")
SampleStats <- read.delim("../Data/SampleStats.txt")
m <- match(Samples, SampleStats$SampleName)
SampleStats <- SampleStats[m,]
nd <- SampleStats$CellNum
nfd <- SampleStats$CellNum - SampleStats$CellNumAfter
names(nd) <- Samples

source("fig_label.R")

pdf("../Figures/QC.pdf", height=13, width=13)
layout(matrix(c(1,2,3,4,5), 5, 1, byrow = TRUE),
   heights=c(3, 0.2, 3, 1, 3))
par(mar=c(9,5,2,2))
boxplot(seqDepth_vec ~ Labels, varwidth=FALSE, xlab="", ylab="",log="y", las=2, range=0, medlwd=1, 
    col=col, main="Library size (before filtering)", frame=FALSE)

par(mar=c(0,0,0,0))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
fig_label("b.", cex=2.5)

par(mar=c(8,5,2,2))
boxplot(NGenes_vec ~ Labels, varwidth=FALSE, xlab="", ylab="",log="y", las=2, range=0, medlwd=1, 
    col=col, main="Detected genes (before filtering)", frame=FALSE)
fig_label("a.", cex=2.5, pos="topleft")

par(mar=c(0,0,0,0))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend = Type, pch=15, pt.cex=3, cex=1.5, bty='n',  col = col.pLight[1:length(Type)], ncol=3)

### CellNumber
yupper <- 2.3e4
ylim <- c(0, yupper)
par(mar=c(7.5,5,2,2))
bp <- barplot(rbind(nd-nfd, nfd), space=0.3, beside=FALSE, col=c("lightblue", "salmon"), 
    main="Number of cells", legend.text=c("Cells retained", "Cell filtered"), 
    args.legend=list(x="topright",bty="n"), las=2, border=NA, 
    axes=FALSE, ylab="", ylim=ylim)
axis(side=2, at=5e3*(0:4), las=2)
text(x=bp, y=nd + yupper/20, labels=paste0(round(100*(nfd/nd),digits=1),"%"), xpd=TRUE, cex=0.7, srt=90)
fig_label("c.", cex=2.5)
dev.off()

