#cp -r ~/Project/Bhupinder/10X_Human/_10X-Normal-Epi-Combine/Data/ /wehisan/general/academic/lab_smyth/cheny/Project/BPal/10X_Human/ScientificData/NormEpi/

smythlab
cd cheny/Project/BPal/10X_Human/ScientificData/NormEpi
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
#SamplesComb <- c("N1105","N280","N1B","NE","NF","MH0023Epi","MH0064Epi",
#    "PM0095Epi","PM0342Epi","PM0372Epi","MH275Epi")
Samples <- c("N-1105-epi","N-0280-epi","N-0230.16-epi","N-0408-epi","N-1469-epi","N-0123-epi","N-0064-epi",
    "N-0093-epi","N-0342-epi","N-0372-epi","N-0275-epi")
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

pdf("QC-NormEpi-LibSize.pdf", height=6, width=10)
p <- ggplot(dat, aes(x=Sample, y=log10(seqDepth), fill=Sample, color=Sample)) + 
	geom_violin(trim=TRUE, scale="width", show.legend=FALSE, adjust=1) + 
	#geom_jitter(shape=16, alpha=0.5, position=position_jitter(0.3), size=0.8) + 
	scale_fill_manual(values=rep(col.pLight,2)[1:NS]) + 
	scale_color_manual(values=rep(col.pDark,2)[1:NS]) + 
	theme(legend.position="none") + 
	labs(title="Library size", x="Sample", y="log10(lib.size)")
p
dev.off()

pdf("QC-NormEpi-LibSize.pdf", height=5, width=7)
par(mar=c(7,5,2,1))
boxplot(seqDepth_vec ~ Labels, varwidth=FALSE, xlab="", ylab="Library size\n",log="y", las=2, range=0, 
    col="grey70", main="Normal epi samples", frame=FALSE)
dev.off()

pdf("QC-NormEpi-NGenes.pdf", height=5, width=7)
par(mar=c(7,5,2,1))
boxplot(NGenes_vec ~ Labels, varwidth=FALSE, xlab="", ylab="Detected genes",log="y", las=2, range=0, 
    col="grey70", main="Normal epi samples", frame=FALSE)
dev.off()

save(seqDepth_vec, NGenes_vec, Labels, file="QC-NormEpi.RData")

