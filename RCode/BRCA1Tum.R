smythlab
cd cheny/Project/BPal/10X_Human/ScientificData/BRCA1Tum
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


#####################################
# Fig 4 D-I, Fig 5 A-H

### Samples to be combined
Samples <- c("B1-0894", "B1-0033", "B1-0023", "B1-0090", 
    "TN-B1-0131", "TN-B1-0554", "TN-B1-4031", "TN-B1-0177")
SamplesComb <- gsub("-","_",Samples)
SamplesPreB1 <- SamplesComb[1:4]
SamplesTum <- SamplesComb[5:8]

### Readin the data
DGE <- paste0("dge_", SamplesComb)
DD <- paste0("dd_", SamplesComb)
for(i in 1:length(SamplesComb)) {
    eval(call("<-", as.name(DGE[i]), 
        read10X(path=paste0("../Data/",Samples[i]), DGEList=TRUE)))
}

### Gene annotation
for(i in 1:length(SamplesComb)){
    ann <- alias2SymbolUsingNCBI(get(DGE[i])$genes$Symbol, required.columns=c("GeneID","Symbol"), 
        gene.info.file="../Data/180808_Homo_sapiens.gene_info.gz")
    Genes <- cbind(get(DGE[i])$genes, Official=ann$Symbol, GeneID=ann$GeneID)
    eval(parse(text=paste0(DGE[i],"$genes <- Genes")))
}

### Quality control
# Calculate the percentage of reads from mitochondrial genes and the number of expressed genes (with at least one read) in each cell.
for(i in 1:length(SamplesComb)) {
    mito <- grep("^MT-", get(DGE[i])$genes$Symbol) 
    percent.mito <- colSums(get(DGE[i])$counts[mito, ]) / get(DGE[i])$samples$lib.size
    sp <- cbind( get(DGE[i])$samples, percent.mito=percent.mito, 
        nGenes=colSums(get(DGE[i])$counts!=0) )
    eval(parse(text=paste0(DGE[i],"$samples <- sp")))
}

# Cell filtering (potential doublets)
SampleStats <- read.delim("../Data/SampleStats.txt")
m <- match(Samples, SampleStats$SampleName)
SampleStats <- SampleStats[m,]

mito_upper <- SampleStats$Mito
nGenes_lower <- SampleStats$GeneLower
nGenes_upper <- SampleStats$GeneUpper
lib_upper <- SampleStats$LibSize

for(i in 1:length(SamplesComb)){
    y <- get(DGE[i])
    keep.mito <- y$samples$percent.mito < mito_upper[i]
    keep.nGenes <- y$samples$nGenes > nGenes_lower[i] & y$samples$nGenes < nGenes_upper[i]
    keep.nUMIs <- y$samples$lib.size < lib_upper[i]
    keep <- keep.mito & keep.nGenes & keep.nUMIs
    eval( parse(text=paste0(DGE[i],"<-",DGE[i],"[,keep]")) )
}

# Gene filtering (expressed in fewer than 1% of the cells, invalid symbol, duplicated genes)
for(i in 1:length(SamplesComb)) {
    y <- get(DGE[i])
    colnames(y) <- paste(SamplesComb[i], y$samples$Barcode, sep="_")
    o <- order(rowSums(y$counts), decreasing=TRUE)
    y <- y[o, ]
    keep1 <- rowSums(y$counts > 0) >= ncol(y)*0.01
    keep2 <- !is.na(y$genes$Official)
    keep3 <- !duplicated(y$genes$Official)
    yall <- y[keep2 & keep3, , keep=FALSE]
    rownames(yall) <- yall$genes$Official
    eval( parse(text=paste0(DD[i],"<- yall")) )
    y <- y[keep1 & keep2 & keep3, , keep=FALSE]
    rownames(y) <- y$genes$Official
    eval( parse(text=paste0(DGE[i],"<- y")) )
}


### Individual Seurat object
for(i in 1:length(SamplesComb)){
    # Readin
    y <- get(DGE[i])
    so <- CreateSeuratObject(y$counts, project=SamplesComb[i])
    so <- NormalizeData(so)
    # HVG
    so <- FindVariableFeatures(so, selection.method="vst", nfeatures=1500)
    # Row-scale 
    so <- ScaleData(so)
    so@meta.data$group <- SamplesComb[i]
    eval( parse(text=paste0(SamplesComb[i],"<- so")) )
}

### Integration analysis
# List
CombSeurat <- list()
for(i in 1:length(SamplesComb)) CombSeurat[[i]] <- get(SamplesComb[i])
names(CombSeurat) <- SamplesComb

# Anchors
dimUsed <- 30
Anchors <- FindIntegrationAnchors(object.list=CombSeurat, dims=1:dimUsed,
    anchor.features=1000, scale=TRUE, k.anchor=5, k.filter=30,
    k.score=20, max.features=100)

# Integration
BRCA1Tum <- IntegrateData(anchorset=Anchors, dims=1:dimUsed, k.weight=100)

# Dim reduction
DefaultAssay(BRCA1Tum) <- "integrated"
BRCA1Tum <- ScaleData(BRCA1Tum, verbose=FALSE)
BRCA1Tum <- RunPCA(BRCA1Tum, npcs=dimUsed, verbose=FALSE)
BRCA1Tum <- RunTSNE(BRCA1Tum, dims=1:dimUsed, seed.use=2018)
tSNE <- BRCA1Tum@reductions$tsne@cell.embeddings

# Fig 4D
Group <- factor(BRCA1Tum@meta.data$group, levels=SamplesComb)
col <- col.p[Group]
plotOrd <- sample(ncol(BRCA1Tum))
pdf("Fig4D.pdf", height=9, width=9)
plot(tSNE[plotOrd,], pch=16, col=col[plotOrd], cex=0.7, xlab="tSNE-1", ylab="tSNE-2",
    main="t-SNE - By sample")
dev.off()

# Condition
col.p2 <- c("darkblue", "gold2")
Group2 <- c("Preneo", "Tumour")[Group %in% SamplesTum + 1L]
Group2 <- factor(Group2, levels=c("Preneo", "Tumour"))

# Cell clustering
resolution <- 0.15
BRCA1Tum <- FindNeighbors(BRCA1Tum, dims=1:dimUsed, verbose=FALSE)
BRCA1Tum <- FindClusters(BRCA1Tum, resolution=resolution, verbose=FALSE)
Cluster <- as.integer(BRCA1Tum@meta.data$seurat_clusters)
ncls <- length(table(Cluster))

# Fig 4E
col <- col.p[Cluster]
pdf("Fig4E.pdf", height=9, width=9)
plot(tSNE, pch=16, col=col, cex=0.7, xlab="tSNE-1", ylab="tSNE-2",
    main="t-SNE - By cluster")
dev.off()

# Fig 4H
col <- col.p2[Group2]
pdf("Fig4H.pdf", height=9, width=9)
plot(tSNE[plotOrd,], pch=16, col=col[plotOrd], cex=0.7, xlab="tSNE-1", ylab="tSNE-2",
    main="t-SNE - By condition")
dev.off()

# Combine raw data
allGenes <- rownames(get(DD[1]))
allGenesFilter <- c()
for(i in DGE) allGenesFilter <- c(allGenesFilter, rownames(get(i)))
allGenesFilter <- names(table(allGenesFilter))[table(allGenesFilter) >= 2]
u <- get(DD[1])
d <- u[allGenes, ]
y <- u[allGenesFilter, ]
for(i in 2:length(SamplesComb)) {
    u <- get(DD[i])
    d <- cbind(d, u[allGenes, ])
    y <- cbind(y, u[allGenesFilter, ])
}
d$samples$group <- y$samples$group <- Group

# log-CPM
prior.count <- 1
z <- edgeR::cpm(d, log=TRUE, prior.count=prior.count)

load("../Data/Human-PosSigGenes.RData")
z_Basal <- z[rownames(z) %in% Basal, ]
z_LP <- z[rownames(z) %in% LP, ]
z_ML <- z[rownames(z) %in% ML, ]
z_Str <- z[rownames(z) %in% Str, ]
score_Basal <- colMeans(z_Basal)
score_LP <- colMeans(z_LP)
score_ML <- colMeans(z_ML)
score_Str <- colMeans(z_Str)
titles <- paste0(c("Basal", "LP", "ML", "Stroma"), " signatures")
dat <- data.frame( cbind(score_Basal, score_LP, score_ML, score_Str) )

# Fig 4I
pdf("Fig4I.pdf", height=7, width=10)
par(mfrow=c(2,2))
for(i in 1:4)
    boxplot(dat[,i] ~ Cluster, col=col.p[1:ncls], range=0, main=titles[i], medlwd = 1)
dev.off()

#saveRDS(BRCA1Tum, file="SeuratObject_BRCA1Tum.rds")
#save(z, d, y, SamplesComb, SamplesPreB1, SamplesTum, Group, Group2, ncls, Cluster, col.p, col.p2, tSNE, file="BRCA1Tum.RData")



####################################################################################################
### Micro-environment
Sub <- !(Cluster %in% c(1,6,9,10))

cellNamesSub <- rownames(BRCA1Tum@meta.data)[Sub]
DGESub <- paste0("dge_sub_", SamplesComb)
DDSub <- paste0("dd_sub_", SamplesComb)
for(i in 1:length(SamplesComb)) {
    d <- get(DD[i])
    d <- d[, colnames(d) %in% cellNamesSub]
    keep1 <- rowSums(d$counts > 0) >= ncol(d)*0.01
    eval( parse(text=paste0(DDSub[i],"<- d")) )
    d <- d[keep1, , keep=FALSE]
    eval( parse(text=paste0(DGESub[i],"<- d")) )
}

# HVGs for each subsetted sample
SamplesCombSub <- paste0(SamplesComb, "_Sub")
CombSeuratSub <- list()
for(i in 1:length(SamplesCombSub)) {
    d <- get(DGESub[i])
    so <- CreateSeuratObject(d$counts, project=SamplesCombSub[i])
    so <- NormalizeData(so)
    so <- FindVariableFeatures(so, selection.method="vst", nfeatures=1500)
    so <- ScaleData(so)
    so@meta.data$group <- SamplesComb[i]
    CombSeuratSub[[i]] <- so
    eval( parse(text=paste0(SamplesCombSub[i],"<- so")) )
}
names(CombSeuratSub) <- SamplesCombSub

# Integration
dimUsed <- 30
AnchorsSub <- FindIntegrationAnchors(object.list=CombSeuratSub, dims=1:dimUsed,
    anchor.features=1000, scale=TRUE, k.anchor=5, k.filter=30,
    k.score=20, max.features=100)
BRCA1TumSub <- IntegrateData(anchorset=AnchorsSub, dims=1:dimUsed, k.weight=100)

# Dimensional reduction
dimUsed <- 30
DefaultAssay(BRCA1TumSub) <- "integrated"
BRCA1TumSub <- ScaleData(BRCA1TumSub, verbose=FALSE)
BRCA1TumSub <- RunPCA(BRCA1TumSub, npcs=dimUsed, verbose=FALSE)
BRCA1TumSub <- RunTSNE(BRCA1TumSub, dims=1:dimUsed, seed.use=2018)
tSNE <- BRCA1TumSub@reductions$tsne@cell.embeddings

GroupSub <- factor(BRCA1TumSub@meta.data$group, levels=SamplesComb)
plotOrd <- sample(ncol(BRCA1TumSub))

col.p2 <- c("darkblue", "gold2")
GroupSub2 <- c("Preneo", "Tumour")[GroupSub %in% SamplesTum + 1L]
GroupSub2 <- factor(GroupSub2, levels=c("Preneo", "Tumour"))

# Fig 5B
pdf("Fig5B.pdf", height=9, width=9)
plot(tSNE[plotOrd,], pch=16, col=col.p2[GroupSub2][plotOrd], cex=0.7, xlab="tSNE-1", ylab="tSNE-2",
    main="t-SNE - By condition")
dev.off()

# clustering
resolution <- 0.1
BRCA1TumSub <- FindNeighbors(BRCA1TumSub, dims=1:dimUsed, verbose=FALSE)
BRCA1TumSub <- FindClusters(BRCA1TumSub, resolution=resolution, verbose=FALSE)
ClusterSub <- as.integer(BRCA1TumSub@meta.data$seurat_clusters)
ncls <- length(table(ClusterSub))

# Fig 5A
pdf("Fig5A.pdf", height=9, width=9)
col <- col.p[ClusterSub]
plot(tSNE, pch=16, col=col, cex=0.7, xlab="tSNE-1", ylab="tSNE-2", main="t-SNE - By cluster")
dev.off()

# Fig 5E-left
pdf("Fig5E-left.pdf", height=7, width=14)
par(mfrow=c(1,2))
plot(tSNE[GroupSub2==levels(GroupSub2)[1],], pch=16, col=col.p2[1], cex=0.7, xlab="tSNE-1", ylab="tSNE-2",
    main="Preneo")
plot(tSNE[GroupSub2==levels(GroupSub2)[2],], pch=16, col=col.p2[2], cex=0.7, xlab="tSNE-1", ylab="tSNE-2",
    main="Tumour")
dev.off()

# Combine raw data
allGenes <- rownames(get(DDSub[1]))
allGenesFilter <- c()
for(i in DGESub) allGenesFilter <- c(allGenesFilter, rownames(get(i)))
allGenesFilter <- names(table(allGenesFilter))[table(allGenesFilter) >= 2]
u <- get(DDSub[1])
d <- u[allGenes, ]
y <- u[allGenesFilter, ]
for(i in 2:length(SamplesComb)) {
    u <- get(DDSub[i])
    d <- cbind(d, u[allGenes, ])
    y <- cbind(y, u[allGenesFilter, ])
}
d$samples$group <- y$samples$group <- GroupSub

# log-CPM
prior.count <- 1
z <- edgeR::cpm(d, log=TRUE, prior.count=prior.count)

#save(z, d, y, SamplesCombSub, GroupSub, GroupSub2, ncls, ClusterSub, col.p, col.p2, tSNE, file="BRCA1TumSub.RData")
#saveRDS(BRCA1TumSub, file="SeuratObject_BRCA1TumSub.rds")

cellNum <- table(ClusterSub, GroupSub2)
cellProp <- t( t(cellNum) / rowSums(t(cellNum)) )

# Fig 5E-right
pdf("Fig5E-right.pdf", height=6, width=8)
barplot(t(cellProp*100), beside=TRUE, main="", xlab="Cluster (Pre | Post)", ylab="Cell number Percentage", 
    col=rep(col.p2, ncls), names=1:ncls, space=c(0,0.3))
dev.off()

# Marker genes
markers <- list()
for(i in 1:ncls){
    markers[[i]] <- FindMarkers(BRCA1TumSub, ident.1=i-1, min.pct=0.1, only.pos=TRUE)
}
names(markers) <- paste0("Cluster", 1:ncls)

top <- 20
nTop <- pmin(top, sapply(markers, nrow))
topMarkers <- c()
for(i in 1:ncls) topMarkers <- c(topMarkers, rownames(markers[[i]][1:nTop[i],]))
topMarkers <- topMarkers[!duplicated(topMarkers)]
topMarkers

# Fig 5H
pdf("Fig5H.pdf", height=12, width=10)
mat <- t(scale(t(z[topMarkers, ])))
ord <- order(ClusterSub, GroupSub2)
annot <- data.frame(Condition=GroupSub2[ord], Cluster=paste0("Cluster ", ClusterSub[ord]))
rownames(annot) <- colnames(mat)[ord]
ann_colors <- list(Cluster=col.p[1:ncls], Condition=col.p2)
names(ann_colors$Cluster) <- paste0("Cluster ", 1:ncls)
names(ann_colors$Condition) <- levels(GroupSub2)
library(pheatmap)
pheatmap(mat[,ord], color=colorRampPalette(c("blue","white","red"))(100), 
    breaks=seq(-2,2,length.out=101), cluster_cols=FALSE, scale="none", 
    labels_col=rep("",ncol(mat)), fontsize_row=7, 
    clustering_method="ward.D2", annotation_col=annot, annotation_colors=ann_colors)
dev.off()


