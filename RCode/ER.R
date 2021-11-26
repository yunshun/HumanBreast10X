### R v3.6.1
library(Seurat)
library(edgeR)
library(vcd)
library(ggplot2)
library(pheatmap)

col.pMedium <- scater:::.get_palette("tableau10medium")
col.pDark <- scater:::.get_palette("tableau20")[2*(1:10)-1]
col.pLight <- scater:::.get_palette("tableau20")[2*(1:10)]
col.p <- c(col.pDark, col.pLight)


#####################################
# Fig EV3, 6

### Samples to be combined
Samples <- c("ER-0114-T3", "ER-0360", "ER-0167-T", "ER-0151", "ER-0032", "ER-0125",
    "ER-0043-T", "ER-0025", "ER-0001", "ER-0042", "ER-0319", "ER-0040-T", "ER-0163")
SamplesComb <- gsub("-","_",Samples)

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
mito_upper[7] <- 0.2
nGenes_upper[7] <- 4500

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
ERTotal <- IntegrateData(anchorset=Anchors, dims=1:dimUsed, k.weight=100)

# Dim reduction
DefaultAssay(ERTotal) <- "integrated"
ERTotal <- ScaleData(ERTotal, verbose=FALSE)
ERTotal <- RunPCA(ERTotal, npcs=dimUsed, verbose=FALSE)
ERTotal <- RunTSNE(ERTotal, dims=1:dimUsed, seed.use=2018)
tSNE <- ERTotal@reductions$tsne@cell.embeddings

# Fig EV 2D
Group <- factor(ERTotal@meta.data$group, levels=SamplesComb)
col <- col.p[Group]
plotOrd <- sample(ncol(ERTotal))
pdf("FigEV2D.pdf", height=9, width=9)
plot(tSNE[plotOrd,], pch=16, col=col[plotOrd], cex=0.7, xlab="tSNE-1", ylab="tSNE-2",
    main="t-SNE - By sample")
dev.off()

# Cell clustering
resolution <- 0.1
ERTotal <- FindNeighbors(ERTotal, dims=1:dimUsed, verbose=FALSE)
ERTotal <- FindClusters(ERTotal, resolution=resolution, verbose=FALSE)
Cluster <- as.integer(ERTotal@meta.data$seurat_clusters)
ncls <- length(table(Cluster))

# Fig 6C
col <- col.p[Cluster]
pdf("Fig6C.pdf", height=9, width=9)
plot(tSNE, pch=16, col=col, cex=0.7, xlab="tSNE-1", ylab="tSNE-2",
    main="t-SNE - By cluster")
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
z <- edgeR::cpm(y, log=TRUE, prior.count=prior.count)

load("../Data/Human-PosSigGenes.RData")
z_Basal <- z[rownames(z) %in% Basal, ]
z_LP <- z[rownames(z) %in% LP, ]
z_ML <- z[rownames(z) %in% ML, ]
score_Basal <- colMeans(z_Basal)
score_LP <- colMeans(z_LP)
score_ML <- colMeans(z_ML)
titles <- c("Basal", "LP", "ML")
dat <- data.frame(score_Basal, score_LP, score_ML)

# Boxplots Fig EV3 A
pdf("FigEV3A-ER.pdf", height=9, width=4)
par(mfrow=c(3,1))
for(j in 1:3)
boxplot(dat[,j] ~ Cluster, varwidth=FALSE, xlab="Cluster", ylab="AveLogCPM", range=0, medlwd = 1, 
    col=col.p[1:ncls], main=titles[j])
dev.off()

# log-CPM
prior.count <- 1
z <- edgeR::cpm(d, log=TRUE, prior.count=prior.count)

#saveRDS(ERTotal, file="SeuratObject_ERTotal.rds")
#save(z, d, y, SamplesComb, Group, ncls, Cluster, col.p, tSNE, file="ERTotal.RData")


####################################################################################################
### Micro-environment
Sub <- !(Cluster %in% c(1,5,6,7))

cellNamesSub <- rownames(ERTotal@meta.data)[Sub]
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
ERTotalSub <- IntegrateData(anchorset=AnchorsSub, dims=1:dimUsed, k.weight=100)

# Dimensional reduction
dimUsed <- 30
DefaultAssay(ERTotalSub) <- "integrated"
ERTotalSub <- ScaleData(ERTotalSub, verbose=FALSE)
ERTotalSub <- RunPCA(ERTotalSub, npcs=dimUsed, verbose=FALSE)
ERTotalSub <- RunTSNE(ERTotalSub, dims=1:dimUsed, seed.use=2018)
tSNE <- ERTotalSub@reductions$tsne@cell.embeddings

GroupSub <- factor(ERTotalSub@meta.data$group, levels=SamplesComb)
plotOrd <- sample(ncol(ERTotalSub))

# clustering
resolution <- 0.1
ERTotalSub <- FindNeighbors(ERTotalSub, dims=1:dimUsed, verbose=FALSE)
ERTotalSub <- FindClusters(ERTotalSub, resolution=resolution, verbose=FALSE)
ClusterSub <- as.integer(ERTotalSub@meta.data$seurat_clusters)
ncls <- length(table(ClusterSub))

col.p3 <- col.p[c(1,2,4,9,18,8,17,6,10,15,3,12,13)]

# Fig 7C
pdf("Fig7C.pdf", height=9, width=9)
col <- col.p3[ClusterSub]
plot(tSNE, pch=16, col=col, cex=0.7, xlab="tSNE-1", ylab="tSNE-2", main="t-SNE - By cluster")
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

#save(z, d, y, SamplesCombSub, GroupSub, ncls, ClusterSub, col.p, col.p3, tSNE, file="ERTotalSub.RData")
#saveRDS(ERTotalSub, file="SeuratObject_ERTotalSub.rds")


############
# Pseudo-bulk DE analysis
sampClust <- paste(y$samples$group, ClusterSub, sep="_Clst")
counts2 <- t(rowsum(t(y$counts), group=sampClust))
yComb <- DGEList(counts2)
yComb$samples$Patient <- gsub("_Clst.*$","",colnames(yComb))
yComb$samples$Cluster <- as.numeric(gsub("^.*_Clst","",colnames(yComb)))
yComb$samples$group <- yComb$samples$Cluster

# Focus on all clusters
N <- 1:ncls
yClstSub <- yComb[, yComb$samples$Cluster %in% N]

# Filter out some pseudo samples with low library size.
sel <- yClstSub$samples$lib.size > 3e4
yClstSub2 <- yClstSub[,sel]

# Gene filtering
keep <- filterByExpr(yClstSub2, min.count=10, min.total.count=20)
yClstSub2 <- yClstSub2[keep,,keep=FALSE]

# TMM
yClstSub2 <- calcNormFactors(yClstSub2)
sn <- paste0(rep(SamplesComb, each=length(N)), "_Clst", N)
sel <- sn %in% colnames(yClstSub2)
sn <- sn[sel]
m <- match(sn, colnames(yClstSub2))
yClstSub2 <- yClstSub2[,m]

# design matrix
Cls <- as.factor(yClstSub2$samples$Cluster)
Pat <- factor(yClstSub2$samples$Patient, levels=SamplesComb)
design <- model.matrix(~ Cls + Pat)
colnames(design)

# Estimate NB dispersion and QL dispersion.
yClstSub2 <- estimateDisp(yClstSub2, design=design)
qfit2 <- glmQLFit(yClstSub2, design)

# log-CPM
prior.count <- 1
zClstSub2 <- edgeR::cpm(yClstSub2, log=TRUE, prior.count=prior.count)

# Contrast
contr <- rbind( matrix(1/(1-ncls), ncls, ncls),
                matrix(0, ncol(design)-ncls, ncls) )
diag(contr) <- 1
contr[1,] <- 0
rownames(contr) <- colnames(design)
colnames(contr) <- paste0("Cls", 1:ncls)
ctest <- list()
for(i in 1:ncls) ctest[[i]] <- glmQLFTest(qfit2, contrast=contr[,i])

top <- 30
pseudoMakers <- list()
for(i in 1:ncls) {
    ord <- order(ctest[[i]]$table$PValue, decreasing=FALSE)
    upreg <- ctest[[i]]$table$logFC > 0
    pseudoMakers[[i]] <- rownames(yClstSub2)[ord[upreg][1:top]]
}
Markers <- unlist(pseudoMakers)
Markers <- Markers[!duplicated(Markers)]

annot2 <- data.frame(Cluster=paste0("Cluster ", Cls), Patient=Pat)
rownames(annot2) <- colnames(zClstSub2)
ann_colors2 <- list(Cluster=col.p[N], Patient=col.p[1:length(SamplesComb)])
names(ann_colors2$Cluster) <- paste0("Cluster ", N)
names(ann_colors2$Patient) <- SamplesComb

# Fig 7C right
pdf("Fig7C-right.pdf", height=15, width=9)
mat4 <- t(scale(t(zClstSub2[Markers, ])))
pheatmap(mat4, color=colorRampPalette(c("blue","white","red"))(100), border_color="NA",
    breaks=seq(-2,2,length.out=101), cluster_cols=TRUE, scale="none", fontsize_row=4,
    show_colnames=FALSE, treeheight_row=70, treeheight_col=70, cutree_cols=1,
    clustering_method="ward.D2", annotation_col=annot2, annotation_colors=ann_colors2)
dev.off()



#####################################
### T-cell population
# Fig 
N <- c(1,8)
TC <- ClusterSub %in% N

cellNamesTC <- rownames(ERTotalSub@meta.data)[TC]
DGETC <- paste0("dge_tc_", SamplesComb)
DDTC <- paste0("dd_tc_", SamplesComb)
for(i in 1:length(SamplesComb)) {
    d <- get(DD[i])
    d <- d[, colnames(d) %in% cellNamesTC]
    keep1 <- rowSums(d$counts > 0) >= ncol(d)*0.01
    eval( parse(text=paste0(DDTC[i],"<- d")) )
    d <- d[keep1, , keep=FALSE]
    eval( parse(text=paste0(DGETC[i],"<- d")) )
}

# HVGs
SamplesCombTC <- paste0(SamplesComb, "_TC")
CombSeuratTC <- list()
for(i in 1:length(SamplesCombTC)) {
    d <- get(DGETC[i])
    so <- CreateSeuratObject(d$counts, project=SamplesCombTC[i])
    so <- NormalizeData(so)
    so <- FindVariableFeatures(so, selection.method="vst", nfeatures=1500)
    so <- ScaleData(so)
    so@meta.data$group <- SamplesComb[i]
    CombSeuratTC[[i]] <- so
    eval( parse(text=paste0(SamplesCombTC[i],"<- so")) )
}
names(CombSeuratTC) <- SamplesCombTC

# keep samples with more than 200 T-cells
keep.samples <- which(sapply(CombSeuratTC, ncol) >=200)
keep.samples
CombSeuratTC <- CombSeuratTC[keep.samples]

# Integration
dimUsed <- 30
AnchorsTC <- FindIntegrationAnchors(object.list=CombSeuratTC, dims=1:dimUsed,
    anchor.features=1000, scale=TRUE, k.anchor=5, k.filter=30,
    k.score=20, max.features=100)
ERTotalTC <- IntegrateData(anchorset=AnchorsTC, dims=1:dimUsed, k.weight=100)

# Dimensional reduction
dimUsed <- 20
DefaultAssay(ERTotalTC) <- "integrated"
ERTotalTC <- ScaleData(ERTotalTC, verbose=FALSE)
ERTotalTC <- RunPCA(ERTotalTC, npcs=dimUsed, verbose=FALSE)
ERTotalTC <- RunTSNE(ERTotalTC, dims=1:dimUsed, seed.use=2018)
tSNE <- ERTotalTC@reductions$tsne@cell.embeddings

GroupTC <- factor(ERTotalTC@meta.data$group, levels=SamplesComb)
plotOrd <- sample(ncol(ERTotalTC))

# Clustering
resolution <- 0.2
ERTotalTC <- FindNeighbors(ERTotalTC, dims=1:dimUsed, verbose=FALSE)
ERTotalTC <- FindClusters(ERTotalTC, resolution=resolution, verbose=FALSE)
ClusterTC <- as.integer(ERTotalTC@meta.data$seurat_clusters)
ncls <- length(table(ClusterTC))

# Fig EV 4A
pdf("FigEV4A-ERTotal.pdf", height=9, width=9)
col <- col.p[ClusterTC]
plot(tSNE, pch=16, col=col, cex=0.7, xlab="tSNE-1", ylab="tSNE-2", main="t-SNE - By cluster")
dev.off()

# Combine raw data
allGenes <- rownames(get(DDTC[keep.samples][1]))
allGenesFilter <- c()
for(i in DGETC[keep.samples]) allGenesFilter <- c(allGenesFilter, rownames(get(i)))
allGenesFilter <- names(table(allGenesFilter))[table(allGenesFilter) >= 2]
u <- get(DDTC[keep.samples][1])
d <- u[allGenes, ]
y <- u[allGenesFilter, ]
for(i in 2:length(SamplesComb[keep.samples])) {
    u <- get(DDTC[keep.samples][i])
    d <- cbind(d, u[allGenes, ])
    y <- cbind(y, u[allGenesFilter, ])
}
d$samples$group <- y$samples$group <- GroupTC

# log-CPM
prior.count <- 1
z <- edgeR::cpm(d, log=TRUE, prior.count=prior.count)

#save(z, d, y, SamplesCombTC, GroupTC, ncls, ClusterTC, col.p, tSNE, file="ERTotalTC.RData")
#saveRDS(ERTotalTC, file="SeuratObject_ERTotalTC.rds")


############
# Pseudo-bulk DE analysis
sampClust <- paste(y$samples$group, ClusterTC, sep="_Clst")
counts2 <- t(rowsum(t(y$counts), group=sampClust))
yComb <- DGEList(counts2)
yComb$samples$Patient <- gsub("_Clst.*$","",colnames(yComb))
yComb$samples$Cluster <- as.numeric(gsub("^.*_Clst","",colnames(yComb)))
yComb$samples$group <- yComb$samples$Cluster

N <- 1:ncls
yClstSub2 <- yComb[, yComb$samples$Cluster %in% N]

# Filter out some pseudo samples with low library size.
sel <- yClstSub2$samples$lib.size > 1e4
table(sel)
yClstSub2 <- yClstSub2[,sel]

# Filter out lowly expressed genes.
keep <- filterByExpr(yClstSub2, min.count=7, min.total.count=15)
table(keep)
yClstSub2 <- yClstSub2[keep,,keep=FALSE]

# TMM
yClstSub2 <- calcNormFactors(yClstSub2)
sn <- paste0(rep(SamplesComb, each=length(N)), "_Clst", N)
sel <- sn %in% colnames(yClstSub2)
sn <- sn[sel]
m <- match(sn, colnames(yClstSub2))
yClstSub2 <- yClstSub2[,m]
yClstSub2$samples

# design
Cls <- as.factor(yClstSub2$samples$Cluster)
Pat <- factor(yClstSub2$samples$Patient, levels=SamplesComb[keep.samples])
design <- model.matrix(~ Cls + Pat)
colnames(design)

# dispersion
yClstSub2 <- estimateDisp(yClstSub2, design=design)
qfit2 <- glmQLFit(yClstSub2, design)

# Test
qftest <- glmQLFTest(qfit2, coef=2:ncls)
top <- 100
ord <- order(qftest$table$PValue, decreasing=FALSE)
topDE <- rownames(yClstSub2)[ ord[1:top] ]

prior.count <- 1
zClstSub2 <- edgeR::cpm(yClstSub2, log=TRUE, prior.count=prior.count)

annot2 <- data.frame(Cluster=paste0("Cluster ", Cls), Patient=Pat)
rownames(annot2) <- colnames(zClstSub2)
ann_colors2 <- list(Cluster=col.p[N], Patient=col.p[1:length(SamplesComb)][keep.samples])
names(ann_colors2$Cluster) <- paste0("Cluster ", N)
names(ann_colors2$Patient) <- SamplesComb[keep.samples]
mat3 <- t(scale(t(zClstSub2[topDE, ])))

pdf("FigEV4C-ERTotal.pdf", height=14, width=8)
pheatmap(mat3, color=colorRampPalette(c("blue","white","red"))(100), border_color="NA",
    breaks=seq(-2,2,length.out=101), cluster_cols=TRUE, scale="none", fontsize_row=6, 
    show_colnames=FALSE, treeheight_row=70, treeheight_col=70, cutree_cols=1,
    clustering_method="ward.D2", annotation_col=annot2, annotation_colors=ann_colors2)
dev.off()



#####################################
### Tumour population
# Fig 
Tum <- Cluster %in% c(1,5)

cellNamesTum <- rownames(ERTotal@meta.data)[Tum]
DGETum <- paste0("dge_tum_", SamplesComb)
DDTum <- paste0("dd_tum_", SamplesComb)
for(i in 1:length(SamplesComb)) {
    d <- get(DD[i])
    d <- d[, colnames(d) %in% cellNamesTum]
    keep1 <- rowSums(d$counts > 0) >= ncol(d)*0.01
    eval( parse(text=paste0(DDTum[i],"<- d")) )
    d <- d[keep1, , keep=FALSE]
    eval( parse(text=paste0(DGETum[i],"<- d")) )
}

# HVGs
SamplesCombTum <- paste0(SamplesComb, "_Tum")
CombSeuratTum <- list()
for(i in 1:length(SamplesCombTum)) {
    d <- get(DGETum[i])
    so <- CreateSeuratObject(d$counts, project=SamplesCombTum[i])
    so <- NormalizeData(so)
    so <- FindVariableFeatures(so, selection.method="vst", nfeatures=1500)
    so <- ScaleData(so)
    so@meta.data$group <- SamplesComb[i]
    CombSeuratTum[[i]] <- so
    eval( parse(text=paste0(SamplesCombTum[i],"<- so")) )
}
names(CombSeuratTum) <- SamplesCombTum

# Integration
dimUsed <- 30
AnchorsTum <- FindIntegrationAnchors(object.list=CombSeuratTum, dims=1:dimUsed,
    anchor.features=1000, scale=TRUE, k.anchor=5, k.filter=30,
    k.score=20, max.features=100)
ERTotalTum <- IntegrateData(anchorset=AnchorsTum, dims=1:dimUsed, k.weight=100)

# Dimensional reduction
dimUsed <- 30
DefaultAssay(ERTotalTum) <- "integrated"
ERTotalTum <- ScaleData(ERTotalTum, verbose=FALSE)
ERTotalTum <- RunPCA(ERTotalTum, npcs=dimUsed, verbose=FALSE)
ERTotalTum <- RunTSNE(ERTotalTum, dims=1:dimUsed, seed.use=2018)
tSNE <- ERTotalTum@reductions$tsne@cell.embeddings

GroupTum <- factor(ERTotalTum@meta.data$group, levels=SamplesComb)
plotOrd <- sample(ncol(ERTotalTum))

# Clustering
resolution <- 0.05
ERTotalTum <- FindNeighbors(ERTotalTum, dims=1:dimUsed, verbose=FALSE)
ERTotalTum <- FindClusters(ERTotalTum, resolution=resolution, verbose=FALSE)
ClusterTum <- as.integer(ERTotalTum@meta.data$seurat_clusters)
ncls <- length(table(ClusterTum))

# Fig 6E
pdf("Fig6E.pdf", height=9, width=9)
col <- col.p[ClusterTum]
plot(tSNE, pch=16, col=col, cex=0.7, xlab="tSNE-1", ylab="tSNE-2", main="t-SNE - By cluster")
dev.off()

# Combine raw data
allGenes <- rownames(get(DDTum[1]))
allGenesFilter <- c()
for(i in DGETum) allGenesFilter <- c(allGenesFilter, rownames(get(i)))
allGenesFilter <- names(table(allGenesFilter))[table(allGenesFilter) >= 2]
u <- get(DDTum[1])
d <- u[allGenes, ]
y <- u[allGenesFilter, ]
for(i in 2:length(SamplesComb)) {
    u <- get(DDTum[i])
    d <- cbind(d, u[allGenes, ])
    y <- cbind(y, u[allGenesFilter, ])
}
d$samples$group <- y$samples$group <- GroupTum

# log-CPM
prior.count <- 1
z <- edgeR::cpm(d, log=TRUE, prior.count=prior.count)

#save(z, d, y, SamplesCombTum, GroupTum, ncls, ClusterTum, col.p, tSNE, file="ERTotalTum.RData")
#saveRDS(ERTotalTum, file="SeuratObject_ERTotalTum.rds")


############
# Pseudo-bulk DE analysis

sampClust <- paste(y$samples$group, ClusterTum, sep="_Clst")
counts2 <- t(rowsum(t(y$counts), group=sampClust))
yComb <- DGEList(counts2)
yComb$samples$Patient <- gsub("_Clst.*$","",colnames(yComb))
yComb$samples$Cluster <- as.numeric(gsub("^.*_Clst","",colnames(yComb)))
yComb$samples$group <- yComb$samples$Cluster
N <- 1:ncls
yClstSub2 <- yComb[, yComb$samples$Cluster %in% N]

# Filter out some pseudo samples with low library size.
sel <- yClstSub2$samples$lib.size > 3e4
table(sel)
yClstSub2 <- yClstSub2[,sel]

# Filter out lowly expressed genes.
keep <- filterByExpr(yClstSub2, min.count=10, min.total.count=20)
yClstSub2 <- yClstSub2[keep,,keep=FALSE]

yClstSub2 <- calcNormFactors(yClstSub2)
sn <- paste0(rep(SamplesComb, each=length(N)), "_Clst", N)
sel <- sn %in% colnames(yClstSub2)
sn <- sn[sel]
m <- match(sn, colnames(yClstSub2))
yClstSub2 <- yClstSub2[,m]

# design
Cls <- as.factor(yClstSub2$samples$Cluster)
Pat <- factor(yClstSub2$samples$Patient, levels=SamplesComb)
design <- model.matrix(~ Cls + Pat)
colnames(design)

# dispersion
yClstSub2 <- estimateDisp(yClstSub2, design=design)
qfit2 <- glmQLFit(yClstSub2, design)

# Contrast
contr <- makeContrasts(
    Cls1 = -(Cls2+Cls3)/2,
    Cls2 = Cls2-Cls3/2,
    Cls3 = Cls3-Cls2/2, levels=design)
ctest1 <- glmQLFTest(qfit2, contrast=contr[,1])
ctest2 <- glmQLFTest(qfit2, contrast=contr[,2])
ctest3 <- glmQLFTest(qfit2, contrast=contr[,3])

# KEGG
library(org.Hs.eg.db)
geneid <- mapIds(org.Hs.eg.db, keys=rownames(qfit2), column="ENTREZID", 
    keytype="SYMBOL")
keg1 <- kegga(ctest1, geneid=geneid, species="Hs")
keg2 <- kegga(ctest2, geneid=geneid, species="Hs")
keg3 <- kegga(ctest3, geneid=geneid, species="Hs")
topKEGG(keg1, sort="up", n=20)
topKEGG(keg2, sort="up", n=20)
topKEGG(keg3, sort="up", n=20)


prior.count <- 1
zClstSub2 <- edgeR::cpm(yClstSub2, log=TRUE, prior.count=prior.count)

annot2 <- data.frame(Cluster=paste0("Cluster ", Cls), Patient=Pat)
rownames(annot2) <- colnames(zClstSub2)
ann_colors2 <- list(Cluster=col.p[N], Patient=col.p[1:length(SamplesComb)])
names(ann_colors2$Cluster) <- paste0("Cluster ", N)
names(ann_colors2$Patient) <- SamplesComb

PAM50 <- read.delim("../Data/PAM50.txt", header=TRUE, stringsAsFactors=FALSE)
keep <- PAM50$Gene %in% rownames(zClstSub2)
PAM50 <- PAM50[keep, ]
mat5 <- t(scale(t(zClstSub2[PAM50$Gene, ])))

# Fig 6F PAM50 Heat map
pdf("Fig6F.pdf", height=8, width=7)
pheatmap(mat5, color=colorRampPalette(c("blue","white","red"))(100), border_color="NA",
    breaks=seq(-2,2,length.out=101), cluster_cols=TRUE, scale="none", fontsize_row=7, 
    show_colnames=FALSE, treeheight_row=70, treeheight_col=70, cutree_cols=1,
    clustering_method="ward.D2", annotation_col=annot2, annotation_colors=ann_colors2)
dev.off()

# We obtained TCGA signatures from the DE analyses of the TCGA data.
# Boxplots were made for visulization.
load("../Data/TCGA_Signatures.RData")
z_Basal <- zClstSub2[rownames(zClstSub2) %in% TCGA_Basal, ]
z_Her2 <- zClstSub2[rownames(zClstSub2) %in% TCGA_Her2, ]
z_LumA <- zClstSub2[rownames(zClstSub2) %in% TCGA_LumA, ]
z_LumB <- zClstSub2[rownames(zClstSub2) %in% TCGA_LumB, ]
score_Basal <- colMeans(z_Basal)
score_Her2 <- colMeans(z_Her2)
score_LumA <- colMeans(z_LumA)
score_LumB <- colMeans(z_LumB)

# Fig EV 3E
pdf("FigEV3E.pdf", height=7, width=7)
par(mfrow=c(2,2))
boxplot(score_Basal ~ Cls, xlab="", ylab="AveLogCPM", range=0, medlwd = 1, 
    col=col.p[N], main="Basal", frame=FALSE)
boxplot(score_Her2 ~ Cls, xlab="", ylab="AveLogCPM", range=0, medlwd = 1, 
    col=col.p[N], main="Her2", frame=FALSE)
boxplot(score_LumA ~ Cls, xlab="", ylab="AveLogCPM", range=0, medlwd = 1, 
    col=col.p[N], main="LumA", frame=FALSE)
boxplot(score_LumB ~ Cls, xlab="", ylab="AveLogCPM", range=0, medlwd = 1, 
    col=col.p[N], main="LumB", frame=FALSE)
dev.off()


