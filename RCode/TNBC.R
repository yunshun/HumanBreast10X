#cp -r ~/Project/Bhupinder/10X_Human/_10X-Normal-Total-Combine/Data/ /wehisan/general/academic/lab_smyth/cheny/Project/BPal/10X_Human/ScientificData/NormTotal/

smythlab
cd cheny/Project/BPal/10X_Human/ScientificData/TNBC
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
# Fig 6 A, 7 A, EV2 A

### Samples to be combined
Samples <- c("TN-B1-0554", "TN-B1-0177", "TN-0135", "TN-B1-4031", 
    "TN-B1-0131", "TN-0126", "TN-0106", "TN-0114-T2")
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
TNBC <- IntegrateData(anchorset=Anchors, dims=1:dimUsed, k.weight=100)

# Dim reduction
DefaultAssay(TNBC) <- "integrated"
TNBC <- ScaleData(TNBC, verbose=FALSE)
TNBC <- RunPCA(TNBC, npcs=dimUsed, verbose=FALSE)
TNBC <- RunTSNE(TNBC, dims=1:dimUsed, seed.use=2018)
tSNE <- TNBC@reductions$tsne@cell.embeddings


# Fig EV 2B
Group <- factor(TNBC@meta.data$group, levels=SamplesComb)
col <- col.p[Group]
plotOrd <- sample(ncol(TNBC))
pdf("FigEV2B.pdf", height=9, width=9)
plot(tSNE[plotOrd,], pch=16, col=col[plotOrd], cex=0.7, xlab="tSNE-1", ylab="tSNE-2",
    main="t-SNE - By sample")
dev.off()

# Cell clustering
resolution <- 0.1
TNBC <- FindNeighbors(TNBC, dims=1:dimUsed, verbose=FALSE)
TNBC <- FindClusters(TNBC, resolution=resolution, verbose=FALSE)
Cluster <- as.integer(TNBC@meta.data$seurat_clusters)
ncls <- length(table(Cluster))

# Fig 6A
col <- col.p[Cluster]
pdf("Fig6A.pdf", height=9, width=9)
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

# Fig EV3 A
load("../Data/Human-PosSigGenes.RData")
z_Basal <- z[rownames(z) %in% Basal, ]
z_LP <- z[rownames(z) %in% LP, ]
z_ML <- z[rownames(z) %in% ML, ]
score_Basal <- colMeans(z_Basal)
score_LP <- colMeans(z_LP)
score_ML <- colMeans(z_ML)
titles <- c("Basal", "LP", "ML")
dat <- data.frame(score_Basal, score_LP, score_ML)

pdf("FigEV3A-TNBC.pdf", height=9, width=4)
par(mfrow=c(3,1))
for(j in 1:3)
boxplot(dat[,j] ~ Cluster, varwidth=FALSE, xlab="Cluster", ylab="AveLogCPM", range=0, medlwd = 1, 
    col=col.p[1:ncls], main=titles[j])
dev.off()

# log-CPM
prior.count <- 1
z <- edgeR::cpm(d, log=TRUE, prior.count=prior.count)

#saveRDS(TNBC, file="SeuratObject_TNBC.rds")
#save(z, d, y, SamplesComb, Group, ncls, Cluster, col.p, tSNE, file="TNBC.RData")



####################################################################################################
### Micro-environment
Sub <- !(Cluster %in% c(1,3))

cellNamesSub <- rownames(TNBC@meta.data)[Sub]
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
TNBCSub <- IntegrateData(anchorset=AnchorsSub, dims=1:dimUsed, k.weight=100)

# Dimensional reduction
dimUsed <- 30
DefaultAssay(TNBCSub) <- "integrated"
TNBCSub <- ScaleData(TNBCSub, verbose=FALSE)
TNBCSub <- RunPCA(TNBCSub, npcs=dimUsed, verbose=FALSE)
TNBCSub <- RunTSNE(TNBCSub, dims=1:dimUsed, seed.use=2018)
tSNE <- TNBCSub@reductions$tsne@cell.embeddings

GroupSub <- factor(TNBCSub@meta.data$group, levels=SamplesComb)
plotOrd <- sample(ncol(TNBCSub))

# clustering
resolution <- 0.136
TNBCSub <- FindNeighbors(TNBCSub, dims=1:dimUsed, verbose=FALSE)
TNBCSub <- FindClusters(TNBCSub, resolution=resolution, verbose=FALSE)
ClusterSub <- as.integer(TNBCSub@meta.data$seurat_clusters)
ncls <- length(table(ClusterSub))

# Fig 7A left
pdf("Fig7A-left.pdf", height=9, width=9)
col <- col.p[ClusterSub]
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

#save(z, d, y, SamplesCombSub, GroupSub, ncls, ClusterSub, col.p, tSNE, file="TNBCSub.RData")
#saveRDS(TNBCSub, file="SeuratObject_TNBCSub.rds")


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
sel <- yClstSub$samples$lib.size > 1e4
yClstSub2 <- yClstSub[,sel]
yClstSub2$samples

# Gene filtering
keep <- filterByExpr(yClstSub2, min.count=7, min.total.count=15)
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

# Fig 7A right
pdf("Fig7A-right.pdf", height=15, width=9)
mat4 <- t(scale(t(zClstSub2[Markers, ])))
pheatmap(mat4, color=colorRampPalette(c("blue","white","red"))(100), border_color="NA",
    breaks=seq(-2,2,length.out=101), cluster_cols=TRUE, scale="none", fontsize_row=4,
    show_colnames=FALSE, treeheight_row=70, treeheight_col=70, cutree_cols=1,
    clustering_method="ward.D2", annotation_col=annot2, annotation_colors=ann_colors2)
dev.off()



#####################################
### T-cell population
# Fig 
N <- c(1,5)
TC <- ClusterSub %in% N

cellNamesTC <- rownames(TNBCSub@meta.data)[TC]
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
TNBCTC <- IntegrateData(anchorset=AnchorsTC, dims=1:dimUsed, k.weight=100)

# Dimensional reduction
dimUsed <- 20
DefaultAssay(TNBCTC) <- "integrated"
TNBCTC <- ScaleData(TNBCTC, verbose=FALSE)
TNBCTC <- RunPCA(TNBCTC, npcs=dimUsed, verbose=FALSE)
TNBCTC <- RunTSNE(TNBCTC, dims=1:dimUsed, seed.use=2018)
tSNE <- TNBCTC@reductions$tsne@cell.embeddings

GroupTC <- factor(TNBCTC@meta.data$group, levels=SamplesComb)
plotOrd <- sample(ncol(TNBCTC))

# Clustering
resolution <- 0.2
TNBCTC <- FindNeighbors(TNBCTC, dims=1:dimUsed, verbose=FALSE)
TNBCTC <- FindClusters(TNBCTC, resolution=resolution, verbose=FALSE)
ClusterTC <- as.integer(TNBCTC@meta.data$seurat_clusters)
ncls <- length(table(ClusterTC))

# Fig EV 4A
pdf("FigEV4A.pdf", height=9, width=9)
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

#save(z, d, y, SamplesCombTC, GroupTC, ncls, ClusterTC, col.p, tSNE, file="TNBCTC.RData")
#saveRDS(TNBCTC, file="SeuratObject_TNBCTC.rds")



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

pdf("FigEV4C-TNBC.pdf", height=14, width=8)
pheatmap(mat3, color=colorRampPalette(c("blue","white","red"))(100), border_color="NA",
    breaks=seq(-2,2,length.out=101), cluster_cols=TRUE, scale="none", fontsize_row=6, 
    show_colnames=FALSE, treeheight_row=70, treeheight_col=70, cutree_cols=1,
    clustering_method="ward.D2", annotation_col=annot2, annotation_colors=ann_colors2)
dev.off()




#####################################
### Tumour population
# Fig 
Tum <- Cluster %in% c(1,3)

cellNamesTum <- rownames(TNBC@meta.data)[Tum]
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
TNBCTum <- IntegrateData(anchorset=AnchorsTum, dims=1:dimUsed, k.weight=100)

# Dimensional reduction
dimUsed <- 30
DefaultAssay(TNBCTum) <- "integrated"
TNBCTum <- ScaleData(TNBCTum, verbose=FALSE)
TNBCTum <- RunPCA(TNBCTum, npcs=dimUsed, verbose=FALSE)
TNBCTum <- RunTSNE(TNBCTum, dims=1:dimUsed, seed.use=2018)
tSNE <- TNBCTum@reductions$tsne@cell.embeddings

GroupTum <- factor(TNBCTum@meta.data$group, levels=SamplesComb)
plotOrd <- sample(ncol(TNBCTum))

# Clustering
resolution <- 0.1
TNBCTum <- FindNeighbors(TNBCTum, dims=1:dimUsed, verbose=FALSE)
TNBCTum <- FindClusters(TNBCTum, resolution=resolution, verbose=FALSE)
ClusterTum <- as.integer(TNBCTum@meta.data$seurat_clusters)
ncls <- length(table(ClusterTum))

# Fig EV 3B
pdf("FigEV3B.pdf", height=9, width=9)
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

#save(z, d, y, SamplesCombTum, GroupTum, ncls, ClusterTum, col.p, tSNE, file="TNBCTum.RData")
#saveRDS(TNBCTum, file="SeuratObject_TNBCTum.rds")


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

# QL F-test
qftest <- glmQLFTest(qfit2, coef=2)
summary(decideTestsDGE(qftest))
topTags(qftest, n=30)

# KEGG
library(org.Hs.eg.db)
geneid <- mapIds(org.Hs.eg.db, keys=rownames(qfit2), column="ENTREZID", 
    keytype="SYMBOL")
keg <- kegga(qftest, geneid=geneid, species="Hs")
res <- topKEGG(keg, sort="up", n=20)
sel <- c(1:12)[-c(4,10)]

# Fig EV 3D
pdf("FigEV3D.pdf", height=7, width=7)
par(mar=c(5, 12, 2, 2))
barplot(-log10(res$P.Up[sel]), horiz=TRUE, col="grey70", names=res$Pathway[sel], xlab="-Log10(P-Value)", ylab="", las=2)
dev.off()


