#cp -r ~/Project/Bhupinder/10X_Human/_10X-Normal-Total-Combine/Data/ /wehisan/general/academic/lab_smyth/cheny/Project/BPal/10X_Human/ScientificData/NormTotal/

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


#####################################
# Fig 2 A-I, Fig EV1 F-H

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

### Gene annotation
for(i in 1:length(SamplesComb)){
    ann <- alias2SymbolUsingNCBI(get(DGE[i])$genes$Symbol, required.columns=c("GeneID","Symbol"), 
        gene.info.file="~/Annotation/180808_Homo_sapiens.gene_info.gz")
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
mito_upper <- rep(0.2, length(SamplesComb))
mito_upper[6] <- 0.3
nGenes_lower <- rep(500, length(SamplesComb))
nGenes_upper <- 1e3 * c(7.5, 9, 6.5, 8, 7, 7, 6.5, 5, 7, 8, 2.7, 6.5, 6)
lib_upper <- 1e4 * c(8, 12, 7, 12, 6, 6, 6, 3, 6, 8, 1, 6, 5)
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
Anchors <- FindIntegrationAnchors(object.list=CombSeurat, dims=1:dimUsed)

# Integration
NormTotal <- IntegrateData(anchorset=Anchors, dims=1:dimUsed, k.weight=100)

# Dim reduction
DefaultAssay(NormTotal) <- "integrated"
NormTotal <- ScaleData(NormTotal, verbose=FALSE)
NormTotal <- RunPCA(NormTotal, npcs=dimUsed, verbose=FALSE)
NormTotal <- RunTSNE(NormTotal, dims=1:dimUsed, seed.use=2018)
tSNE <- NormTotal@reductions$tsne@cell.embeddings


# Fig 2A
Group <- factor(NormTotal@meta.data$group, levels=SamplesComb)
col <- col.p[Group]
plotOrd <- sample(ncol(NormTotal))
pdf("Fig2A.pdf", height=9, width=9)
plot(tSNE[plotOrd,], pch=16, col=col[plotOrd], cex=0.7, xlab="tSNE-1", ylab="tSNE-2",
    main="t-SNE - By sample")
dev.off()

# Menopausal condition
col.p2 <- c("darkblue", "gold2")
Group2 <- c("Pre-meno", "Post-meno")[Group %in% c("N_0342_total", "N_0372_total","N_0021_total","N_0275_total","N_0288_total") + 1L]
Group2 <- factor(Group2, levels=c("Pre-meno", "Post-meno"))

# Cell clustering
resolution <- 0.05
NormTotal <- FindNeighbors(NormTotal, dims=1:dimUsed, verbose=FALSE)
NormTotal <- FindClusters(NormTotal, resolution=resolution, verbose=FALSE)
Cluster <- as.integer(NormTotal@meta.data$seurat_clusters)
ncls <- length(table(Cluster))

# Fig 2B
col <- col.p[Cluster]
pdf("Fig2B.pdf", height=9, width=9)
par(mar=c(5.1, 4.1, 4.1, 9.5), xpd=TRUE)
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

# Fig 2C
load("../Data/Human-PosSigGenes.RData")
Basal <- intersect(Basal, y$genes$Official)
LP <- intersect(LP, y$genes$Official)
ML <- intersect(ML, y$genes$Official)
IN2 <- matrix(0L, ncol(y), 3L)
colnames(IN2) <- c("Basal", "LP", "ML")
IN2[, "Basal"] <- colSums(y$counts[Basal,] > 0L)
IN2[, "LP"] <- colSums(y$counts[LP,] > 0L)
IN2[, "ML"] <- colSums(y$counts[ML,] > 0L)
pdf("Fig2C.pdf", height=8, width=8)
ternaryplot(IN2[plotOrd,c(2,3,1)], cex=0.3, pch=16, col=col.p[Cluster][plotOrd], 
    main="10X human", grid=TRUE)
dev.off()

# log-CPM
prior.count <- 1
z <- edgeR::cpm(d, log=TRUE, prior.count=prior.count)

#saveRDS(NormTotal, file="SeuratObject_NormTotal.rds")
#save(z, d, y, SamplesComb, Group, Group2, ncls, Cluster, col.p, col.p2, tSNE, file="NormTotal.RData")



####################################################################################################
### Micro-environment
Sub <- !(Cluster %in% c(1,3,4))

cellNamesSub <- rownames(NormTotal@meta.data)[Sub]
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
AnchorsSub <- FindIntegrationAnchors(object.list=CombSeuratSub, dims=1:dimUsed)
NormTotalSub <- IntegrateData(anchorset=AnchorsSub, dims=1:dimUsed, k.weight=100)

# Dimensional reduction
dimUsed <- 30
DefaultAssay(NormTotalSub) <- "integrated"
NormTotalSub <- ScaleData(NormTotalSub, verbose=FALSE)
NormTotalSub <- RunPCA(NormTotalSub, npcs=dimUsed, verbose=FALSE)
NormTotalSub <- RunTSNE(NormTotalSub, dims=1:dimUsed, seed.use=2018)
tSNE <- NormTotalSub@reductions$tsne@cell.embeddings

GroupSub <- factor(NormTotalSub@meta.data$group, levels=SamplesComb)
plotOrd <- sample(ncol(NormTotalSub))

col.p2 <- c("darkblue", "gold2")
GroupSub2 <- c("Pre-meno", "Post-meno")[GroupSub %in% c("N_0342_total", "N_0372_total","N_0021_total","N_0275_total","N_0288_total") + 1L]
GroupSub2 <- factor(GroupSub2, levels=c("Pre-meno", "Post-meno"))

# Fig 2G
pdf("Fig2G.pdf", height=9, width=9)
plot(tSNE[plotOrd,], pch=16, col=col.p2[GroupSub2][plotOrd], cex=0.7, xlab="tSNE-1", ylab="tSNE-2",
    main="t-SNE - By menopausal type")
dev.off()

# clustering
resolution <- 0.05
NormTotalSub <- FindNeighbors(NormTotalSub, dims=1:dimUsed, verbose=FALSE)
NormTotalSub <- FindClusters(NormTotalSub, resolution=resolution, verbose=FALSE)
ClusterSub <- as.integer(NormTotalSub@meta.data$seurat_clusters)
ncls <- length(table(ClusterSub))

# Fig 2D
pdf("Fig2D.pdf", height=9, width=9)
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

#save(z, d, y, SamplesCombSub, GroupSub, GroupSub2, ncls, ClusterSub, col.p, col.p2, tSNE, file="NormTotalSub.RData")
#saveRDS(NormTotalSub, file="SeuratObject_NormTotalSub.rds")


############
# Pseudo-bulk DE analysis
sampClust <- paste(y$samples$group, ClusterSub, sep="_Clst")
counts2 <- t(rowsum(t(y$counts), group=sampClust))
yComb <- DGEList(counts2)
yComb$samples$Patient <- gsub("_Clst.*$","",colnames(yComb))
yComb$samples$Cluster <- as.numeric(gsub("^.*_Clst","",colnames(yComb)))
yComb$samples$group <- yComb$samples$Cluster

# Focus on the 7 major clusters
N <- 1:7
yClstSub <- yComb[, yComb$samples$Cluster %in% N]
keep <- filterByExpr(yClstSub, min.count=10, min.total.count=20)
yClstSub <- yClstSub[keep,,keep=FALSE]
yClstSub <- calcNormFactors(yClstSub)
sn <- paste0(rep(SamplesComb, each=length(N)), "_Clst", N)
sel <- sn %in% colnames(yClstSub)
sn <- sn[sel]
m <- match(sn, colnames(yClstSub))
yClstSub <- yClstSub[,m]
Meno <- c("Pre", "Post")[SamplesComb %in% c("N_0342_total", "N_0372_total","N_0021_total","N_0275_total","N_0288_total") + 1L]
Meno <- rep(Meno, each=length(N))[sel]
Meno <- factor(Meno, levels=c("Pre", "Post"))
yClstSub$samples <- cbind(yClstSub$samples, Meno=Meno)

# filter out some pseudo samples with low library size.
sel <- yClstSub$samples$lib.size > 3.5e4
yClstSub2 <- yClstSub[,sel]
yClstSub2$samples

# Fig 2E
mds <- plotMDS(yClstSub2, plot=FALSE)
pdf("Fig2E.pdf", height=6, width=6)
plotMDS(mds, pch=16, col=col.p[yClstSub2$samples$Cluster], main="")
legend("bottomright", legend=paste0("Cluster", N), pch=16, col=col.p[N])
dev.off()

# design matrix
Cls <- as.factor(yClstSub2$samples$Cluster)
Pat <- factor(yClstSub2$samples$Patient, levels=SamplesComb)
Meno <- yClstSub2$samples$Meno
design <- model.matrix(~ Cls + Pat)

# Estimate NB dispersion and QL dispersion.
yClstSub2 <- estimateDisp(yClstSub2, design=design)
qfit2 <- glmQLFit(yClstSub2, design)

# log-CPM
prior.count <- 1
zClstSub2 <- edgeR::cpm(yClstSub2, log=TRUE, prior.count=prior.count)

annot2 <- data.frame(Cluster=paste0("Cluster ", Cls), Patient=Pat, Menopause=Meno)
rownames(annot2) <- colnames(zClstSub2)
ann_colors2 <- list(Cluster=col.p[N], Patient=col.p[1:length(SamplesComb)], Menopause=col.p2[1:2])
names(ann_colors2$Cluster) <- paste0("Cluster ", N)
names(ann_colors2$Patient) <- SamplesComb
names(ann_colors2$Menopause) <- c("Pre","Post")

contr <- makeContrasts(
    Cls1 = -(Cls2+Cls3+Cls4+Cls5+Cls6+Cls7)/6,
    Cls2 = Cls2-(Cls3+Cls4+Cls5+Cls6+Cls7)/6,
    Cls3 = Cls3-(Cls2+Cls4+Cls5+Cls6+Cls7)/6,
    Cls4 = Cls4-(Cls2+Cls3+Cls5+Cls6+Cls7)/6,
    Cls5 = Cls5-(Cls2+Cls3+Cls4+Cls6+Cls7)/6,
    Cls6 = Cls6-(Cls2+Cls3+Cls4+Cls5+Cls7)/6,
    Cls7 = Cls7-(Cls2+Cls3+Cls4+Cls5+Cls6)/6, levels=design)
ctest <- list()
for(i in 1:ncls) ctest[[i]] <- glmQLFTest(qfit2, contrast=contr[,i])
summary(decideTestsDGE(ctest[[1]]))
summary(decideTestsDGE(ctest[[2]]))
summary(decideTestsDGE(ctest[[3]]))
summary(decideTestsDGE(ctest[[4]]))
summary(decideTestsDGE(ctest[[5]]))
summary(decideTestsDGE(ctest[[6]]))
summary(decideTestsDGE(ctest[[7]]))

top <- 20
pseudoMakers <- list()
for(i in 1:ncls) {
    ord <- order(ctest[[i]]$table$PValue, decreasing=FALSE)
    upreg <- ctest[[i]]$table$logFC > 0
    pseudoMakers[[i]] <- rownames(yClstSub2)[ord[upreg][1:top]]
}
Markers <- unlist(pseudoMakers)
Markers <- Markers[!duplicated(Markers)]

# Fig 2F
pdf("Fig2F.pdf",height=11, width=9)
mat4 <- t(scale(t(zClstSub2[Markers, ])))
pheatmap(mat4, color=colorRampPalette(c("blue","white","red"))(100), border_color="NA",
    breaks=seq(-2,2,length.out=101), cluster_cols=TRUE, scale="none", fontsize_row=5, 
    show_colnames=FALSE, treeheight_row=70, treeheight_col=70, cutree_cols=7,
    clustering_method="ward.D2", annotation_col=annot2, annotation_colors=ann_colors2)
dev.off()

# Fig EV1 H
pdf("FigEV1H.pdf",height=10, width=9)
SigImm <- read.delim("../Data/ImmuneMarkers2.txt", header=TRUE, stringsAsFactors=FALSE)
Sig <- SigImm[SigImm$Signatures %in% rownames(yClstSub2),]
Sig <- Sig[!duplicated(Sig$Signatures),]
mat2 <- t(scale(t(zClstSub2[Sig$Signatures, ])))
pheatmap(mat2, color=colorRampPalette(c("blue","white","red"))(100), border_color="NA",
    breaks=seq(-2,2,length.out=101), cluster_cols=TRUE, scale="none", fontsize_row=10, 
    show_colnames=FALSE, treeheight_row=70, treeheight_col=70, cutree_cols=7,
    clustering_method="ward.D2", annotation_col=annot2, annotation_colors=ann_colors2)
dev.off()


###################################
### Quasi multinomial test
cellNumMeno <- table(ClusterSub, GroupSub2)
cellPropMeno <- t( t(cellNumMeno) / rowSums(t(cellNumMeno)) )

pdf("Fig2G-right.pdf",height=5, width=6)
barplot(t(cellPropMeno*100), beside=TRUE, main="", xlab="Cluster (Pre | Post)", ylab="Cell number Percentage", 
    col=rep(col.p2, ncls), names=1:ncls, space=c(0,0.3))
dev.off()

# Quasi-Poisson GLM
cellNum <- table(ClusterSub, GroupSub)
grp <- rep(1:2, c(8,5))
y <- DGEList(counts=cellNum, group=grp)
pat <- gl(13,7)
PostMeno <- factor(rep(y$sample$group, each=7))
clust <- gl(7,1,13*7)

v <- as.vector(cellNum)
fit <- glm(v ~ clust + pat + PostMeno:clust, family=quasipoisson())
summary(fit)
anova(fit, test="F")

colnames(cellNum) <- Samples
rownames(cellNum) <- paste("Cluster", rownames(cellNum))
write.csv(cellNum, file="NormTotalSub-CellCounts.csv")


# Quasi-Poisson GLM: cluster 1
cellNum1 <- rowsum(cellNum, c(1,2,2,2,2,2,2))
grp1 <- rep(1:2, c(8,5))
y1 <- DGEList(counts=cellNum1, group=grp1)
pat1 <- gl(13,2)
PostMeno1 <- factor(rep(y1$sample$group, each=2))
clust1 <- gl(2,1,13*2)

v1 <- as.vector(cellNum1)
fit1 <- glm(v1 ~ clust1 + pat1 + PostMeno1:clust1, family=quasipoisson())
summary(fit1)
anova(fit1, test="F")

# Quasi-Poisson GLM: cluster 2
cellNum2 <- rowsum(cellNum, c(2,1,2,2,2,2,2))
grp2 <- rep(1:2, c(8,5))
y2 <- DGEList(counts=cellNum2, group=grp2)
pat2 <- gl(13,2)
PostMeno2 <- factor(rep(y2$sample$group, each=2))
clust2 <- gl(2,1,13*2)

v2 <- as.vector(cellNum2)
fit2 <- glm(v2 ~ clust2 + pat2 + PostMeno2:clust2, family=quasipoisson())
summary(fit2)
anova(fit2, test="F")



#####################################
### Fibroblast
# Fig 3 C-F
Fib <- ClusterSub == 1

cellNamesFib <- rownames(NormTotalSub@meta.data)[Fib]
DGEFib <- paste0("dge_fib_", SamplesComb)
DDFib <- paste0("dd_fib_", SamplesComb)
for(i in 1:length(SamplesComb)) {
    d <- get(DD[i])
    d <- d[, colnames(d) %in% cellNamesFib]
    keep1 <- rowSums(d$counts > 0) >= ncol(d)*0.01
    eval( parse(text=paste0(DDFib[i],"<- d")) )
    d <- d[keep1, , keep=FALSE]
    eval( parse(text=paste0(DGEFib[i],"<- d")) )
}

# HVGs
SamplesCombFib <- paste0(SamplesComb, "_Fib")
CombSeuratFib <- list()
for(i in 1:length(SamplesCombFib)) {
    d <- get(DGEFib[i])
    so <- CreateSeuratObject(d$counts, project=SamplesCombFib[i])
    so <- NormalizeData(so)
    so <- FindVariableFeatures(so, selection.method="vst", nfeatures=1500)
    so <- ScaleData(so)
    so@meta.data$group <- SamplesComb[i]
    CombSeuratFib[[i]] <- so
    eval( parse(text=paste0(SamplesCombFib[i],"<- so")) )
}
names(CombSeuratFib) <- SamplesCombFib

# Integration
dimUsed <- 30
AnchorsFib <- FindIntegrationAnchors(object.list=CombSeuratFib, dims=1:dimUsed,
    anchor.features=1000, scale=TRUE, k.anchor=5, k.filter=30,
    k.score=20, max.features=100)
NormTotalFib <- IntegrateData(anchorset=AnchorsFib, dims=1:dimUsed, k.weight=100)

# Dimensional reduction
dimUsed <- 30
DefaultAssay(NormTotalFib) <- "integrated"
NormTotalFib <- ScaleData(NormTotalFib, verbose=FALSE)
NormTotalFib <- RunPCA(NormTotalFib, npcs=dimUsed, verbose=FALSE)
NormTotalFib <- RunTSNE(NormTotalFib, dims=1:dimUsed, seed.use=2018)
tSNE <- NormTotalFib@reductions$tsne@cell.embeddings

GroupFib <- factor(NormTotalFib@meta.data$group, levels=SamplesComb)
plotOrd <- sample(ncol(NormTotalFib))

col.p2 <- c("darkblue", "gold2")
GroupFib2 <- c("Pre-meno", "Post-meno")[GroupFib %in% c("N_0342_total", "N_0372_total","N_0021_total","N_0275_total","N_0288_total") + 1L]
GroupFib2 <- factor(GroupFib2, levels=c("Pre-meno", "Post-meno"))

# Clustering
resolution <- 0.1
NormTotalFib <- FindNeighbors(NormTotalFib, dims=1:dimUsed, verbose=FALSE)
NormTotalFib <- FindClusters(NormTotalFib, resolution=resolution, verbose=FALSE)
ClusterFib <- as.integer(NormTotalFib@meta.data$seurat_clusters)
ncls <- length(table(ClusterFib))

# Fig 3C
pdf("Fig3C.pdf", height=9, width=9)
col <- col.p[GroupFib]
plot(tSNE[plotOrd,], pch=16, col=col[plotOrd], cex=0.7, xlab="tSNE-1", ylab="tSNE-2",
    main="t-SNE - By sample")
dev.off()

# Fig 3D
pdf("Fig3D.pdf", height=9, width=9)
col <- col.p[ClusterFib]
plot(tSNE, pch=16, col=col, cex=0.7, xlab="tSNE-1", ylab="tSNE-2", main="t-SNE - By cluster")
dev.off()


# Combine raw data
allGenes <- rownames(get(DDFib[1]))
allGenesFilter <- c()
for(i in DGEFib) allGenesFilter <- c(allGenesFilter, rownames(get(i)))
allGenesFilter <- names(table(allGenesFilter))[table(allGenesFilter) >= 2]
u <- get(DDFib[1])
d <- u[allGenes, ]
y <- u[allGenesFilter, ]
for(i in 2:length(SamplesComb)) {
    u <- get(DDFib[i])
    d <- cbind(d, u[allGenes, ])
    y <- cbind(y, u[allGenesFilter, ])
}
d$samples$group <- y$samples$group <- GroupFib

# log-CPM
prior.count <- 1
z <- edgeR::cpm(d, log=TRUE, prior.count=prior.count)

#save(z, d, y, SamplesCombFib, GroupFib, GroupFib2, ncls, ClusterFib, col.p, col.p2, tSNE, file="NormTotalFib.RData")
#saveRDS(NormTotalFib, file="SeuratObject_NormTotalFib.rds")


# Marker genes
markers <- list()
for(i in 1:ncls){
    markers[[i]] <- FindMarkers(NormTotalFib, ident.1=i-1, min.pct=0.1, only.pos=TRUE)
}
names(markers) <- paste0("Cluster", 1:ncls)

top <- 20
nTop <- pmin(top, sapply(markers, nrow))
topMarkers <- c()
for(i in 1:ncls) topMarkers <- c(topMarkers, rownames(markers[[i]][1:nTop[i],]))
topMarkers <- topMarkers[!duplicated(topMarkers)]
topMarkers

# Fig 3F
pdf("Fig3F.pdf", height=9, width=7)
mat <- t(scale(t(z[topMarkers, ])))
ord <- order(ClusterFib, GroupFib2)
annot <- data.frame(Condition=GroupFib2[ord], Cluster=paste0("Cluster ", ClusterFib[ord]))
rownames(annot) <- colnames(mat)[ord]
ann_colors <- list(Cluster=col.p[1:ncls], Condition=col.p2)
names(ann_colors$Cluster) <- paste0("Cluster ", 1:ncls)
names(ann_colors$Condition) <- levels(GroupFib2)
library(pheatmap)
pheatmap(mat[,ord], color=colorRampPalette(c("blue","white","red"))(100), 
    breaks=seq(-2,2,length.out=101), cluster_cols=FALSE, scale="none", 
    labels_col=rep("",ncol(mat)), fontsize_row=7, 
    clustering_method="ward.D2", annotation_col=annot, annotation_colors=ann_colors)
dev.off()

# KEGG analysis
mk1vs2 <- FindMarkers(NormTotalFib, ident.1=0, ident.2=1, min.pct=0.1, only.pos=TRUE)
mk1vs2 <- mk1vs2[mk1vs2$p_val_adj < 0.05, ]

anno <- read.delim(file="../Data/200702_Homo_sapiens.gene_info.gz", 
    header=TRUE)[,-1]
Universe <- rownames(NormTotalFib@assays$RNA)
m <- match(Universe, anno$Symbol)
UniverseID <- anno$GeneID[m]
UniverseID <- UniverseID[!is.na(UniverseID)]
m <- match(rownames(mk1vs2), anno$Symbol)
mk1v2ID <- anno$GeneID[m]
mk1v2ID <- mk1v2ID[!is.na(mk1v2ID)]
keg1v2 <- kegga(mk1v2ID, universe=UniverseID, species="Hs")
topKEGG(keg1v2, truncate=43)
topKEGG(keg1v2)

res <- topKEGG(keg1v2)
sel <- c(2:5,8,6,11,13,15)

# Fig 3E
pdf("Fig3E.pdf", height=9, width=6)
par(mar=c(24, 5, 2, 2))
barplot(-log10(res$P.DE[sel]), col="darkblue", names=res$Pathway[sel], ylab="-Log10(P-Value)", xlab="", las=2)
dev.off()



