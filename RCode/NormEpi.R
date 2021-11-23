### R v3.6.1
library(Seurat)
library(edgeR)
library(vcd)
library(pheatmap)

col.pMedium <- scater:::.get_palette("tableau10medium")
col.pDark <- scater:::.get_palette("tableau20")[2*(1:10)-1]
col.pLight <- scater:::.get_palette("tableau20")[2*(1:10)]
col.p <- c(col.pDark, col.pLight)


#dge <- read10X(path="../Data/N280/", DGEList=TRUE)

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
mito_upper[8] <- 0.3
nGenes_lower <- rep(500, length(SamplesComb))
nGenes_upper <- 1e3 * c(5.5, 4.5, 3, 4, 3.3, 4.5, 2.5, 5, 5.5, 5, 4.5)
lib_upper <- 1e4 * c(3.5, 2.5, 1, 2.5, 1.5, 2, 0.9, 3.2, 3.5, 2.8, 3)
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
NormEpi <- IntegrateData(anchorset=Anchors, dims=1:dimUsed, k.weight=100)

# Dim reduction
DefaultAssay(NormEpi) <- "integrated"
NormEpi <- ScaleData(NormEpi, verbose=FALSE)
NormEpi <- RunPCA(NormEpi, npcs=dimUsed, verbose=FALSE)
seed <- 200
NormEpi <- RunTSNE(NormEpi, dims=1:dimUsed, seed.use=2018)
tSNE <- NormEpi@reductions$tsne@cell.embeddings

# Fig EV1 C (left)
Group <- factor(NormEpi@meta.data$group, levels=SamplesComb)
col <- col.p[Group]
plotOrd <- sample(ncol(NormEpi))
pdf("tSNE-NormEpi.pdf", height=10, width=11)
par(mar=c(5.1, 4.1, 4.1, 9.5), xpd=TRUE)
plot(tSNE[plotOrd,], pch=16, col=col[plotOrd], cex=0.7, xlab="tSNE-1", ylab="tSNE-2",
    main="t-SNE - By sample")
legend("topright", inset=c(-0.2,0), legend=paste0(SamplesComb, " - ", table(Group)), 
    pch=16, col=col.p[1:length(SamplesComb)], title="Sample", title.adj=0.1, bty="n")
dev.off()

# Menopausal condition
col.p2 <- c("darkblue", "gold2")
Group2 <- c("Pre-meno", "Post-meno")[Group %in% c("N_0342_epi", "N_0372_epi", "N_0275_epi")+1L]
Group2 <- factor(Group2, levels=c("Pre-meno", "Post-meno"))

# Cell clustering
resolution <- 0.015
NormEpi <- FindNeighbors(NormEpi, dims=1:dimUsed, verbose=FALSE)
NormEpi <- FindClusters(NormEpi, resolution=resolution, verbose=FALSE)
Cluster <- as.integer(NormEpi@meta.data$seurat_clusters)
ncls <- length(table(Cluster))

# Fig EV1 C (right)
col <- col.p[Cluster]
pdf("tSNE-NormEpi-Cluster.pdf", height=10, width=11)
par(mar=c(5.1, 4.1, 4.1, 9.5), xpd=TRUE)
plot(tSNE, pch=16, col=col, cex=0.7, xlab="tSNE-1", ylab="tSNE-2",
    main="t-SNE - By cluster")
legend("topright", inset=c(-0.2,0), legend=paste0("Cluster", 1:ncls, " - ", table(Cluster)), 
    pch=16, col=col.p[1:ncls], title="Cluster", title.adj=0.1, bty="n")
dev.off()

#saveRDS(NormEpi, file="SeuratObject_NormEpi.rds")

allGenes <- rownames(get(DD[1]))
allGenesFilter <- c()
for(i in DGE) allGenesFilter <- c(allGenesFilter, rownames(get(i)))
allGenesFilter <- names(table(allGenesFilter))[table(allGenesFilter) >= 2]
u <- get(DD[1])
#d <- u[allGenes, ]
y <- u[allGenesFilter, ]
for(i in 2:length(SamplesComb)) {
    u <- get(DD[i])
#    d <- cbind(d, u[allGenes, ])
    y <- cbind(y, u[allGenesFilter, ])
}
#d$samples$group <- y$samples$group <- Group
y$samples$group <- Group

# log-CPM
prior.count <- 1
z <- edgeR::cpm(y, log=TRUE, prior.count=prior.count)

# Boxplots Fig EV1 D
load("../Data/Human-PosSigGenes.RData")
z_Basal <- z[rownames(z) %in% Basal, ]
z_LP <- z[rownames(z) %in% LP, ]
z_ML <- z[rownames(z) %in% ML, ]
z_Str <- z[rownames(z) %in% Str, ]
score_Basal <- colMeans(z_Basal)
score_LP <- colMeans(z_LP)
score_ML <- colMeans(z_ML)
score_Str <- colMeans(z_Str)
titles <- c("Basal", "LP", "ML", "Stroma")
dat <- data.frame(score_Basal, score_LP, score_ML, score_Str)
pdf("Boxplots-NormEpi.pdf", height=2.5, width=10)
par(mfrow=c(1,4))
for(j in 1:4)
boxplot(dat[,j] ~ Cluster, varwidth=FALSE, xlab="Cluster", ylab="AveLogCPM", range=0, medlwd = 1, 
    col=col.p[1:ncls], main=titles[j])
dev.off()

#save(z, y, SamplesComb, Group, Group2, ncls, Cluster, col.p, col.p2, tSNE, file="NormEpi.RData")


####################################################################################################
### Remove Stroma
Sub <- Cluster %in% c(1,2,3)
cellNamesSub <- rownames(NormEpi@meta.data)[Sub]
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
NormEpiSub <- IntegrateData(anchorset=AnchorsSub, dims=1:dimUsed, k.weight=100)

# Dimensional reduction
dimUsed <- 30
DefaultAssay(NormEpiSub) <- "integrated"
NormEpiSub <- ScaleData(NormEpiSub, verbose=FALSE)
NormEpiSub <- RunPCA(NormEpiSub, npcs=dimUsed, verbose=FALSE)
NormEpiSub <- RunTSNE(NormEpiSub, dims=1:dimUsed, seed.use=2018)
tSNE <- NormEpiSub@reductions$tsne@cell.embeddings

# Fig 1C
GroupSub <- factor(NormEpiSub@meta.data$group, levels=SamplesComb)
col <- col.p[GroupSub]
plotOrd <- sample(ncol(NormEpiSub))
pdf("Fig1C.pdf", height=9, width=9)
plot(tSNE[plotOrd,], pch=16, col=col[plotOrd], cex=0.7, xlab="tSNE-1", ylab="tSNE-2",
    main="t-SNE - By sample")
dev.off()

# Menopausal condition
col.p2 <- c("darkblue", "gold2")
GroupSub2 <- c("Pre-meno", "Post-meno")[GroupSub %in% c("N_0342_epi", "N_0372_epi", "N_0275_epi")+ 1L]
GroupSub2 <- factor(GroupSub2, levels=c("Pre-meno", "Post-meno"))

# Fig 1D
pdf("Fig1D.pdf", height=9, width=18)
par(mfrow=c(1,2))
plot(tSNE[GroupSub2==levels(GroupSub2)[1],], pch=16, col=col.p2[1], cex=0.7, xlab="tSNE-1", ylab="tSNE-2",
    main="Pre-menopausal")
plot(tSNE[GroupSub2==levels(GroupSub2)[2],], pch=16, col=col.p2[2], cex=0.7, xlab="tSNE-1", ylab="tSNE-2",
    main="Post-menopausal")
dev.off()

# Clustering
resolution <- 0.015
NormEpiSub <- FindNeighbors(NormEpiSub, dims=1:dimUsed, verbose=FALSE)
NormEpiSub <- FindClusters(NormEpiSub, resolution=resolution, verbose=FALSE)
ClusterSub <- as.integer(NormEpiSub@meta.data$seurat_clusters)
ncls <- length(table(ClusterSub))

# Fig 1E
pdf("Fig1E.pdf", height=9, width=9)
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

# Ternary plot (Fig 1F)
load("../Data/Human-PosSigGenes.RData")
Basal <- intersect(Basal, y$genes$Official)
LP <- intersect(LP, y$genes$Official)
ML <- intersect(ML, y$genes$Official)
IN2 <- matrix(0L, ncol(y), 3L)
colnames(IN2) <- c("Basal", "LP", "ML")
IN2[, "Basal"] <- colSums(y$counts[Basal,] > 0L)
IN2[, "LP"] <- colSums(y$counts[LP,] > 0L)
IN2[, "ML"] <- colSums(y$counts[ML,] > 0L)
pdf("Fig1F.pdf", height=8, width=8)
ternaryplot(IN2[plotOrd,c(2,3,1)], cex=0.3, pch=16, col=col.p[ClusterSub][plotOrd], 
    main="10X human", grid=TRUE)
dev.off()

# log-CPM
prior.count <- 1
z <- edgeR::cpm(d, log=TRUE, prior.count=prior.count)

#save(z, d, y, SamplesCombSub, GroupSub, GroupSub2, ncls, ClusterSub, col.p, col.p2, tSNE, file="NormEpiSub.RData")
#saveRDS(NormEpiSub, file="SeuratObject_NormEpiSub.rds")


############
# Pseudo-bulk DE analysis
sampClust <- paste(y$samples$group, ClusterSub, sep="_Clst")
counts2 <- t(rowsum(t(y$counts), group=sampClust))
yComb <- DGEList(counts2)
yComb$samples$Patient <- gsub("_Clst.*$","",colnames(yComb))
yComb$samples$Cluster <- as.numeric(gsub("^.*_Clst","",colnames(yComb)))
yComb$samples$group <- yComb$samples$Cluster

# Focus on the three major clusters
N <- 1:3
yClstSub <- yComb[, yComb$samples$Cluster %in% N]
keep <- filterByExpr(yClstSub, min.count=15, min.total.count=30)
yClstSub <- yClstSub[keep,,keep=FALSE]
yClstSub <- calcNormFactors(yClstSub)
m <- match(paste0(rep(SamplesComb,each=3), "_Clst", 1:3), colnames(yClstSub))
yClstSub <- yClstSub[,m]
Meno <- c("Pre", "Post")[SamplesComb %in% c("N_0342_epi", "N_0372_epi", "N_0275_epi") + 1L]
Meno <- rep(Meno, each=3)
Meno <- factor(Meno, levels=c("Pre", "Post"))
yClstSub$samples <- cbind(yClstSub$samples, Meno=Meno)

# design matrix
Cls <- as.factor(yClstSub$samples$Cluster)
Pat <- factor(yClstSub$samples$Patient, levels=SamplesComb)
design <- model.matrix(~ Cls + Pat)

# Estimate NB dispersion and QL dispersion.
yClstSub <- estimateDisp(yClstSub, design=design)
qfit <- glmQLFit(yClstSub, design)

# Contrast
contr <- makeContrasts(
    Cls1 = -(Cls2+Cls3)/2,
    Cls2 = Cls2-Cls3/2,
    Cls3 = Cls3-Cls2/2, levels=design) 
ctest1 <- glmQLFTest(qfit, contrast=contr[,1])
ctest2 <- glmQLFTest(qfit, contrast=contr[,2])
ctest3 <- glmQLFTest(qfit, contrast=contr[,3])

# Top 50 markers for heat map
top <- 50
ord <- order(ctest1$table$PValue, decreasing=FALSE)
upreg <- ctest1$table$logFC > 0
Markers1 <- rownames(yClstSub)[ord[upreg][1:top]]
ord <- order(ctest2$table$PValue, decreasing=FALSE)
upreg <- ctest2$table$logFC > 0
Markers2 <- rownames(yClstSub)[ord[upreg][1:top]]
ord <- order(ctest3$table$PValue, decreasing=FALSE)
upreg <- ctest3$table$logFC > 0
Markers3 <- rownames(yClstSub)[ord[upreg][1:top]]
Markers <- c(Markers1,Markers2,Markers3)
Markers <- Markers[!duplicated(Markers)]

# log-CPM
prior.count <- 1
zClstSub <- edgeR::cpm(yClstSub, log=TRUE, prior.count=prior.count)

# Fig EV1 E
pdf("FigEV1E.pdf",height=12, width=9)
annot2 <- data.frame(Cluster=paste0("Cluster ", Cls), Patient=Pat, Menopause=Meno)
rownames(annot2) <- colnames(zClstSub)
ann_colors2 <- list(Cluster=col.p[N], Patient=col.p[1:length(SamplesComb)], Menopause=col.p2[1:2])
names(ann_colors2$Cluster) <- paste0("Cluster ", N)
names(ann_colors2$Patient) <- SamplesComb
names(ann_colors2$Menopause) <- c("Pre","Post")
mat2 <- t(scale(t(zClstSub[Markers, ])))
pheatmap(mat2, color=colorRampPalette(c("blue","white","red"))(100), border_color="NA",
    breaks=seq(-2,2,length.out=101), cluster_cols=TRUE, scale="none", fontsize_row=5, 
    show_colnames=FALSE, treeheight_row=70, treeheight_col=70, cutree_cols=3,
    clustering_method="ward.D2", annotation_col=annot2, annotation_colors=ann_colors2)
dev.off()



###################################
### Quasi multinomial test
cellNumMeno <- table(ClusterSub, GroupSub2)
cellPropMeno <- t( t(cellNumMeno) / rowSums(t(cellNumMeno)) )

# Quasi-Poisson GLM
cellNum <- table(ClusterSub, GroupSub)
grp <- rep(1:2, c(8,3))
y <- DGEList(counts=cellNum, group=grp)
pat <- gl(11,4)
PostMeno <- factor(rep(y$sample$group, each=4))
clust <- gl(4,1,11*4)

v <- as.vector(cellNum)
fit <- glm(v ~ clust + pat + PostMeno:clust, family=quasipoisson())
summary(fit)
anova(fit, test="F")

colnames(cellNum) <- Samples
rownames(cellNum) <- paste("Cluster", rownames(cellNum))
write.csv(cellNum, file="NormEpiSub-CellCounts.csv")




################
library(destiny)
setwd("DestinySV33Sub/")
knit("DestinySV33Sub.Rnw")
setwd("../")

##### Primed cells on diffusion plot
setwd("~/Project/Bhupinder_SingleCell/Human/_10X-Normal-Epi-Combine/SeuratV3")
load("NormEpiSub_tSNE.RData")
source("../../pipFromLocator.R")





