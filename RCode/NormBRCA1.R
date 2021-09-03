
smythlab
cd cheny/Project/BPal/10X_Human/ScientificData/NormBRCA1
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
# Fig 4 A-I, Fig App S1 A-C

### Samples to be combined
#SamplesComb <- c("PM0019","PM0233","PM0092","PM0230","PM0095Total",
#    "MH0023Total","MH0064Total","MH0169","KCF0894","MH0033Total",
#    "MH0023","MH0090")
Samples <- c("N-0019-total","N-0233-total","N-0092-total","N-0230.17-total","N-0093-total",
    "N-0123-total","N-0064-total","N-0169-total","B1-0894","B1-0033","B1-0023","B1-0090")
SamplesComb <- gsub("-","_",Samples)
SamplesNormal <- SamplesComb[1:8]
SamplesPreB1 <- SamplesComb[9:12]

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
NormB1Total <- IntegrateData(anchorset=Anchors, dims=1:dimUsed, k.weight=100)

# Dim reduction
DefaultAssay(NormB1Total) <- "integrated"
NormB1Total <- ScaleData(NormB1Total, verbose=FALSE)
NormB1Total <- RunPCA(NormB1Total, npcs=dimUsed, verbose=FALSE)
NormB1Total <- RunTSNE(NormB1Total, dims=1:dimUsed, seed.use=2018)
tSNE <- NormB1Total@reductions$tsne@cell.embeddings

# Fig 4A
Group <- factor(NormB1Total@meta.data$group, levels=SamplesComb)
col <- col.p[Group]
plotOrd <- sample(ncol(NormB1Total))
pdf("Fig4A.pdf", height=9, width=9)
plot(tSNE[plotOrd,], pch=16, col=col[plotOrd], cex=0.7, xlab="tSNE-1", ylab="tSNE-2",
    main="t-SNE - By sample")
dev.off()

# Condition
col.p2 <- c("darkblue", "gold2")
Group2 <- c("Normal", "BRCA1")[Group %in% SamplesPreB1 + 1L]
Group2 <- factor(Group2, levels=c("Normal", "BRCA1"))

# Fig App S1A-left
pdf("FigS1A-left.pdf", height=7, width=14)
par(mfrow=c(1,2))
plot(tSNE[Group2==levels(Group2)[1],], pch=16, col=col.p2[1], cex=0.7, xlab="tSNE-1", ylab="tSNE-2",
    main="Normal")
plot(tSNE[Group2==levels(Group2)[2],], pch=16, col=col.p2[2], cex=0.7, xlab="tSNE-1", ylab="tSNE-2",
    main="BRCA1")
dev.off()


# Cell clustering
resolution <- 0.12
NormB1Total <- FindNeighbors(NormB1Total, dims=1:dimUsed, verbose=FALSE)
NormB1Total <- FindClusters(NormB1Total, resolution=resolution, verbose=FALSE)
Cluster <- as.integer(NormB1Total@meta.data$seurat_clusters)
ncls <- length(table(Cluster))

# Fig 4B
col <- col.p[Cluster]
pdf("Fig4B.pdf", height=9, width=9)
plot(tSNE, pch=16, col=col, cex=0.7, xlab="tSNE-1", ylab="tSNE-2", main="t-SNE - By cluster")
dev.off()

# Fig App S1A-right
cellNum <- table(Cluster, Group2)
cellProp <- t( t(cellNum) / rowSums(t(cellNum)) )
pdf("FigS1A-right.pdf", height=6, width=9)
barplot(t(cellProp*100), beside=TRUE, main="", xlab="Cluster (Pre | Post)", ylab="Cell number Percentage", 
    col=rep(col.p2, ncls), names=1:ncls, space=c(0,0.3))
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

#saveRDS(NormB1Total, file="SeuratObject_NormB1Total.rds")
#save(z, d, y, SamplesComb, Group, Group2, ncls, Cluster, col.p, col.p2, tSNE, file="NormB1Total.RData")


####################################################################################################
### Micro-environment
Sub <- !(Cluster %in% c(2,4,5,8))

cellNamesSub <- rownames(NormB1Total@meta.data)[Sub]
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
NormB1TotalSub <- IntegrateData(anchorset=AnchorsSub, dims=1:dimUsed, k.weight=100)

# Dimensional reduction
dimUsed <- 30
DefaultAssay(NormB1TotalSub) <- "integrated"
NormB1TotalSub <- ScaleData(NormB1TotalSub, verbose=FALSE)
NormB1TotalSub <- RunPCA(NormB1TotalSub, npcs=dimUsed, verbose=FALSE)
NormB1TotalSub <- RunTSNE(NormB1TotalSub, dims=1:dimUsed, seed.use=2018)
tSNE <- NormB1TotalSub@reductions$tsne@cell.embeddings

GroupSub <- factor(NormB1TotalSub@meta.data$group, levels=SamplesComb)
plotOrd <- sample(ncol(NormB1TotalSub))

col.p2 <- c("darkblue", "gold2")
GroupSub2 <- c("Normal", "BRCA1")[GroupSub %in% SamplesPreB1 + 1L]
GroupSub2 <- factor(GroupSub2, levels=c("Normal", "BRCA1"))

# clustering
resolution <- 0.08
NormB1TotalSub <- FindNeighbors(NormB1TotalSub, dims=1:dimUsed, verbose=FALSE)
NormB1TotalSub <- FindClusters(NormB1TotalSub, resolution=resolution, verbose=FALSE)
ClusterSub <- as.integer(NormB1TotalSub@meta.data$seurat_clusters)
ncls <- length(table(ClusterSub))

# Fig 4C
pdf("Fig4C.pdf", height=9, width=9)
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

#save(z, d, y, SamplesCombSub, GroupSub, GroupSub2, ncls, ClusterSub, col.p, col.p2, tSNE, file="NormB1TotalSub.RData")
#saveRDS(NormB1TotalSub, file="SeuratObject_NormB1TotalSub.rds")


############
# Pseudo-bulk DE analysis
sampClust <- paste(y$samples$group, ClusterSub, sep="_Clst")
counts2 <- t(rowsum(t(y$counts), group=sampClust))
yComb <- DGEList(counts2)
yComb$samples$Patient <- gsub("_Clst.*$","",colnames(yComb))
yComb$samples$Cluster <- as.numeric(gsub("^.*_Clst","",colnames(yComb)))
yComb$samples$group <- yComb$samples$Cluster

# Focus on all clusters
N <- 1:9
yClstSub <- yComb[, yComb$samples$Cluster %in% N]
keep <- filterByExpr(yClstSub, min.count=7, min.total.count=15)
yClstSub <- yClstSub[keep,,keep=FALSE]
yClstSub <- calcNormFactors(yClstSub)
sn <- paste0(rep(SamplesComb, each=length(N)), "_Clst", N)
sel <- sn %in% colnames(yClstSub)
sn <- sn[sel]
m <- match(sn, colnames(yClstSub))
yClstSub <- yClstSub[,m]
Population <- yClstSub$samples$Patient
Population[Population %in% SamplesNormal] <- "Normal"
Population[Population %in% SamplesPreB1] <- "BRCA1-Pre"
Population <- factor(Population, levels=c("Normal", "BRCA1-Pre"))
yClstSub$samples <- cbind(yClstSub$samples, Population=Population)
yClstSub$samples

# filter out some pseudo samples with low library size.
sel <- yClstSub$samples$lib.size > 4e4
yClstSub2 <- yClstSub[,sel]
yClstSub2$samples

# design matrix
Cls <- as.factor(yClstSub2$samples$Cluster)
Pat <- factor(yClstSub2$samples$Patient, levels=SamplesComb)
Pop <- yClstSub2$samples$Population
design <- model.matrix(~ Cls + Pat)

# Estimate NB dispersion and QL dispersion.
yClstSub2 <- estimateDisp(yClstSub2, design=design)
qfit2 <- glmQLFit(yClstSub2, design)

# log-CPM
prior.count <- 1
zClstSub2 <- edgeR::cpm(yClstSub2, log=TRUE, prior.count=prior.count)

annot2 <- data.frame(Cluster=paste0("Cluster ", Cls), Patient=Pat, Population=Pop)
rownames(annot2) <- colnames(zClstSub2)
ann_colors2 <- list(Cluster=col.p[N], Patient=col.p[1:length(SamplesComb)], 
    Population=col.p2[1:2])
names(ann_colors2$Cluster) <- paste0("Cluster ", N)
names(ann_colors2$Patient) <- SamplesComb
names(ann_colors2$Population) <- c("Normal", "BRCA1-Pre")

contr <- rbind( matrix(1/(1-ncls), ncls, ncls),
                matrix(0, ncol(design)-ncls, ncls) )
diag(contr) <- 1
contr[1,] <- 0
rownames(contr) <- colnames(design)
colnames(contr) <- paste0("Cls", 1:ncls)
ctest <- list()
for(i in 1:ncls) ctest[[i]] <- glmQLFTest(qfit2, contrast=contr[,i])

top <- 15
pseudoMakers <- list()
for(i in 1:ncls) {
    ord <- order(ctest[[i]]$table$PValue, decreasing=FALSE)
    upreg <- ctest[[i]]$table$logFC > 0
    pseudoMakers[[i]] <- rownames(yClstSub2)[ord[upreg][1:top]]
}
Markers <- unlist(pseudoMakers)
Markers <- Markers[!duplicated(Markers)]

# Fig App S1C
pdf("FigS1C.pdf",height=12, width=9)
mat4 <- t(scale(t(zClstSub2[Markers, ])))
pheatmap(mat4, color=colorRampPalette(c("blue","white","red"))(100), border_color="NA",
    breaks=seq(-2,2,length.out=101), cluster_cols=TRUE, scale="none", fontsize_row=5, 
    show_colnames=FALSE, treeheight_row=70, treeheight_col=70, 
    clustering_method="ward.D2", annotation_col=annot2, annotation_colors=ann_colors2)
dev.off()


