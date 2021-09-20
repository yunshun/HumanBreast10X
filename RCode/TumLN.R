smythlab
cd cheny/Project/BPal/10X_Human/ScientificData/TumLN
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
# Fig 6 B, 7 B, EV2 B

### Samples to be combined
Samples <- c("ER-0040-T", "ER-0040-LN", "ER-0043-T", "ER-0043-LN", 
    "ER-0056-T", "ER-0056-LN", "ER-0064-T", "ER-0064-LN", "ER-0167-T", 
    "ER-0167-LN", "ER-0173-T", "ER-0173-LN", "mER-0068-T", "mER-0068-LN")
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

# Anchors + Integration
dimUsed <- 30
TumLN <- list()
for(i in 1:6){
    Anchors <- FindIntegrationAnchors(object.list=CombSeurat[(2*i-1):(2*i)], dims=1:dimUsed,
        anchor.features=1000, scale=TRUE, k.anchor=5, k.filter=30,
        k.score=20, max.features=100)
    so <- IntegrateData(anchorset=Anchors, dims=1:dimUsed, k.weight=100)
    DefaultAssay(so) <- "integrated"
    so <- ScaleData(so, verbose=FALSE)
    so <- RunPCA(so, npcs=dimUsed, verbose=FALSE)
    so <- RunTSNE(so, dims=1:dimUsed, seed.use=2018)
    TumLN[[i]] <- so
}

# mER-0068 Tum+LN
Anchors <- FindIntegrationAnchors(object.list=CombSeurat[13:14], dims=1:dimUsed)
so <- IntegrateData(anchorset=Anchors, dims=1:dimUsed, k.weight=100)
DefaultAssay(so) <- "integrated"
so <- ScaleData(so, verbose=FALSE)
so <- RunPCA(so, npcs=dimUsed, verbose=FALSE)
so <- RunTSNE(so, dims=1:dimUsed, seed.use=2018)
TumLN[[7]] <- so

# Names for Seurat list
SN <- SamplesComb[2*(1:7)]
SN <- gsub("_LN","",SN)
SN <- gsub("_","-",SN)
names(TumLN) <- SN

# Cell clustering
resolution <- c(0.1, 0.05, 0.1, 0.05, 0.1, 0.1, 0.1)
ncls <- c()
Cluster <- tSNE <- Group <- list()
for(i in 1:7){
    tSNE[[i]] <- TumLN[[i]]@reductions$tsne@cell.embeddings
    TumLN[[i]] <- FindNeighbors(TumLN[[i]], dims=1:dimUsed, verbose=FALSE)
    TumLN[[i]] <- FindClusters(TumLN[[i]], resolution=resolution[i], verbose=FALSE)
    Cluster[[i]] <- as.integer(TumLN[[i]]@meta.data$seurat_clusters)
    ncls[i] <- length(table(Cluster[[i]]))
    Group[[i]] <- factor(TumLN[[i]]@meta.data$group, levels=SamplesComb[(2*i-1):(2*i)])
}
names(Cluster) <- names(tSNE) <- names(Group) <- names(ncls) <- SN

col.p2 <- c("darkblue", "gold2")

# Fig 9A
pdf("Fig9A.pdf", height=10, width=35)
par(mfcol=c(2,7))
for(i in 1:7){
    col <- col.p2[Group[[i]]]
    plotOrd <- sample(ncol(TumLN[[i]]))
    plot(tSNE[[i]][plotOrd,], pch=16, col=col[plotOrd], cex=0.7, xlab="tSNE-1", ylab="tSNE-2", main=SN[i])
    col <- col.p[Cluster[[i]]]
    plot(tSNE[[i]], pch=16, col=col, cex=0.7, xlab="tSNE-1", ylab="tSNE-2", main="")
}
dev.off()


# Combine raw data
d <- y <- z <- list()
prior.count <- 1
for(j in 1:7){
    allGenes <- rownames(get(DD[2*j-1]))
    allGenesFilter <- c()
    for(i in DGE[(2*j-1):(2*j)]) allGenesFilter <- c(allGenesFilter, rownames(get(i)))
    allGenesFilter <- names(table(allGenesFilter))[table(allGenesFilter) >= 1]
    u <- get(DD[2*j-1])
    d0 <- u[allGenes, ]
    y0 <- u[allGenesFilter, ]
    u <- get(DD[2*j])
    d0 <- cbind(d0, u[allGenes, ])
    y0 <- cbind(y0, u[allGenesFilter, ])
    d0$samples$group <- y0$samples$group <- Group[[j]]
    d[[j]] <- d0
    y[[j]] <- y0
    z[[j]] <- edgeR::cpm(y0, log=TRUE, prior.count=prior.count)
}

names(d) <- names(y) <- names(z) <- SN

saveRDS(TumLN, file="SeuratObject_TumLN.rds")
save(z, d, y, SamplesComb, Group, ncls, Cluster, col.p, col.p2, tSNE, file="TumLN.RData")




