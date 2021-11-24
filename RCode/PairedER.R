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
# Fig 6 B, 7 B, EV2 B

### Samples to be combined
Samples <- c("ER-0029-7C", "ER-0029-9C")
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
Anchors <- FindIntegrationAnchors(object.list=CombSeurat, dims=1:dimUsed)

# Integration
PairedER <- IntegrateData(anchorset=Anchors, dims=1:dimUsed, k.weight=100)

# Dim reduction
DefaultAssay(PairedER) <- "integrated"
PairedER <- ScaleData(PairedER, verbose=FALSE)
PairedER <- RunPCA(PairedER, npcs=dimUsed, verbose=FALSE)
PairedER <- RunTSNE(PairedER, dims=1:dimUsed, seed.use=2018)
tSNE <- PairedER@reductions$tsne@cell.embeddings

col.p2 <- c("darkblue", "gold2")

# Fig 6H left
Group <- factor(PairedER@meta.data$group, levels=SamplesComb)
col <- col.p2[Group]
plotOrd <- sample(ncol(PairedER))
pdf("Fig6H-left.pdf", height=9, width=9)
plot(tSNE[plotOrd,], pch=16, col=col[plotOrd], cex=0.7, xlab="tSNE-1", ylab="tSNE-2",
    main="t-SNE - By sample")
dev.off()

# Cell clustering
resolution <- 0.1
PairedER <- FindNeighbors(PairedER, dims=1:dimUsed, verbose=FALSE)
PairedER <- FindClusters(PairedER, resolution=resolution, verbose=FALSE)
Cluster <- as.integer(PairedER@meta.data$seurat_clusters)
ncls <- length(table(Cluster))

# Fig 6H right
col <- col.p[Cluster]
pdf("Fig6H-right.pdf", height=9, width=9)
plot(tSNE, pch=16, col=col, cex=0.7, xlab="tSNE-1", ylab="tSNE-2",
    main="t-SNE - By cluster")
dev.off()

# Combine raw data
allGenes <- rownames(get(DD[1]))
allGenesFilter <- c()
for(i in DGE) allGenesFilter <- c(allGenesFilter, rownames(get(i)))
allGenesFilter <- names(table(allGenesFilter))[table(allGenesFilter) >= 1]
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

#saveRDS(PairedER, file="SeuratObject_PairedER.rds")
#save(z, d, y, SamplesComb, Group, ncls, Cluster, col.p, col.p2, tSNE, file="PairedER.RData")


