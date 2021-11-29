### R code from vignette source
col.pMedium <- scater:::.get_palette("tableau10medium")
col.pDark <- scater:::.get_palette("tableau20")[2*(1:10)-1]
col.pLight <- scater:::.get_palette("tableau20")[2*(1:10)]
col.p <- c(col.pDark, col.pLight)

counts <- read.delim("../Data/GeneCountMatrix.txt.gz", stringsAsFactors=FALSE)
genes <- counts[,1:2]
counts <- counts[,-c(1:2)]
group <- gsub("[.].*$","",colnames(counts))
group <- factor(group)

library(edgeR)
y <- DGEList(counts, genes=genes)
y$samples$group <- group

patient <- colnames(y)
patient <- gsub("[.][R,L,1,2]$","",patient)
patient <- gsub("^.*[.]","",patient)
patient <- factor(paste0("P", patient))

### Gene annotation
anno <- read.delim(file="../Data/170929_Homo_sapiens.gene_info.gz",
    header=FALSE, skip=1)[,-1]
colnames(anno) <- c("GeneID", "Symbol", "LocusTag", "Synonyms", "dbXrefs", "Chr",
    "map_location", "description", "Type", "Symbol_from_nomenclature_authority",
    "Full_name_from_nomenclature_authority", "Nomenclature_status", 
    "Other_designations", "Modification_date")
m <- match(genes$GeneID, anno$GeneID)

### Remove unwanted genes
genes <- cbind(genes, anno[m,c(2,8),drop=FALSE])[,c(1,3,2,4)]
y$genes <- genes
y <- y[!is.na(y$genes$Symbol),,keep.lib.size=FALSE]
ig1 <- grep("^IG", y$genes$Symbol)
ig2 <- grep("immunoglob", y$genes$description)
ig <- ig1[ig1 %in% ig2]
y <- y[-ig,,keep.lib.size=FALSE]
y$genes <- y$genes[,-4]
dim(y)

### Filtering and normalization
isexpr <- rowSums(cpm(y)>0.3) >=3
table(isexpr)
y <- y[isexpr,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

### MDS
#col.p <- c("black","red","dodgerblue","purple","orange","limegreen")
mds <- plotMDS(y, plot=FALSE)
pdf("../Figures/MDS.pdf", height=5, width=5)
par(mar=c(5,5,2,2))
plotMDS(mds, pch=16, col=col.p[group], main="", )
legend("right", legend=levels(group), pch=16, 
       text.col=col.p[1:nlevels(group)], col=col.p[1:nlevels(group)])
dev.off()

### Design
design <- model.matrix(~ 0 + group)
colnames(design) <- gsub("group", "", colnames(design))
colSums(design)

### BCV
y <- estimateDisp(y, design, robust=TRUE)
summary(y$prior.df)
plotBCV(y)

### voom with random effect
v <- voom(y, design, plot=TRUE)
dupcor <- duplicateCorrelation(v, design=design, block=patient)
dupcor$consensus
vfit <- lmFit(v, design, block=patient, correlation=dupcor$consensus)

### Contrast
contr <- makeContrasts(
             Basal_vs_LP = Basal - LP,
             Basal_vs_ML = Basal - ML,
             Basal_vs_Str = Basal - Stroma,
             LP_vs_ML = LP - ML,
             LP_vs_Str = LP - Stroma,
             ML_vs_Str = ML - Stroma, levels=design)

### DE analysis with TREAT
vfit <- contrasts.fit(vfit, contr)
trt <- treat(vfit, lfc=log2(1.5))
dt <- decideTests(trt)
summary(dt)

### Signature genes
sel_Basal <- dt[,1]==1 & dt[,2]==1 & dt[,3]==1
sel_LP <- dt[,1]==-1 & dt[,4]==1 & dt[,5]==1
sel_ML <- dt[,2]==-1 & dt[,4]==-1 & dt[,6]==1
sel_Str <- dt[,3]==-1 & dt[,5]==-1 & dt[,6]==-1
Basal <- as.character( vfit$genes$Symbol[sel_Basal] )
LP <- as.character( vfit$genes$Symbol[sel_LP] )
ML <- as.character( vfit$genes$Symbol[sel_ML] )
Str <- as.character( vfit$genes$Symbol[sel_Str] )
length(Basal)
length(LP)
length(ML)
length(Str)



##################################
### QC

dat <- read.csv("../Data/BulkQC.csv")
nd <- dat$NumMapped
nfd <- dat$NumTotal - dat$NumMapped
names(nd) <- dat$Sample

pdf("Mapping.pdf", height=5, width=15)
yupper <- 8.5e7
ylim <- c(0, yupper)
par(mar=c(7.5,5,2,1))
bp <- barplot(rbind(nd-nfd, nfd), space=0.3, beside=FALSE, col=c("grey70", "grey30"), 
    main="Number of read pairs", legend.text=c("Mapped", "Unmapped"), 
    args.legend=list(x="topright",bty="n"), las=2, border=NA, 
    axes=FALSE, ylab="", ylim=ylim)
axis(side=2, at=2e7*(0:4), las=2)
#text(x=bp, y=nd + yupper/20, labels=paste0(round(100*(nfd/nd),digits=1),"%"), xpd=TRUE, cex=0.7, srt=90)
dev.off()


source("fig_label.R")

pdf("../Figures/BulkRNAseq.pdf", height=8, width=9)
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE),
   heights=c(1.2, 1))
### Mapping statistics
par(mar=c(9,5,2,2))
yupper <- 8.5e7
ylim <- c(0, yupper)
bp <- barplot(rbind(nd-nfd, nfd), space=0.3, beside=FALSE, col=c("grey70", "grey30"), 
    main="Number of read pairs", legend.text=c("Mapped", "Unmapped"), 
    args.legend=list(x="topright",bty="n"), las=2, border=NA, 
    axes=FALSE, ylab="", ylim=ylim)
axis(side=2, at=2e7*(0:4), las=2)
fig_label("a.", cex=2)
### MDS
par(mar=c(6,6,2,3))
plotMDS(mds, pch=16, col=col.p[group], main="", )
legend("right", legend=levels(group), pch=16, 
       text.col=col.p[1:nlevels(group)], col=col.p[1:nlevels(group)])
fig_label("b.", cex=2)
### BCV
par(mar=c(6,6,2,3))
plotBCV(y)
fig_label("c.", cex=2)
dev.off()


