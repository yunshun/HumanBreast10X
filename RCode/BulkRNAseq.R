### R code 
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
col.p <- c("black","red","dodgerblue","purple","orange","limegreen")
mds <- plotMDS(y, plot=FALSE)
plotMDS(mds, pch=16, col=col.p[group], main="MDS")
legend("right", legend=levels(group), pch=16, 
       text.col=col.p[1:nlevels(group)], col=col.p[1:nlevels(group)])

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


