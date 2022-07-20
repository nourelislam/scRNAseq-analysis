##### Normalization ###
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scater")
BiocManager::install("scran")

library(scater)
library(scran)
library(Seurat)
############# QC by SCE ###
alldata_sce = SingleCellExperiment(assays = list(counts = as.matrix(v3.1k)))
alldata_sce = addPerCellQC(alldata_sce, subset = list(MT = grep("^MT-", rownames(alldata_sce))))
alldata_sce = addPerFeatureQC(alldata_sce)
head(rowData(alldata_sce))
ncol(alldata_sce)
alldata_sce = alldata_sce[, alldata_sce$subsets_MT_percent < 20]
ncol(alldata_sce)
alldata_sce = alldata_sce[, alldata_sce$detected > 1000 & alldata_sce$detected < 4100]
ncol(alldata_sce)

assay(alldata_sce,"logcounts_raw") <- log2(counts(alldata_sce) + 1)
head(alldata_sce)
plotRLE(alldata_sce[,1:50], exprs_values = "logcounts_raw", style = "full")

raw.sce = runPCA(alldata_sce, exprs_values = "logcounts_raw")
p1 = scater::plotPCA(raw.sce, colour_by = "total")
p2 = plotReducedDim(raw.sce, dimred = "PCA", by_exprs_values = "logcounts_raw",
                    colour_by = "GNLY") ## color could be customized to gene or col.id #
p1+p2

####create suerat object for further normalization ##
pbmc.seu = CreateSeuratObject(counts(alldata_sce, project = "PBMC"))
pbmc.seu = NormalizeData(pbmc.seu)
pbmc.seu.sce = as.SingleCellExperiment(pbmc.seu)
pbmc.seu.sce = addPerCellQC(pbmc.seu.sce)
colData(pbmc.seu.sce)

pbmc.seu.sce = runPCA(pbmc.seu.sce)
p1 = scater::plotPCA(pbmc.seu.sce, colour_by = "LYZ")
p2 = plotReducedDim(pbmc.seu.sce, dimred = "PCA", colour_by = "GNLY")
p1+p2

### Normalization by scran ## 
qclust = scran::quickCluster(alldata_sce)
pbmc.sce = scran::computeSumFactors(alldata_sce, clusters = qclust)
colData(pbmc.sce)
summary(head(sizeFactors(pbmc.sce)))
pbmc.sce = logNormCounts(pbmc.sce)
plotRLE(pbmc.sce[,1:50], exprs_values = "logcounts", exprs_logged = FALSE, 
        style = "full")

pbmc.sce = runPCA(pbmc.sce)
p1 = scater::plotPCA(pbmc.sce, colour_by = "total")
p2 = plotReducedDim(pbmc.sce, dimred = "PCA", colour_by = "GNLY")
p1+p2


pbmc.seu = FindVariableFeatures(pbmc.seu, selection.method = "vst")
top10 = head(VariableFeatures(pbmc.seu), 10)
vplot = VariableFeaturePlot(pbmc.seu)
LabelPoints(plot = vplot, points = top10, repel = T)

head(pbmc.seu[["RNA"]][[]])
saveRDS(pbmc.seu, file = "pbmc3k.rds")









