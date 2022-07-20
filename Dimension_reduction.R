library(scater)
library(scran)
library(Seurat)
library(ggplot2)

pbmc = readRDS(file = "pbmc3k.rds")
labels = read.delim(file = "celltype_labels.tsv", row.names = 1)
pbmc = AddMetaData(pbmc, metadata = labels)
pbmc = ScaleData(pbmc, verbose = F)
pbmc = RunPCA(pbmc, npcs = 50, verbose= F)
DimPlot(pbmc, reduction = "pca", group.by = "celltype")
DimHeatmap(pbmc, dims = 1:6, cells = 500)
ElbowPlot(pbmc, ndims = 50)

##### tSNE ###
pbmc = RunTSNE(pbmc)
DimPlot(pbmc, reduction  = "tsne", group.by = "celltype")
pbmc = RunTSNE(pbmc, dims = 1:20, reduction = "pca", perplexity = 50, max_iter = 1000) ##dim number based on the PCs observed in the previous SD with PCs #
DimPlot(pbmc, group.by = "celltype", reduction = "tsne" ) + ggtitle("20_PCs, perp50")
### UMAP ####
pbmc = RunUMAP(pbmc, dims = 1:30, verbose = F, n.neighbors = 30,
               min.dist = 0.5, n.epochs = 500, metric = "cosine")
DimPlot(pbmc, group.by = "celltype", label = T, pt.size = 0.6,
        repel = T) + NoLegend()

topfeatures_pc1 = TopFeatures(pbmc, nfeatures = 2, dim = 1)
topfeatures_pc2 = TopFeatures(pbmc, nfeatures = 2, dim = 2)
topfeatures_pcs = c(topfeatures_pc1, topfeatures_pc2)

FeaturePlot(pbmc, features = topfeatures_pcs, reduction = "pca")
FeaturePlot(pbmc, features = topfeatures_pcs, reduction = "tsne")
FeaturePlot(pbmc, features = topfeatures_pc2[80:90], reduction = "tsne")
interactivePlot = FeaturePlot(pbmc, reduction = "tsne", features = "CD3D")
HoverLocator(plot = interactivePlot, information = FetchData(pbmc,
                                                             vars = c("celltype", topfeatures_pcs)))

saveRDS(pbmc, file = "pbmc3k.rds")






















