library(Seurat)
library(ggplot2)

split_Seurat = readRDS('data_created/split_Seurat')
integ_features = SelectIntegrationFeatures(split_Seurat, nfeatures = 3000)
split_Seurat = PrepSCTIntegration(split_Seurat, anchor.features = integ_features)
integ_anchors = FindIntegrationAnchors(split_Seurat, anchor.features = integ_features,
                                      normalization.method = 'SCT')
integ_Seurat = IntegrateData(anchorset = integ_anchors, normalization.method = 'SCT')

### visualization ##
integ_Seurat = RunPCA(integ_Seurat)
ElbowPlot(integ_Seurat, ndims = 50)
DimPlot(integ_Seurat, group.by = 'sample', split.by = 'sample')

integ_Seurat = RunUMAP(integ_Seurat, reduction = 'pca', dims = 1:40)
p1 = DimPlot(integ_Seurat, reduction = 'umap', group.by = 'sample') + ggtitle("SCT_integration")
p2 = DimPlot(merged.data, reduction = 'umap', group.by = 'sample') + ggtitle("Unintegration")
p1+p2
saveRDS(integ_Seurat, 'data_created/integ_Seurat')













