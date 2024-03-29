install.packages("tidyverse")
install.packages('BiocManager')
BiocManager::install('multtest')
install.packages('metap')

library(multtest)
library(metap)
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
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

## visualize the most variable PCs ##
integ_Seurat = readRDS('data_created/integ_Seurat')
ElbowPlot(integ_Seurat, ndims = 15) ### utilizing the 1st 8 PCs is the best option for clustering ##
DimHeatmap(integ_Seurat, dims = 1:9, cells = 500)
print(integ_Seurat[['pca']], dims = 1:10)
### starting clustering by determine the number of Neighbors ##
integ_Seurat = FindNeighbors(integ_Seurat, dims = 1:40)
integ_Seurat = FindClusters(integ_Seurat, resolution = c(0.4,0.6,0.8,
                                                         1.0,1.2,1.4))
# head(integ_Seurat@meta.data)
Idents(integ_Seurat) = "integrated_snn_res.0.8" 
DimPlot(integ_Seurat, reduction = 'umap', label = T) +ggtitle('Clustering_0.8_resolution/Neighbor_Dim20')
### change the identity ##
# Idents(integ_Seurat) = 'integrated_snn_res.0.4'
# DimPlot(integ_Seurat, label = T, reduction = 'umap') + ggtitle("Clustering_0.4_resolution")

#### QC of clustering ##
colnames(integ_Seurat@meta.data)
n_cells = FetchData(integ_Seurat, vars = c('orig.ident', 'ident')) %>%
  dplyr::count(ident, orig.ident) %>% tidyr::spread(ident, n)
DimPlot(integ_Seurat, reduction = 'umap', split.by = 'sample', label = T, group.by = 'ident')
### QC by cell cycle features ##
DimPlot(integ_Seurat, group.by = 'Phase', split.by = 'Phase') +ggtitle('QC by cellCycle features')
## QC Mito ##
DimPlot(integ_Seurat, split.by = 'mit.frc', group.by = "mit.frc")
features = c('nFeature_RNA','nCount_RNA','G2M.Score','S.Score')
FeaturePlot(integ_Seurat, features = features, label = T)

columns = c(paste0('PC_',1:16), 'ident','UMAP_1','UMAP_2')
PC_data = FetchData(integ_Seurat, vars = columns)
head(PC_data)
### the following code this is where the reduction score stored ##
integ_Seurat@reductions$pca@cell.embeddings[1:5, 1:2]
print(integ_Seurat[['pca']], dim = 1:5, nfeatures = 5)
DefaultAssay(integ_Seurat) = 'RNA'
integ_Seurat = NormalizeData(integ_Seurat, verbose = F)
### visualise the most variable gene and marker genes ##
top_variable_genes = FeaturePlot(integ_Seurat, features = c('FTL', 'GNLY', 'CCL2', 'HBB', 'HBA2', 'CXCL10', 'HBA1', 'CCL7', 'CCL3', 'CCL4'),
            reduction = 'umap', label = T)
FeaturePlot(integ_Seurat, reduction = 'umap', features = c('HBB', 'HBA2'), label = T)
##############
integ_Seurat = readRDS('data_created/integ_Seurat')
### make the default assay teh normalized values ##
DefaultAssay(integ_Seurat) = 'RNA'
cluster0_conserved_marker = FindConservedMarkers(integ_Seurat, ident.1 = 0,
                                                 grouping.var = 'sample', logfc.threshold = 0.25)
get_conserved = function(cluster_number){
  cluster0_conserved_marker = FindConservedMarkers(integ_Seurat, ident.1 = 0,
                                                   grouping.var = 'sample', logfc.threshold = 0.25)
}







