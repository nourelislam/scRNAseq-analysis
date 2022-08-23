install.packages("tidyverse")

library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
library(ggplot2)

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

  
  
  





