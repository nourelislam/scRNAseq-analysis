install.packages('BiocManager')
BiocManager::install('multtest')
install.packages('metap')

library(multtest)
library(metap)
library(Seurat)
integ_Seurat = readRDS('data_created/integ_Seurat')
### make the default assay teh normalized values ##
DefaultAssay(integ_Seurat) = 'RNA'
cluster0_conserved_marker = FindConservedMarkers(integ_Seurat, ident.1 = 0,
                                                 grouping.var = 'sample', logfc.threshold = 0.25)
get_conserved = function(cluster_number){
  cluster0_conserved_marker = FindConservedMarkers(integ_Seurat, ident.1 = 0,
                                                   grouping.var = 'sample', logfc.threshold = 0.25)
}

















