
##### Clustering ###
suppressMessages(require(Seurat))

pbmc = readRDS("pbmc3k.rds")
scaled_pbmc = pbmc@assays$RNA@scale.data
dis_euclidean = dist(t(scaled_pbmc))
ward_hclust = hclust(dis_euclidean, method = "ward.D2")
plot(ward_hclust, main = "dist = eucledian, Ward linkage", labels=FALSE)

cluster_hclust = cutree(ward_hclust, k = 7)
#table(pbmc@meta.data$celltype)
pbmc@meta.data$cluster_hclust = factor(cluster_hclust)
p1 = DimPlot(pbmc, reduction = "tsne", group.by = "cluster_hclust")
p2 = DimPlot(pbmc, reduction = "tsne", group.by = "celltype")
p1 + p2

FeaturePlot(pbmc, reduction = "tsne", features = c( 'GNLY','NKG7'))
head(pbmc@assays$RNA@var.features)

#### clustering by seurat ###
pbmc = FindNeighbors(pbmc, dims = 1:10)
pbmc = FindClusters(pbmc, resolution = 0.22) ##resolution is determine the distance between clusters #
# head(pbmc@meta.data$seurat_clusters)
p1 = DimPlot(pbmc, group.by = "seurat_clusters", reduction = 'tsne') +NoLegend()
p2 = DimPlot(pbmc, group.by = "celltype", reduction = 'tsne',
             label = T) + NoLegend()
p1+p2

FeaturePlot(pbmc, reduction = 'umap', features = c('CD19', 'CD3D', 'CD14', 'NKG7'))

new_ids = c('B cells','T cells and NK cells','Monocyte') ## which is which ????????###
names(new_ids) = levels(pbmc)
pbmc = RenameIdents(pbmc, new_ids)
DimPlot(pbmc, reduction = 'tsne', label = T)+ NoLegend()

saveRDS(pbmc, file = "pbmc3k.rds")







