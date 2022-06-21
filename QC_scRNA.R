if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scater")

install.packages('Seurat')
install.packages("hdf5r")

library(Seurat)
library(hdf5r)
suppressMessages(require(Seurat))
suppressMessages(require(scater))
suppressMessages(require(Matrix))

## read the data ##
v3.1k <- Read10X_h5("pbmc_1k_v3_filtered_feature_bc_matrix.h5", use.names = T)
v2.1k <- Read10X_h5("pbmc_1k_v2_filtered_feature_bc_matrix.h5", use.names = T)
p3.1k <- Read10X_h5("pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5", use.names = T)
colnames(v2.1k)

# select only gene expression data from the CITE-seq data.
p3.1k <- p3.1k$`Gene Expression`
head(p3.1k)
##First, create Seurat objects for each of the datasets, and then merge into one large seurat object##
sdata.v2.1k <- CreateSeuratObject(v2.1k, project = "v2.1k")
sdata.v3.1k <- CreateSeuratObject(v3.1k, project = "v3.1k")
sdata.p3.1k <- CreateSeuratObject(p3.1k, project = "p3.1k")

# merge into one single seurat object. Add cell ids just in case you have overlapping barcodes between the datasets.
alldata <- merge(sdata.v2.1k, c(sdata.v3.1k,sdata.p3.1k), add.cell.ids=c("v2.1k","v3.1k","p3.1k"))
# also add in a metadata column that indicates v2 vs v3 chemistry
chemistry <- rep("v3",ncol(alldata))
chemistry[Idents(alldata) == "v2.1k"] <- "v2"
alldata <- AddMetaData(alldata, chemistry, col.name = "Chemistry")
head(alldata$orig.idents)
table(chemistry)
# check number of cells from each sample, is stored in the orig.ident slot of metadata and is autmatically set as active ident.
table(Idents(alldata))
head(alldata@meta.data)

##We will manually calculate the proportion of mitochondrial reads and add to the metadata table ##
mito.percent = PercentageFeatureSet(alldata, pattern = "^MT-")
alldata = AddMetaData(alldata, mito.percent, col.name = "mito.percent")
head(alldata@meta.data)
#Calculate ribosomal proportion
ribo_percent = PercentageFeatureSet(alldata, pattern = "^RP[SL]")
alldata = AddMetaData(alldata, ribo_percent, col.name = "ribo_percent")
head(alldata@meta.data)
#another way of calculate and adding metadata ##
# mt.genes <- rownames(alldata)[grep("^MT-",rownames(alldata))]
# C<-GetAssayData(object = alldata, slot = "counts")
# 
# percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
# alldata <- AddMetaData(alldata, percent.mito, col.name = "percent.mito")

#########
VlnPlot(object = alldata, features = c("nCount_RNA", "nFeature_RNA", "mito.percent",
                                       "ribo_percent"), ncol = 2,)+ NoLegend()
#QC-measures as scatter plots
p1 = FeatureScatter(object = alldata, feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA") + NoLegend()
p2 = FeatureScatter(object = alldata, feature1 = "mito.percent", 
                    feature2 = "nFeature_RNA") + NoLegend()
p3 = FeatureScatter(object = alldata, feature1 = "ribo_percent", 
                    feature2 = "nFeature_RNA") + NoLegend()
P1 + p2 + p3

# visualize specfic sample
FeatureScatter(alldata, feature1 = "nFeature_RNA", 
               feature2 = "mito.percent", cells = WhichCells(alldata,
                                                             expression = Chemistry == "v2"))
#filtering the low quality cells###
idx = which(alldata$mito.percent < 25)
cell_filtered = WhichCells(alldata, cells = idx)
length(cell_filtered)
head(cell_filtered)
# length(WhichCells(alldata))
selected_cells = subset(alldata, cells = cell_filtered)
head(selected_cells)
VlnPlot(object = selected_cells, features = c("nFeature_RNA", 
                                       "mito.percent"))
high_v3_exp = WhichCells(object = selected_cells, expression = nFeature_RNA > 4100)
high_V2_exp = WhichCells(object = selected_cells, expression = nFeature_RNA > 2000 &
                         orig.ident == "v2.1k")
selected_cells = subset(selected_cells, cells = c(high_V2_exp, high_v3_exp), invert = T)
ncol(selected_cells)

VlnPlot(object = selected_cells, features = "nFeature_RNA")
# filtered cells have low expression genes #
low_v3_exp = WhichCells(object = selected_cells, expression = nFeature_RNA < 1000 
                        & orig.ident != "v2.1k")
low_v2_exp = WhichCells(object = selected_cells, expression = nFeature_RNA < 500)
selected_cells = subset(selected_cells, cells = c(low_v2_exp, low_v3_exp), invert = T)
ncol(selected_cells)
VlnPlot(object = selected_cells, features = "nFeature_RNA")

table(Idents(alldata))
table(Idents(selected_cells))

# calculating cell cycle scores ##
selected_cells = CellCycleScoring(object = selected_cells, s.features = cc.genes$s.genes,
                                  g2m.features = cc.genes$g2m.genes)
VlnPlot(selected_cells, features = c("S.Score", "G2M.Score"))









