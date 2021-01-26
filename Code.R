## Code has been modified from the Seurat Guided Clustering Tutorial that uses the Peripheral Blood Mononuclear Cells (pbmc) dataset obtained from 10x Genomics ##

# Expression matrix of mouse retinal progenitor cells across 10 timepoints of retinal development from 10x Genomics Chromium Single Cell system was downloaded using GEO Accession ID GSE118614 along with its metadata #

setwd('../../Desktop/UAlberta-Data-Analysis/Dr.Graf-RNA-Seq-Data/scRNA_retina/')

install.packages("Matrix")
install.packages("Seurat")
install.packages("dplyr")
install.packages("umap")
library(Matrix)
library(dplyr)
library(Seurat)

############# retina #############

retina.data <- Read10X(data.dir = "retina_github/")
retina.data

retina.metadata <- read.csv("retina_github/metadata_ex.csv", header = TRUE, row.names = NULL, sep = ',')
retina.metadata
class(retina.data)

metadata <- as.data.frame(retina.metadata)
head(metadata)
metadata$barcode
rownames(metadata) <- metadata$barcode

summary(colSums(retina.data))

# 27,998 genes and 1,07,052 single cells
dim(retina.data)

retina <- CreateSeuratObject(counts = retina.data, meta.data = metadata, project = "mouse_retina")
retina
table(Idents(retina))

#dense.size <- object.size(as.matrix(retina.data))
#dense.size

slotNames(retina)

CellsMetaIden = retina@active.ident
head(CellsMetaIden)

CellsMeta = retina@meta.data
head(CellsMeta)

sparse.size <- object.size(retina.data)
sparse.size

#dense.size/sparse.size

head(retina@meta.data)

retina[["percent.mt"]] <- PercentageFeatureSet(retina, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(retina, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(retina, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(retina, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#CombinePlots(plots = list(plot1, plot2))

retina <- subset(retina, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

retina <- NormalizeData(retina, normalization.method = "LogNormalize", scale.factor = 10000)

retina <- NormalizeData(retina)

retina <- FindVariableFeatures(retina, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(retina), 10)
top10
list(top10)

head(retina@var.genes)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(retina)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
#CombinePlots(plots = list(plot1, plot2))

all.genes <- rownames(retina)
#retina <- ScaleData(retina, features = all.genes) #computationally intensive
retina <- ScaleData(retina)

retina <- RunPCA(retina, features = VariableFeatures(object = retina))
print(retina[["pca"]], dims = 1:5, nfeatures = 5)

# NOTE: This process can take a long time for big datasets, comment out for expediency. More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
retina <- JackStraw(retina, num.replicate = 100)
retina <- ScoreJackStraw(retina, dims = 1:20)

JackStrawPlot(retina, dims = 1:15)

ElbowPlot(retina)

retina <- FindNeighbors(retina, dims = 1:10)
retina <- FindClusters(retina, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(retina), 5)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages ='umap-learn')#
retina <- RunUMAP(retina, dims = 1:10)

TSNEPlot(object = retina, group.by = "umap_CellType")

# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters

DimPlot(retina, reduction = "umap", label = TRUE, group.by = "umap_CellType")
DimPlot(retina, reduction = "umap", label = TRUE, group.by = "age")

saveRDS(retina, file = "retina_github/retina_ids_11.rds")

### Other Analysis ###

retina <- readRDS("retina_ids_11.rds")

# find all markers of cluster 1
cluster1.markers <- FindMarkers(retina, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(retina, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
retina.markers <- FindAllMarkers(retina, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
retina.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

cluster1.markers <- FindMarkers(retina, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(retina, features = c("Opn1sw", "Pax6"))
#pdf(file = 'volcanoplot_markers.pdf'), v)
# you can plot raw counts as well
VlnPlot(retina, features = c("Opn1sw", "Pax6"), slot = "counts", log = TRUE)
FeaturePlot(retina, features = c("Bmp7", "Calb1", "Nefh", "Gfap", "Pax6", "Rho", "Rcvrn", "Opn1sw", "Opn1mw", "Opn1lw"))
#single gene features needed!
FeaturePlot(retina, features = "Bmp7")
FeaturePlot(retina, features = "Calb1")
FeaturePlot(retina, features = "Nefh")
FeaturePlot(retina, features = "Gfap")
FeaturePlot(retina, features = "Pax6")
FeaturePlot(retina, features = "Rho")
FeaturePlot(retina, features = "Rcvrn")
FeaturePlot(retina, features = "Opn1sw")
FeaturePlot(retina, features = "Opn1mw")

top10 <- retina.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10
DoHeatmap(retina, features = top10$gene) + NoLegend()

#Idents(object = retina) <- "age"
#table(Idents(retina))

#E11   E12   E14   E16   E18    P0   P14    P2    P5    P8 
#7806   401 23760  4747 19786  8762  3458 15370  5352 10155 

Idents(object = retina) <- "umap_CellType"
table(Idents(retina))

#levels(retina)
#names(new.cluster.ids) <- levels(retina)
retina <- RenameIdents(retina, new.cluster.ids)

DimPlot(retina, reduction = "umap", label = TRUE, group.by = "umap_CellType" , pt.size = 0.5) + NoLegend()

saveRDS(retina, file = "retina_ids_11_final.rds")

### rds to AnnData ###

install.packages("devtools")
library(devtools)
devtools::install_github("cellgeni/sceasy")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("LoomExperiment", "SingleCellExperiment"))

install.packages('reticulate')

Sys.setenv(RETICULATE_PYTHON="C:/Users/priya/AppData/Local/r-miniconda/envs/r-reticulate/python")
conda_install(envname = "r-reticulate", packages="loompy") ### Fixes loompy issues ###
conda_install(envname = "r-reticulate", packages="anndata")
library(sceasy)
library(reticulate)
use_condaenv('r-reticulate')
loompy <- reticulate::import('loompy')
anndata <- import('anndata')

### Conversion ###

sceasy::convertFormat(retina, from="seurat", to="anndata",outFile= 'retina_ids_11_final.h5ad')
