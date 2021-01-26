# Retina

The downstream analysis of the Single Cell RNA-Seq data from GEO Accession ID GSE118614 was performed using Seurat. The results were then visualized using non-linear embedding methods including tSNE and UMAP and a sparse matrix R data file ( .rds) file was generated. 
Sceasy was used to convert the file format from seurat’s .rds object to annData’s .h5ad  object

The final .h5ad file can be visualized in the interactive single cell transcriptome data explorer, cellxgene using the command: 

cellxgene launch retina_cellxgene.h5ad --open
