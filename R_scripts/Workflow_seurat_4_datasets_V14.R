# Note: if a library does not load, then use this below
# install.packages('PACKAGE_NAME)
# And if that does not work, then:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("PACKAGE_NAME")
#if (!require("pacman")) install.packages("pacman")
pacman::p_load(here,
               tidyverse,
               janitor, # Cleaning column names
               scales, # Transform axis scales
               ggrepel)

# loading in all important libraries
require("devtools")
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tidyr)
library(mvoutlier)
library(limma)
library(knitr)
library(SingleCellExperiment)
library(scater)
library(RColorBrewer)
library(plot3D)
library(stringr)
library(SAVER)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(ggpubr)
library(circlize)
library(cowplot)
library(clustree)
library(grid)
library(gridExtra)
library(ape)
library(ggplot2)
library(DO.db)
library(clusterProfiler)
library(BiocParallel)
library(scANANSESeurat)
library(viridis)
library(DESeq2)
#install.packages("ggalluvial")
library(ggalluvial)
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
library(fgsea)

###########################################################################
# Storing results in directories of your choice
workdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/"

# setting up results directory
dateoftoday <- gsub("-", "", as.character(Sys.Date()))
resultsdir <- paste0(workdir, dateoftoday)
system(paste("mkdir -p ", resultsdir))

# Loading in all functions
source("modify_vlnplot_v1.R")
source("stacked_vlnplot_v1.R")
source("QC_seurobj_v1.R")
source("PCA_seurobj_v1.R")
source("QC_seurobj_clustering_v1.R")
source("clustres_seurobj_v1.R")
source("clustres_seurobj_int_v1.R")

###########################################################################
# Importing the separate datasets and removing the doublets
reference <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_Collin_2021/20220113 CollinRNAannotated.rds")
query <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_Catala2021/20220114 CatalaRNAannotated.rds")
query2 <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_Li2021/20220114 LiRNAannotated.rds")
query3 <-readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_Gautam2022/...")

reference <- AddMetaData(reference, metadata="Co", col.name="Condition")

query <-AddMetaData(query, metadata="Ca", col.name="Condition")

query2 <- AddMetaData(query2, metadata="Li", col.name="Condition")

query3 <- AddMetaData(query3, metadata="Ga", col.name="Condition")

###########################################################################
# Remove doublets for all datasets

# Collin
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
seu_Collin <- reference
sweep.res.list_Collin <- paramSweep_v3(seu_Collin, PCs = 1:30, sct = FALSE)
sweep.stats_Collin <- summarizeSweep(sweep.res.list_Collin, GT = FALSE)
bcmvn_Collin <- find.pK(sweep.stats_Collin)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(seu_Collin@meta.data$costum_clustering)           ## ex: annotations <- seu_Collin@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(seu_Collin@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------

# Choose pK 0.06 based on the highest BC metrix
seu_Collin <- doubletFinder_v3(seu_Collin, PCs = 1:30, pN = 0.25, pK = 0.06, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
#seu_Collin <- doubletFinder_v3(seu_Collin, PCs = 1:30, pN = 0.25, pK = 0.06, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

DimPlot(seu_Collin, label=F,label.size = 8, group.by = 'DF.classifications_0.25_0.06_1806') + ggtitle("Doublet estimation")

reference <- subset(x=seu_Collin, subset = (DF.classifications_0.25_0.06_1806 == "Singlet"))
saveRDS(reference, file = paste0(resultsdir,"/Collin_singlets.rds"))

# Catala
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
seu_Catala <- query
sweep.res.list_Catala <- paramSweep_v3(seu_Catala, PCs = 1:18, sct = FALSE)
sweep.stats_Catala <- summarizeSweep(sweep.res.list_Catala, GT = FALSE)
bcmvn_Catala <- find.pK(sweep.stats_Catala)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(seu_Catala@meta.data$costum_clustering)           ## ex: annotations <- seu_Catala@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(seu_Catala@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------

# Choose pK 0.005 based on the highest BC metrix
seu_Catala <- doubletFinder_v3(seu_Catala, PCs = 1:18, pN = 0.25, pK = 0.005, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
#seu_Catala <- doubletFinder_v3(seu_Catala, PCs = 1:18, pN = 0.25, pK = 0.06, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

DimPlot(seu_Catala, label=F,label.size = 8, group.by = 'DF.classifications_0.25_0.005_1289') + ggtitle("Doublet estimation")

query <- subset(x=seu_Catala, subset = (DF.classifications_0.25_0.005_1289 == "Singlet"))
saveRDS(query, file = paste0(resultsdir,"/Catala_singlets.rds"))

# Li
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
seu_Li <- query2
sweep.res.list_Li <- paramSweep_v3(seu_Li, PCs = 1:16, sct = FALSE)
sweep.stats_Li <- summarizeSweep(sweep.res.list_Li, GT = FALSE)
bcmvn_Li <- find.pK(sweep.stats_Li)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(seu_Li@meta.data$costum_clustering)           ## ex: annotations <- seu_Li@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(seu_Li@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------

# Choose pK 0.18 based on the highest BC metrix
seu_Li <- doubletFinder_v3(seu_Li, PCs = 1:16, pN = 0.25, pK = 0.18, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
#seu_Li <- doubletFinder_v3(seu_Li, PCs = 1:16, pN = 0.25, pK = 0.18, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

DimPlot(seu_Li, label=F,label.size = 8, group.by = 'DF.classifications_0.25_0.18_1347') + ggtitle("Doublet estimation")

query2 <- subset(x=seu_Li, subset = (DF.classifications_0.25_0.18_1347 == "Singlet"))
saveRDS(query2, file = paste0(resultsdir,"/Li_singlets.rds"))

# Gautam
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
seu_Gautam <- query3
sweep.res.list_Gautam <- paramSweep_v3(seu_Gautam, PCs = 1:24, sct = FALSE)
sweep.stats_Gautam <- summarizeSweep(sweep.res.list_Gautam, GT = FALSE)
bcmvn_Gautam <- find.pK(sweep.stats_Gautam)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(seu_Gautam@meta.data$costum_clustering)           ## ex: annotations <- seu_Gautam@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(seu_Gautam@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------

# Choose pK 0.18 based on the highest BC metrix
seu_Gautam <- doubletFinder_v3(seu_Gautam, PCs = 1:24, pN = 0.25, pK = 0.14, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
#seu_Gautam <- doubletFinder_v3(seu_Gautam, PCs = 1:24, pN = 0.25, pK = 0.14, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

DimPlot(seu_Gautam, label=F,label.size = 8, group.by = 'DF.classifications_0.25_0.14_831') + ggtitle("Doublet estimation")

query3 <- subset(x=seu_Gautam, subset = (DF.classifications_0.25_0.14_831 == "Singlet"))
saveRDS(query3, file = paste0(resultsdir,"/Gautam_singlets.rds"))

############################################################
# Reload datasets with only Singlets
reference <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220207/Collin_singlets.rds")
query <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220207/Catala_singlets.rds")
query2 <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220207/Li_singlets.rds")
query3 <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220211/Gautam_singlets.rds")

fullcornea_epi <- merge(reference, y = c(query,query2,query3), add.cell.ids = c("Collin", "Catala","Li","Gautam"), project = "fullcornea_epi")
unname(fullcornea_epi$Condition)

# Setting the correct parameters for the anndata conversion
fullcornea_epi$batch <-unname(fullcornea_epi$Condition)
fullcornea_epi@meta.data$batch <- unname(fullcornea_epi$Condition)
# generate training data for the machine learning model for every group

# Save it as a separate object raw integration of the four datasets
saveRDS(fullcornea_epi, file = paste0(resultsdir,"/fullcornea_epi.rds"))

############################################################
# generate a contingency table for 4 dataset integration
# import the datasets
scVI_4 <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220214/fullcornea_epi.rds")
scVI_labels_4 <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/scVI_leiden_clusters_fc_epi_0.35_15022022.tsv', sep= '\t', header = F, row.names=1,comment.char = "")

# check if the rownames are the same after processing
sum(names(scVI_4@active.ident) == rownames(scVI_labels_4))

# add the IDs if correct
scVI_4$scVI_label <- as.factor(scVI_labels_4$V2)

# Re-name the clustering of the UMAP
cluster_order <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)

# Generate new cluster labeling:
current.cluster.ids <- c(6,3,2,8,0,9,1,12,18,4,7,10,11,5,13,15,16,14,17,20,19)

new.cluster.ids <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)

scVI_4$scVI_label <- plyr::mapvalues(x = as.factor(scVI_4$scVI_label), from = current.cluster.ids, to = new.cluster.ids)
scVI_4$scVI_label <- as.factor(scVI_4$scVI_label)

# Generate the probability table
cont_table_4 <- prop.table(table(scVI_4$scVI_label,scVI_4$costum_clustering),margin=2)*100

# Order the rows and columns to make it more clear
pdf(paste(resultsdir,'heatmap_cont_table_4_dataset.pdf',sep="/") ,width=10,height=5,paper='special')
f1 = colorRamp2(c(0, 100), c("white", "darkred"), space = "RGB")

# Order the rows and columns to make it more clear
cont_df <- as.data.frame.matrix(cont_table_4)
vec <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)
vec_num<- vec
cont_df<-cont_df[order(match(rownames(cont_df), vec)), , drop = FALSE]
vec2 <- c("Ca_C10_CLC","Co_C6_LPC","Li_C6_PC","Ga_C5_Cj","Ga_C7_Cj","Li_C5_DC","Ga_C3_THCEC","Ga_C9_THCEC","Ca_C8_LSC","Co_C5_LNCP","Li_C1_DC","Li_C3_DC","Ca_C7_WSEL","Li_C8_TALSC","Li_C4_PC","Li_C2_Cj","Ga_C2_Cj",
          "Co_C3_CjB","Co_C2_CjS","Ga_C6_Cj","Ca_C5_BCE","Co_C1_SCE","Ca_C6_LCE","Ga_C4_EHCEC","Co_C10_BCE",
          "Co_C4_CSSC","Co_C7_LSt","Ca_C1_ASK","Ca_C4_TSK","Ca_C3_GSK","Ca_C2_GSK","Ga_C1_CF","Co_C8_LF","Ga_C10_CF","Ca_C11_CES",
          "Ca_C9_CEM","Li_C9_LSC","Co_C14_LV","Co_C9_BV","Ga_C8_Mel","Li_C7_Mel",
          "Ca_C12_UN","Co_C11_MEC","Ca_C13_ESD","Li_C10_UN","Ga_C12_CTC","Co_C13_IC","Ga_C11_Mon","Co_C12_FCEC")

cont_df<-cont_df[vec2]
matz <- as.matrix(cont_df)

vec3 <- as.character(vec)

ht1 = Heatmap(matz, col = f1, cluster_columns = F,cluster_rows = F, name = "Percentage_found"
              ,right_annotation  = rowAnnotation(foo = anno_text(vec3, location = 0.5, just = "center")))
htlist1 <- ht1

v = c("Limbal stem cell","Limbal stem cell","Limbal epithelial","Limbal epithelial","Conjunctival","Conjunctival"
      ,"Corneal epithelial","Corneal epithelial","Corneal epithelial","Stromal","Stromal","Stromal",
      "Stromal","Stromal","Stromal","Endothelial","Vessels","Melanocytes","Immune cells","Unknown","Unknown")
f2 <-viridis(10,option = "H",alpha = 0.8)
z = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)

# Calculating the number of cells
num_cells <- c()
for (i in rownames(matz)){
  print(i)
  num_cells <- c(num_cells,sum(scVI_4@meta.data$scVI_label == i))
}

# Determining the colors of the subsets
show_col(viridis_pal(option = "H",alpha = 0.8)(10))
vec0<- viridis(10,option = "H",alpha = 0.8)
show_col(viridis_pal(option = "H",alpha = 0.8)(40))
vec <- viridis(40,option = "H",alpha = 0.8)
print(vec)

newcol <- c(vec[19:20],vec[15:16],vec[2:3],vec[5:7],vec[26:31],vec0[3],vec0[10],vec0[7],vec0[4],"#FF0000",vec[37])

ha = rowAnnotation(foo = anno_text(as.character(num_cells), location = 0.5, just = "center",
                                   gp = gpar(fill = "black", col = "white", border = "black"),
                                   width = max_text_width(as.character(num_cells))*1.2))

f3 <- newcol
small_df <- data.frame(num_cells,row.names = vec_num)
small_mat <- as.matrix(small_df)

ht_list = htlist1 + Heatmap(v, col = f2,name = "Cell origin", width = unit(0.3, "cm"), heatmap_legend_param = list(
  at = unique(v),
  labels = unique(v),
  title = "Cell origin",
  col = f2,
  legend_height = unit(4, "cm"))) +Heatmap(z, name = "# of cells", col = f3,cluster_rows = F, width = unit(1.5, "cm"), show_heatmap_legend = FALSE,
                                                                                              cell_fun = function(j, i, x, y, width, height, fill) {
                                                                                                grid.text(sprintf("%.1f", small_mat[i, j]), x, y, gp = gpar(fontsize = 11))
                                                                                              })
print(draw(ht_list))
dev.off()

z2 <- z-1
colmat <- matrix(c(z2, newcol), ncol = 2)
# Write the colors as a table
write.table(colmat, file = paste0(resultsdir,'/unannotated_meta_colors.tsv'), sep = '\t')

marker_gene_file <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/LSC_Marker_Genes.csv"
marker_genes_df <- read.table(marker_gene_file, header = TRUE, sep = ',',comment.char = "#", stringsAsFactors = F)
seur_obj <- scVI_4
DefaultAssay(seur_obj) <- "RNA"

seur_obj@active.ident <- seur_obj$scVI_label
seur_obj@active.ident <- factor(seur_obj@active.ident, 
                                levels=vec_num)

fullmarkers2<-NULL

# Make dotplots for all datasets with their cell populations and marker genes
for (paper in unique(marker_genes_df$abreviation)){
  marker_dir <- paste0(paste0(resultsdir, '/marker_genes_scVI_4_int/'),paper)
  print(marker_dir)
  system(paste("mkdir -p ", marker_dir))
  sub_marker_df <- marker_genes_df[marker_genes_df$abreviation == paper,]
  
  fullmarkers <- NULL
  for (cell_type in unique(sub_marker_df$name_plot)){
    marker_genes <- sub_marker_df[sub_marker_df$name_plot == cell_type,]$marker_genes
    marker_genes <- unlist(strsplit(marker_genes, ','))
    fullmarkers <- c(fullmarkers,marker_genes)
    m <- length(fullmarkers)
    n <- length(unique(sub_marker_df$name_plot))
    if (m<30){
      w <- 10
    }else{
      w <- m/5
    }
  }
  fullmarkers <- unique(fullmarkers)
  fullmarkers2<-c(fullmarkers2,fullmarkers)
  pdf(paste0(paste(paste(marker_dir, paper, sep = '/'), sep = '_'),'full_markers.pdf'), width =w, height = 5)
  print(DotPlot(object = seur_obj, features = fullmarkers,cluster.idents = T,)+ RotatedAxis())
  dev.off()
}

fullmarkers2 <- unique(fullmarkers2)
m <- length(fullmarkers2)

#cluster markergenes based on expression by complex heatmap
# Generate the pseudobulk-table for correlation between datasets
sce_qc <- SingleCellExperiment(assays = list(counts = GetAssayData(object = seur_obj, slot = "data")))
pseudobulk_df <- NULL
pseudobulk_df <- as.data.frame(row.names(counts(sce_qc)))
for (i in unique(seur_obj@active.ident)){
pseudobulk_df[[i]] <- as.vector(rowSums(counts(sce_qc)[,seur_obj@active.ident == i]))
}
pseudobulk_df <- pseudobulk_df[pseudobulk_df$`row.names(counts(sce_qc))`%in%fullmarkers2,]
rownames(pseudobulk_df) <- pseudobulk_df$`row.names(counts(sce_qc))`
pseudobulk_df<-pseudobulk_df[,-1]

matz <- as.matrix(pseudobulk_df)
ht1 = Heatmap(matz, col = f1, cluster_columns = T,cluster_rows = F)
colord <- row_order(ht1)
matz <- matz[colord,]
new_order <- rownames(matz)

pdf(paste0(resultsdir,'/full_markers_scVI4.pdf'), width =m/5, height = 5)
print(DotPlot(object = seur_obj, features = new_order,cluster.idents = F)+ RotatedAxis())
dev.off()

neworder2<- c("ACTA2","CDH19","CCL3","MITF","KIT","ACKR1","COL4A3","COL1A1","COL3A1","COL4A1","COL5A1","FBLN1","KERA","LUM","CD34","MMP2","MMP3","MMP1","PAX6","KRT3","S100A9","KRT7","MUC1","TP63","KRT6A","KRT14","CPVL","GPHA2")

seur_obj@active.ident <- seur_obj$scVI_label
vec4 <- rev(vec_num)
seur_obj@active.ident <- factor(seur_obj@active.ident, 
                                levels=vec4)
neworder3<-rev(neworder2)

pdf(paste0(resultsdir,'/sub_markers_scVI4.pdf'), width =7, height = 5)
print(DotPlot(object = seur_obj, features = neworder3,cluster.idents = F)+ RotatedAxis())
dev.off()

# mesenchymal stromal cells (cd73 = NT5E; thy1 = CD90; CD105 = END)

neworder4<- c("KRT19","PROM1","CD34","ALDH1A1","ALDH3A1","ACTA2","NT5E","THY1","ENG","MMP2",
              "MMP3","MMP1","MMP12","MMP","MMP14","COL1A1","COL3A1","COL4A1","COL5A1","B3GNT7",
              "CHST6","FN1","EDA","ABCG2","BMI1","ALCAM","NOTCH1","SIX2","CXADR","PTGDS","PDK4","LY6D","LY6E","LY6H","LY6K","LYPD2")
pdf(paste0(resultsdir,'/new_markers_scVI4.pdf'), width =8, height = 5)
print(DotPlot(object = seur_obj, features = neworder4,cluster.idents = F)+ RotatedAxis())
dev.off()

neworder4<- c("MUC1","KRT7","KRT8","MMP2","MMP9","MMP14")
pdf(paste0(resultsdir,'/new_markers_Cj_scVI4.pdf'), width =8, height = 5)
print(DotPlot(object = seur_obj, features = neworder4,cluster.idents = F)+ RotatedAxis())
dev.off()


neworder4<-c("KRT15","SOX9","ACTN1","FZD7","KRT17","ATF3","IFITM3","CD63","MT1A","SOCS3")
pdf(paste0(resultsdir,'/new_markers_LESC_scVI4.pdf'), width =8, height = 5)
print(DotPlot(object = seur_obj, features = neworder4,cluster.idents = F)+ RotatedAxis())
dev.off()


neworder4<-c("MMP10","KRT24","KRT15","GPHA2","KRT12","KRT3","KRT4","CEACAM7","PHLDA1","LAMA5","LAMA2","LAMA4","LAMB2","LAMB3",
             "LAMC2","FLG","IVL","GJB2","GJB6","GJA1","HSPG2","ITGA3","ITGB4","CAV1","CXCL14","CKS2","MOXD1","LAMA3","TNFRSF21")
pdf(paste0(resultsdir,'/new_markers_staining_scVI4.pdf'), width =8, height = 5)
print(DotPlot(object = seur_obj, features = neworder4,cluster.idents = F)+ RotatedAxis())
dev.off()

neworder2<- c("ACTA2","SCN7A","NGFR","CDH19","SOX10","CCL3","MITF",
              "LYVE1","ACKR1","CA3","COL5A1","COL1A1","FBLN1",
              "MMP3","MMP2","MMP1","CD34","KERA",
              "LUM","AREG","KRT24","KRT12","KRT3",
              "S100A9","KRT7","MUC1","CXCL17",
              "CPVL","GPHA2","KRT15","TP63","S100A2","KRT14","PAX6")

neworder3<-rev(neworder2)

pdf(paste0(resultsdir,'/final_markers_scVI4.pdf'), width =9, height = 5)
print(DotPlot(object = seur_obj, features = neworder3,cluster.idents = F)+ RotatedAxis())
dev.off()

neworder2<- c("IGFBP3","KRT78","LAMA3","GJB6","GJB2","ABCB5","VIM","KRT19","ABCG2","ACTA2","SCN7A","NGFR","CDH19","SOX10","CCL3","MITF",
              "LYVE1","ACKR1","CA3","COL5A1","COL1A1","FBLN1",
              "MMP3","MMP2","MMP1","CD34","KERA",
              "LUM","AREG","KRT24","KRT12","KRT3",
              "S100A9","KRT7","MUC1","CXCL17",
              "CPVL","GPHA2","KRT15","TP63","S100A2","KRT14","PAX6")

neworder3<-rev(neworder2)

pdf(paste0(resultsdir,'/extended_markers_scVI4.pdf'), width =12, height = 5)
print(DotPlot(object = seur_obj, features = neworder3,cluster.idents = F)+ RotatedAxis())
dev.off()

# These are the last markers remove above ...
neworder2<-c("PAX6", "KRT14", "S100A2", "CXCL14", "TP63", "KRT15", "GPHA2", "CPVL", "MKI67", "CXCL17", 
  "MUC1", "KRT7", "S100A8", "S100A9", "KRT13", "AQP5", "KRT3", "KRT12", "KRT24", 
  "AREG", "GJB6", "VIM", "LUM", "KERA", "CD34", "AQP1", "MMP1", "MMP2", "MMP3", "THY1", 
  "NT5E", "FBLN1", "PDGFRA", "COL8A2", "CA3", "SLC4A11", "ACKR1", "PECAM1", "LYVE1", "MITF", 
  "PMEL", "CCL3", "SOX10", "CDH19", "NGFR", "SCN7A", "TAGLN", "ACTA2", "NOTCH3")

new.cluster.ids<- c(
  "1: LSC", "2: LESC", "3: LE", "4: LE", "5: Cj", "6: Cj", "7: CE", "8: CE", "9: CE",
  "10: qSK", "11: SK", "12: SK", "13: TSK", "14: CF", "15: CF", "16: EC", "17: Ves",
  "18: Mel", "19: IC", "20: nm-cSC", "21: MC"
)

cluster_order<-new.cluster.ids

current.cluster.ids <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)

scVI_4$scVI_label <- plyr::mapvalues(x = as.factor(scVI_4$scVI_label), from = current.cluster.ids, to = new.cluster.ids)
scVI_4$scVI_label <- as.factor(scVI_4$scVI_label)

seur_obj <- scVI_4
DefaultAssay(seur_obj) <- "RNA"

seur_obj@active.ident <- seur_obj$scVI_label
seur_obj@active.ident <- factor(seur_obj@active.ident, 
                                levels=rev(cluster_order))

pdf(paste0(resultsdir,'/extended_markers_final_scVI4.pdf'), width =13, height = 4.5)
print(DotPlot(object = seur_obj, features = neworder2,cluster.idents = F)+ RotatedAxis())
dev.off()
################################################################################
# Correlation analysis of cells determine if the clusters can be joined or not
# Generate the pseudobulk-table for correlation between datasets
sce_qc <- SingleCellExperiment(assays = list(counts = GetAssayData(object = seur_obj, slot = "data")))

sce_qc <- SingleCellExperiment(assays = list(counts = GetAssayData(object = seur_obj, slot = "data")))

sce_qc$sample <- seur_obj@active.ident
sce_qc$sample

pseudobulk_df <- NULL

pseudobulk_df <- as.data.frame(row.names(counts(sce_qc)))
pseudobulk_df[["C3"]] <- as.vector(rowSums(counts(sce_qc)[,sce_qc$sample == 3]))
pseudobulk_df[["C4"]] <- as.vector(rowSums(counts(sce_qc)[,sce_qc$sample == 4]))

pseudobulk_df[["C5"]] <- as.vector(rowSums(counts(sce_qc)[,sce_qc$sample == 5]))
pseudobulk_df[["C6"]] <- as.vector(rowSums(counts(sce_qc)[,sce_qc$sample == 6]))

pseudobulk_df[["C7"]] <- as.vector(rowSums(counts(sce_qc)[,sce_qc$sample == 7]))
pseudobulk_df[["C8"]] <- as.vector(rowSums(counts(sce_qc)[,sce_qc$sample ==8]))
pseudobulk_df[["C9"]] <- as.vector(rowSums(counts(sce_qc)[,sce_qc$sample == 9]))

pseudobulk_df[["C11"]] <- as.vector(rowSums(counts(sce_qc)[,sce_qc$sample == 11]))
pseudobulk_df[["C12"]] <- as.vector(rowSums(counts(sce_qc)[,sce_qc$sample == 12]))

pseudobulk_df[["C14"]] <- as.vector(rowSums(counts(sce_qc)[,sce_qc$sample == 14]))
pseudobulk_df[["C15"]] <- as.vector(rowSums(counts(sce_qc)[,sce_qc$sample == 15]))

row.names(pseudobulk_df) <- pseudobulk_df$`row.names(counts(sce_qc))`
pseudobulk_df$`row.names(counts(sce_qc))` <- NULL

# Generate a function to perform spearman correlation
var <- FindVariableFeatures(seur_obj)

pseudobulk_df <- pseudobulk_df[rownames(pseudobulk_df) %in%var@assays$RNA@var.features,]

library(Kendall)

noms_x <- c()
noms_y  <- c()
vals <- c()
p.val <- c()
for (i in unique(colnames(pseudobulk_df))){
  for (j in unique(colnames(pseudobulk_df))){
    K <- Kendall(pseudobulk_df[[i]], pseudobulk_df[[j]])
    p.val<-c(p.val,K$sl[1])
    vals<-c(vals,K$tau[1])
    noms_x <- c(noms_x,i)
    noms_y <- c(noms_y,j)
    
  }
}

#colSums(matrix(vals, nrow=4))

df <- as.data.frame(x=vals)
df$cond1 <- noms_x
df$cond2 <- noms_y
df$cond3 <- noms_x==noms_y
# split the columns
df2 <-  df[df$cond3==T,]
df3 <-  df[df$cond3==F,]
df3 <- df3[match(unique(df3$vals), df3$vals),]
df3 <- df3[order(df3$cond1, df3$cond2,decreasing = T), ]
#df <- rbind(df2,df3)

M <- matrix(1, length(unique(noms_x)), length(unique(noms_x)))

transform(df3, Freq = ave(df3$f, df3$cond1, FUN = length))
df3$Freq = unname(table(df3$cond1)[df3$cond1])


# filling the matrix based on number of conditions
df3 <- df3[order(-ave(df3$Freq, df3$cond1, FUN = max), -df3$Freq), ]

df3$Freq2 = unname(table(df3$cond2)[df3$cond2])

df3 <- df3 %>% 
  arrange(desc(Freq),Freq2)

M[lower.tri(M)] <- df3[[1]] 


noms<-c(df3[1,2],unique(df3$cond2))

rownames(M)<-noms
colnames(M)<-noms

df3 <- df3 %>% 
  arrange(Freq,desc(Freq2))

M[upper.tri(M)]<- df3[[1]] 




pdf(paste(resultsdir,'Correlation_datasets_3_4.pdf',sep="/") ,width=5,height=5,paper='special')

# test for normality with qq-plots
ggqqplot(pseudobulk_df$C3, ylab = "3")
ggqqplot(pseudobulk_df$C4, ylab = "4")

#library(ConsRank)

#X<- as.matrix(pseudobulk_df)
#Tau_X(X, Y=NULL)

# data not normally distributed thus spearman correlation
spearmancortest <- cor.test(pseudobulk_df$C3, pseudobulk_df$C4,  method = "spearman",exact = F, alternative = "greater")
spearmancorval <- round(as.numeric(as.character(unname(spearmancortest['estimate'][[1]]))),2)
spearmanpval <- round(as.numeric(as.character(unname(spearmancortest['p.value'][[1]]))),2)

cor3_4 <- ggplot(pseudobulk_df, aes(x=C3, y=C4)) +
        geom_point()+
        geom_smooth(method=lm)+
        scale_color_manual(aesthetics = "color")+
        geom_text(x= max(pseudobulk_df$C3,na.rm=T)/10*1,y = max(pseudobulk_df$C4,na.rm = T)/10*9, label = paste0("Rho: ",spearmancorval), parse = TRUE,check_overlap = T)
print(cor3_4)
dev.off()

pdf(paste(resultsdir,'Correlation_datasets_5_6.pdf',sep="/") ,width=5,height=5,paper='special')

# test for normality with qq-plots
ggqqplot(pseudobulk_df$C5, ylab = "5")
ggqqplot(pseudobulk_df$C6, ylab = "6")

# data not normally distributed thus spearman correlation
spearmancortest <- cor.test(pseudobulk_df$C5, pseudobulk_df$C6,  method = "spearman", use = "complete.obs")
spearmancorval <- round(as.numeric(as.character(unname(spearmancortest['estimate'][[1]]))),2)

cor5_6 <- ggplot(pseudobulk_df, aes(x=C5, y=C6)) +
  geom_point()+
  geom_smooth(method=lm)+
  scale_color_manual(aesthetics = "color")+
  geom_text(x= max(pseudobulk_df$C5,na.rm=T)/10*1,y = max(pseudobulk_df$C6,na.rm = T)/10*9, label = paste0("Rho: ",spearmancorval), parse = TRUE,check_overlap = T)
print(cor5_6)

dev.off()

pdf(paste(resultsdir,'Correlation_datasets_7_8.pdf',sep="/") ,width=5,height=5,paper='special')

# test for normality with qq-plots
ggqqplot(pseudobulk_df$C7, ylab = "7")
ggqqplot(pseudobulk_df$C8, ylab = "8")

# data not normally distributed thus spearman correlation
spearmancortest <- cor.test(pseudobulk_df$C7, pseudobulk_df$C8,  method = "spearman", use = "complete.obs")
spearmancorval <- round(as.numeric(as.character(unname(spearmancortest['estimate'][[1]]))),2)

cor7_8 <- ggplot(pseudobulk_df, aes(x=C7, y=C8)) +
  geom_point()+
  geom_smooth(method=lm)+
  scale_color_manual(aesthetics = "color")+
  geom_text(x= max(pseudobulk_df$C7,na.rm=T)/10*1,y = max(pseudobulk_df$C8,na.rm = T)/10*9, label = paste0("Rho: ",spearmancorval), parse = TRUE,check_overlap = T)
print(cor7_8)
dev.off()

pdf(paste(resultsdir,'Correlation_datasets_7_9.pdf',sep="/") ,width=5,height=5,paper='special')

# test for normality with qq-plots
ggqqplot(pseudobulk_df$C7, ylab = "7")
ggqqplot(pseudobulk_df$C9, ylab = "9")

# data not normally distributed thus spearman correlation
spearmancortest <- cor.test(pseudobulk_df$C7, pseudobulk_df$C9,  method = "spearman", use = "complete.obs")
spearmancorval <- round(as.numeric(as.character(unname(spearmancortest['estimate'][[1]]))),2)

cor7_9 <- ggplot(pseudobulk_df, aes(x=C7, y=C9)) +
  geom_point()+
  geom_smooth(method=lm)+
  scale_color_manual(aesthetics = "color")+
  geom_text(x= max(pseudobulk_df$C7,na.rm=T)/10*1,y = max(pseudobulk_df$C9,na.rm = T)/10*9, label = paste0("Rho: ",spearmancorval), parse = TRUE,check_overlap = T)
print(cor7_9)
dev.off()

pdf(paste(resultsdir,'Correlation_datasets_8_9.pdf',sep="/") ,width=5,height=5,paper='special')

# test for normality with qq-plots
ggqqplot(pseudobulk_df$C9, ylab = "9")
ggqqplot(pseudobulk_df$C8, ylab = "8")

# data not normally distributed thus spearman correlation
spearmancortest <- cor.test(pseudobulk_df$C8, pseudobulk_df$C9,  method = "spearman", use = "complete.obs")
spearmancorval <- round(as.numeric(as.character(unname(spearmancortest['estimate'][[1]]))),2)

cor8_9 <- ggplot(pseudobulk_df, aes(x=C8, y=C9)) +
  geom_point()+
  geom_smooth(method=lm)+
  scale_color_manual(aesthetics = "color")+
  geom_text(x= max(pseudobulk_df$C8,na.rm=T)/10*1,y = max(pseudobulk_df$C9,na.rm = T)/10*9, label = paste0("Rho: ",spearmancorval), parse = TRUE,check_overlap = T)
print(cor8_9)
dev.off()

pdf(paste(resultsdir,'Correlation_datasets_11_12.pdf',sep="/") ,width=5,height=5,paper='special')

# test for normality with qq-plots
ggqqplot(pseudobulk_df$C11, ylab = "11")
ggqqplot(pseudobulk_df$C12, ylab = "12")

# data not normally distributed thus spearman correlation
spearmancortest <- cor.test(pseudobulk_df$C11, pseudobulk_df$C12,  method = "spearman", use = "complete.obs")
spearmancorval <- round(as.numeric(as.character(unname(spearmancortest['estimate'][[1]]))),2)

cor11_12 <- ggplot(pseudobulk_df, aes(x=C11, y=C12)) +
  geom_point()+
  geom_smooth(method=lm)+
  scale_color_manual(aesthetics = "color")+
  geom_text(x= max(pseudobulk_df$C11,na.rm=T)/10*1,y = max(pseudobulk_df$C12,na.rm = T)/10*9, label = paste0("Rho: ",spearmancorval), parse = TRUE,check_overlap = T)
print(cor11_12)
dev.off()

pdf(paste(resultsdir,'Correlation_datasets_14_15.pdf',sep="/") ,width=5,height=5,paper='special')

# test for normality with qq-plots
ggqqplot(pseudobulk_df$C14, ylab = "14")
ggqqplot(pseudobulk_df$C15, ylab = "15")

# data not normally distributed thus spearman correlation
spearmancortest <- cor.test(pseudobulk_df$C14, pseudobulk_df$C15,  method = "spearman", use = "complete.obs")
spearmancorval <- round(as.numeric(as.character(unname(spearmancortest['estimate'][[1]]))),2)
#spearmancorstat <- round(as.numeric(as.character(unname(spearmancortest['statistic'][[1]]))),2)
cor14_15 <- ggplot(pseudobulk_df, aes(x=C14, y=C15)) +
  geom_point()+
  geom_smooth(method=lm)+
  scale_color_manual(aesthetics = "color")+
  geom_text(x= max(pseudobulk_df$C14,na.rm=T)/10*1,y = max(pseudobulk_df$C15,na.rm = T)/10*9, label = paste0("Rho: ",spearmancorval), parse = TRUE,check_overlap = T)
print(cor14_15)

dev.off()

# Joined correlation plots
library(patchwork)
design_cor <- layout <- '
ABCD
EGH#
'

pdf(paste0(resultsdir,'/cor_wrap.pdf') ,width=12,height=6,paper='special')
wrap_plots(list(A=cor3_4,B=cor5_6,C=cor7_8,D=cor7_9,E=cor8_9,G=cor11_12,H=cor14_15),design=design_cor)
dev.off()

################################################################################
# Continue with 4 datasets integration
seur_split <- seur_obj

# Re-name the clustering of the UMAP
cluster_order <- c('LSC','LESC','LE','LE','Cj','Cj','CE',
                   'CE','CE','CSSC','SK','SK','TSK','CF',
                   'CF','EC','Ves','Mel','IC','nm-cSC','MF')

# Generate new cluster labeling:
current.cluster.ids <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)

new.cluster.ids <- c('LSC','LESC','LE','LE','Cj','Cj','CE',
                     'CE','CE','CSSC','SK','SK','TSK','CF',
                     'CF','EC','Ves','Mel','IC','nm-cSC','MF')

seur_obj$scVI_label <- plyr::mapvalues(x = as.factor(seur_obj$scVI_label), from = current.cluster.ids, to = new.cluster.ids)
seur_obj$scVI_label <- as.factor(seur_obj$scVI_label)

# Export labeling for visualization on the UMAP in Python
write.csv(seur_obj$scVI_label,paste0(resultsdir,"/labels_4datasets.csv"),row.names=F)

# Check clusters on  the first two PCs
seur_obj <- FindVariableFeatures(seur_obj)
seur_obj <- ScaleData(seur_obj)
seur_obj <- RunPCA(seur_obj, verbose = FALSE, npcs = 2)

mat <- Seurat::GetAssayData(seur_obj, assay = "RNA", slot = "scale.data")
pca <- seur_obj[["pca"]]

# Get the total variance:
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance

# Checking the first PC
Stdev(object = seur_obj[["pca"]])[1]

seur_obj@active.ident <- seur_obj$scVI_label
pdf(paste0(resultsdir,'/scVI_4datasets.pdf') ,width=6,height=5,paper='special')
pc1_viz <-2
pc2_viz <-1
y_label = "PC2"#paste0(paste0(paste0("PC",pc2_viz), ' stdev:  '),round(Stdev(object = seur_obj[["pca"]])[pc2_viz],3))
x_label = "PC1"#paste0(paste0(paste0("PC",pc1_viz), ' stdev:  '),round(Stdev(object = seur_obj[["pca"]])[pc1_viz],3))
PC_dimred <- DimPlot(seur_obj, reduction = "pca",label = F, dims = c(pc1_viz,pc2_viz))+labs(y= y_label, x = x_label)
print(PC_dimred)
dev.off()

saveRDS(seur_obj,paste0(resultsdir,"/4datasets_annotated_joined.rds"))

#Lets find marker genes for each cluster still split
seur_obj <- seur_split

# Re-name the clustering of the UMAP
cluster_order <- c('1: LSC','2: LESC','3: LE','4: LE','5: Cj','6: Cj',
                   '7: CE','8: CE','9: CE','10: CSSC','11: SK','12: SK','13: TSK',
                   '14: CF','15: CF','16: EC','17: Ves','18: Mel','19: IC','20: nm-cSC','21: MF')

# Generate new cluster labeling:
current.cluster.ids <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)

new.cluster.ids <- c('1: LSC','2: LESC','3: LE','4: LE','5: Cj','6: Cj',
                     '7: CE','8: CE','9: CE','10: CSSC','11: SK','12: SK','13: TSK',
                     '14: CF','15: CF','16: EC','17: Ves','18: Mel','19: IC','20: nm-cSC','21: MF')

seur_obj$scVI_label <- plyr::mapvalues(x = as.factor(seur_obj$scVI_label), from = current.cluster.ids, to = new.cluster.ids)
seur_obj$scVI_label <- as.factor(seur_obj$scVI_label)
seur_obj@active.ident <- seur_obj$scVI_label

cluster.markers <- FindAllMarkers(seur_obj, only.pos = TRUE)

#filter on only sigi
cluster.markers$gene_name_shorter <- str_sub(cluster.markers$gene,end = -1)
cluster.markers_fc <- cluster.markers[cluster.markers$p_val_adj < 0.01,]

cluster.markers_fc <- cluster.markers_fc[(cluster.markers_fc$avg_log2FC > 0.58) ,]

# avg_logFC...
table(cluster.markers_fc$cluster) # number of significantly expressed genes with a high log2FC (>0.58)

# Write the markers as a table
write.table(cluster.markers_fc, file = paste0(resultsdir,'/cluster_markers_all_split.csv'), sep = ',')

# Now for all joined clusters
seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220829/4datasets_annotated_joined.rds")

#Lets find marker genes for each cluster
cluster.markers <- FindAllMarkers(seur_obj, only.pos = TRUE)

# making a nice heatmap of expression in each cluster not significant yet
n_genes <- 100
heatmap.markers <- cluster.markers %>% group_by(cluster) %>% top_n(n_genes, avg_log2FC)
DoHeatmap(seur_obj, features = heatmap.markers$gene) + NoLegend()

#filter on only sigi
cluster.markers$gene_name_shorter <- str_sub(cluster.markers$gene,end = -1)
cluster.markers_fc <- cluster.markers[cluster.markers$p_val_adj < 0.01,]

cluster.markers_fc <- cluster.markers_fc[(cluster.markers_fc$avg_log2FC > 0.58) ,]
# avg_logFC...
table(cluster.markers_fc$cluster) # number of significantly expressed genes with a high log2FC (>0.58)

# GO terms
plot_list <- list()
cluster_GOs <- c()
#expressed_genes <- row.names(counts(sce_qc)[rowSums(counts(sce_qc))>20,])
expressed_genes <- rownames(cluster.markers_fc)

unique(expressed_genes)

# Write the markers as a table
write.table(cluster.markers_fc, file = paste0(resultsdir,'/cluster_markers_all.csv'), sep = ',')
write.table(heatmap.markers, file = paste0(resultsdir,'/heatmap_markers_all.csv'), sep = ',')

#saveRDS(seur_obj, file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAunannotated2.rds") # 10 PC
#saveRDS(seur_obj, file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAunannotated_4clus.rds")
#seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAunannotated2.rds") # 10 PC
#seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAunannotated_4clus.rds")

###########################################################################
# Appending the GO terms of the unannotated clusters
plot_list <- list()
cluster_GOs <- c()
#expressed_genes <- rownames(cluster.markers_fc)

cluster.markers <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220829/cluster_markers_all.csv', sep = ',', header = TRUE, row.names = 1)

df <- NULL
df_go <- NULL

for (cluster in unique(cluster.markers$cluster)){
  print(cluster)
  df_subset <- cluster.markers[cluster.markers$cluster == cluster,]
  cluster_genes <- unique(str_sub(df_subset$gene,end = -1))
  PC1_ego <- enrichGO(gene = cluster_genes,
                      #universe = expressed_genes,
                      OrgDb         = 'org.Hs.eg.db',
                      keyType       = "SYMBOL",
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      #pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)
  if(is.null(PC1_ego)){
    go_plot1 <- rectGrob(gp=gpar(col=NA))
  }else{
    go_plot1 <- barplot(PC1_ego, showCategory=10)
  }
  PC1_ego <- simplify(PC1_ego, cutoff = 0.8, by = "p.adjust", select_fun = min)
  df<-PC1_ego@result
  df$condition  <-rep((cluster),times=nrow(df))
  df_go <- rbind(df_go,df)
  plot_list <- c(plot_list, list(go_plot1))
  cluster_GOs <- append(cluster_GOs, (paste('cluster of cell type ', cluster, sep = '')))
}
dev.off()

# plotting the GO-terms
grob2 <- ggarrange(plotlist= plot_list, nrow = 21, ncol = 1, labels = cluster_GOs) #, align = 'v')
pdf(paste(resultsdir,'/all_DEGS_heatmap.pdf',sep="") ,width=10,height=40,paper='special')
print(ggarrange(grob2, ncol =1 , nrow = 18, widths= c(6, 4)))
print(grob2)
dev.off()

# Heatmap cluster marking genes:
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

heatmap.markers <- cluster.markers_fc
cluster.markers_fc$log2FC <- cluster.markers_fc$avg_log2FC
top10 <- cluster.markers_fc %>%  group_by(cluster) %>%  top_n(n = 5)
top5_genes <- DoHeatmap(slot = "data",seur_obj, features = top10$gene) + NoLegend()
#top10 <- cluster.markers_fc %>% group_by(cluster) %>% top_n(n = 10, wt = 'avg_log2FC')
#marker_top5 <- DoHeatmap(slot= "data",seur_obj, features = top10$gene) + NoLegend()

n_genes <- 100
heatmap.markers <- cluster.markers_fc %>% group_by(cluster) #%>% top_n(n_genes, avg_logFC)

cell_type_plot <- DimPlot(seur_obj, label=TRUE,pt.size = 2) + ggtitle("Timepoint")
cluster_plot <- DimPlot(seur_obj, label=TRUE, group.by = 'costum_clustering', pt.size = 2) + ggtitle("Louvain Clustering")
grob_clustering = ggarrange(plotlist = list(cell_type_plot, cluster_plot) , ncol =1) #labels = cluster_GOs) #align = 'v')
nclust<- length(unique(as.numeric(cluster.markers_fc$cluster)))
RGB_colours_ggplot <- as.list(as.character(gg_color_hue(nclust)))
names <- as.numeric(unique(heatmap.markers$cluster))
col_fun = colorRamp2(names, RGB_colours_ggplot)#make sure the rows correspond to ggplot colour mapping

seur_obj2 <- seur_obj
seur_obj@assays$SCT <- NULL
seur_obj@assays$integrated <- NULL

mat <- as.matrix(seur_obj@assays$RNA[heatmap.markers$gene,]) # removed scale data
mat <- rbind(mat,seur_obj@active.ident)
mat <- mat[,order(mat[nrow(mat),])]

column_ha <- HeatmapAnnotation(cluster = mat[nrow(mat),], col =list(cluster = col_fun))
breaks <- mat[nrow(mat),]
mat <- mat[-nrow(mat),]
row_ha = rowAnnotation(adj_p_vallue = heatmap.markers$p_val_adj)
clust_heatmap <- grid.grabExpr(draw(Heatmap(mat, column_split = breaks,row_split = as.numeric(heatmap.markers$cluster), cluster_columns = F, cluster_rows = F,show_row_names = F,show_column_names = F,row_names_gp = gpar(fontsize = 6), top_annotation = column_ha, left_annotation = row_ha, row_names_rot = -40)))

###########################################################################
# Plotting the GO terms of the un-annotated clusters and the heatmap

#expressed_genes <- rownames(cluster.markers_fc)

pdf(paste0(resultsdir,'/GO_terms_clusters_epi_sub.pdf') ,width=12,height=120,paper='special')
#cowplot::plot_grid(plotlist =grob2, ncol = 1, nrow = 1)#, rel_heights = c(1, 4,2), rel_widths = c(3,2,2))
print(grob2)
dev.off()

pdf(paste0(resultsdir,'/complex_heatmap_epi_sub.pdf') ,width=40,height=24,paper='special')
cowplot::plot_grid(top5_genes,grob_clustering, ncol = 3, nrow = 1, rel_heights = c(1, 4,2), rel_widths = c(3,2,2))
dev.off()

# Write the markers as a table if not done already
#write.table(cluster.markers, file = paste0(resultsdir,'/cluster_markers.csv'), sep = ',')


cluster.markers <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220829/cluster_markers_all.csv', sep = ',', header = TRUE, row.names = 1)
df <- NULL
df_go <- NULL
plot_list <- list()
cluster_GOs <- c()

# Subselecting the top n cluster markers:
#n_genes <- 150
#cluster.markers <- cluster.markers %>% group_by(cluster) %>% top_n(n_genes, avg_log2FC)

for (cluster in unique(cluster.markers$cluster)){
  print(cluster)
  df_subset <- cluster.markers[cluster.markers$cluster == cluster,]
  cluster_genes <- unique(str_sub(df_subset$gene,end = -1))
  PC1_ego <- enrichGO(gene = cluster_genes,
                      #universe = expressed_genes,
                      OrgDb         = 'org.Hs.eg.db',
                      keyType       = "SYMBOL",
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      #pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)
  if(is.null(PC1_ego)){
    go_plot1 <- rectGrob(gp=gpar(col=NA))
  }else{
    go_plot1 <- barplot(PC1_ego, showCategory=10)
  }
  PC1_ego <- simplify(PC1_ego, cutoff = 0.5, by = "p.adjust", select_fun = min)
  df<-PC1_ego@result
  df$condition  <-rep((cluster),times=nrow(df))
  df_go <- rbind(df_go,df)
  plot_list <- c(plot_list, list(go_plot1))
  cluster_GOs <- append(cluster_GOs, (paste('cluster of cell type ', cluster, sep = '')))
}
dev.off()

# Go term dotplot
write.table(df_go, file = paste0(resultsdir,'/go_cluster_markers_0.5_all.csv'), sep = ',',row.names = F)

#### TRY GENE RATIO >0.1 or 0.11 -> then it will make more sense

#df_go$gene.ratio <- foo$num$new
vec <- c()
for (i in 1:length(df_go$GeneRatio)){
  vec <- c(vec,strsplit(df_go$GeneRatio,split = "/")[[i]][2])
}
df_go$total <- as.numeric(vec)
df_go$GR <- round(df_go$Count/df_go$total,digits = 1)

vec <- c()
for (i in 1:length(df_go$BgRatio)){
  vec <- c(vec,strsplit(df_go$BgRatio,split = "/")[[i]][1])
}
df_go$BGcount <- as.numeric(vec)

df_go2 <- df_go[df_go$GR>=0.1 &df_go$BGcount>=200 &df_go$p.adjust<0.01,]
df_go2 <- df_go2 %>% group_by(condition)  %>% slice_head(n=3)#%>% top_n(n = 10, wt = 'Count')

pdf(paste0(resultsdir,'/go_dot_epi_sub_sub.pdf') ,width=14,height=7,paper='special')

vec_cond<-c("LSC","LESC","LE","CE","Cj","Mel","IC","CSSC","SK","TSK","CF","MF","EC","Ves","nm-cSC")
vec_cond <- rev(vec_cond)

df_go2<-df_go2[order(match(as.factor(df_go2$condition), vec_cond)), , drop = FALSE]

df_go2 <- df_go2[as.factor(str_order(df_go2$Description,decreasing = T)),]

ggplot(data = df_go2, aes(x = factor(condition,level=vec_cond), y = Description, 
                          color = `p.adjust`, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GO enrichment analysis") + coord_flip()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

pdf(paste0(resultsdir,'/go_dot_epi.pdf') ,width=10,height=5,paper='special')

vec_cond<-c("LSC","LESC","LE","CE","Cj")
df_go2 <- df_go[df_go$GR>=0.1 &df_go$BGcount>=200&df_go$p.adjust<0.01,] #df_go$Count>12&
df_go2 <- df_go2[df_go2$condition %in% vec_cond,]

df_go2 <- df_go2 %>% group_by(condition) %>% arrange(desc(Count),.by_group = T)

df_go2 <- df_go2 %>% group_by(condition)  %>% slice_head(n=3)#%>% top_n(n = 10, wt = 'Count')
# df_go2<-df_go2[order(match(as.factor(df_go2$condition), vec_cond)), , drop = FALSE]
# 
# df_go2 <- df_go2[as.factor(str_order(df_go2$Description,decreasing = T)),]
#mid<-1.0e−06
p1 <-ggplot(data = df_go2, aes(x = factor(condition,level=vec_cond), y = Description, 
                          color = -log10(`p.adjust`), size = Count, show.legend = FALSE)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GO enrichment analysis")+
  theme(axis.text.x = element_text(angle = 45, hjust=1),legend.position="bottom")#+
  #scale_size_area()
p1
dev.off()

#df_go$gene.ratio <- foo$num$new
df_go3 <- df_go[df_go$Count>5 &df_go$p.adjust<0.05,]
#df_go3 <- df_go3 %>% group_by(condition)  %>% slice_head(n=3)#%>% top_n(n = 10, wt = 'Count')

pdf(paste0(resultsdir,'/go_dot_non_epi.pdf') ,width=10,height=5,paper='special')

vec_cond2<-c("CSSC","SK","TSK","CF","EC","Ves","Mel","IC","MF","nm-cSC")
#df_go2 <- df_go[df_go$Count>10 &df_go$p.adjust<0.0005,]

df_go3 <- df_go[df_go$GR>=0.1 &df_go$BGcount>=200&df_go$p.adjust<0.01,] #df_go$Count>12&
df_go3<- df_go3[df_go3$condition %in% vec_cond2,]

df_go3 <- df_go3 %>% group_by(condition) %>% arrange(desc(Count),.by_group = T)

df_go3 <- df_go3 %>% group_by(condition)  %>% slice_head(n=3)#%>% top_n(n = 10, wt = 'Count')


p2<-ggplot(data = df_go3, aes(x = factor(condition,level=vec_cond2), y = Description, 
                          color = -log10(`p.adjust`), size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GO enrichment analysis")+
  theme(axis.text.x = element_text(angle = 45, hjust=1),legend.position="bottom")#+
  #scale_size_area()
p2
dev.off()


# Joined go plots
library(patchwork)
design_go <- "AB"

pdf(paste0(resultsdir,'/go_wrap.pdf') ,width=10,height=5,paper='special')
wrap_plots(list(A=p1,B=p2),design=design_go)
dev.off()
write.table(df_go, file = paste0(resultsdir,'/go_cluster_markers.csv'), sep = ',',row.names = F)

###############################################################################
#Pahway analysis with progeny
#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("progeny")

library(progeny)
seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220829/4datasets_annotated_joined.rds")

seur_obj <- progeny(seur_obj, scale=TRUE, organism="Human", top=500, perm=1, 
                return_assay = TRUE)

## We create a data frame with the specification of the cells that belong to 
## each cluster to match with the Progeny scores.
CellsClusters <- data.frame(Cell = names(Idents(seur_obj)), 
                            CellType = as.character(Idents(seur_obj)),
                            stringsAsFactors = FALSE)

## We can now directly apply Seurat functions in our Progeny scores. 
## For instance, we scale the pathway activity scores. 
seur_obj <- Seurat::ScaleData(seur_obj, assay = "progeny") 

## We transform Progeny scores into a data frame to better handling the results
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(seur_obj, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 

## We match Progeny scores with the cell clusters.
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

mydata_full <- progeny_scores_df[,c("Pathway","Activity","CellType")]
for (i in unique(seur_obj$scVI_label)){
  
  pdf(paste0(resultsdir,"/",i,'Progeny.pdf') ,width=8,height=6,paper='special')
  mydata <- mydata_full[mydata_full$CellType==i,]
  mydata <- mydata[,c("Pathway","Activity")]
  names(mydata) <- c("group", "value")

# function for computing mean, DS, max and min values
min.mean.sd.max <- function(x) {
  r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

p1 <- ggplot(aes(y = value, x = factor(group)), data = mydata)
p1 <- p1 + stat_summary(fun.data = min.mean.sd.max, geom = "boxplot") #+ geom_jitter(position=position_jitter(width=.2), size=3) + ggtitle("Boxplot con media, 95%CI, valore min. e max.") + xlab("Gruppi") + ylab("Valori")
print(p1)
dev.off()
}

pdf(paste0(resultsdir,"/heatmap_progeny_raw.pdf") ,width=5,height=4,paper='special')
mat <- t(summarized_progeny_scores_df)
mat <- mat[,c("LSC","LESC","LE","Cj","CE","CSSC","SK","TSK","CF","MF","EC","Ves","Mel","IC","nm-cSC")]

f1 = colorRamp2(c(-1.5, 0, 1.5), c("blue", "#EEEEEE", "red"), space = "RGB")

ht1 = Heatmap(mat, col = f1, cluster_columns = F,cluster_rows = T, name = "progeny_raw_score")
print(ht1)
dev.off()

###############################################################################
# Generate barplot overview first (added to the matrix later)
newdf <- NULL

for (xx in unique(seur_obj$scVI_label)) {
  print(xx)
  for (j in unique(seur_obj$Condition)){
    val<-sum(seur_obj$scVI_label==xx & seur_obj$Condition==j)
    #col<- dfcol[rownames(dfcol)%in% xx,]
    df <- data.frame(values=val,Dataset=j,cell=xx)
    newdf <- rbind(newdf, df)
  }
}

pdf(paste(resultsdir,'Cell_dist_datasets_epi_full.pdf',sep="/") ,width=10,height=8,paper='special')
# Stacked barchart showing the composition
print(ggplot(newdf, aes(fill=Dataset, y=values, x=cell)) +
        geom_bar(position="stack", stat="identity") + labs(x="Cell type",y="Number of cells"))
dev.off()

# Re-order that the colors make sense for the bar_plot
newdf <- newdf[order(newdf$cell),]
#dfcol2 <- dfcol[order(rownames(dfcol)),]

pdf(paste(resultsdir,'Cell_dist_datasets_stacked_epi_full.pdf',sep="/") ,width=5,height=5,paper='special')
print(ggplot(newdf, aes(fill=cell, y=values, x=Dataset)) +
        geom_bar(position="fill", stat="identity", width = 0.3) + labs(x="Dataset",y="Cell proportions"))#+scale_fill_manual(values = dfcol2)
dev.off()
################################################################################
# Generate a z-score table to determine relative expression in-between cell 
# populations
# Splitting the pseudobulk datasets unbiased into two replicates for Z-score calculation and over pseudotime
seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220229/4datasets_annotated_joined.rds")
sce_qc <- SingleCellExperiment(assays = list(counts = GetAssayData(object = seur_obj, slot = "data")))

pseudobulk_df <- as.data.frame(row.names(counts(sce_qc)))

sce_qc <- SingleCellExperiment(assays = list(counts = GetAssayData(object = seur_obj, slot = "data")))

sce_qc$sample <- seur_obj$scVI_label

sce_qc$sample <- factor(sce_qc$sample)
sce_qc$sample

sce_qc$rep <- seur_obj@meta.data$orig.ident
sce_qc$rep

unique(seur_obj$scVI_label)

for (sample in unique(sce_qc$sample)){
  pseudobulk_df[ , paste0(sample)] <- rowSums(counts(sce_qc)[,sce_qc$sample == sample])
}

row.names(pseudobulk_df) <- pseudobulk_df$`row.names(counts(sce_qc))`
pseudobulk_df$`row.names(counts(sce_qc))` <- NULL

write.table(data.frame("ID"=rownames(pseudobulk_df),pseudobulk_df), file = paste0(resultsdir,'/pseudobulkdf_full.tsv'), sep = '\t',quote = F, row.names = F)

pseudobulk_df <- NULL
pseudobulk_df <- as.data.frame(row.names(counts(sce_qc)))

# putting values of artificial rep2 and rep1 into the pseudobulk dataframe columns
for (sample in unique(sce_qc$sample)){
      sce_qc$reps <- "rep1"
      
      # use rep2 for Catala and one subject from Li
      sce_qc$reps[sce_qc$rep == "GautamSRR11470710"| sce_qc$rep =="GautamSRR14742510"| sce_qc$rep =="GautamSRR14742511"] <- "rep2"
      sce_qc$reps[sce_qc$rep == "CatalaGSM5651509" | sce_qc$rep == "CatalaGSM5651511" | sce_qc$rep == "CatalaGSM5651513" | sce_qc$rep == "CatalaGSM56511515"| sce_qc$rep == "CatalaGSM5651117"| sce_qc$rep == "CatalaGSM5651119"] <- "rep2"
      sce_qc$reps[sce_qc$rep == "CatalaGSM5651510" | sce_qc$rep == "CatalaGSM5651512" | sce_qc$rep == "CatalaGSM5651514" | sce_qc$rep == "CatalaGSM56511516"| sce_qc$rep == "CatalaGSM5651118"|sce_qc$rep == "CatalaGSM5651120"] <- "rep2"
      
      pseudobulk_df[ , paste0(sample,"_1")] <- rowSums(counts(sce_qc)[,sce_qc$sample == sample & sce_qc$reps == "rep1"])
      pseudobulk_df[ , paste0(sample,"_2")] <- rowSums(counts(sce_qc)[,sce_qc$sample == sample & sce_qc$reps == "rep2"])
}

row.names(pseudobulk_df) <- pseudobulk_df$`row.names(counts(sce_qc))`
pseudobulk_df$`row.names(counts(sce_qc))` <- NULL

write.table(data.frame("ID"=rownames(pseudobulk_df),pseudobulk_df), file = paste0(resultsdir,'/pseudobulk_reps_DE_datasets_markers_split2.tsv'), sep = '\t',quote = F, row.names = F)

pseudobulk_df<- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220830/pseudobulk_reps_DE_datasets_markers_split2.tsv',header = T,row.names = 1)
pseudobulk_df <- pseudobulk_df[,1:28]

# For the deseq2 matrix vs ESC
lakocountfile <- pseudobulk_df

lakovst <- lakocountfile

# Generate coldata dataframe
coldata <- NULL

# conditions
j <- unlist(strsplit(colnames(pseudobulk_df), split = "_"))
j <- j[seq(1,length(j),2)]

# Change _1 & _2 to condition 1 and _3 & _4 to condition 2
c <- c("1","2")
d <- paste(j,c,sep="_")

# reps
cols <- colnames(pseudobulk_df)

# type
xx <- rep("paired-end",length(cols))

coldata <- data.frame(cols=cols,condition=d,condition2=j,type=xx)

rownames(coldata)<- coldata$cols
coldata <- coldata[,c("condition","condition2","type")]

# setting the correct columns for the coldata
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

# checking the data
head(lakovst,2)

coldata

# setting the rowdata in coldata similar to the coldata in the count matrix
rownames(coldata) <- sub("fb", "", rownames(coldata))
all(rownames(coldata) %in% colnames(lakovst))

all(rownames(coldata) == colnames(lakovst))

lakovst <- lakovst[, rownames(coldata)]
all(rownames(coldata) == colnames(lakovst))

dds <- DESeqDataSetFromMatrix(countData = round(lakovst),
                              colData = coldata,
                              design = ~ condition)
dds
###################################################################
# Complex heatmap Z-score table

# transpose and scale the matrix per row Z-score per gene per sample is scaled upon
vsd <- assay(vst(dds,blind = T))
Z <- t(scale(t(vsd)))
Z
Z_score <- as.data.frame(Z)

# joining the columns on the means sequential
n <- 2
Z_joined <- t(rowMeans(t(Z), as.integer(gl(ncol(Z_score), n, ncol(Z_score)))) / n)

#  generate list of factors
vec <- rownames(coldata)
x <- split(vec, ceiling(seq_along(vec)/2))

scoretable <- as.data.frame(do.call(cbind, lapply(x, function(i) rowMeans(Z_score[, i]))), row.names = rownames(lakovst))
colnames(scoretable) <- unique(coldata$condition2)
write.table(data.frame("ID"=rownames(scoretable),scoretable), file = paste0(resultsdir,'/Zscoretable_markers_split2.tsv'), sep = '\t',quote = F, row.names = F)
scoretable1 <- scoretable

################################################################################
# Add the markers top 3 for each cell fate
markers<- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220829/cluster_markers_all.csv', sep = ',', header = TRUE, row.names = 1)
data_new2 <- markers %>%                                      # Top N highest values by group
  arrange(desc(avg_log2FC)) %>% 
  group_by(cluster) %>%
  slice(1:3) 
data_new2    

markers_vec <- unique(data_new2$gene_name_shorter)

################################################################################
# Generate the complex heatmap of the z-score
scoretable1 <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220830/Zscoretable_markers_split2.tsv", sep="\t",row.names=1)

# nm.cSC to nm-cSC in the colnames
colnames(scoretable1) <- c("Cj","CSSC","LESC","CE","LSC","Mel","CF","Ves","MF","SK","IC","nm-cSC","LE","EC")

pdf(paste(resultsdir,'heatmap_z-score_table_4_dataset.pdf',sep="/") ,width=5,height=11,paper='special')
f1 = colorRamp2(c(-2,0, 2), c("blue","white", "red"), space = "RGB")

newdf <- NULL

for (xx in unique(seur_obj$scVI_label)) {
  print(xx)
  for (j in unique(seur_obj$Condition)){
    val<-sum(seur_obj$scVI_label==xx & seur_obj$Condition==j)
    #col<- dfcol[rownames(dfcol)%in% xx,]
    df <- data.frame(values=val,Dataset=j,cell=xx)
    newdf <- rbind(newdf, df)
  }
}

# Order the rows and columns to make it more clear
cont_df <- as.data.frame.matrix(scoretable1)

# Re-order and leave out populations where only one single dataset is present
cont_df <- cont_df[,colnames(cont_df)%in%c("LSC","LESC","LE","CE","Cj","Mel","CSSC","CF","IC","Ves","nm-cSC")]
cont_df <- cont_df[c("LSC","LESC","LE","CE","Cj","Mel","CSSC","CF","IC","Ves","nm-cSC")]

cont_df<-na.omit(cont_df)

matz <- as.matrix(cont_df)
rownames(matz)<-rownames(cont_df)
vec3 <- markers_vec

# Defining multiple vectors to show relative contributions of datasets
newdf <- newdf[newdf$cell%in%colnames(cont_df),]

m3 <- newdf[newdf$Dataset == "Co",]$values
co <- m3
m3 <- newdf[newdf$Dataset == "Ca",]$values
ca <- m3
m3 <- newdf[newdf$Dataset == "Ga",]$values
ga <- m3
m3 <- newdf[newdf$Dataset == "Li",]$values
li <- m3

# Creating matrix
m <- cbind(co, ca, ga, li)
rownames(m) <- unique(newdf$cell)

m2 <- m/rowSums(m)

m2 <- m2[c("LSC","LESC","LE","CE","Cj","Mel","CSSC","CF","IC","Ves","nm-cSC"),,drop=FALSE]

anno = anno_barplot(m2, gp = gpar(fill = 2:5), bar_width = 0.5, height = unit(1, "cm"))

ha = HeatmapAnnotation(Composition=anno,
                       show_legend = c("dataset" = TRUE))

# Barplots don't generate legends so need to do it manually
lgd_list = list(
  Legend(labels = c("Collin", "Catala","Gautam","Li"), title = "Dataset annotation", 
         legend_gp = gpar(fill = 2:5)))

# Simplify the dendrogram
hc = hclust(dist(matz))
group = cutree(hc, k = 10)

ht1 = Heatmap(matz, col = f1, cluster_columns = F, name = "Relative expression",
              top_annotation  = ha, cluster_rows = cluster_within_group(t(matz), group), 
              row_split = 10, border = TRUE, row_names_gp = gpar(fontsize = 10))

# Annotate rownames
vec1 <- which(rownames(matz) %in% vec3, arr.ind = T)
vec3 <- rownames(matz[rownames(matz) %in% vec3 ,])

htlist1 <- ht1 + rowAnnotation(link = anno_mark(at =  vec1,labels = vec3))
draw(htlist1, 
     heatmap_legend_list = lgd_list)
dev.off()

################################################################################
# Gene set enrichment analysis with GSEA and epifactors database
epig_factors <- read.table(file = "genes.txt",sep = '\t',header = T)
library(plyr)

epig_factors <- epig_factors[epig_factors$Function!="-",]

 table <- NULL
 for(i in unique(epig_factors$Function)){
   print(i)
   epi_sel <- epig_factors[epig_factors$Function==i,]
   col <- epi_sel$HGNC.approved.symbol
   print(col)
   table[[i]] <- col
 }

max.length <- max(sapply(table, length))

l <- table
## Add NA values to list elements
l <- lapply(l, function(v) { c(v, rep(NA, max.length-length(v)))})
## Rbind
table2 <-do.call(rbind, l)

pseudobulk_df<- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220830/pseudobulk_reps_DE_datasets_markers_split2.tsv',header = T,row.names = 1)
pseudobulk_df <- pseudobulk_df[,1:28]

# For the deseq2 matrix vs ESC
countfile <- pseudobulk_df

lakovst <- countfile

# Generate coldata dataframe
coldata <- NULL

# conditions
j <- unlist(strsplit(colnames(pseudobulk_df), split = "_"))
j <- j[seq(1,length(j),2)]

# Change _1 & _2 to condition 1 and _3 & _4 to condition 2
c <- c("1","2")
d <- paste(j,c,sep="_")

# reps
cols <- colnames(pseudobulk_df)

# type
xx <- rep("paired-end",length(cols))

coldata <- data.frame(cols=cols,condition=d,condition2=j,type=xx)

rownames(coldata)<- coldata$cols
coldata <- coldata[,c("condition","condition2","type")]

# setting the correct columns for the coldata
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

# checking the data
head(lakovst,2)

coldata

# setting the rowdata in coldata similar to the coldata in the count matrix
rownames(coldata) <- sub("fb", "", rownames(coldata))
all(rownames(coldata) %in% colnames(lakovst))

all(rownames(coldata) == colnames(lakovst))

lakovst <- lakovst[, rownames(coldata)]
all(rownames(coldata) == colnames(lakovst))

# importing the pathways table
pathways.hallmark <- gmtPathways("msigdb.v7.2.symbols.gmt")

# show a few lines from the pathways file
head(pathways.hallmark)

pathways.hallmark <- pathways.hallmark[1:1000]

pathways.epig <- table

pathways <- pathways.epig

df_gsea <- NULL
LESC_gsea <- NULL
LE_gsea <- NULL
CSSC_gsea <- NULL
nm.cSC_gsea <- NULL
res_up_LESC <- NULL
res_down_LESC <- NULL
ranks_LESC <- NULL
res_up_LE <- NULL
res_down_LE <- NULL
ranks_LE <- NULL
res_up_CSSC <- NULL
res_down_CSSC <- NULL
ranks_CSSC <- NULL
res_up_nm.cSC <- NULL
res_down_nm.cSC <- NULL
ranks_nm.cSC <- NULL

for(i in unique(coldata$condition2)){
  print(i)
  coldata2 <- coldata
  coldata2$condition2
  
  coldata2$condition2[which(coldata2$condition2!=i)] <- "others"
  
  dds <- DESeqDataSetFromMatrix(countData = round(lakovst),
                                colData = coldata2,
                                design = ~ condition2)
  dds
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  rm(keep)

  #________________DE_analysis_____________#
  dds <- DESeq(dds) #This would take some time
  res <- results(dds, alpha=0.05)
  summary(res)
  
  #_________________GSEA___________________#
  # Steps toward doing gene set enrichment analysis (GSEA):
  
  # 1- obtaining stats for ranking genes in your experiment,
  # 2- creating a named vector out of the DESeq2 result
  # 3- Obtaining a gene set from mysigbd
  # 4- doing analysis
  
  
  # already we performed DESeq2 analysis and have statistics for working on it
  res$SYMBOL <- rownames(res)
  
  # creating  a named vector [ranked genes]
  ranks <- res$stat
  names(ranks) <- res$SYMBOL
  
  library(snow)
  SnowParam(workers = 1)
  
  #Running fgsea algorithm:
  fgseaRes <- fgseaMultilevel(pathways=pathways, stats=ranks)
  
  # Tidy the results:
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES)) # order by normalized enrichment score (NES)
  
  # add gene_ratio
  vec <- lengths(pathways)
  
  df <- data.frame(lapply(vec, type.convert), stringsAsFactors=FALSE)
  df2 <- as.data.frame(t(df))
  
  names_pathway <- as.character(names(pathways))
  names(names_pathway)<-NULL
  
  df2$pathway <- names_pathway
  names(df2) <- c("total","pathway")
  
  fgseaResTidy_ratio <- merge(fgseaResTidy, df2, by = "pathway")
  fgseaResTidy_ratio$ratio <- fgseaResTidy_ratio$size/fgseaResTidy_ratio$total
  
  pathways.hallmark[2]
  
  fgseaResTidy_ratio$condition  <-rep((i),times=nrow(fgseaResTidy_ratio))
  res_up <- res[res$log2FoldChange>0,]
  res_down <- res[res$log2FoldChange<0,]
  if (i=="LESC"){
    LESC_gsea <- fgseaResTidy
    res_up_LESC <- res_up
    res_down_LESC <- res_down
    ranks_LESC <- ranks} else if ((i=="LE")){
      LE_gsea <- fgseaResTidy
      res_up_LE <- res_up
      res_down_LE <- res_down
      ranks_LE <- ranks} else if (i=="CSSC"){
        CSSC_gsea <- fgseaResTidy
        res_up_CSSC <- res_up
        res_down_CSSC <- res_down
        ranks_CSSC <- ranks} else if (i=="nm.cSC"){
          nm.cSC_gsea <- fgseaResTidy
          res_up_nm.cSC <- res_up
          res_down_nm.cSC <- res_down
          ranks_nm.cSC <- ranks}
  df_gsea <- rbind(df_gsea,fgseaResTidy_ratio)
}

# Epigen term dotplot
str(df_gsea)
df_gsea$leadingEdge <- as.character(df_gsea$leadingEdge)
write.table(df_gsea, file = paste0(resultsdir,'/gsea_DEG_one_v_all.csv'), sep = ',',row.names = T)

df_gsea2 <- df_gsea[df_gsea$size>2 &df_gsea$padj<0.05,]
df_gsea2 <- df_gsea2 %>% group_by(condition)  %>% slice_head(n=10)#%>% top_n(n = 10, wt = 'Count')

vec_cond<-c("LESC","LE","CSSC","nm.cSC")
vec_cond <- rev(vec_cond)

df_gsea2<-df_gsea2[order(match(as.factor(df_gsea2$condition), vec_cond)), , drop = FALSE]
df_gsea2 <- df_gsea2[!is.na(df_gsea2$pathway),]
df_gsea2 <- df_gsea2[df_gsea2$condition%in%vec_cond,]

pdf(paste0(resultsdir,'/gsea_dot.pdf') ,width=6,height=3,paper='special')
ggplot(data = df_gsea2, aes(x = factor(condition,level=rev(vec_cond)), y = pathway, 
                            color = `padj`, size = size)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("fgsea enrichment analysis") #+ coord_flip()+
#theme(axis.text.x = element_text(angle = 270, hjust=0))
dev.off()

df_gsea3 <- df_gsea2[,c("pathway","leadingEdge","condition")]

#########################################################################
# Look at the up and down regulation of associated factors

###
# LESC 
# Histone modification write
ranks <- ranks_LESC
pdf(paste0(resultsdir,'/gsea_hmodwr_LESC_es.pdf') ,width=7,height=4,paper='special')
plotGseaTable(pathways["Histone modification write"], ranks, LESC_gsea, 
              gseaParam=0.5,)
dev.off()

LESC_up <- res_up_LESC$SYMBOL[which(res_up_LESC$SYMBOL %in% unlist(unname(pathways["Histone modification write"])))]
LESC_down <-res_down_LESC$SYMBOL[which(res_down_LESC$SYMBOL %in% unlist(unname(pathways["Histone modification write"])))]
LESC_genes <- c(LESC_up,LESC_down)

LESC_expr_up <-data.frame(res_up_LESC@listData$log2FoldChange,res_up_LESC@listData$padj,row.names = res_up_LESC@listData$SYMBOL)
colnames(LESC_expr_up) <- c("log2FC","padj")

LESC_expr_down <-data.frame(res_down_LESC@listData$log2FoldChange,res_down_LESC@listData$padj,row.names = res_down_LESC@listData$SYMBOL)
colnames(LESC_expr_down) <- c("log2FC","padj")

LESC_df<-rbind(LESC_expr_up,LESC_expr_down)

LESC_df <- LESC_df[rownames(LESC_df)%in%LESC_genes,]

PcG_up <- epig_factors[epig_factors$HGNC.approved.symbol%in%LESC_up,c("HGNC.approved.symbol","Protein.complex","Target.entity","Product")]
PcG_down <- epig_factors[epig_factors$HGNC.approved.symbol%in%LESC_down,c("HGNC.approved.symbol","Protein.complex","Target.entity","Product")]

###
# LE
term <- "Histone modification read"
ranks <- ranks_LE
pdf(paste0(resultsdir,'/gsea_H_READ_LE_es.pdf') ,width=7,height=4,paper='special')
plotGseaTable(pathways[term], ranks, LE_gsea, 
              gseaParam=0.5,)
dev.off()

LE_up <-res_up_LE$SYMBOL[which(res_up_LE$SYMBOL %in% unlist(unname(pathways[term])))]
LE_down <-res_down_LE$SYMBOL[which(res_down_LE$SYMBOL %in% unlist(unname(pathways[term])))]

LE_up <- res_up_LE$SYMBOL[which(res_up_LE$SYMBOL %in% unlist(unname(pathways["Histone modification read"])))]
LE_down <-res_down_LE$SYMBOL[which(res_down_LE$SYMBOL %in% unlist(unname(pathways["Histone modification read"])))]
LE_genes <- c(LE_up,LE_down)

LE_expr_up <-data.frame(res_up_LE@listData$log2FoldChange,res_up_LE@listData$padj,row.names = res_up_LE@listData$SYMBOL)
colnames(LE_expr_up) <- c("log2FC","padj")

LE_expr_down <-data.frame(res_down_LE@listData$log2FoldChange,res_down_LE@listData$padj,row.names = res_down_LE@listData$SYMBOL)
colnames(LE_expr_down) <- c("log2FC","padj")

LE_df<-rbind(LE_expr_up,LE_expr_down)

LE_df <- LE_df[rownames(LE_df)%in%LE_genes,]
LE_df_read <- LE_df

term <- "Histone modification erase cofactor"
pdf(paste0(resultsdir,'/gsea_H_ERCO_LE_es.pdf') ,width=7,height=4,paper='special')
plotGseaTable(pathways[term], ranks, LE_gsea, 
              gseaParam=0.5,)
dev.off()

LE_up <-res_up_LE$SYMBOL[which(res_up_LE$SYMBOL %in% unlist(unname(pathways[term])))]
LE_down <-res_down_LE$SYMBOL[which(res_down_LE$SYMBOL %in% unlist(unname(pathways[term])))]

LE_up <- res_up_LE$SYMBOL[which(res_up_LE$SYMBOL %in% unlist(unname(pathways["Histone modification erase cofactor"])))]
LE_down <-res_down_LE$SYMBOL[which(res_down_LE$SYMBOL %in% unlist(unname(pathways["Histone modification erase cofactor"])))]
LE_genes <- c(LE_up,LE_down)

LE_df<-rbind(LE_expr_up,LE_expr_down)
LE_df <- LE_df[rownames(LE_df)%in%LE_genes,]
LE_genes[LE_genes%in%markers$gene]

ERCO_LE_up <- epig_factors[epig_factors$HGNC.approved.symbol%in%LE_up,c("HGNC.approved.symbol","Protein.complex","Target.entity","Product")]
ERCO_LE_down <- epig_factors[epig_factors$HGNC.approved.symbol%in%LE_down,c("HGNC.approved.symbol","Protein.complex","Target.entity","Product")]

###
# CSSC
ranks <- ranks_CSSC
term <- "DNA modification, RNA modification"
pdf(paste0(resultsdir,'/gsea_DNA_RNA_CSSC_es.pdf') ,width=7,height=4,paper='special')
plotGseaTable(pathways[term], ranks, CSSC_gsea, 
              gseaParam=0.5,)
dev.off()

CSSC_up <-res_up_CSSC$SYMBOL[which(res_up_CSSC$SYMBOL %in% unlist(unname(pathways[term])))]
CSSC_down <-res_down_CSSC$SYMBOL[which(res_down_CSSC$SYMBOL %in% unlist(unname(pathways[term])))]

CSSC_up <- res_up_CSSC$SYMBOL[which(res_up_CSSC$SYMBOL %in% unlist(unname(pathways["DNA modification, RNA modification"])))]
CSSC_down <-res_down_CSSC$SYMBOL[which(res_down_CSSC$SYMBOL %in% unlist(unname(pathways["DNA modification, RNA modification"])))]
CSSC_genes <- c(CSSC_up,CSSC_down)

CSSC_expr_up <-data.frame(res_up_CSSC@listData$log2FoldChange,res_up_CSSC@listData$padj,row.names = res_up_CSSC@listData$SYMBOL)
colnames(CSSC_expr_up) <- c("log2FC","padj")

CSSC_expr_down <-data.frame(res_down_CSSC@listData$log2FoldChange,res_down_CSSC@listData$padj,row.names = res_down_CSSC@listData$SYMBOL)
colnames(CSSC_expr_down) <- c("log2FC","padj")

CSSC_df<-rbind(CSSC_expr_up,CSSC_expr_down)


CSSC_df <- CSSC_df[rownames(CSSC_df)%in%CSSC_genes,]
CSSC_genes[CSSC_genes%in%markers$gene]

DR_CSSC_up <- epig_factors[epig_factors$HGNC.approved.symbol%in%CSSC_up,c("HGNC.approved.symbol","Protein.complex","Target.entity","Product")]
DR_CSSC_down <- epig_factors[epig_factors$HGNC.approved.symbol%in%CSSC_down,c("HGNC.approved.symbol","Protein.complex","Target.entity","Product")]

###
# nm.cSC
ranks <- ranks_nm.cSC
term <- "Scaffold protein, RNA modification"
pdf(paste0(resultsdir,'/gsea_DNA_RNA_nm.cSC_es.pdf') ,width=7,height=4,paper='special')
plotGseaTable(pathways[term], ranks, nm.cSC_gsea, 
              gseaParam=0.5,)
dev.off()

nm.cSC_up <-res_up_nm.cSC$SYMBOL[which(res_up_nm.cSC$SYMBOL %in% unlist(unname(pathways[term])))]
nm.cSC_down <-res_down_nm.cSC$SYMBOL[which(res_down_nm.cSC$SYMBOL %in% unlist(unname(pathways[term])))]

nm.cSC_up <- res_up_nm.cSC$SYMBOL[which(res_up_nm.cSC$SYMBOL %in% unlist(unname(pathways["Scaffold protein, RNA modification"])))]
nm.cSC_down <-res_down_nm.cSC$SYMBOL[which(res_down_nm.cSC$SYMBOL %in% unlist(unname(pathways["Scaffold protein, RNA modification"])))]
nm.cSC_genes <- c(nm.cSC_up,nm.cSC_down)

nm.cSC_expr_up <-data.frame(res_up_nm.cSC@listData$log2FoldChange,res_up_nm.cSC@listData$padj,row.names = res_up_nm.cSC@listData$SYMBOL)
colnames(nm.cSC_expr_up) <- c("log2FC","padj")

nm.cSC_expr_down <-data.frame(res_down_nm.cSC@listData$log2FoldChange,res_down_nm.cSC@listData$padj,row.names = res_down_nm.cSC@listData$SYMBOL)
colnames(nm.cSC_expr_down) <- c("log2FC","padj")

nm.cSC_df<-rbind(nm.cSC_expr_up,nm.cSC_expr_down)

nm.cSC_df <- nm.cSC_df[rownames(nm.cSC_df)%in%nm.cSC_genes,]
nm.cSC_genes[nm.cSC_genes%in%markers$gene]