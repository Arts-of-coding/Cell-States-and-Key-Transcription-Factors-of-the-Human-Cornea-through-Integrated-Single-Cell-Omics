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

###########################################################################
# Storing results in directories of your choice
workdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_Gautam2022/"

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
# load in all datasets
path2 <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/Gautam2021/"

list.files (path = path2, pattern = "SRR")
nam3 <- NULL
v2 <- c()
for (i in list.files (path = path2, pattern = "SRR")){
  print(i)
  v <- paste0(path2,i,"/outs/filtered_feature_bc_matrix/")
  nam <- paste("A", i, sep = "")
  nam <- Read10X(data.dir = v)
  nam2 <- paste("B", i, sep = "")
  nam2 <- CreateSeuratObject(counts = nam, project = paste0("Gautam",i), min.cells = 3, min.features = 200)
  assign(paste0(i), nam2)
  v2 <- c(v2,i)
}

Combined <- SRR11470710
remove(nam)
remove(nam2)
remove(SRR11470710)

# merge the first with the rest of the datasets
seur_obj_all <- merge(Combined, y = c(SRR14742510,SRR14742511), add.cell.ids = v2, project = "Gautam2022all")

# Remove all separate objects
rm(list = ls()[grep("SRR", ls())])

# Perform quality control
seur_obj <- QCSeurObj(SeuratObject = seur_obj_all,pct_MT = 30)
rm(seur_obj_all)

seur_obj <- PCASeurObj(SeuratObject = seur_obj)

seur_obj <- QCSeurObjClustering(SeuratObject = seur_obj,dimensions=24)

ClusteringResSeurObj(SeuratObject = seur_obj,dimensions=24)

chosen_res <- 0.1

# Final chosen resolution
pdf(paste(resultsdir,'8b.clusdendrogram_0.1_24dims.pdf',sep="/") ,width=8,height=8,paper='special')
seur_obj <- FindNeighbors(seur_obj, dims = 1:24)
res <- chosen_res
cluster_variable_name <- paste0("RNA_snn_res.", res)
seur_obj <- FindClusters(seur_obj, verbose = FALSE, resolution = res, graph.name = "RNA_snn")
seur_obj <- BuildClusterTree(seur_obj)
PlotClusterTree(seur_obj, direction = "downwards")
DimPlot(seur_obj, label=TRUE, label.size = 6, group.by = 'seurat_clusters') + ggtitle("Louvain Clustering") + ggtitle(paste0("cluster resolution ", res))
dev.off()

############################################################################
# Rename your cluster numbers to cellnumbers
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11)
new.cluster.ids <- c('cell1','cell2','cell3','cell4','cell5','cell6','cell7','cell8','cell9','cell10','cell11','cell12')
seur_obj$costum_clustering <- plyr::mapvalues(x = seur_obj$RNA_snn_res.0.1, from = current.cluster.ids, to = new.cluster.ids)
cluster_order <- c('cell1','cell2','cell3','cell4','cell5','cell6','cell7','cell8','cell9','cell10','cell11','cell12')

# Re-level object@ident
seur_obj$costum_clustering <- factor(x = seur_obj$costum_clustering, levels = cluster_order)

# Add custom labels to the phylogenetic tree
data.tree <- Tool(object = seur_obj, slot = "BuildClusterTree")
data.tree$tip.label <- levels(seur_obj$costum_clustering)

pdf(paste(resultsdir, "8c.clustering_costum_tree.pdf", sep = '/'), width = 8, height = 8)
#p1 <- DimPlot(seur_obj, label=TRUE, group.by = 'cell_type')
p1 <- DimPlot(seur_obj, label=TRUE,label.size = 8)
p2 <- DimPlot(seur_obj, label=TRUE,label.size = 8, group.by = 'costum_clustering') + ggtitle("Louvain Clustering")
p3 <- DimPlot(seur_obj, label=F,label.size = 5, group.by = 'orig.ident') + ggtitle("Louvain Clustering")
ape::plot.phylo(x = data.tree, direction = "downwards",)
#p1
p2
p3
dev.off()

pdf(paste(resultsdir, "8d.clustering_costum_tree_orig.pdf", sep = '/'), width = 16, height = 8)
p3 <- DimPlot(seur_obj, label=F,label.size = 5, group.by = 'orig.ident') + ggtitle("Louvain Clustering")
p3
dev.off()

saveRDS(seur_obj, file = (paste(resultsdir,'/gautam_unannotated.rds',sep="/")))

seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_Gautam2022/20220210/gautam_unannotated.rds")
marker_gene_file <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/LSC_Marker_Genes_Gautam.csv"
marker_genes_df <- read.table(marker_gene_file, header = TRUE, sep = ',',comment.char = "#", stringsAsFactors = F)

# Setting two cores as max on the system
SnowParam(workers = 2)

# This vector will contain all marker genes across datsets
fullmarkers2<-NULL

# Make dotplots for all datasets with their cell populations and marker genes
for (paper in unique(marker_genes_df$abreviation)){
  seur_obj@active.ident <- seur_obj$costum_clustering
  marker_dir <- paste0(paste0(resultsdir, '/marker_genes_int/'),paper)
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

cluster_order <- c('Ga_C1_CF','Ga_C2_Cj',
                   'Ga_C3_THCEC','Ga_C4_EHCEC','Ga_C5_Cj',
                   'Ga_C6_Cj','Ga_C7_Cj','Ga_C8_Mel',
                   'Ga_C9_THCEC','Ga_C10_CF','Ga_C11_Mon',
                   'Ga_C12_CTC')

# Generate new cluster labeling:
current.cluster.ids <- c('cell1','cell2',
                         'cell3','cell4','cell5',
                         'cell6','cell7','cell8',
                         'cell9','cell10','cell11',
                         'cell12')


new.cluster.ids <- c('Ga_C1_CF','Ga_C2_Cj',
                     'Ga_C3_THCEC','Ga_C4_EHCEC','Ga_C5_Cj',
                     'Ga_C6_Cj','Ga_C7_Cj','Ga_C8_Mel',
                     'Ga_C9_THCEC','Ga_C10_CF','Ga_C11_Mon',
                     'Ga_C12_CTC')

seur_obj$costum_clustering <- plyr::mapvalues(x = seur_obj$costum_clustering, from = current.cluster.ids, to = new.cluster.ids)

# Re-level object@ident
seur_obj$costum_clustering <- factor(x = seur_obj$costum_clustering, levels = cluster_order)
seur_obj$costum_clustering

# Add custom labels to the phylogenetic tree
data.tree <- Tool(object = seur_obj, slot = "BuildClusterTree")
data.tree$tip.label <- levels(seur_obj$costum_clustering)

# Plot the tree and the final clustering
pdf(paste(resultsdir, "13.clustering_costum_final.pdf", sep = '/'), width = 6, height = 5)
#p1 <- DimPlot(seur_obj, label=TRUE, group.by = 'cell_type')
p1 <- DimPlot(seur_obj, label=F,label.size = 8)
p2 <- DimPlot(seur_obj, label=F,label.size = 8, group.by = 'costum_clustering') + ggtitle("Louvain Clustering")
p3 <- DimPlot(seur_obj, label=F,label.size = 5, group.by = 'orig.ident') + ggtitle("Louvain Clustering")
ape::plot.phylo(x = data.tree, direction = "downwards",)
#
p1
p2
p3
#print(multiplot(p1, p2, p3, cols = 3))
dev.off()

marker_gene_file <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/LSC_Marker_Genes_Gautam.csv"
marker_genes_df <- read.table(marker_gene_file, header = TRUE, sep = ',',comment.char = "#", stringsAsFactors = F)

# Setting two cores as max on the system
SnowParam(workers = 2)

# This vector will contain all marker genes across datsets
fullmarkers2<-NULL

# Make dotplots for all datasets with their cell populations and marker genes
for (paper in unique(marker_genes_df$abreviation)){
  seur_obj@active.ident <- seur_obj$costum_clustering
  marker_dir <- paste0(paste0(resultsdir, '/marker_genes_int/'),paper)
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
  pdf(paste0(paste(paste(marker_dir, paper, sep = '/'), sep = '_'),'full_markers_annotated.pdf'), width =7, height = 5)
  print(DotPlot(object = seur_obj, features = fullmarkers,cluster.idents = T,)+ RotatedAxis())
  dev.off()
}

###########################################################################
# Saving the annotated file
saveRDS(seur_obj, file = paste(resultsdir,"GautamRNAannotated.rds"))





# BELOW IS EXTRA
###########################################################################
# Annotating your cell clusters based on marker gene files
Idents(seur_obj) <- "costum_clustering"

#marker_gene_file <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/LSC_Marker_Genes.csv" #lapointe
marker_gene_file <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/OM_Marker_Genes_V1.csv"
marker_genes_df <- read.table(marker_gene_file, header = TRUE, sep = ',',comment.char = "#", stringsAsFactors = F)

# Setting two cores as max on the system
SnowParam(workers = 2)

for (paper in unique(marker_genes_df$abreviation)){
  marker_dir <- paste0(paste0(resultsdir, '/marker_genes_0.05/'),paper)
  print(marker_dir)
  system(paste("mkdir -p ", marker_dir))
  sub_marker_df <- marker_genes_df[marker_genes_df$abreviation == paper,]

  for (cell_type in unique(sub_marker_df$name_plot)){
    #print(cell_type)
    marker_genes <- sub_marker_df[sub_marker_df$name_plot == cell_type,]$marker_genes
    marker_genes <- unlist(strsplit(marker_genes, ','))

    m <- length(unique(sub_marker_df$name_plot))
    n <- length(marker_genes)
    k <- 10 ## your LEN
    marker_gene_sets <- split(marker_genes, rep(1:ceiling(n/k), each=k)[1:n])

    pdf(paste0(paste(paste(marker_dir, cell_type, sep = '/'), sep = '_'),'stacked_markers.pdf'), width = 10, height = n*2)
    print(StackedVlnPlot(obj = seur_obj, features = unlist(marker_gene_sets)))
    dev.off()

    for (count_type in c('raw_counts', 'seurat_norm_counts','SAVER_imputed_counts')){        #loop  over the different type of normalized data
    #for (count_type in c('seurat_norm_counts')){        #loop  over the different type of normalized data

      print(count_type)

      pdf(paste0(paste(paste(marker_dir, cell_type, sep = '/'), count_type, sep = '_'),'.pdf'), width = 20, height = 15)

      for (genes in  marker_gene_sets){#loop over subsets of max 12 genes to vizualize
        print(genes)
        plot_list1 <- DimPlot(seur_obj, label=TRUE, combine = F, pt.size = 0.1)
        #plot_list1 <- append(plot_list1, DimPlot(seur_obj, label=TRUE, group.by = 'custom_clustering', combine = F, pt.size = 0.1))
        count_plots <- plot_list1
        count_plots_clusters <- plot_list1

        if (count_type == 'raw_counts'){
          count_plots <- append(count_plots, try(VlnPlot(seur_obj, combine=F, features = genes, ncol = 4, group.by = 'cell_type', assay = 'RNA', slot = 'counts')))
          count_plots_clusters <- append(count_plots_clusters, try(VlnPlot(seur_obj, combine=F, features = genes, ncol = 4, group.by = 'costum_clustering', assay = 'RNA', slot = 'counts')))}
        if(count_type == 'seurat_norm_counts'){
          count_plots <- append(count_plots, try(VlnPlot(seur_obj, combine=F, features = genes, ncol = 4, group.by = 'cell_type', assay = 'RNA')))
          count_plots_clusters <- append(count_plots_clusters, try(VlnPlot(seur_obj, combine=F, features = genes, ncol = 4, group.by = 'costum_clustering', assay = 'RNA')))}

        if(count_type == 'SAVER_imputed_counts'){
          count_plots <- append(count_plots, try(VlnPlot(seur_obj, combine=F, features = genes, ncol = 4, group.by = 'cell_type', assay = 'RNA_SAVER')))
          count_plots_clusters <- append(count_plots_clusters, try(VlnPlot(seur_obj, combine=F, features = genes, ncol = 4, group.by = 'costum_clustering', assay = 'RNA_SAVER')))}


        plot_list1 <- append(plot_list1,FeaturePlot(seur_obj, features = genes, pt.size = 0.2, ncol = 4, combine=F))
        print(cowplot::plot_grid(plotlist = plot_list1))
        #print(FeaturePlot(seur_obj, features = genes, pt.size = 0.2, ncol = 4))

        p1 <- cowplot::plot_grid(plotlist = count_plots)
        title <- ggdraw() + draw_label(paste0(count_type,"of marker genes per timepoint"), fontface = 'bold')
        print(cowplot::plot_grid(title, p1, ncol = 1, rel_heights = c(0.1, 1)))
        p1 <- cowplot::plot_grid(plotlist = count_plots_clusters)
        title <- ggdraw() + draw_label(paste0(count_type,"of marker genes per cluster"), fontface = 'bold')
        print(cowplot::plot_grid(title, p1, ncol = 1, rel_heights = c(0.1, 1)))
      }
      dev.off()}
  }
}

###########################################################################
# Lets find marker genes unbiased for each cell cluster
cluster.markers <- FindAllMarkers(seur_obj, only.pos = TRUE)

# making a nice heatmap of expression in each cluster not significant yet
n_genes <- 40
heatmap.markers <- cluster.markers %>% group_by(cluster) %>% top_n(n_genes, avg_log2FC)
DoHeatmap(seur_obj, features = heatmap.markers$gene) + NoLegend()

#filter on only significance adjusted value
cluster.markers$gene_name_shorter <- str_sub(cluster.markers$gene,end = -1)
cluster.markers_fc <- cluster.markers[cluster.markers$p_val_adj < 0.01,]

cluster.markers_fc <- cluster.markers_fc[(cluster.markers_fc$avg_log2FC > 0.58) ,]
table(cluster.markers_fc$cluster) # number of significantly expressed genes with a high log2FC (>0.58)

###########################################################################
# Saving the file to make sure no progress is lost
#saveRDS(seur_obj, file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAunannotated2.rds")
#seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAunannotated2.rds")

###########################################################################
# Appending the GO terms of the unannotated clusters
plot_list <- list()
cluster_GOs <- c()
expressed_genes <- rownames(cluster.markers_fc)

for (cluster in unique(cluster.markers_fc$cluster)){
  print(cluster)
  df_subset <- cluster.markers_fc[cluster.markers_fc$cluster == cluster,]
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
  plot_list <- c(plot_list, list(go_plot1))
  cluster_GOs <- append(cluster_GOs, (paste('cluster of cell type ', cluster, sep = '')))
}
dev.off()

# plotting the GO-terms
grob2 <- ggarrange(plotlist= plot_list, nrow = 18, ncol = 1, labels = cluster_GOs) #, align = 'v')
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
top5_genes <- DoHeatmap(seur_obj, features = top10$gene) + NoLegend()
top10 <- cluster.markers_fc %>% group_by(cluster) %>% top_n(n = 10, wt = 'avg_log2FC')
marker_top5 <- DoHeatmap(seur_obj, features = top10$gene) + NoLegend()

n_genes <- 100
heatmap.markers <- cluster.markers_fc %>% group_by(cluster) #%>% top_n(n_genes, avg_logFC)

cell_type_plot <- DimPlot(seur_obj, label=TRUE,pt.size = 2) + ggtitle("Timepoint")
cluster_plot <- DimPlot(seur_obj, label=TRUE, group.by = 'costum_clustering', pt.size = 2) + ggtitle("Louvain Clustering")
grob_clustering = ggarrange(plotlist = list(cell_type_plot, cluster_plot) , ncol =1) #labels = cluster_GOs) #align = 'v')
nclust<- length(unique(as.numeric(cluster.markers_fc$cluster)))
RGB_colours_ggplot <- as.list(as.character(gg_color_hue(nclust)))
names <- as.numeric(unique(heatmap.markers$cluster))
col_fun = colorRamp2(names, RGB_colours_ggplot)#make sure the rows correspond to ggplot colour mapping

mat <- as.matrix(seur_obj@assays$RNA@scale.data[heatmap.markers$gene,])
mat <- rbind(mat,seur_obj@meta.data$costum_clustering)
mat <- mat[,order(mat[nrow(mat),])]

column_ha <- HeatmapAnnotation(cluster = mat[nrow(mat),], col =list(cluster = col_fun))
breaks <- mat[nrow(mat),]
mat <- mat[-nrow(mat),]
row_ha = rowAnnotation(adj_p_vallue = heatmap.markers$p_val_adj)
clust_heatmap <- grid.grabExpr(draw(Heatmap(mat, column_split = breaks,row_split = as.numeric(heatmap.markers$cluster), cluster_columns = F, cluster_rows = F,show_row_names = F,show_column_names = F,row_names_gp = gpar(fontsize = 6), top_annotation = column_ha, left_annotation = row_ha, row_names_rot = -40)))

###########################################################################
# Plotting the GO terms of the un-annotated clusters and the heatmap

pdf(paste0(resultsdir,'/9.GO_terms_clusters.pdf') ,width=120,height=4,paper='special')
cowplot::plot_grid(plotlist =grob2, ncol = 1, nrow = 1)#, rel_heights = c(1, 4,2), rel_widths = c(3,2,2))
dev.off()

pdf(paste0(resultsdir,'/9.complex_heatmap.pdf') ,width=40,height=24,paper='special')
cowplot::plot_grid(marker_top5,grob_clustering, ncol = 3, nrow = 1, rel_heights = c(1, 4,2), rel_widths = c(3,2,2))
dev.off()

# Write the markers as a table
write.table(cluster.markers_fc, file = paste0(resultsdir,'/9b.cluster_markers.csv'), sep = ',')

###########################################################################
# Determining the cell types based upon markers
# Annotation notes: see github file regarding annotation

cluster_order <- c('StCSC1','nsKer1',
                   'StCSC2','nsKer2','tdCE',
                   'LSCs','LSE','BCE',
                   'sKer','LSC2','EC1',
                   'cell12','EC2','Mel')

# Generate new cluster labeling:
current.cluster.ids <- c('cell1','cell2',
                         'cell3','cell4','cell5',
                         'cell6','cell7','cell8',
                         'cell9','cell10','cell11',
                         'cell12','cell13','cell14')

# sKer = stromal keratocytes
# nsKer = non stromal keratocytes
# BCE = basal corneal epithelium
# SCE = suprabasal corneal epithelium
# tdCE = terminally differentiated corneal epithelium
# EC = endothelial cells
# StCSC = stromal cells and stromal stem cells
# LSE = limbal superficial epithelium; this is migratory
# LSC = limbal stem cells

new.cluster.ids <- c('StCSC1','nsKer1',
                     'StCSC2','nsKer2','tdCE',
                     'LSCs','LSE','BCE',
                     'sKer','LSC2','EC1',
                     'cell12','EC2','Mel')

seur_obj$costum_clustering <- plyr::mapvalues(x = seur_obj$costum_clustering, from = current.cluster.ids, to = new.cluster.ids)

# Re-level object@ident
seur_obj$costum_clustering <- factor(x = seur_obj$costum_clustering, levels = cluster_order)
seur_obj$costum_clustering

# Add custom labels to the phylogenetic tree
data.tree <- Tool(object = seur_obj, slot = "BuildClusterTree")
data.tree$tip.label <- levels(seur_obj$costum_clustering)

# Plot the tree and the final clustering
pdf(paste(resultsdir, "13.clustering_costum_final.pdf", sep = '/'), width = 8, height = 8)
#p1 <- DimPlot(seur_obj, label=TRUE, group.by = 'cell_type')
p1 <- DimPlot(seur_obj, label=TRUE,label.size = 8)
p2 <- DimPlot(reference, label=TRUE,label.size = 8, group.by = 'costum_clustering') + ggtitle("Louvain Clustering")
p3 <- DimPlot(seur_obj, label=F,label.size = 5, group.by = 'orig.ident') + ggtitle("Louvain Clustering")
ape::plot.phylo(x = data.tree, direction = "downwards",)
#
p1
p2
p3
#print(multiplot(p1, p2, p3, cols = 3))
dev.off()

###########################################################################
# Saving the annotated file
saveRDS(seur_obj, file = paste(resultsdir,"laPointeRNAannotated.rds"))
#seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAannotated.rds")
#seur_obj<- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_lapointe2021/20211201 laPointeRNAannotated.rds")
############################################################################
# Generating metadata file for the cell counts
RNA_metadata <- seur_obj@meta.data

# pre-allocate space
metnames <- c()
metcounts <- character()

for(i in unique(RNA_metadata$costum_clustering)){
  metnames <- append(metnames,i)
  metcounts <- append(metcounts, nrow(RNA_metadata[RNA_metadata$costum_clustering == i,]))
}
popdf <- data.frame(metnames,metcounts)
names(popdf) <- c("population","number_of_cells")
write.table(popdf, file = paste0(resultsdir,'/','cellcounts.tsv'), sep = '\t',quote = F,row.names = F)

###########################################################################
# Make a pseudobulk table
seur_obj@meta.data$costum_clustering

sce_qc <- SingleCellExperiment(assays = list(counts = GetAssayData(object = seur_obj, slot = "data")))

pseudobulk_df <- as.data.frame(row.names(counts(sce_qc)))

sce_qc <- SingleCellExperiment(assays = list(counts = GetAssayData(object = seur_obj, slot = "data")))

pseudobulk_df <- as.data.frame(row.names(counts(sce_qc)))
#write.table(pseudobulk_df, file = paste0(resultsdir,'/pseudobulk.tsv'), sep = '\t') # this is only if you do not want the replicates

###########################################################################
# Splitting the pseudobulk datasets unbiased into two replicates for DEG
sce_qc$sample <- seur_obj@meta.data$costum_clustering
sce_qc$sample

sce_qc$rep <- seur_obj@meta.data$orig.ident
sce_qc$rep

sce_qc$reps <- "rep1"
sce_qc$reps[sce_qc$rep == "lako1" | sce_qc$rep == "lako3" | sce_qc$rep == "lako5" | sce_qc$rep == "lako7"] <- "rep2"

# putting values of artificial rep2 and rep1 into the pseudobulkdataframecolumns
for (sample in unique(sce_qc$sample)){
  pseudobulk_df[[sample]] <- rowSums(counts(sce_qc)[,sce_qc$sample == sample & sce_qc$reps == "rep2"])
  pseudobulk_df[ , paste0(sample,"1")] <- rowSums(counts(sce_qc)[,sce_qc$sample == sample & sce_qc$reps == "rep1"])
}

row.names(pseudobulk_df) <- pseudobulk_df$`row.names(counts(sce_qc))`
pseudobulk_df$`row.names(counts(sce_qc))` <- NULL
write.table(data.frame("ID"=rownames(pseudobulk_df),pseudobulk_df), file = paste0(resultsdir,'/pseudobulk_4cells.tsv'), sep = '\t',quote = F, row.names = F)
###########################################################################
# Additional: interesting genes
Markergenes <- c("IRF1", "IRF5", "FOXC1", "NR2F2", "PAX6", "ELF3", "TNF", "CXCL1", "CXCL2", "CXCL3", "PTGS2", "NFKBIA", "OTX1", "ELK3", "HES4")
pdf(paste(resultsdir,'interesting.pdf'), width = 20, height = 15)
FeaturePlot(seur_obj, features = Markergenes)
dev.off()

###########################################################################
# Integrating Lako and Lapointe -> better/correct annotation
query <- seur_obj
reference <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAannotated.rds")

# SnowParam(workers = 2)
# anchors <- FindTransferAnchors(reference = reference, query = query,reduction = "cca")
#
# # transfer labels
# SnowParam(workers = 2)
# predictions <- TransferData(
#   anchorset = anchors,
#   refdata = reference$costum_clustering,
#   weight.reduction = 'cca'
# )
#
# SnowParam(workers = 2)
# query <- AddMetaData(object = query, metadata = predictions)

reference <- AddMetaData(reference, metadata="lako", col.name="Condition")

query <-AddMetaData(query, metadata="lapointe", col.name="Condition")

your_list <- list(reference, query)



SnowParam(workers = 2)
your_list <- lapply(X = your_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

SnowParam(workers = 2)
immune.anchors <- FindIntegrationAnchors(object.list = your_list, dims = 1:30,reduction = 'cca')

SnowParam(workers = 2)
cornea.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:30)

# Running the UMAP

# original unmodified data still resides in the 'RNA' assay
DefaultAssay(cornea.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
cornea.combined <- ScaleData(cornea.combined, verbose = FALSE)
cornea.combined <- RunPCA(cornea.combined, npcs = 10, verbose = FALSE)
cornea.combined <- RunUMAP(cornea.combined, reduction = "pca", dims = 1:10)

# cornea.combined <- QCSeurObj(cornea.combined)
# cornea.combined <- RunPCA(cornea.combined,npcs=10)
# cornea.combined <- QCSeurObjClustering(cornea.combined,dimensions = 10)

# Plot the tree and the final clustering
pdf(paste(resultsdir, "14.clustering_costum_transfered_cca.pdf", sep = '/'), width = 16, height = 16)
#p1 <- DimPlot(seur_obj, label=TRUE, group.by = 'cell_type')
p1 <- DimPlot(cornea.combined, reduction = 'umap', label=TRUE,label.size = 8)
p2 <- DimPlot(cornea.combined, label=TRUE,label.size = 8, group.by = 'costum_clustering') + ggtitle("Louvain Clustering")
p3 <- DimPlot(cornea.combined, label=F,label.size = 5, group.by = 'orig.ident') + ggtitle("Louvain Clustering")
p4 <- DimPlot(cornea.combined, label=F,label.size = 5, group.by = 'Condition') + ggtitle("Louvain Clustering")
#ape::plot.phylo(x = data.tree, direction = "downwards",)
#
p1
p2
p3
p4
#print(multiplot(p1, p2, p3, cols = 3))
dev.off()

###########################################################################
# Additional: interesting genes
Markergenes <- c("IRF1", "IRF5", "FOXC1", "NR2F2", "PAX6", "ELF3", "TNF", "CXCL1", "CXCL2", "CXCL3", "PTGS2", "NFKBIA", "OTX1", "ELK3", "HES4","KRT14")
pdf(paste(resultsdir,'interesting_cca.pdf'), width = 20, height = 15)
FeaturePlot(cornea.combined, features = Markergenes)
dev.off()

saveRDS(cornea.combined, file = paste(resultsdir,"cornea_combined_cca.rds"))

cornea.combined <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_lapointe2021/20211202 cornea_combined_cca.rds")

cornea.combined <- QCSeurObjClustering(cornea.combined,dimensions = 10)
ClusteringResSeurObjIntegrated(SeuratObject = cornea.combined,dimensions=10) # 10 based on lako data

chosen_res <- 0.15

# Final chosen resolution
pdf(paste(resultsdir,'20b.clusdendrogram_combined_cca.pdf',sep="/") ,width=8,height=8,paper='special')
cornea.combined <- FindNeighbors(cornea.combined, dims = 1:10)
res <- chosen_res
#cluster_variable_name <- paste0("RNA_snn_res.", res)
cornea.combined <- FindClusters(cornea.combined, verbose = FALSE, resolution = res)
cornea.combined <- BuildClusterTree(cornea.combined)
PlotClusterTree(cornea.combined, direction = "downwards")
DimPlot(cornea.combined, label=TRUE, label.size = 6, group.by = 'integrated_snn_res.0.15') + ggtitle("Louvain Clustering") + ggtitle(paste0("cluster resolution ", res))
dev.off()
##########################################################################
# rPCA

reference <- AddMetaData(reference, metadata="lako", col.name="Condition")

query <-AddMetaData(query, metadata="lapointe", col.name="Condition")

your_list <- list(reference, query)

SnowParam(workers = 2)
your_list <- lapply(X = your_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

SnowParam(workers = 2)
immune.anchors <- FindIntegrationAnchors(object.list = your_list, dims = 1:30,reduction = 'rpca')

SnowParam(workers = 2)
cornea.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:30)

# Running the UMAP

# original unmodified data still resides in the 'RNA' assay
DefaultAssay(cornea.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
cornea.combined <- ScaleData(cornea.combined, verbose = FALSE)
cornea.combined <- RunPCA(cornea.combined, npcs = 10, verbose = FALSE)
cornea.combined <- RunUMAP(cornea.combined, reduction = "pca", dims = 1:10)

# cornea.combined <- QCSeurObj(cornea.combined)
# cornea.combined <- RunPCA(cornea.combined,npcs=10)
# cornea.combined <- QCSeurObjClustering(cornea.combined,dimensions = 10)

# Plot the tree and the final clustering
pdf(paste(resultsdir, "14.clustering_costum_transfered_rpca.pdf", sep = '/'), width = 16, height = 16)
#p1 <- DimPlot(seur_obj, label=TRUE, group.by = 'cell_type')
p1 <- DimPlot(cornea.combined, reduction = 'umap', label=TRUE,label.size = 8)
p2 <- DimPlot(cornea.combined, label=TRUE,label.size = 8, group.by = 'costum_clustering') + ggtitle("Louvain Clustering")
p3 <- DimPlot(cornea.combined, label=F,label.size = 5, group.by = 'orig.ident') + ggtitle("Louvain Clustering")
p4 <- DimPlot(cornea.combined, label=F,label.size = 5, group.by = 'Condition') + ggtitle("Louvain Clustering")
#ape::plot.phylo(x = data.tree, direction = "downwards",)
#
p1
p2
p3
p4
#print(multiplot(p1, p2, p3, cols = 3))
dev.off()

#saveRDS(cornea.combined, file = paste(resultsdir,"cornea_combined_rpca.rds"))
cornea.combined<- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_lapointe2021/20211202 cornea_combined_rpca.rds")
cornea.combined@assays$integrated<-NULL


###########################################################################
# Additional: interesting genes
Markergenes <- c("IRF1", "IRF5", "FOXC1", "NR2F2", "PAX6", "ELF3", "TNF", "CXCL1", "CXCL2", "CXCL3", "PTGS2", "NFKBIA", "OTX1", "ELK3", "HES4","KRT14")
pdf(paste(resultsdir,'interesting_rpca.pdf'), width = 20, height = 15)
FeaturePlot(cornea.combined, features = Markergenes)
dev.off()

##########################################################################
# Rerunning clustering & annotation

cornea.combined <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_lapointe2021/20211202 cornea_combined_rpca.rds")

cornea.combined <- QCSeurObjClustering(cornea.combined,dimensions = 10)

DefaultAssay(cornea.combined) <- "RNA"

#cornea.combined <- FindNeighbors(cornea.combined, dims = 1:10,assay = seur_obj@assays)

ClusteringResSeurObj(SeuratObject = cornea.combined,dimensions=10) # 10 based on lako data

chosen_res <- 0.23

# Final chosen resolution
pdf(paste(resultsdir,'20b.clusdendrogram_combined_rpca.pdf',sep="/") ,width=8,height=8,paper='special')

res <- chosen_res
#cluster_variable_name <- paste0("RNA_snn_res.", res)
cornea.combined <- FindClusters(cornea.combined, verbose = FALSE, resolution = res)
cornea.combined <- BuildClusterTree(cornea.combined)
PlotClusterTree(cornea.combined, direction = "downwards")
DimPlot(cornea.combined, label=TRUE, label.size = 6, group.by = 'integrated_snn_res.0.23') + ggtitle("Louvain Clustering") + ggtitle(paste0("cluster resolution ", res))
dev.off()

current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11,12,13,14)

new.cluster.ids <- c('cell1','cell2','cell3','cell4','cell5','cell6','cell7','cell8','cell9','cell10','cell11','cell12','cell13','cell14','cell15')#,'cell16','cell17')
cornea.combined$costum_clustering <- plyr::mapvalues(x = cornea.combined$integrated_snn_res.0.23, from = current.cluster.ids, to = new.cluster.ids)
cluster_order <- c('cell1','cell2','cell3','cell4','cell5','cell6','cell7','cell8','cell9','cell10','cell11','cell12','cell13','cell14','cell15')#,'cell16')'cell17') # change if nessesary

# Re-level object@ident
cornea.combined$costum_clustering <- factor(x = cornea.combined$costum_clustering, levels = cluster_order)

# Add custom labels to the phylogenetic tree
data.tree <- Tool(object = cornea.combined, slot = "BuildClusterTree")
data.tree$tip.label <- levels(cornea.combined$costum_clustering)

# Generating a color matrix of all cell_types
cells <- unique(levels(cornea.combined$costum_clustering))
col_mat <- viridis(length(cells),alpha = 0.5,option = "H",begin = 0.1)
dfcol <- as.data.frame(col_mat,cells)

pdf(paste(resultsdir, "20c.clustering_costum_tree.pdf", sep = '/'), width = 8, height = 8)
#p1 <- DimPlot(cornea.combined, label=TRUE, group.by = 'cell_type')
p1 <- DimPlot(cornea.combined, label=TRUE,label.size = 8)
p2 <- DimPlot(cornea.combined, label=TRUE,label.size = 8, group.by = 'costum_clustering',cols = dfcol$col_mat) + ggtitle("Louvain Clustering")
p3 <- DimPlot(cornea.combined, label=F,label.size = 5, group.by = 'orig.ident') + ggtitle("Louvain Clustering")
ape::plot.phylo(x = data.tree, direction = "downwards",)
#p1
p2
p3
dev.off()

Idents(cornea.combined) <- "costum_clustering"

#marker_gene_file <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/LSC_Marker_Genes.csv" #lapointe
marker_gene_file <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/LSC_Marker_Genes.csv"
marker_genes_df <- read.table(marker_gene_file, header = TRUE, sep = ',',comment.char = "#", stringsAsFactors = F)

# Setting two cores as max on the system
SnowParam(workers = 2)

for (paper in unique(marker_genes_df$abreviation)){
  marker_dir <- paste0(paste0(resultsdir, '/marker_genes_rpca_int/'),paper)
  print(marker_dir)
  system(paste("mkdir -p ", marker_dir))
  sub_marker_df <- marker_genes_df[marker_genes_df$abreviation == paper,]

  for (cell_type in unique(sub_marker_df$name_plot)){
    #print(cell_type)
    marker_genes <- sub_marker_df[sub_marker_df$name_plot == cell_type,]$marker_genes
    marker_genes <- unlist(strsplit(marker_genes, ','))

    m <- length(unique(sub_marker_df$name_plot))
    n <- length(marker_genes)
    k <- 10 ## your LEN
    marker_gene_sets <- split(marker_genes, rep(1:ceiling(n/k), each=k)[1:n])

    pdf(paste0(paste(paste(marker_dir, cell_type, sep = '/'), sep = '_'),'stacked_markers.pdf'), width = 10, height = n*2)
    print(StackedVlnPlot(obj = cornea.combined, cols2 = dfcol$col_mat,features = unlist(marker_gene_sets)))
    dev.off()

    for (count_type in c('raw_counts', 'seurat_norm_counts','SAVER_imputed_counts')){        #loop  over the different type of normalized data
      #for (count_type in c('seurat_norm_counts')){        #loop  over the different type of normalized data

      print(count_type)

      pdf(paste0(paste(paste(marker_dir, cell_type, sep = '/'), count_type, sep = '_'),'.pdf'), width = 20, height = 15)

      for (genes in  marker_gene_sets){#loop over subsets of max 12 genes to vizualize
        print(genes)
        plot_list1 <- DimPlot(cornea.combined, label=TRUE, combine = F, pt.size = 0.1)
        #plot_list1 <- append(plot_list1, DimPlot(cornea.combined, label=TRUE, group.by = 'custom_clustering', combine = F, pt.size = 0.1))
        count_plots <- plot_list1
        count_plots_clusters <- plot_list1

        if (count_type == 'raw_counts'){
          count_plots <- append(count_plots, try(VlnPlot(cornea.combined, combine=F, features = genes, ncol = 4, group.by = 'cell_type', assay = 'RNA', slot = 'counts')))
          count_plots_clusters <- append(count_plots_clusters, try(VlnPlot(cornea.combined, combine=F, features = genes, ncol = 4, group.by = 'costum_clustering', assay = 'RNA', slot = 'counts')))}
        if(count_type == 'seurat_norm_counts'){
          count_plots <- append(count_plots, try(VlnPlot(cornea.combined, combine=F, features = genes, ncol = 4, group.by = 'cell_type', assay = 'RNA')))
          count_plots_clusters <- append(count_plots_clusters, try(VlnPlot(cornea.combined, combine=F, features = genes, ncol = 4, group.by = 'costum_clustering', assay = 'RNA')))}

        if(count_type == 'SAVER_imputed_counts'){
          count_plots <- append(count_plots, try(VlnPlot(cornea.combined, combine=F, features = genes, ncol = 4, group.by = 'cell_type', assay = 'RNA_SAVER')))
          count_plots_clusters <- append(count_plots_clusters, try(VlnPlot(cornea.combined, combine=F, features = genes, ncol = 4, group.by = 'costum_clustering', assay = 'RNA_SAVER')))}


        plot_list1 <- append(plot_list1,FeaturePlot(cornea.combined, features = genes, pt.size = 0.2, ncol = 4, combine=F))
        print(cowplot::plot_grid(plotlist = plot_list1))
        #print(FeaturePlot(cornea.combined, features = genes, pt.size = 0.2, ncol = 4))

        p1 <- cowplot::plot_grid(plotlist = count_plots)
        title <- ggdraw() + draw_label(paste0(count_type,"of marker genes per timepoint"), fontface = 'bold')
        print(cowplot::plot_grid(title, p1, ncol = 1, rel_heights = c(0.1, 1)))
        p1 <- cowplot::plot_grid(plotlist = count_plots_clusters)
        title <- ggdraw() + draw_label(paste0(count_type,"of marker genes per cluster"), fontface = 'bold')
        print(cowplot::plot_grid(title, p1, ncol = 1, rel_heights = c(0.1, 1)))
      }
      dev.off()}
  }
}

cornea.combined <- PCASeurObj(SeuratObject = cornea.combined,pcs = 10)
saveRDS(cornea.combined, file = paste(resultsdir,"cornea_combined_rpca_PCA_UMAP.rds"))
cornea.combined <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_lapointe2021/20211210 cornea_combined_rpca_PCA_UMAP.rds")
