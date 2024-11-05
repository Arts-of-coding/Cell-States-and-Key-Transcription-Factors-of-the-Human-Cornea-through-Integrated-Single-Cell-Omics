#' Alter the standard violin plot in Seurat
#'
#' @param obj A seurat object of interest
#' @param feature A gene of interest
#' @param pt.size The size of the points, default is 0
#' @param plot.margin The margins of the plot standard -0.05, 0, -0.05, 0
#' @param ...
#'
#' @return A modified violin plot is generated
#' @export
#'
#' @examples modify_vlnplot(Seur_obj, "TP63")
#' @examples modify_vlnplot(obj = Seur_obj, feature = "TP63")
modify_vlnplot<- function(obj,
                          feature,
                          pt.size = 0,
                          plot.margin = unit(c(-0.05, 0, -0.05, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature,pt.size = pt.size, ... )  +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = rel(1), angle = 0),
          axis.text.y = element_text(size = rel(1)),
          plot.margin = plot.margin )
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,cols2,
                          pt.size = 0,
                          plot.margin = unit(c(-0.05, 0, -0.05, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x,cols=cols2, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x +
                            scale_y_continuous(breaks = c(y)) +
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

QCSeurObj <- function(SeuratObject,pct_MT=30) {
  ###########################################################################
  ## Filtering of the datasets
  seur_obj_all <- SeuratObject
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  seur_obj_all[["percent.mt"]] <- PercentageFeatureSet(seur_obj_all, pattern = "^MT-")
  
  # QC tresholds of counts, features and mt percentage
  seur_obj <- subset(seur_obj_all, subset = nCount_RNA > 2000 & nFeature_RNA > 1000 & percent.mt < pct_MT)
  
  ###########################################################################
  # Plotting percentage MTs before and after
  if (length(unique(unlist(seur_obj$orig.ident)))<6){
    width2 <- 6
  } else {
    width2 <- length(unique(unlist(seur_obj$orig.ident)))/2
  }
  
  pdf(paste(resultsdir,'3.QC_depth_ERCC_MT.pdf',sep="/") ,width=width2,height=6,paper='special')
  print(VlnPlot(object = seur_obj_all, features = ("nCount_RNA"), assay = 'RNA') + scale_y_continuous(limits = c(0,20000)) + ggtitle('total UMI counts all cells'))
  print(VlnPlot(object = seur_obj, features = c("nCount_RNA"), assay = 'RNA') + scale_y_continuous(limits = c(0,20000)) + ggtitle('total UMI counts filtered cells'))
  print(VlnPlot(object = seur_obj_all, features = c("nFeature_RNA"), assay = 'RNA') + scale_y_continuous(limits = c(0,5000)) + ggtitle('genes measured all cells'))
  print(VlnPlot(object = seur_obj, features = c("nFeature_RNA"), assay = 'RNA')  + scale_y_continuous(limits = c(0,5000)) + ggtitle('genes measured filtered cells'))
  print(VlnPlot(object = seur_obj_all, features = c("percent.mt"))+ scale_y_continuous(limits = c(0,40)) + ggtitle('% mitochondrial reads all cells'))
  print(VlnPlot(object = seur_obj, features = c("percent.mt"))+ scale_y_continuous(limits = c(0,40)) + ggtitle('% mitochondrial reads filtered cells'))
  print(FeatureScatter(seur_obj_all, feature1 = ("nCount_RNA"), feature2 = ("nFeature_RNA")) + ggtitle('% counts vs genes measured all cells'))
  print(FeatureScatter(seur_obj, feature1 = ("nCount_RNA"), feature2 = ("nFeature_RNA")) + ggtitle('% counts vs genes measured filtered cells'))
  dev.off()
  
  ###########################################################################
  # normalization of the data
  Idents(seur_obj) <- "cell_type"
  Idents(seur_obj_all) <- "cell_type"
  
  seur_obj <- NormalizeData(
    object = seur_obj, assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
  
  ###########################################################################
  # finding variable features within the datasets and determine cell_cycle influence
  seur_obj <- FindVariableFeatures(seur_obj, selection.method = "vst", nfeatures = 2000)
  seur_obj <- ScaleData(seur_obj, features = rownames(seur_obj), assay = 'RNA')
  
  seur_obj <- CellCycleScoring(seur_obj, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE)
  seur_obj <- RunPCA(seur_obj, features = c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes))
  
  pdf(paste(resultsdir,'5.cell_cycle_markers.pdf',sep="/") ,width=12,height=6,paper='special')
  Idents(seur_obj) <- "Phase"
  print(PCAPlot(seur_obj))
  print(RidgePlot(seur_obj, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), assay= 'RNA', ncol = 2, slot = "data"))
  print(RidgePlot(seur_obj, features = c("S.Score", "G2M.Score")))
  dev.off()
  
  return(seur_obj)
}

PCASeurObj <- function (SeuratObject,
                        pcs=50) {
  seur_obj <- SeuratObject
  ###########################################################################
  # PCA plots for all genes
  seur_obj <- RunPCA(seur_obj, verbose = FALSE, npcs = pcs)
  mat <- Seurat::GetAssayData(seur_obj, assay = "RNA", slot = "scale.data")
  pca <- seur_obj[["pca"]]
  
  # Get the total variance:
  total_variance <- sum(matrixStats::rowVars(mat))
  eigValues = (pca@stdev)^2  ## EigenValues
  varExplained = eigValues / total_variance
  
  # Checking the first PC
  Stdev(object = seur_obj[["pca"]])[1]
  
  Idents(seur_obj) <- "cell_type"
  pdf(paste0(resultsdir,'/6.Principle_components.pdf') ,width=12,height=6,paper='special')
  print(ElbowPlot(seur_obj))
  for (pcs in c(1:15)){
    pc1_viz <- pcs*2-1
    pc2_viz <- pcs*2
    y_label = paste0(paste0(paste0("PC",pc2_viz), ' stdev:  '),round(Stdev(object = seur_obj[["pca"]])[pc2_viz],3))
    x_label = paste0(paste0(paste0("PC",pc1_viz), ' stdev:  '),round(Stdev(object = seur_obj[["pca"]])[pc1_viz],3))
    PC_dimred <- DimPlot(seur_obj, reduction = "pca", dims = c(pc1_viz,pc2_viz))+labs(y= y_label, x = x_label)
    PC1_genes <- DimHeatmap(seur_obj, dims = c(pc1_viz), fast = FALSE)
    PC2_genes <- DimHeatmap(seur_obj, dims = c(pc2_viz), fast = FALSE)
    final_plot <- grid.arrange(as_grob(PC_dimred),as_grob(PC2_genes),as_grob(PC1_genes),
                               ncol=2,
                               as.table=TRUE)
    #heights=c(3,1))
    print(final_plot)}
  dev.off()
  return(seur_obj)
}

QCSeurObjClustering <- function (SeuratObject,dimensions=30) {
  seur_obj <- SeuratObject
  seur_obj <- RunUMAP(seur_obj, dims = 1:dimensions)
  
  # plotting the quality control upon the UMAP to see if clustering is not driven by a specific factor
  pdf(paste(resultsdir,'7a.norm_umap.pdf',sep="/") ,width=8,height=8,paper='special')
  print(DimPlot(seur_obj, label=FALSE)+ ggtitle("original identity"))
  print(FeaturePlot(seur_obj, features = 'nCount_RNA', cols =c('white','purple')))
  print(FeaturePlot(seur_obj, features = 'nFeature_RNA', cols =c('white','purple')))
  print(FeaturePlot(seur_obj, features = 'percent.mt', cols =c('white','purple')))
  print(DimPlot(seur_obj, group.by = 'Phase'))
  print(FeaturePlot(seur_obj, features = 'S.Score', cols =c('white','dodgerblue3')))
  print(FeaturePlot(seur_obj, features = 'G2M.Score', cols =c('white','green4')))
  dev.off()
  return(seur_obj)
}

ClusteringResSeurObj<- function(SeuratObject,dimensions=24) {
  chosen_dims <- dimensions
  ###########################################################################
  # Perform clustering with resolution, to see which resolution fits the best
  pdf(paste(resultsdir, "8.clustering_resolution.pdf", sep = '/'), width = 15, height = 8)
  seur_obj <- SeuratObject
  seur_obj <- FindNeighbors(seur_obj, dims = 1:chosen_dims)
  for (i in 1:20){
    res <- i/20
    #cluster_variable_name <- paste0("RNA_snn_res.", res)
    cluster_variable_name <- paste0("RNA_snn_res.", res)
    seur_obj <- FindClusters(seur_obj, verbose = FALSE, resolution = res, graph.name = "RNA_snn")
    seur_obj <- BuildClusterTree(seur_obj)
    
    if (i != 1){
      
      #
      # seur_obj <- FindClusters(seur_obj, verbose = FALSE, resolution = res)
      # seur_obj <- BuildClusterTree(seur_obj)
      #p1 <- DimPlot(seur_obj,label=TRUE, group.by = 'seurat_clusters') + ggtitle("Louvain Clustering") + ggtitle(paste0("cluster resolution ", res))
      p2 <- ggplot(seur_obj@meta.data, aes(eval(parse(text= i))))+geom_bar(stat="count")
      print(clustree(
        seur_obj,
        prefix = "RNA_snn_res.",
        exprs = c("data", "counts", "scale.data"),
        assay = NULL,
        node_colour = "sc3_stability"
      ))
      print(p2)
    }
    #print(p1)
    
  }
  dev.off()
}


ClusteringResSeurObjIntegrated<- function(SeuratObject,dimensions=24) {
  chosen_dims <- dimensions
  ###########################################################################
  # Perform clustering with resolution, to see which resolution fits the best
  pdf(paste(resultsdir, "8.clustering_resolution.pdf", sep = '/'), width = 15, height = 8)
  seur_obj <- SeuratObject
  seur_obj <- FindNeighbors(seur_obj, dims = 1:chosen_dims)
  
  for (i in 1:20){
    res <- i/20
    #cluster_variable_name <- paste0("RNA_snn_res.", res)
    cluster_variable_name <- paste0("integrated_snn_res.", res)
    seur_obj <- FindClusters(seur_obj, verbose = FALSE, resolution = res, graph.name = "integrated_snn")
    seur_obj <- BuildClusterTree(seur_obj)
    
    if (i != 1){
      
      #
      # seur_obj <- FindClusters(seur_obj, verbose = FALSE, resolution = res)
      # seur_obj <- BuildClusterTree(seur_obj)
      #p1 <- DimPlot(seur_obj,label=TRUE, group.by = 'seurat_clusters') + ggtitle("Louvain Clustering") + ggtitle(paste0("cluster resolution ", res))
      p2 <- ggplot(seur_obj@meta.data, aes(eval(parse(text= i))))+geom_bar(stat="count")
      print(clustree(
        seur_obj,
        prefix = "integrated_snn_res.",
        exprs = c("data", "counts", "scale.data"),
        assay = NULL,
        node_colour = "sc3_stability"
      ))
      print(p2)
    }
    #print(p1)
    
  }
  dev.off()
}
