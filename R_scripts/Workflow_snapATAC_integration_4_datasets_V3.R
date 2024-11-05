###########################################################################
# SnapATAC workflow

# setting up libraries for installing SnapATAC package
library(doSNOW)
library(devtools)
#install.packages("plot3D")
library(plot3D)
library(rhdf5)
library(Rhdf5lib)
library(BiocParallel)
library(GenomicRanges)
#devtools::install_github("r3fang/SnapATAC") #skip all updates
library(SnapATAC)
library(plyr)
library(readr)
library(leiden)
library(umap)
library(ggplot2)

###########################################################################
# Storing results
workdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scATAC/"

# setting up results directory
dateoftoday <- gsub("-", "", as.character(Sys.Date()))
resultsdir <- paste0(workdir, dateoftoday)
system(paste("mkdir -p ", resultsdir))

###########################################################################
# loading in and joining the snap files data
x.sp1 = createSnap(file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/snap_files/ataclako4.snap", sample="atac_lako4",num.cores=1)
summarySnap(x.sp1)
x.sp2 = createSnap(file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/snap_files/ataclako3.snap", sample="atac_lako3",num.cores=1)
summarySnap(x.sp2)
x.sp3 = createSnap(file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/snap_files/ataclako1.snap", sample="atac_lako1",num.cores=1)
summarySnap(x.sp3)
x.sp0 = createSnap(file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/snap_files/ataclako2.snap", sample="atac_lako2",num.cores=1)
summarySnap(x.sp0)

x.sp = snapRbind(x.sp1,x.sp2)
x.sp = snapRbind(x.sp,x.sp3)
x.sp = snapRbind(x.sp,x.sp0)

###########################################################################
# loading in and joining the barcode files
mydir = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/single/"
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
myfiles

barcodes = ldply(myfiles, read_csv)
barcodes = barcodes[2:nrow(barcodes),]
head(barcodes)

barcodes$logUMI = log10(barcodes$passed_filters + 1)
head(barcodes$logUMI)

###########################################################################
# QC of the snap files

pdf(paste(resultsdir,'barcodeQC1_all.pdf',sep=""),width=6,height=6,paper='special')
for (name in unique(x.sp@sample)){
  plotBarcode(
    main.title = name,
    obj=x.sp[x.sp@sample==name,],
    pdf.file.name=NULL,
    pdf.width=7,
    pdf.height=7,
    col="grey",
    border="grey",
    breaks=50)
  mtext(name, outer=TRUE,  cex=1, line=-1)
}
dev.off()

# Figure with barcodes before filtering
pdf(paste(resultsdir,'barcodeQC_UMI1_all.pdf',sep=""),width=6,height=6,paper='special')
hist(barcodes$logUMI, main= "frequency of unique counts barcodes",xlab = "log10 unique counts")
dev.off()

# Filter the barcodes
idd = which(barcodes$logUMI>=3.75 & barcodes$logUMI <= 5)
barcodes <- barcodes[idd,]

# Figure with barcodes after filtering
pdf(paste(resultsdir,'barcodeQC_UMI2_all.pdf',sep=""),width=6,height=6,paper='special')
hist(barcodes$logUMI, main= "frequency of unique counts barcodes",xlab = "log10 unique counts")
dev.off()

###########################################################################
# QC of fragment number, mitochondrial ratio and duplicate ratio
x.sp = filterCells(
  obj=x.sp,
  subset.names=c("fragment.num"),
  low.thresholds=c(5000),
  high.thresholds=c(500000)
)

x.sp = filterCells(
  obj=x.sp,
  subset.names=c('mito.ratio'),
  low.thresholds=c(0),
  high.thresholds=c(0.2)
)

x.sp = filterCells(
  obj=x.sp,
  subset.names=c('dup.ratio'),
  low.thresholds=c(0),
  high.thresholds=c(0.8)
)
summary(x.sp)

###########################################################################
# Plotting the quality control
pdf(paste(resultsdir,'barcodeQC2_all.pdf',sep=""),width=6,height=6,paper='special')
for (name in unique(x.sp@sample)){
  plotBarcode(
    main.title = name,
    obj=x.sp[x.sp@sample==name,],
    pdf.file.name=NULL,
    pdf.width=7,
    pdf.height=7,
    col="grey",
    border="grey",
    breaks=50)
  mtext(name, outer=TRUE,  cex=1, line=-1)
}
plotBarcode(
  main.title = "all",
  obj=x.sp,
  pdf.file.name=NULL,
  pdf.width=7,
  pdf.height=7,
  col="grey",
  border="grey",
  breaks=50)
mtext(name, outer=TRUE,  cex=1, line=-1)
dev.off()

###########################################################################
# Cross-reference the filtered barcodes with the snap object
barcode.shared = intersect(x.sp@barcode, barcodes$barcode);
x.sp = x.sp[match(barcode.shared, x.sp@barcode),];
barcodes = barcodes[match(barcode.shared, barcodes$barcode),];
x.sp@metaData = barcodes;
summary(x.sp)

# Matrix binarisation
x.sp = addBmatToSnap(x.sp, bin.size=5000)
x.sp = makeBinary(x.sp, mat="bmat")

# Adding the promoter ratio to the barcodes file
promoter.df = read.table("promoter.bed")
promoter.gr = GRanges(promoter.df[,1], IRanges(promoter.df[,2], promoter.df[,3]))
ov = findOverlaps(x.sp@feature, promoter.gr)
x.sp@feature

# Determine promoter ratio and exclude overlapping bins
idy = queryHits(ov)
idy = unique(idy)

log_cov = log10(SnapATAC::rowSums(x.sp, mat="bmat")+1)

x.sp@metaData$log_cov = log_cov
promoter_ratio = Matrix::rowSums(x.sp@bmat[,idy]) / Matrix::rowSums(x.sp@bmat);

# Plot the log count vs promoter ratio before filtering
pdf(paste(resultsdir,'barcodeQC3_promcount_all.pdf',sep=""),width=6,height=6,paper='special')
plot(log_cov, promoter_ratio, cex=0.5, col="grey", xlab="log(count)", ylab="FIP Ratio", ylim=c(0,1 ))
dev.off()

idz = which(promoter_ratio > 0.2 & promoter_ratio < 0.8 & log_cov > 3.6)
x.sp = x.sp[idz,]
x.sp

promoter_ratio_filt = Matrix::rowSums(x.sp@bmat[,idy]) / Matrix::rowSums(x.sp@bmat);

# Plot the log count vs promoter ratio after filtering
pdf(paste(resultsdir,'barcodeQC4_promcount_all.pdf',sep=""),width=6,height=6,paper='special')
plot(x.sp@metaData$log_cov, promoter_ratio_filt, cex=0.5, col="grey", xlab="log(count)", ylab="FIP Ratio", ylim=c(0,1 ))
dev.off()

###########################################################################
# Bin filtering and blacklist exclusion for unwanted regions
system("wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz");
library(GenomicRanges);
black_list = read.table("hg38.blacklist.bed.gz")
black_list.gr = GRanges(
  black_list[,1],
  IRanges(black_list[,2], black_list[,3])
);
idy = queryHits(findOverlaps(x.sp@feature, black_list.gr));
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]};
x.sp

# Remove unwanted chromosomes
chr.exclude = seqlevels(x.sp@feature)[grep("random|chrM", seqlevels(x.sp@feature))];
idy = grep(paste(chr.exclude, collapse="|"), x.sp@feature);
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]};
x.sp

# Remove invariant features
bin.cov = log10(Matrix::colSums(x.sp@bmat)+1);
binhist <- hist(
  bin.cov[bin.cov > 0],
  xlab="log10(bin cov)",
  main="log10(Bin Cov)",
  col="lightblue",
  xlim=c(0, 5)
);

# Plotting the bin distribution
pdf(paste(resultsdir,'bins_all.pdf',sep=""),width=6,height=6,paper='special')
plot(binhist)
dev.off()

# Filter on the rowsums being > 1000
idx = which(Matrix::rowSums(x.sp@bmat) > 1000);
x.sp = x.sp[idx,];
x.sp

bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95);
idy = which(bin.cov <= bin.cutoff & bin.cov > 0);
x.sp = x.sp[, idy, mat="bmat"];
x.sp

# Removing low mapq cells based on downstream analysis
x.sp = x.sp[x.sp@metaData$lowmapq<5000,]
x.sp

#sum(x.sp@sample == "atac_lako1")

###########################################################################
# Determine the number of PCAs to use for generating the UMAP

# Setting core usage for the server below is intensive
SnowParam(workers = 2)

# Dimensionality reduction
x.sp = runDiffusionMaps(
  obj=x.sp,
  input.mat="bmat",
  num.eigs=30
)

pdf(paste(resultsdir,'/numPCAs_all.pdf',sep=""),width=6,height=6,paper='special')
plotDimReductPW(
  obj=x.sp,
  eigs.dims=1:30,
  point.size=0.3,
  point.color="grey",
  point.shape=19,
  point.alpha=0.6,
  down.sample=500,
  pdf.file.name=NULL,
  pdf.height=7,
  pdf.width=7
)
plotDimReductElbow(x.sp)
dev.off()

# Selecting a specific number of dimensions based on the elbow plot and PCAs
x.sp = runKNN(
  obj=x.sp,
  eigs.dims=1:22,
  k=10
)

###########################################################################
# Run clustering for the snap object

###########################################################################
# Use this block of code below if the Leiden clusttering package does not install correctly
#
# install.packages("leiden")
#  library("leiden")
#  library(reticulate)
#  library(igraph)
#  path_to_python <- "/vol/mbconda/julian/envs/kb_snapatac/bin/python"
#  use_python(path_to_python, required = T)
#  py_config()
###########################################################################

# Running the clustering with a specific resolution (unbiased)
x.sp@cluster
library(SnapATAC)
library("leiden")

x.sp=runCluster(
  obj=x.sp,
  resolution=0.3,
  tmp.folder=tempdir(),
  louvain.lib="leiden",
  seed.use=10
)

x.sp@metaData$cluster = x.sp@cluster
x.sp@cluster

###########################################################################
# Visualization upon the UMAP
x.sp = runViz(
  obj=x.sp,
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:10,
  method="umap",
  seed.use=10
)

pdf(paste(resultsdir,'/visualisation_UMAP_all.pdf',sep=""),width=6,height=6,paper='special')
plotViz(
  obj=x.sp,
  method='umap',
  point.size=0.1,
  point.shape=19,
  point.alpha=0.8,
  point.color=x.sp@cluster,
  text.add=TRUE,
  text.size=1.2,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=10000,
  pdf.file.name=NULL,
  pdf.width=7,
  pdf.height=7
)

plotViz(
  obj= x.sp,
  method="umap",
  main="Cluster",
  point.color=x.sp@cluster,
  point.size=0.1,
  point.shape=19,
  text.add=TRUE,
  text.size=0.5,
  text.color="black",
  down.sample=10000,
  legend.add=FALSE
)

plotViz(
  obj= x.sp,
  method="umap",
  main="Sample",
  point.size=0.1,
  point.shape=19,
  point.color=x.sp@sample,
  text.add=FALSE,
  text.size=0.5,
  text.color="black",
  down.sample=10000,
  legend.add=TRUE
)
dev.off()

# Checking the contribution of each replicate to the snap object
count(x.sp@sample == "atac_lako1")
count(x.sp@sample == "atac_lako2")
count(x.sp@sample == "atac_lako3")
count(x.sp@sample == "atac_lako4")

# Calculate the ensemble signals for each cluster
ensemble.ls = lapply(split(seq(length(x.sp@cluster)), x.sp@cluster), function(x){
  SnapATAC::colMeans(x.sp[x,], mat="bmat");
})

# Cluster using 1-cor as distance and generating the clustree
hc = hclust(as.dist(1 - cor(t(do.call(rbind, ensemble.ls)))), method="ward.D2")
plot(hc, hang=-1, xlab="")

pdf(paste(resultsdir,'cluster_tree_all_unannotated.pdf',sep=""),width=6,height=6,paper='special')
plot(hc, hang=-1, xlab="")
dev.off()

##########################################################################################################3
# Annotating genes for each cluster

# Importing the gene annotation file
gene_df <- read.table("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/gencode.v30.gene.annotation.gtf")
genes.gr = GRanges(gene_df[,1],
                   IRanges(gene_df[,2], gene_df[,3]),
                   name=gene_df[,4]
)
genes.gr2  <- genes.gr[!duplicated(genes.gr$name),]

# For a limited amount of markers put them in this vector
marker.genes = c(
  "CPVL","PAX6","TP63","ACKR1","POSTN",'KERA')

genes.sel.gr2 <- genes.gr2[which(genes.gr2$name %in% marker.genes)]

x.sp4 = createGmatFromMat(
  obj=x.sp,
  genes= genes.sel.gr2,
  do.par=TRUE
)

x.sp4 = scaleCountMatrix(
  obj=x.sp4,
  cov=as.integer(SnapATAC::rowSums(x.sp, mat="bmat")),
  mat="gmat",
  method="logRPM")

x.sp4 = runMagic(
   obj=x.sp4,
   input.mat="gmat",
   step.size=3
 )

pdf(paste(resultsdir,'Marker_genes_sample_all.pdf',sep=""),width=6,height=6,paper='special')
plotViz(
  obj=x.sp4,
  method="umap",
  main="lako_atac_all",
  point.color=x.sp4@cluster,
  point.size=0.5,
  point.shape=19,
  point.alpha=0.8,
  text.add=FALSE,
  text.size=1.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=10000,
  legend.add=TRUE
)

for(i in 1:length(marker.genes)){
  plotFeatureSingle(
    obj=x.sp4,
    feature.value=x.sp4@gmat[, marker.genes[i]],
    method="umap",
    main=marker.genes[i],
    point.size=0.1,
    point.shape=19,
    down.sample=10000,
    quantiles=c(0, 1)
  )}
dev.off()

##############################################################################
# Annotation for all marker genes of interest

# Loading in the marker gene .csv file
marker_gene_file <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/LSC_Marker_Genes3.csv"
marker_genes_df <- read.table(marker_gene_file, header = TRUE, sep = ',',comment.char = "#", stringsAsFactors = F)

# Add an empty list where all markers can be stored
mylist <- c()

for (paper in unique(marker_genes_df$abreviation)){
  marker_dir <- paste0(paste0(resultsdir, '/marker_genes/'),paper)
  print(marker_dir)
  system(paste("mkdir -p ", marker_dir))
  sub_marker_df <- marker_genes_df[marker_genes_df$abreviation == paper,]

  for (cell_type in unique(sub_marker_df$name_plot)){
    marker_genes <- sub_marker_df[sub_marker_df$name_plot == cell_type,]$marker_genes
    marker_genes <- unlist(strsplit(marker_genes, ','))

    n <- length(marker_genes)
    k <- 10 # specifying a maximal length of target genes for each celltype
    marker_gene_sets <- split(marker_genes, rep(1:ceiling(n/k), each=k)[1:n])

    for (genes in  marker_gene_sets){
      print(genes)

      for (z in genes){
        print(z)
        mylist <- c(mylist, z)
      }
    }
  }
}

# Remove duplicates from the list
mylist <- unique(mylist)

# Adding the gene matrix for all genes of interest to the snap object
genes.sel.gr2 <- genes.gr2[which(genes.gr2$name %in% mylist)]
x.sp5 = createGmatFromMat(
  obj=x.sp,
  genes= genes.sel.gr2,
  do.par=TRUE
)

x.sp5 = scaleCountMatrix(
  obj=x.sp5,
  cov=as.integer(SnapATAC::rowSums(x.sp, mat="bmat")),
  mat="gmat",
  method="logRPM")

# RunMagic for imputation of data in sparce surrounding datapoints
x.sp5 = runMagic(
  obj=x.sp5,
  input.mat="gmat",
  step.size=3
)

##############################################################################
# Plotting all marker genes for each cell upon the UMAP

# Specifying a low number of cores due to high computational demand
SnowParam(workers = 2)

for (paper in unique(marker_genes_df$abreviation)){
  marker_dir <- paste0(paste0(resultsdir, '/marker_genes_magic/'),paper)
  print(marker_dir)
  system(paste("mkdir -p ", marker_dir))
  sub_marker_df <- marker_genes_df[marker_genes_df$abreviation == paper,]

  for (cell_type in unique(sub_marker_df$name_plot)){
    print(cell_type)
    marker_genes <- sub_marker_df[sub_marker_df$name_plot == cell_type,]$marker_genes
    marker_genes <- unlist(strsplit(marker_genes, ','))

    n <- length(marker_genes)
    k <- 10
    marker_gene_sets <- split(marker_genes, rep(1:ceiling(n/k), each=k)[1:n])


    pdf(paste0(paste(paste(marker_dir, cell_type, sep = '/'), sep = '_'),'.pdf'), width = 10, height = 6)

    for (genes in  marker_gene_sets){
      print(genes)

      for (z in genes){
        print(z)
        plotFeatureSingle(
           obj=x.sp5,
           feature.value=x.sp5@gmat[,z],
           method="umap",
           main=z,
           point.size=0.5,
           point.shape=19,
           down.sample=10000,
           quantiles=c(0.01, 0.99)
         )
        boxPlotFeature(x.sp5,feature = x.sp5@gmat[,z],main = z,ylab = "logRPM norm accesibility")
      }
    }
    dev.off()}
}

##############################################################################
# Save snap object between
saveRDS(x.sp,file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/lakoatacxspall.rds")
x.sp <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/lakoatacxspall.rds")

##############################################################################
# Plotting QC value of each cluster to check if no unwanted clustering is occurring

pdf(paste(resultsdir,'QC clusters.pdf',sep=""),width=6,height=6,paper='special')
data <- data.frame(name= as.factor(x.sp@metaData$cluster), value= x.sp@metaData$log_cov)
ggplot(data, aes(x=as.factor(data$name), y= data$value)) +
geom_boxplot() +
labs(x= "cluster",y="log10 coverage")

data <- data.frame(name= as.factor(x.sp@metaData$cluster), value= x.sp@metaData$logUMI)
ggplot(data, aes(x=as.factor(data$name), y= data$value)) +
geom_boxplot()+
labs(x= "cluster",y="log10 UMI count")

data <- data.frame(name= as.factor(x.sp@metaData$cluster), value= x.sp@metaData$mitochondrial)
ggplot(data, aes(x=as.factor(data$name), y= data$value)) +
geom_boxplot()+
labs(x= "cluster",y="counts mitochondrial")

data <- data.frame(name= as.factor(x.sp@metaData$cluster), value= x.sp@metaData$passed_filters)
ggplot(data, aes(x=as.factor(data$name), y= data$value)) +
geom_boxplot()+
labs(x= "cluster",y="value passed filters")
dev.off()

data <- data.frame(name= as.factor(x.sp@metaData$cluster), value= x.sp@metaData$duplicate)
ggplot(data, aes(x=as.factor(data$name), y= data$value)) +
  geom_boxplot()+
  labs(x= "cluster",y="duplicate value")

plotFeatureSingle(
  obj=x.sp,
  feature.value=x.sp@metaData$lowmapq,
  method="umap",
  main="QC low mapq",
  point.size=0.5,
  point.shape=19,
  down.sample=10000,
  quantiles=c(0.01, 0.99)
)

plotFeatureSingle(
  obj=x.sp,
  feature.value=x.sp@metaData$duplicate,
  method="umap",
  main="QC duplicate",
  point.size=0.5,
  point.shape=19,
  down.sample=10000,
  quantiles=c(0.01, 0.99)
)

plotFeatureSingle(
  obj=x.sp,
  feature.value=x.sp@metaData$passed_filters,
  method="umap",
  main="QC passed filters",
  point.size=0.5,
  point.shape=19,
  down.sample=10000,
  quantiles=c(0.01, 0.99)
)

plotFeatureSingle(
  obj=x.sp,
  feature.value=x.sp@metaData$chimeric,
  method="umap",
  main="QC chimeric",
  point.size=0.5,
  point.shape=19,
  down.sample=10000,
  quantiles=c(0.01, 0.99)
)

plotFeatureSingle(
  obj=x.sp,
  feature.value=x.sp@metaData$unmapped,
  method="umap",
  main="QC unmapped",
  point.size=0.5,
  point.shape=19,
  down.sample=10000,
  quantiles=c(0.01, 0.99)
)

plotFeatureSingle(
  obj=x.sp,
  feature.value=x.sp@metaData$logUMI,
  method="umap",
  main="QC log unique counts",
  point.size=0.5,
  point.shape=19,
  down.sample=10000,
  quantiles=c(0.01, 0.99)
)

plotFeatureSingle(
  obj=x.sp,
  feature.value=x.sp@metaData$log_cov,
  method="umap",
  main="QC log coverage",
  point.size=0.5,
  point.shape=19,
  down.sample=10000,
  quantiles=c(0.01, 0.99)
)
dev.off()

##############################################################################
# Recluster according to accessibility of marker genes
current.cluster.ids <-  c("1","2","3","4","5","6")
new.cluster.ids <- c( "cell1","cell2","cell3","cell4","cell5","cell6")

x.sp@cluster <- plyr::mapvalues(
   x = x.sp@cluster,
   from = current.cluster.ids,
   to = new.cluster.ids
)

x.sp@cluster

# Own annotation of the clusters
current.cluster.ids <-  c("cell1","cell2","cell3","cell4","cell5","cell6")
new.cluster.ids <- c( "cell1","cornealcells","cell3","cell4","stromalcells","LPCs")
#
x.sp@cluster <- plyr::mapvalues(
  x = x.sp@cluster,
  from = current.cluster.ids,
  to = new.cluster.ids
)

x.sp@cluster

pdf(paste(resultsdir,'visualisation_annotated_UMAP_all.pdf',sep=""),width=6,height=6,paper='special')
plotViz(
  obj= x.sp,
  method="umap",
  main="Cluster",
  point.color=x.sp@cluster,
  point.size=0.1,
  point.shape=19,
  text.add=TRUE,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=FALSE
)
dev.off()

##############################################################################
# Using MACS2 for peak calling within SnapATAC (credits to Rebecca and Jelle)

# Setting cluster dirs
dir.create(paste0(resultsdir, "clusters/"))
setwd(paste0(resultsdir, "clusters/"))

# Get amount of cells per cluster:
summary(x.sp@cluster)

# Call peaks for all cluster with more than 100 cells
clusters.sel = names(table(x.sp@cluster))[which(table(x.sp@cluster) > 100)]

# Checking the OS path of snaptools and MACS2
system("which snaptools")
system("which macs2")

# Setting the cores
SnowParam(workers = 2)

peaks.ls = mclapply(seq(clusters.sel), function(i){
  print(clusters.sel[i])
  peaks = runMACS(
    obj=x.sp[which(x.sp@cluster==clusters.sel[i]),],
    output.prefix=paste0("LAKO.", gsub(" ", "_", clusters.sel)[i]),
    path.to.snaptools="/vol/mbconda/julian/envs/kb_snapatac/bin/snaptools",
    path.to.macs="/vol/mbconda/julian/envs/kb_snapatac/bin/macs2",
    gsize="hs", # mm, hs, etc
    buffer.size=500,
    num.cores=1,
    macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
    tmp.folder="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/R_peaks_scATAC"
  )
  peaks
})

Sys.setenv(bdg_location=paste0("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/peaks_scATAC_20210603/"))

# note: MACS2 options are based on the discordance between the tn5 insertion site and the reads themselves
# note 2: I personally do not recommend doing this step in R, after splitting the bamfiles based on barcode
# you can do this all in bash (see Github documentation) at a later stage

##############################################################################
# Go to bash, see the documentation on Github
# for the correct commands to generate bigwigs

##############################################################################
# Extract binarized counts
binarized_counts <- Matrix::colSums(x.sp@bmat)

# Scale scoring for bed file using the maximum n of reads in a bin and create a name for the bins
scaled_binarized_counts <- binarized_counts * (1000/max(binarized_counts))
names_binarized_counts <- paste0(binarized_counts, "of", length(x.sp@barcode))

# Create bed datafame with regions of bins extracted from x.sp@feature and scores based on sum of all cells (colSums of x.sp@bmat)
bed_df <- data.frame(seqnames=seqnames(x.sp@feature),
                     starts=start(x.sp@feature)-1,
                     ends=end(x.sp@feature),
                     names=names_binarized_counts,
                     scores=scaled_binarized_counts,
                     strands=strand(x.sp@feature))

# Rename strands
levels(bed_df$strands)[levels(bed_df$strands)=="*"] <- "."

# Write out table
write.table(bed_df, file=paste0("LAKO.", gsub(" ", "_", clusters.sel)[i], "_binary_matrix.bed"), quote=FALSE, sep= "\t", row.names=FALSE, col.names=FALSE)
#```

##############################################################################
# Make bed files with bins for separate clusters

clusters.sel = names(table(x.sp@cluster))

# Create bed_df for every cluster
for (cluster_name in seq(clusters.sel)) {
  binarized_counts_cluster <- Matrix::colSums(x.sp@bmat[which(x.sp@cluster==clusters.sel[cluster_name]),])

  # Scale scoring for bed file using the maximum n of reads in a bin and create a name for the bins
  scaled_binarized_counts <- binarized_counts_cluster * (1000/max(binarized_counts))
  names_binarized_counts <- paste0(binarized_counts_cluster, "of", length(x.sp@barcode))

  # Create bed datafame with regions of bins extracted from x.sp@feature and scores based on sum of all cells (colSums of x.sp@bmat)
  bed_df <- data.frame(seqnames=seqnames(x.sp@feature),
                       starts=start(x.sp@feature)-1,
                       ends=end(x.sp@feature),
                       names=names_binarized_counts,
                       scores=scaled_binarized_counts,
                       strands=strand(x.sp@feature))
  # Rename strands
  levels(bed_df$strands)[levels(bed_df$strands)=="*"] <- "."

  # Write out table
  write.table(format(bed_df, scientific=FALSE), file=paste0(resultsdir, "/", "_cluster", cluster_name, "_binary_matrix.bed"), quote=FALSE, sep= "\t", row.names=FALSE, col.names=FALSE)
}

##############################################################################
# Peak matrix generation

#Create GRanges object with all peak locations from clusters
peaks.names = system("ls | grep narrowPeak", intern = TRUE)
peak.gr.ls = lapply(peaks.names, function(x){
  peak.df = read.table(x)
  GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
})
peak.gr = suppressWarnings(reduce(Reduce(c, peak.gr.ls)))

# Create a cell-by-peak matrix
peaks.df = as.data.frame(peak.gr)[,1:3]

write.table(peaks.df, file=paste0(resultsdir,"/", "_peaks_combined.bed"), append = FALSE, quote=FALSE, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")

# Delete weird characters
peaks.df$seqnames<-gsub("b'","",as.character(peaks.df$seqnames))
peaks.df
peaks.df$seqnames<-gsub("'","",as.character(peaks.df$seqnames))
peaks.df

write.table(peaks.df, file=paste0(resultsdir,"/", "_peaks_combined_nochar.bed"), append = FALSE, quote=FALSE, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")

##############################################################################
# Save the current snap object into an RDS
saveRDS(x.sp, file=paste0(resultsdir, ".snap.rds"))

# Add pmat to snaps in bash
# See Github documentation

##############################################################################
# Add pmat to snap within R
x.sp = readRDS(x.sp, file= "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/20210518.snap.rds")
x.sp = addPmatToSnap(x.sp,do.par = T,num.cores = 10)
x.sp@pmat

# Create a directory and save unbinarized cell by peak matrix for scRNA integration (part3)
scRNA_int_dir <- paste0(resultsdir, "scRNA_integration/")
dir.create(scRNA_int_dir)

# Get peak matrix and add peak locations as column names, transpose
peak_matrix <- x.sp@pmat
peak_matrix@Dimnames[2] =list(x.sp@peak@elementMetadata@listData$name)
peak_matrix <- t(peak_matrix)

# Save for later use with Seurat integration with scRNA-seq
saveRDS(peak_matrix, file=paste0(scRNA_int_dir, "cell_peak_matrix.rds"))

# Extract metadata and save for later use
ATAC_metadata <- x.sp@metaData
write.csv(ATAC_metadata, paste0(scRNA_int_dir, "ATAC_metadata.csv"), row.names=FALSE)

# Binarize pmat for further usage in snap object
x.sp = makeBinary(x.sp, mat="pmat")

# Define the differential accessible regions for the clusters
idy.ls = lapply(levels(x.sp@cluster), function(cluster_i){
  DARs = findDAR(
    obj=x.sp,
    input.mat="pmat",
    cluster.pos=cluster_i,
    cluster.neg=NULL,
    cluster.neg.method="knn",
    test.method="exactTest",
    bcv=0.01, #Only one plate, so bcv = 0.01
    seed.use=10
  )
  DARs$FDR = p.adjust(DARs$PValue, method="BH")
  idy = which(DARs$FDR < 0.05 & DARs$logFC > 0)
  #For small clusters, we here select the top 2000 peaks
  if((x=length(idy)) < 2000L){
    PValues = DARs$PValue
    PValues[DARs$logFC < 0] = 1
    idy = order(PValues, decreasing=FALSE)[1:2000]
    rm(PValues); # free memory
  }
  idy
})

covs = Matrix::rowSums(x.sp@pmat)
names(idy.ls) = levels(x.sp@cluster)

# Plotting the differential regions per cluster
pdf(paste(resultsdir,'DAR_clusters.pdf',sep=""),width=6,height=6,paper='special')
for(cluster_i in levels(x.sp@cluster)){
  idy = idy.ls[[cluster_i]]
  vals = Matrix::rowSums(x.sp@pmat[,idy]) / covs
  vals.zscore = (vals - mean(vals)) / sd(vals)
  plotFeatureSingle(
    obj=x.sp,
    feature.value=vals.zscore,
    method="umap",
    main=paste0("Cluster ", cluster_i, " z-score"),
    point.size=0.5,
    point.shape=19,
    down.sample=5000,
    quantiles=c(0.01, 0.99)
  )
}
dev.off()

#Create a list with all the directories and file names (OBSOLETE?)
folder_info <- list("working_dir" = workdir,
                    "results_folder" = resultsdir)

#Save the snap object for further analysis
saveRDS(x.sp, file=paste0(resultsdir, ".snap.rds"))

#Save the folder_info list for further analysis
saveRDS(folder_info, paste0(workdir, "_folder_locations.rds"))

##############################################################################
# Imputation of clusters in scATAC-seq from scRNA-seq
# Note: if the atac object has been generated already with variable features
# and you only want new clusters then proceed to STEP 3 directly

# STEP 1: Go to the kb_scrna_R_seurat4 environment for determining clusters from scRNA-seq
library(Seurat)

meta.rna = readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220228/4datasets_annotated_joined.rds")
meta.rna2 = readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAannotated.rds")

# Finding Variable Features only needs to be done on un-integrated data; here only the lako dataset
DefaultAssay(object = meta.rna2) <- "RNA"
variable.genes = VariableFeatures(object = meta.rna2,assay = "RNA")
saveRDS(variable.genes,"/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/variable.genes.lako.rds")

# Variable genes can be found as features of integrated data meta-atlas
meta.rna <- FindVariableFeatures(meta.rna)
variable.genes <- meta.rna@assays$RNA@var.features
saveRDS(variable.genes,"/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/variable.genes.rds")

##############################################################################
# STEP2: Return to kb_snapatac env
library(SnapATAC)
library(GenomicRanges)

# loading in the snap object and adding all variable features for the gene matrix
variable.genes <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/variable.genes.rds")
x.sp = readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/20210609.snap.rds")

genes.df <- read.table("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/gencode.v30.gene.annotation.gtf") #already ran this
genes.gr3 = GRanges(genes.df[,1], IRanges(genes.df[,2], genes.df[,3]), name=genes.df[,4])
genes.sel.gr3 = genes.gr3[which(genes.gr3$name %in% variable.genes)]

x.sp3 = createGmatFromMat(
  obj=x.sp,
  input.mat="bmat",
  genes=genes.sel.gr3,
  do.par=TRUE#,
  #num.cores=10
)

x.sp3

# Convert snap to seurat object
meta_all.atac <- snapToSeurat(obj=x.sp3,
                             eigs.dims=1:10,
                             norm=TRUE,
                             scale=TRUE
)

# Save the seurat object as a rds file
saveRDS(meta_all.atac,file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/atac_meta_all.rds")

##############################################################################
# Move again to the different environment kb_scrna_R_seurat4
# STEP 3: imputing the clusters

library(Seurat)
#lako.rna = readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAannotated.rds")
#lako.rna = readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAannotated_4cells.rds")

meta.rna = readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220228/4datasets_annotated_joined.rds")

meta_all.atac = readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/atac_meta_all.rds")

meta_all.atac@meta.data$cluster

#devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories()) # update none
#install.packages("ArchR")
#library(ArchR)

# projHeme2 <- addGeneIntegrationMatrix(
#   ArchRProj = projHeme2,
#   useMatrix = "GeneScoreMatrix",
#   matrixName = "GeneIntegrationMatrix",
#   reducedDims = "IterativeLSI",
#   seRNA = meta.rna,
#   addToArrow = FALSE,
#   groupRNA = "BioClassification",
#   nameCell = "predictedCell_Un",
#   nameGroup = "predictedGroup_Un",
#   nameScore = "predictedScore_Un"
# )

# Label these cluster after unconstrained integration
meta.rna$Condition <- meta.rna@meta.data$Condition
as.factor(meta.rna$Condition)

meta.rna2 <-subset(x=meta.rna, subset = (Condition== "Ga" |Condition == "Ca"|Condition == "Li"),invert = TRUE)
DefaultAssay(object = meta.rna2) <- "RNA"

variable.genes <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/variable.genes.rds")

transfer.anchors <- FindTransferAnchors(
  reference = meta.rna,
  query = meta_all.atac,
  features = variable.genes,
  reference.assay = "RNA",
  query.assay = "ACTIVITY",
  reduction = "cca"
)

meta.rna$scVI_label

celltype.predictions <- TransferData(
  anchorset = transfer.anchors,
  refdata = meta.rna$scVI_label,
  weight.reduction = meta_all.atac[["SnapATAC"]],
  dims = 1:10
)

saveRDS(celltype.predictions,file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/celltypes_predicitons_all_meta.rds")

meta_all.atac@meta.data$predict.id = celltype.predictions$predicted.id;
meta_all.atac@meta.data$predict.max.score = apply(celltype.predictions[,-1], 1, max);

#saveRDS(lakoall.atac,file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/lakoatacall.rds")
saveRDS(meta_all.atac,file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/meta_all_atac_predictions.rds")

# refdata <- GetAssayData(
#    object = lako.rna,
#    assay = "RNA",
#    slot = "data"
#  )
# 
# imputation <- TransferData(
#    anchorset = transfer.anchors,
# refdata = refdata,
#    weight.reduction = lakoall.atac[["SnapATAC"]],
#    dims = 1:10
#  )

#saveRDS(imputation,file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/imputation.rds")
#saveRDS(imputation,file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/imputation_4cells.rds")


##############################################################################
# STEP4: Go back to kb_snapatac environment
library(SnapATAC)
library(Seurat)
library(Matrix)
library(GenomicRanges)

# reload the gmat upon the original x.sp (not the seurat converted file)
x.sp = readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/20210609.snap.rds")

variable.genes <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/variable.genes.rds")

genes.df <- read.table("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/gencode.v30.gene.annotation.gtf") # Run again if nessesary
genes.gr3 = GRanges(genes.df[,1], IRanges(genes.df[,2], genes.df[,3]), name=genes.df[,4])
genes.sel.gr3 = genes.gr3[which(genes.gr3$name %in% variable.genes)]

x.sp = createGmatFromMat(
  obj=x.sp,
  input.mat="bmat",
  genes=genes.sel.gr3,
  do.par=TRUE#,
  #num.cores=10
)

#celltype.predictions <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/celltypes_predicitons_all_meta.rds")
#celltype.predictions <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/celltypes_predicitons_all_meta_constrained.rds")
celltype.predictions <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/celltypes_predicitons_all_meta.rds")

x.sp@metaData$predict.id = celltype.predictions$predicted.id;
x.sp@metaData$predict.max.score = apply(celltype.predictions[,-1], 1, max);

###########################################################################
# Storing results
workdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scATAC/"

# setting up results directory
dateoftoday <- gsub("-", "", as.character(Sys.Date()))
resultsdir <- paste0(workdir, dateoftoday)
system(paste("mkdir -p ", resultsdir))
###########################################################################

pdf(paste(resultsdir,'/predictscore_all_lako.pdf',sep=""),width=18,height=18,paper='special')
hist(
  x.sp@metaData$predict.max.score,
  xlab="prediction score",
  col="lightblue",
  xlim=c(0, 1),
  main="Meta 10X"
)
abline(v=0.5, col="red", lwd=2, lty=2)
table(x.sp@metaData$predict.max.score > 0.5)
dev.off()

x.sp2 = x.sp[x.sp@metaData$predict.max.score > 0.5,];
x.sp2@metaData
unique(x.sp2@metaData[,"predict.id"])

# Saving imputation files max score 0.5
#saveRDS(x.sp,file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/lakopredicted_05.rds")
#saveRDS(x.sp2,file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/meta_all_05.rds")
saveRDS(x.sp2,file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/meta_all_05.rds")

#x.sp <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/lakopredicted_05.rds")
#x.sp2 <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/lakopredicted_4cells_05.rds")

x.sp3 = x.sp[x.sp@metaData$predict.max.score > 0.4,]
x.sp3@metaData
x.sp3@metaData[,"predict.id"]

# Saving imputation files max score 0.2
#saveRDS(x.sp,file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/lakopredicted_02.rds")
#saveRDS(x.sp3,file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/meta_03.rds")
saveRDS(x.sp3,file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/meta_04.rds")

#x.sp <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/lakopredicted_02.rds")
#x.sp3 <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/meta_03.rds")

library(umap)

pdf(paste(resultsdir,"/", 'visualisation_imputated_UMAP_final.pdf',sep=""),width=10,height=8,paper='special')
plotViz(obj=x.sp, method="umap",text.size = 1.5, point.color=as.factor(x.sp@metaData$predict.id), legend.add = T,text.add = F)
plotViz(obj=x.sp2, method="umap",text.size = 1.5, point.color=x.sp2@metaData$predict.id, legend.add = T)
plotViz(obj=x.sp3, method="umap",text.size = 1.5, point.color=x.sp3@metaData$predict.id, legend.add = T)
dev.off()

###### colors?

x.sp = readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/meta_04.rds")

cols <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220310/unannotated_meta_colors.tsv",sep="\t")

cluster_order <- c('Cj','CE','LE','LESC','CSSC','CF','LSC',
                   'SK','LE','Cj','SK','TSK','CE','CF',
                   'Mel','EC','Ves','IC','CE','MF','CDH19+')

# Generate new cluster labeling:
current.cluster.ids <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)

new.cluster.ids <- c('Cj','CE','LE','LESC','CSSC','CF','LSC',
                     'SK','LE','Cj','SK','TSK','CE','CF',
                     'Mel','EC','Ves','IC','CE','MF','CDH19+')

cols$V1 <- plyr::mapvalues(x = as.factor(cols$V1), from = current.cluster.ids, to = new.cluster.ids)
cols$V1 <- as.factor(cols$V1)

levels(as.factor(x.sp@metaData$predict.id))

library(SnapATAC)
current.cluster.ids <- x.sp@metaData$predict.id


colVector <- cols$V2
cols<-cols[order(match(cols$V1,levels(as.factor(x.sp@metaData$predict.id)) )), , drop = FALSE]
cols$V1 <- as.vector(cols$V1)
cols <- cols[!duplicated(cols$V1),]
x.sp

x.sp2 <- x.sp[x.sp@metaData$predict.max.score > 0.4 & x.sp@metaData$predict.id%in%clusters.sel,]

# Colors altered in Illustrator; no built in function in snapatac
pdf(paste(resultsdir,"/", 'visualisation_imputated_UMAP_subset.pdf',sep=""),width=10,height=8,paper='special')
plotViz(obj=x.sp, method="umap",text.size = 1.5, point.color=as.factor(x.sp@metaData$predict.id), legend.add = T,text.add = F)
plotViz(obj=x.sp2, method="umap",text.size = 1.5, point.color=as.factor(x.sp2@metaData$predict.id), legend.add = T,text.add = F)
dev.off()

##############################################################################
# generating bigwig tracks from imputated data
dir.create(paste0(resultsdir, "clustersimputated/"))
setwd(paste0(resultsdir, "clustersimputated/"))

# call peaks for all cluster with more than 100 cells & a score higher than 0.4
x.sp <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/meta_04.rds")
clusters.sel = names(table(x.sp@metaData$predict.id))[which(table(x.sp@metaData$predict.id) > 100)]
clusters.sel

# Checking the OS path of snaptools and MACS2
system("which snaptools")
system("which macs2")


######################################################
# PEAK CALLING DONE WITH snakefile
# Setting the cores
SnowParam(workers = 2)

peaks.ls = mclapply(seq(clusters.sel), function(i){
  print(clusters.sel[i])
  peaks = runMACS(
    obj=x.sp3[which(x.sp3@metaData$predict.id==clusters.sel[i]),],
    output.prefix=paste0("LAKOIMP02.", gsub(" ", "_", clusters.sel)[i]),
    path.to.snaptools="/vol/mbconda/julian/envs/kb_snapatac/bin/snaptools",
    path.to.macs="/vol/mbconda/julian/envs/kb_snapatac/bin/macs2",
    gsize="hs", # mm, hs, etc
    buffer.size=500,
    num.cores=1,
    macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
    tmp.folder="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/R_peaks_scATAC"
  )
  peaks
})

Sys.setenv(bdg_location=paste0(resultsdir))

##############################################################################
# See bash for bigwig generation (Github)

# Making the bed dataframe
binarized_counts <- Matrix::colSums(x.sp3@bmat)

#Scale scoring for bed file using the maximum n of reads in a bin and create a name for the bins
scaled_binarized_counts <- binarized_counts * (1000/max(binarized_counts))
names_binarized_counts <- paste0(binarized_counts, "of", length(x.sp@barcode))

#Create bed datafame with regions of bins extracted from x.sp@feature and scores based on sum of all cells (colSums of x.sp@bmat)
bed_df <- data.frame(seqnames=seqnames(x.sp3@feature),
                     starts=start(x.sp3@feature)-1,
                     ends=end(x.sp3@feature),
                     names=names_binarized_counts,
                     scores=scaled_binarized_counts,
                     strands=strand(x.sp3@feature))

#Rename strands
levels(bed_df$strands)[levels(bed_df$strands)=="*"] <- "."

# Write out table
for (i in clusters.sel){
write.table(bed_df, file=paste0("LAKOIMP02.", gsub(" ", "_", clusters.sel[i]), "_binary_matrix.bed"), quote=FALSE, sep= "\t", row.names=FALSE, col.names=FALSE)
}

# Create bed_df for every cluster
for (cluster_name in seq(clusters.sel)) {
  binarized_counts_cluster <- Matrix::colSums(x.sp@bmat[which(x.sp@cluster==clusters.sel[cluster_name]),])

  # Scale scoring for bed file using the maximum n of reads in a bin and create a name for the bins
  scaled_binarized_counts <- binarized_counts_cluster * (1000/max(binarized_counts))
  names_binarized_counts <- paste0(binarized_counts_cluster, "of", length(x.sp@barcode))

  # Create bed datafame with regions of bins extracted from x.sp@feature and scores based on sum of all cells (colSums of x.sp@bmat)
  bed_df <- data.frame(seqnames=seqnames(x.sp3@feature),
                       starts=start(x.sp3@feature)-1,
                       ends=end(x.sp3@feature),
                       names=names_binarized_counts,
                       scores=scaled_binarized_counts,
                       strands=strand(x.sp3@feature))
  # Rename strands
  levels(bed_df$strands)[levels(bed_df$strands)=="*"] <- "."

  # Write out table
  write.table(format(bed_df, scientific=FALSE), file=paste0(resultsdir, "_cluster", cluster_name, "_binary_matrix.bed"), quote=FALSE, sep= "\t", row.names=FALSE, col.names=FALSE)

}

# Create GRanges object with all peak locations from clusters
peaks.names = system("ls | grep narrowPeak", intern = TRUE)
peak.gr.ls = lapply(peaks.names, function(x){
  peak.df = read.table(x)
  GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
})
peak.gr = suppressWarnings(reduce(Reduce(c, peak.gr.ls)))

# Create a cell-by-peak matrix
peaks.df = as.data.frame(peak.gr)[,1:3]

write.table(peaks.df, file=paste0(resultsdir, "_peaks_combined.bed"), append = FALSE, quote=FALSE, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")

# Delete weird characters
peaks.df$seqnames<-gsub("b'","",as.character(peaks.df$seqnames))
peaks.df
peaks.df$seqnames<-gsub("'","",as.character(peaks.df$seqnames))
peaks.df

write.table(peaks.df, file=paste0(resultsdir, "_peaks_combined_nochar.bed"), append = FALSE, quote=FALSE, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")

#Save the current snap object into an RDS
x.sp3@feature # seems right
saveRDS(x.sp3, file=paste0(resultsdir, ".snap.rds"))

# Add Pmat after bigwig generation in bash
SnowParam(workers = 6)
x.sp3 = addPmatToSnap(x.sp3,do.par = T,num.cores = 6)
x.sp3@pmat

#saveRDS(x.sp3,file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/lakopredicted_wpmat.rds")
saveRDS(x.sp3,file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/lakopredicted_wpmat_4cells.rds")

x.sp3 <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/lakopredicted_wpmat.rds")
#x.sp3 <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/lakopredicted_wpmat_4cells.rds")

#Create directory and save unbinarized cell by peak matrix for scRNA integration (part3)
scRNA_int_dir <- paste0(resultsdir, "/scRNA_integration_all")
dir.create(scRNA_int_dir)

tracks_int_dir <- paste0(resultsdir, "/tracks_integration_all")
split_int_dir <- paste0("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/split_integration_all/")
dir.create(tracks_int_dir)

#Get peak matrix and add peak locations as column names, transpose
peak_matrix <- x.sp3@pmat
peak_matrix@Dimnames[2] =list(x.sp3@peak@elementMetadata@listData$name)
peak_matrix <- t(peak_matrix)

##############################################################################
# Binarize pmat for further usage in snap object
x.sp <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/meta_04.rds")
clusters.sel = names(table(x.sp@metaData$predict.id))[which(table(x.sp@metaData$predict.id) > 100)]
clusters.sel

x.sp3 <- x.sp

x.sp3 = makeBinary(x.sp3, mat="pmat")

# Add new labels for correct DAG calculation
x.sp3@metaData$predict.id

# current.cluster.ids <- c('CB','CjB',
#                          'CjS','LPCs','IC','LNPCs','Mel',
#                          'Ves','CSSCs',
#                          'LSC','FCECs')
# new.cluster.ids <- c('CB','CSB',
#                      'CjS','LPCs','IC','LNPCs','MEC',
#                      'Ves','CSSCs',
#                      'StC','FCECs')
# x.sp3@metaData$predict.id
# x.sp3@metaData$predict.id <- plyr::mapvalues(x = x.sp3@metaData$predict.id, from = current.cluster.ids, to = new.cluster.ids)

# For filtering the cells
x.sp3 <- x.sp3[x.sp3@metaData$predict.id %in% clusters.sel,]
x.sp4 <- x.sp3[x.sp3@metaData$predict.id != "CF" & x.sp3@metaData$predict.id != "CSSC",]

x.sp4 <- x.sp3
x.sp4@cluster <- as.factor(x.sp4@metaData$predict.id)
x.sp4@cluster
library(umap)


pdf(paste(resultsdir,"/", 'UMAP_final_100.pdf',sep=""),width=6,height=5,paper='special')
plotViz(obj=x.sp3, method="umap",text.add=FALSE,point.color=x.sp3@metaData$predict.id, legend.add = T)
plotViz(obj=x.sp4, method="umap",text.add=FALSE,point.color=x.sp4@metaData$predict.id, legend.add = T)
dev.off()

# Generating a cluster tree
ensemble.ls = lapply(split(seq(length(x.sp4@cluster)), x.sp4@cluster), function(x){
  SnapATAC::colMeans(x.sp4[x,], mat="bmat");
})

# Cluster using 1-cor as distance
hc = hclust(as.dist(1 - cor(t(do.call(rbind, ensemble.ls)))), method="ward.D2")
plot(hc, hang=-1, xlab="")

pdf(paste(resultsdir,'/cluster_tree_all_imputated_filt.pdf',sep=""),width=6,height=6,paper='special')
plot(hc, hang=-1, xlab="")
dev.off()

##############################################################################
# Save for later use
#saveRDS(peak_matrix, file=paste0(scRNA_int_dir, "cell_peak_matrix.rds"))

# Extract metadata for populations over 100 cells and save for later use
x.sp3 <- x.sp3[x.sp3@metaData$predict.id %in% names(table(x.sp3@metaData$predict.id))[which(table(x.sp3@metaData$predict.id) > 100)],]
ATAC_metadata <- x.sp3@metaData
write.csv(ATAC_metadata, paste0(scRNA_int_dir, "/ATAC_metadata.csv"), row.names=FALSE)

##############################################################################
# Generate barcode files for each cell population to be used in "bamslice" for each sample


x.sp3@sample

ATAC_metadata$sample <- x.sp3@sample
ATAC_metadata$sample

scANANSEdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/scANANSE/"

################################################################################
# write as function?

# Setting the scANANSE dir
path <- paste0(scANANSEdir,"ATAC_peaks")
system(paste("mkdir -p ", path))

# Extracting metadata from snap object
snap_obj <- x.sp3
ATAC_metadata <- snap_obj@metaData
write.csv(ATAC_metadata, paste0(path, "/ATAC_metadata.csv"), row.names=FALSE)

# Splitting the barcodes on sample and cell if samples are present
if (!is.null(unique(ATAC_metadata$sample))){
for (z in unique(ATAC_metadata$sample)){
  for (i in unique(ATAC_metadata$predict.id)){
    df <- ATAC_metadata[i==ATAC_metadata$predict.id & ATAC_metadata$sample == z,]
    write.table(df$barcode, paste0(path, "/", i,"_", z,"_atac_metadata_split.csv"),col.names = F, row.names=F, quote =F)
  }
}

# Generate the csv file with cells and barcodes from split
barcodes_csv<-c()
cells <- c()
for (z in unique(ATAC_metadata$sample)){
  for (i in unique(ATAC_metadata$predict.id)){
    cells <- c(cells,paste0(i,"_",z))
    barcodes_csv <- c(barcodes_csv,paste0(path, "/", i,"_",z,"_atac_metadata_split.csv"))
  }
}

df <- as.data.frame(x=barcodes_csv,row.names = cells)
df$library_id <- cells
df <- rev(df)

write.table(df, paste0(path, "/","ATAC_barcodes_split.csv"),col.names = T, row.names=F, quote =F,sep = ",")

# Split the dataset only on cell for bam file generation
for (i in unique(ATAC_metadata$predict.id)){
  df <- ATAC_metadata[i==ATAC_metadata$predict.id,]
  write.table(df$barcode, paste0(path, "/", i,"_atac_metadata_sample.csv"), row.names=F, quote =F,col.names = F)
}

# Generate the csv file with cells and barcodes
barcodes_csv<-c()
cells <- c()

for (i in unique(ATAC_metadata$predict.id)){
  cells <- c(cells,i)
  barcodes_csv <- c(barcodes_csv,paste0(path, "/", i,"_atac_metadata_sample.csv"))
}

df <- as.data.frame(x=barcodes_csv,row.names = cells)
df$library_id <- cells

write.table(df, paste0(path, "/","ATAC_barcodes.csv"),col.names = T, row.names=F, quote =F,sep = ",")

################################################################################
# Split the dataset into blob
ATAC_metadata <- x.sp4@metaData
unique(ATAC_metadata$predict.id)

# Normalization dataframe for cellcounts
dflakocounts <- NULL
for (i in unique(ATAC_metadata$predict.id)){
  dflakocounts[i] <- sum(ATAC_metadata$predict.id ==i)
}
dfcounts <- data.frame(dflakocounts)

dfcounts["ratio"] <- dfcounts$dflakocounts/median(dfcounts$dflakocounts)

# If ratio is below median, then exclude and use original number of cells
if (dfcounts[,"ratio"] <1){
  dfcounts[,"ratio"] <- 1
}

dfcounts$ratio[dfcounts$ratio < 1] <- 1
dfcounts[rownames(dfcounts) == i,"ratio"]

# Small function for random row selection
randomRows = function(df,n){
  return(df[sample(nrow(df),n),])
}

for (i in unique(ATAC_metadata$predict.id)){
  dflako <- ATAC_metadata[i!=ATAC_metadata$predict.id,]
  numberofcells <- dfcounts[rownames(dfcounts) == i,"dflakocounts"]/dfcounts[rownames(dfcounts) == i,"ratio"]
  dflako <- randomRows(dflako,as.integer(numberofcells))
  write.table(dflako$barcode, paste0(tracks_int_dir,"/", i,"_blob_norm_metadata.csv"), row.names=F,col.names = F,quote = F)
}

# Without the cell normalization
for (i in unique(ATAC_metadata$predict.id)){
  dflako <- ATAC_metadata[i!=ATAC_metadata$predict.id,]
  write.table(dflako$barcode, paste0(tracks_int_dir,"/", i,"_blob_metadata.csv"), row.names=F,col.names = F,quote = F)
}

########################################################################################
# Define the differentially accessible regions for the clusters
idy.ls = lapply(levels(x.sp3@cluster), function(cluster_i){
  DARs = findDAR(
    obj=x.sp3,
    input.mat="pmat",
    cluster.pos=cluster_i,
    cluster.neg=NULL,
    cluster.neg.method="knn",
    test.method="exactTest",
    bcv=0.01, #Only one plate, so bcv = 0.01
    seed.use=10
  )
  DARs$FDR = p.adjust(DARs$PValue, method="BH")
  idy = which(DARs$FDR < 0.05 & DARs$logFC > 0)
  #For small clusters, we here select the top 2000 peaks
  if((x=length(idy)) < 2000L){
    PValues = DARs$PValue
    PValues[DARs$logFC < 0] = 1
    idy = order(PValues, decreasing=FALSE)[1:2000]
    rm(PValues); # free memory
  }
  idy
})

covs = Matrix::rowSums(x.sp3@pmat)
names(idy.ls) = levels(x.sp3@cluster)

pdf(paste(resultsdir,'DAR_clusters.pdf',sep=""),width=6,height=6,paper='special')
for(cluster_i in levels(x.sp3@cluster)){
  idy = idy.ls[[cluster_i]]
  vals = Matrix::rowSums(x.sp3@pmat[,idy]) / covs
  vals.zscore = (vals - mean(vals)) / sd(vals)
  plotFeatureSingle(
    obj=x.sp3,
    feature.value=vals.zscore,
    method="umap",
    main=paste0(cluster_i, " z-score"),
    point.size=0.5,
    point.shape=19,
    down.sample=5000,
    quantiles=c(0.01, 0.99)
  )
}
dev.off()
