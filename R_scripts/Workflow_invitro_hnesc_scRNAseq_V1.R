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

#install.packages("ggalluvial")
library(ggalluvial)
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)

###########################################################################
# Storing results in directories of your choice
workdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_hnesc/"

# setting up results directory
dateoftoday <- gsub("-", "", as.character(Sys.Date()))
resultsdir <- paste0(workdir, dateoftoday)
system(paste("mkdir -p ", resultsdir))

# Loading in all functions
source("QC_seurobj_v1.R")

###########################################################################
# load in all datasets
path2 <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/ESC_vitro_sc/SRA_ESC_scrna/scRNA/"

nam <- Read10X(data.dir = path2)

nam2 <- CreateSeuratObject(counts = nam, project = paste0("hnESC"), min.cells = 3, min.features = 200)

# Perform quality control
seur_obj <- QCSeurObj(SeuratObject = nam2,pct_MT = 30)
rm(nam2)
seur_obj$"predicted.id" <- "hnESC"

###########################################################################
# Remove doublets

# hnESC
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
seu_hnESC <- seur_obj
sweep.res.list_hnESC <- paramSweep_v3(seu_hnESC, PCs = 1:30, sct = FALSE)
sweep.stats_hnESC <- summarizeSweep(sweep.res.list_hnESC, GT = FALSE)
bcmvn_hnESC <- find.pK(sweep.stats_hnESC)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(seu_hnESC$predicted.id)           ## ex: annotations <- seu_hnESC@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(seu_hnESC@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------

# Choose pK 0.06 based on the highest BC metrix
seu_hnESC <- doubletFinder_v3(seu_hnESC, PCs = 1:30, pN = 0.25, pK = 0.05, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
#seu_hnESC <- doubletFinder_v3(seu_hnESC, PCs = 1:30, pN = 0.25, pK = 0.06, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

reference <- subset(x=seu_hnESC, subset = (DF.classifications_0.25_0.05_310 == "Singlet"))
saveRDS(reference, file = paste0(resultsdir,"/hnESC_singlets.rds"))

# Add reference to meta-atlas
reference <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_hnesc/20221110/hnESC_singlets.rds")
seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220829/4datasets_annotated_joined.rds")
seur_obj$"predicted.id" <- seur_obj$scVI_label

comp.combined <- merge(seur_obj, y = reference, add.cell.ids = c("meta", "hnESC"), project = "comp")
unique(comp.combined$predicted.id)
DefaultAssay(object = comp.combined) <- "RNA"

# Set raw counts as the main data
comp.combined@assays$RNA@data <- comp.combined@assays$RNA@counts

comp <- DietSeurat(
  comp.combined,
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE,
  features = NULL,
  assays = NULL,
  dimreducs = NULL,
  graphs = NULL
)

saveRDS(comp, file = paste0(resultsdir,"/meta_hnESC.rds"))

# Convert seur object to h5ad for import into AnanseScanpy
library(SeuratDisk)

SaveH5Seurat(comp, filename = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_hnesc/20221110/meta_hnESC.h5Seurat")
Convert("meta_hnESC.h5Seurat", dest = "h5ad")


# Load in the snapatac environment
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
workdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scATAC_hnESC/"

# setting up results directory
dateoftoday <- gsub("-", "", as.character(Sys.Date()))
resultsdir <- paste0(workdir, dateoftoday)
system(paste("mkdir -p ", resultsdir))


# ATAC preprocessing and join with meta-atlas
meta_all.atac = readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/atac_meta_all.rds")

meta_all.atac@meta.data




x.sp = createSnap(file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/ESC_vitro_sc/multiomics/hnESC/outs/hnESC.snap", sample="hnESC",num.cores=1)
summarySnap(x.sp)

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

# Matrix binarisation
x.sp = addBmatToSnap(x.sp, bin.size=5000)
x.sp = makeBinary(x.sp, mat="bmat")

# Adding the promoter ratio to the barcodes file
promoter.df = read.table("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/promoter.bed")
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

x.meta.sp= readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/meta_04.rds")
levels(as.factor(x.meta.sp@metaData$predict.id))
clusters.sel = names(table(x.meta.sp@metaData$predict.id))[which(table(x.meta.sp@metaData$predict.id) > 100)]
x.meta.sp <- x.meta.sp[x.meta.sp@metaData$predict.max.score > 0.4 & x.meta.sp@metaData$predict.id%in%clusters.sel,]


x.meta.sp@metaData$predicted.id <- x.meta.sp@metaData$predict.id
x.meta.sp@metaData <- x.meta.sp@metaData[,c("barcode","predicted.id")]

x.sp@metaData$predicted.id <- rep("hnESC",length(x.sp@metaData$barcode))
x.sp@metaData <- x.sp@metaData[,c("barcode","predicted.id")]
# x.sp@feature
# x.meta.sp@feature
# library(GenomicRanges)
# grl3 <- c(x.sp@feature,x.meta.sp@feature)
# regroup(grl3, names(grl3))
x.sp3 = x.meta.sp

write.csv(x.sp@metaData, paste0(resultsdir, "/hnESC_metadata.csv"), row.names=FALSE)
write.csv(x.sp3@metaData, paste0(resultsdir, "/meta_metadata.csv"), row.names=FALSE)

meta_all.atac@assays$ATAC

meta_all.atac <- snapToSeurat(obj=x.sp3,
                              eigs.dims = 1:20,
                              norm=F,
                              scale=F
)

# Save the seurat object as a rds file
saveRDS(meta_all.atac,file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scATAC_hnESC/atac_meta_all.rds")

genes.df <- read.table("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/gencode.v30.gene.annotation.gtf") #already ran this
genes.gr3 = GRanges(genes.df[,1], IRanges(genes.df[,2], genes.df[,3]), name=genes.df[,4])
#genes.sel.gr3 = genes.gr3[which(genes.gr3$name %in% variable.genes)]

x.sp = createGmatFromMat(
  obj=x.sp,
  input.mat="bmat",
  genes=genes.gr3,
  do.par=TRUE#,
  #num.cores=10
)


hnESC.atac <- snapToSeurat(obj=x.sp,
                              eigs.dims = 1:20,
                              norm=F,
                              scale=F
)

# Save the seurat object as a rds file
saveRDS(hnESC.atac,file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scATAC_hnESC/atac_hnESC.rds")


# Join the two snapobjects and convert to Seurat for generation of chromatin assay
# Checking the OS path of snaptools and MACS2
system("which snaptools")
system("which macs2")


######################################################
# PEAK CALLING DONE WITH snakefile
# Setting the cores
SnowParam(workers = 2)

mclapply(seq(clusters.sel), function(i){
  print(clusters.sel[i])
  peaks = runMACS(
    obj=x.sp3[which(x.sp3@metaData$predicted.id==clusters.sel[i]),],
    output.prefix=paste0("meta.", gsub(" ", "_", clusters.sel)[i]),
    path.to.snaptools="/vol/mbconda/julian/envs/kb_snapatac/bin/snaptools",
    path.to.macs="/vol/mbconda/julian/envs/kb_snapatac/bin/macs2",
    gsize="hs", # mm, hs, etc
    buffer.size=500,
    num.cores=1,
    macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
    tmp.folder="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/hnESC/peaks_tmp"
  )
  #peaks
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
  write.table(bed_df, file=paste0("meta.", gsub(" ", "_", clusters.sel[i]), "_binary_matrix.bed"), quote=FALSE, sep= "\t", row.names=FALSE, col.names=FALSE)
}

# Create bed_df for every cluster
for (cluster_name in seq(clusters.sel)) {
  binarized_counts_cluster <- Matrix::colSums(x.sp3@bmat[which(x.sp3@cluster==clusters.sel[cluster_name]),])
  
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

write.table(peaks.df, file=paste0(resultsdir, "_peaks_combined_meta.bed"), append = FALSE, quote=FALSE, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")

# Delete weird characters
peaks.df$seqnames<-gsub("b'","",as.character(peaks.df$seqnames))
peaks.df
peaks.df$seqnames<-gsub("'","",as.character(peaks.df$seqnames))
peaks.df

write.table(peaks.df, file=paste0(resultsdir, "_peaks_combined_nochar_meta.bed"), append = FALSE, quote=FALSE, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")


######################################################
# PEAK CALLING DONE WITH snakefile
# Setting the cores
SnowParam(workers = 2)
clusters.sel="hnESC"
mclapply(seq(clusters.sel), function(i){
  print(clusters.sel[i])
  peaks = runMACS(
    obj=x.sp[which(x.sp@metaData$predicted.id==clusters.sel[i]),],
    output.prefix=paste0("meta.", gsub(" ", "_", clusters.sel)[i]),
    path.to.snaptools="/vol/mbconda/julian/envs/kb_snapatac/bin/snaptools",
    path.to.macs="/vol/mbconda/julian/envs/kb_snapatac/bin/macs2",
    gsize="hs", # mm, hs, etc
    buffer.size=500,
    num.cores=1,
    macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
    tmp.folder="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/hnESC/peaks_tmp"
  )
  #peaks
})


##############################################################################
# See bash for bigwig generation (Github)

# Making the bed dataframe
binarized_counts <- Matrix::colSums(x.sp@bmat)

#Scale scoring for bed file using the maximum n of reads in a bin and create a name for the bins
scaled_binarized_counts <- binarized_counts * (1000/max(binarized_counts))
names_binarized_counts <- paste0(binarized_counts, "of", length(x.sp@barcode))

#Create bed datafame with regions of bins extracted from x.sp@feature and scores based on sum of all cells (colSums of x.sp@bmat)
bed_df <- data.frame(seqnames=seqnames(x.sp@feature),
                     starts=start(x.sp@feature)-1,
                     ends=end(x.sp@feature),
                     names=names_binarized_counts,
                     scores=scaled_binarized_counts,
                     strands=strand(x.sp@feature))

#Rename strands
levels(bed_df$strands)[levels(bed_df$strands)=="*"] <- "."

# Write out table
for (i in clusters.sel){
  write.table(bed_df, file=paste0("hnESC.", gsub(" ", "_", clusters.sel[i]), "_binary_matrix.bed"), quote=FALSE, sep= "\t", row.names=FALSE, col.names=FALSE)
}

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

write.table(peaks.df, file=paste0(resultsdir, "_peaks_combined_hnESC.bed"), append = FALSE, quote=FALSE, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")

# Delete weird characters
peaks.df$seqnames<-gsub("b'","",as.character(peaks.df$seqnames))
peaks.df
peaks.df$seqnames<-gsub("'","",as.character(peaks.df$seqnames))
peaks.df

write.table(peaks.df, file=paste0(resultsdir, "_peaks_combined_nochar_hnESC.bed"), append = FALSE, quote=FALSE, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")

# Join the bed files into one
hnESC_bed<-read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scATAC_hnESC/20221115_peaks_combined_nochar_hnESC.bed",sep="\t",header = F)
meta_bed<-read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scATAC_hnESC/20221115_peaks_combined_nochar_meta.bed",sep="\t",header=F)

# Add Pmat to snaps with this bed file
join_bed <- rbind(hnESC_bed,meta_bed)

write.table(join_bed, file=paste0(resultsdir, "/joined.bed"), append = FALSE, quote=FALSE, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")

saveRDS(x.sp,file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/hnESC_atac.rds")
saveRDS(x.sp3,file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/meta_atac.rds")

# Reload the object to add Pmats
x.sp <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/hnESC_atac.rds")
x.sp3 <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/meta_atac.rds")
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# library(BiocManager)
# BiocManager::install("rhdf5",force = TRUE)

# Add Pmat after bigwig generation in bash
SnowParam(workers = 48)

#x.sp = createSnap(file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/ESC_vitro_sc/multiomics/hnESC/outs/hnESC.snap", sample="hnESC",num.cores=1)
#summarySnap(x.sp)

x.sp = addPmatToSnap(x.sp,do.par = T,num.cores = 48)
x.sp = makeBinary(x.sp, mat="pmat")
x.sp@pmat
sum(x.sp@pmat)

# Add Pmat after bigwig generation in bash
SnowParam(workers = 12)
x.sp3 = addPmatToSnap(x.sp3,do.par = T,num.cores = 6)
x.sp3@pmat
sum(x.sp3@pmat)

# Load the peak counts from snapobjects to seurat
# snaptoseurat function only converts gmat to seuart assay, not pmat

pmat_meta <- x.sp3@pmat
cols_meta<-x.sp3@peak@elementMetadata@listData$name
cols_meta<-gsub("b'","",as.character(cols_meta))
cols_meta<-gsub("'","",as.character(cols_meta))
colnames(pmat_meta) <- cols_meta
pmat_meta <- t(pmat_meta)

# Remove unwanted chromosomes
pmat_meta <- pmat_meta[-which(grepl("^GL",rownames(pmat_meta))),]    #remove rows containing GL- from matrix
pmat_meta <- pmat_meta[-which(grepl("^KI",rownames(pmat_meta))),]    #remove rows containing GL- from matrix
pmat_meta <- pmat_meta[-which(grepl("^chrM",rownames(pmat_meta))),]

# Check if the removal worked
unique(vapply(strsplit(rownames(pmat_meta),":"), `[`, 1, FUN.VALUE=character(1)))

library(Seurat)
meta_atac<-CreateSeuratObject(pmat_meta, project = "meta_atac")
meta_atac@meta.data=x.sp3@metaData

pmat_hnESC <- x.sp@pmat
cols_hnESC<-x.sp@peak@elementMetadata@listData$name
cols_hnESC<-gsub("b'","",as.character(cols_hnESC))
cols_hnESC<-gsub("'","",as.character(cols_hnESC))
colnames(pmat_hnESC) <- cols_hnESC
pmat_hnESC <- t(pmat_hnESC)

# Remove unwanted chromosomes
pmat_hnESC <- pmat_hnESC[-which(grepl("^GL",rownames(pmat_hnESC))),]    #remove rows containing GL- from matrix
pmat_hnESC <- pmat_hnESC[-which(grepl("^KI",rownames(pmat_hnESC))),]    #remove rows containing GL- from matrix
pmat_hnESC <- pmat_hnESC[-which(grepl("^chrM",rownames(pmat_hnESC))),]

# Check if the removal worked
unique(vapply(strsplit(rownames(pmat_hnESC),":"), `[`, 1, FUN.VALUE=character(1)))

hnESC_atac<-CreateSeuratObject(pmat_hnESC, project = "hnESC_atac")
hnESC_atac@meta.data=x.sp@metaData

# Join the seurat objects with joined peak counts
seur_comp <- merge(hnESC_atac, y = meta_atac, add.cell.ids = c("hnESC", "meta"), project = "joined")

# Check the predicted.id
unique(seur_comp@meta.data$predicted.id)

# Check the chromosomes
unique(vapply(strsplit(rownames(seur_comp@assays$RNA),":"), `[`, 1, FUN.VALUE=character(1)))

saveRDS(seur_comp,"/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scATAC_hnESC/20221116/seur_comp_atac.rds")

# Convert seur object to h5ad for import into AnanseScanpy
library(SeuratDisk)

SaveH5Seurat(seur_comp_atac, filename = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scATAC_hnESC/20221116/seur_comp_atac.h5Seurat")
Convert("seur_comp_atac.h5Seurat", dest = "h5ad")


# seur_comp@assays$RNA %>% 
#   filter(!grepl('GL', Name))

# gene_df <- read.table("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/gencode.v30.gene.annotation.gtf")
# genes.gr = GRanges(gene_df[,1],
#                    IRanges(gene_df[,2], gene_df[,3]),
#                    name=gene_df[,4]
# )
# genes.gr2  <- genes.gr[!duplicated(genes.gr$name),]
# x.sp = createGmatFromMat(
#   obj=x.sp,
#   genes= genes.gr2,
#   do.par=TRUE
# )
# 
# hnESC_atac <- snapToSeurat(x.sp, eigs.dims = 1:10, norm = F, scale = F)
# 
# 
# 
# fragpath <- paste0(workdir, "data/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz")
# 
# # get gene annotations for hg38
# annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# seqlevelsStyle(annotation) <- "UCSC"
# 
# # create ATAC assay and add it to the object
# pbmc[["ATAC"]] <- CreateChromatinAssay(
#   counts = counts$Peaks,
#   sep = c(":", "-"),
#   fragments = fragpath,
#   annotation = annotation
# )
# 
# DefaultAssay(pbmc) <- "ATAC"
# 
# pbmc <- NucleosomeSignal(pbmc)
# pbmc <- TSSEnrichment(pbmc)
# 
# # VlnPlot(
# #   object = pbmc,
# #   features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
# #   ncol = 4,
# #   pt.size = 0
# # )
# 
# # filter out low quality cells
# pbmc <- subset(
#   x = pbmc,
#   subset = nCount_ATAC < 100000 &
#     nCount_RNA < 25000 &
#     nCount_ATAC > 1000 &
#     nCount_RNA > 1000 &
#     nucleosome_signal < 2 &
#     TSS.enrichment > 1
# )
# pbmc
# 
# # call peaks using MACS2
# peaks <- CallPeaks(pbmc, macs2.path = "/vol/mbconda/julian/envs/kb_snapatac/bin/macs2")
# 
# # remove peaks on nonstandard chromosomes and in genomic blacklist regions
# peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
# peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
# 
# # quantify counts in each peak
# macs2_counts <- FeatureMatrix(
#   fragments = Fragments(pbmc),
#   features = peaks,
#   cells = colnames(pbmc)
# )
# 
# # create a new assay using the MACS2 peak set and add it to the Seurat object
# pbmc[["peaks"]] <- CreateChromatinAssay(
#   counts = macs2_counts,
#   fragments = fragpath,
#   annotation = annotation
# )