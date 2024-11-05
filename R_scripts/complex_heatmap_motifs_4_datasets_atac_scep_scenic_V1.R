# COMPLEX HEATMAP OF MOTIFS AND RNA-seq
# Loading in important packages
library(chorddiag)
library(circlize)
library(viridis)
library(ggvenn)
#install.packages('venn')
library(venn)
library(dplyr)
library(tidyr)
library(tidyselect)
library(ComplexHeatmap)

# Setting up the directory
workdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/complexheatmap_motifs/"

# Setting up results directory
dateoftoday <- gsub("-", "", as.character(Sys.Date()))
resultsdir <- paste0(workdir, dateoftoday)
system(paste("mkdir -p ", resultsdir))

# Read in transcription factor annotation files if they are available
expected <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/ANANSE/Expected_tfs.csv",header = T, comment.char = '#')
exgeneral <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/ANANSE/Expectedgeneral_tfs.csv",header = T, comment.char = '#')

## load in the z-score normalized dataset:
lakorna <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/Z-score_datasets20220225/Zscoretable_markers_split2.tsv", sep="\t",row.names=1)

lakorna <- na.omit(lakorna)

# for filtering
lakorna <- lakorna[, !names(lakorna) %in% c("IC", "Mel", "Ves")]

# Setting the colors of the general transcription factors expected to be found
neurexpected <- unlist(strsplit(exgeneral[exgeneral$cell_type == "neural",]$expected_tfs,","))
eyeexpected <- unlist(strsplit(exgeneral[exgeneral$cell_type == "eye",]$expected_tfs,","))
epiexpected <- unlist(strsplit(exgeneral[exgeneral$cell_type == "epidermal",]$expected_tfs,","))

##############################

lakoqnormmotifs <- read.table(file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/scANANSE_11032022/ATAC_motif/final.out.txt",sep = '\t', header = TRUE, row.names=1)
lakoqnormmotifs<- lakoqnormmotifs[,grep("z.score", names(lakoqnormmotifs), value = TRUE)]

m2f <- read.table(file ="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/scANANSE_11032022/ATAC_motif/gimme.vertebrate.v5.0.motif2factors.txt",sep = '\t', header = TRUE)
m2f <- m2f[,1:2]

lakoqnormmotifs$Motif <- rownames(lakoqnormmotifs)
lakoqnormmotifs <- merge(x = lakoqnormmotifs, y = m2f, by = "Motif", all = TRUE)

lakoqnormmotifs <- na.omit(lakoqnormmotifs)
lakoqnormmotifs$Factor <- toupper(lakoqnormmotifs$Factor)
lakoqnormmotifs <- lakoqnormmotifs[match(unique(lakoqnormmotifs$Factor), lakoqnormmotifs$Factor),]
lakoqnormmotifs2 <- lakoqnormmotifs
lakoqnormmotifs2$Motif <- NULL

# convert df to matrix
y <- lakoqnormmotifs2
rownames(y) <- y$Factor
y$Factor <- NULL

maty <- as.matrix(y)

 vec2 <- NULL
 for (i in colnames(y)) {
 sel <- maty[order(maty[,i],decreasing = T),]
 vec <- rownames(sel)[1:30]
 print(vec)
 vec2 <- c(vec2,vec)
 }
 vec2 <- unique(vec2)

# Comparing the relative expression
matz <- as.matrix(lakorna)
matz <- as.matrix(matz[rownames(matz)%in%rownames(y),])

matz <- as.matrix(matz[rownames(matz)%in%rownames(y),])

# roword <- row_order(ht2)
# matz <- as.matrix(matz[roword,])



# Making the complex heatmap
f1 = colorRamp2(c(-3, 0, 3), c("blue", "#EEEEEE", "red"), space = "RGB")
f2 = colorRamp2(c(-5, 0, 5), c("grey39", "#EEEEEE", "mediumseagreen"), space = "RGB")

# Only if you want the top motifs
maty <- subset(maty, rownames(maty) %in% vec2)

maty <- as.matrix(maty)

maty <- maty[rownames(maty)%in%rownames(matz),]
vec <- NULL
for (i in colnames(maty)){
        vec <- c(vec,unlist(strsplit(i, split='.', fixed=TRUE))[3])
}
colnames(maty) <- vec

matz <- matz[rownames(maty),]

# Show motifs  without TFs only in comparison to SCEPIA?
lakoqnormmotifs_uniq <- lakoqnormmotifs[match(unique(lakoqnormmotifs$Motif), lakoqnormmotifs$Motif),]
lakoqnormmotifs_uniq <- lakoqnormmotifs_uniq
#lakoqnormmotifs_uniq$Motif <- NULL

# convert df to matrix
y2 <- lakoqnormmotifs_uniq
rownames(y2) <- y2$Motif
y2$Factor <- NULL
y2$Motif <- NULL

matmot_atac <- as.matrix(y2)

# Load in scepia inferred motifs
# # importing SCEPIA motif data
# 
# 
# # load in the motifs after SCEPIA

#  SCEPIA_LESC <- read.table(file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/LESCmotif_act.tsv",sep = '\t', header = TRUE)
#  names(SCEPIA_LESC) <- c('motif','SCEP_LESC')
#  SCEPIA_LSC <- read.table(file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/LSCmotif_act.tsv",sep = '\t', header = TRUE)
#  names(SCEPIA_LSC) <- c('motif','SCEP_LSC')
#  SCEPIA_LE <- read.table(file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/LEmotif_act.tsv",sep = '\t', header = TRUE)
#  names(SCEPIA_LE) <- c('motif','SCEP_LE')
#  SCEPIA_Cj <- read.table(file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/Cjmotif_act.tsv",sep = '\t', header = TRUE)
#  names(SCEPIA_Cj) <- c('motif','SCEP_Cj')
#  SCEPIA_CE <- read.table(file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/CEmotif_act.tsv",sep = '\t', header = TRUE)
#  names(SCEPIA_CE) <- c('motif','SCEP_CE')
#  SCEPIA_CSSC <- read.table(file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/CSSCmotif_act.tsv",sep = '\t', header = TRUE)
#  names(SCEPIA_CSSC) <- c('motif','SCEP_CSSC')
#  SCEPIA_CF <- read.table(file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/CFmotif_act.tsv",sep = '\t', header = TRUE)
#  names(SCEPIA_CF) <- c('motif','SCEP_CF')
# # 
# # 
#  SCEPIA_full <- merge(SCEPIA_LESC,SCEPIA_LSC, by = 'motif')
#  SCEPIA_full <- merge(SCEPIA_full, SCEPIA_LE, by='motif',all=TRUE)# SCEPIA_full <- merge(SCEPIA_full, SCEPIA_LNSC, by='motif',all=TRUE)
#  SCEPIA_full <- merge(SCEPIA_full, SCEPIA_Cj, by='motif',all=TRUE)
#  SCEPIA_full <- merge(SCEPIA_full, SCEPIA_CE, by='motif',all=TRUE)
#  SCEPIA_full <- merge(SCEPIA_full, SCEPIA_CSSC, by='motif',all=TRUE)
#  SCEPIA_full <- merge(SCEPIA_full, SCEPIA_CF, by='motif',all=TRUE)
# # 
# # 
#  SCEPIA_full$Motif <-SCEPIA_full$motif
#  SCEPIA_full$motif <- NULL
SCEPIA_full<- read.table(file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/pseudobulk_20022022_meta.tsv",sep = '\t', header = TRUE,row.names = 1)
SCEPIA_full <- as.data.frame(t(SCEPIA_full))
SCEPIA_full$Motif <- rownames(SCEPIA_full)
# 
# 
 rownames(lakoqnormmotifs_uniq) <- lakoqnormmotifs_uniq$Motif
 lakoqnormmotifs_uniq$Motif <- NULL
 lakoqnormmotifs_uniq$Factor <- NULL
 vec <- NULL
 for (i in colnames(lakoqnormmotifs_uniq)){
         vec <- c(vec,unlist(strsplit(i, split='.', fixed=TRUE))[3])
 }
 colnames(lakoqnormmotifs_uniq) <- vec
 lakoqnormmotifs_uniq$Motif <- rownames(lakoqnormmotifs_uniq)
 
combined <- merge(lakoqnormmotifs_uniq, SCEPIA_full, by='Motif',all=TRUE)
combined <- na.omit(combined)
rownames(combined) <-combined$Motif
combined$Motif<- NULL

comb_mat <- as.matrix(combined)

comb_mat2 <- comb_mat

# Select the top 10 fators in scATAC motif prediciton
vec2 <- NULL
for (i in colnames(comb_mat2)[1:7]) {
        sel <- comb_mat2[order(comb_mat2[,i],decreasing = T),]
        vec <- rownames(sel)[1:10]
        print(vec)
        vec2 <- c(vec2,vec)
}
vec2 <- unique(vec2)

comb_mat2 <- comb_mat2[rownames(comb_mat2) %in% vec2,]
colord <- c("LSC","LESC","LE","CE","Cj","CSSC","CF","LSC_scepia","LESC_scepia","LE_scepia","CE_scepia","Cj_scepia","CSSC_scepia","CF_scepia")
comb_mat2 <- comb_mat2[, colord]
# split combined mat in two for the scales
comb_mat3 <- comb_mat2[,colnames(comb_mat2) %in%c("LSC_scepia","LESC_scepia","LE_scepia","CE_scepia","Cj_scepia","CSSC_scepia","CF_scepia")]

comb_mat2 <- comb_mat2[,colnames(comb_mat2) %in%c("LSC","LESC","LE","CE","Cj","CSSC","CF")]

ht1 = Heatmap(comb_mat2, col = f2, cluster_columns = F,cluster_rows = T, name = "Relative motif accesibility")

roword <- row_order(ht1)
comb_mat2 <- comb_mat2[roword,]

comb_mat3 <- comb_mat3[rownames(comb_mat2),]

ht1 = Heatmap(comb_mat2, col = f2, cluster_columns = F,cluster_rows = F, name = "Relative motif accesibility")
f4 = colorRamp2(c(-2e-3, 0, 2e-3), c("grey39", "#EEEEEE", "orange"), space = "RGB")
ht2 = Heatmap(comb_mat3, col = f4, cluster_columns = F,cluster_rows = F, name = "Relative motif accesibility Scepia")

vec2 <- rownames(comb_mat2)
vec1 <- which(rownames(comb_mat2) %in% vec2, arr.ind = T)

# Adding the motif names on the right
v <- rownames(comb_mat2)

# Specify motif type
library(stringr)
v2 <- v
v2 <- v2 %>% str_replace("GM.5.0.", "")
v2 <- gsub("\\..*","",v2)

f3 <-viridis(length(unique(v2)),option = "H",alpha = 0.8)

pdf(paste(resultsdir,'/complexheatmap_motifs_marked_motifs_scepia.pdf',sep="/") ,width=10,height=10,paper='special')
htlist <-ht1 + ht2 + Heatmap(v2, col = f3,name = "Motif type", width = unit(0.5, "cm"), heatmap_legend_param = list(
        at = unique(v2),
        labels = unique(v2),
        title = "Motif type",
        col = f3,
        legend_height = unit(4, "cm"))) + rowAnnotation(link = anno_mark(at =  vec1,labels = vec2))
draw(htlist, column_title = "inferred motifs comparison")
dev.off()

# Color the cell identities & the motifs with the same color palette
# vec <- viridis(40,option = "H",alpha = 0.8)
# print(vec)
# 
# newcol <- c(vec[5],vec[31],vec[2],vec[28],vec[15],vec[20],vec[19])
# f4 <- rbind(newcol)
# 
# # Remove colors already in cell identities
# f3 <- sample (vec, size=length(unique(v2))+length(colnames(matz)), replace =F)
# vec_col <- f3[!f3 %in% newcol]
# f3 <- sample (vec_col, size=length(unique(v2)), replace =F)
# #matscep<- matscep[rownames(maty),colnames(maty)]
#############################################################
# Show the top TF for each motif
colord <- c("LSC","LESC","LE","CE","Cj","CSSC","CF","LSC_scepia","LESC_scepia","LE_scepia","CE_scepia","Cj_scepia","CSSC_scepia","CF_scepia")
comb_mat <- comb_mat[, colord]
comb_mat2 <- comb_mat


# Select the top 10 fators in scATAC motif prediciton
vec2 <- NULL
for (i in colnames(comb_mat2)[1:7]) {
        sel <- comb_mat2[order(comb_mat2[,i],decreasing = T),]
        vec <- rownames(sel)[1:10]
        print(vec)
        vec2 <- c(vec2,vec)
}
vec2 <- unique(vec2)

comb_mat2 <- comb_mat2[rownames(comb_mat2) %in% vec2,]
colord <- c("LSC","LESC","LE","CE","Cj","CSSC","CF","LSC_scepia","LESC_scepia","LE_scepia","CE_scepia","Cj_scepia","CSSC_scepia","CF_scepia")
comb_mat2 <- comb_mat2[, colord]
# split combined mat in two for the scales
comb_mat3 <- comb_mat2[,colnames(comb_mat2) %in%c("LSC_scepia","LESC_scepia","LE_scepia","CE_scepia","Cj_scepia","CSSC_scepia","CF_scepia")]

comb_mat2 <- comb_mat2[,colnames(comb_mat2) %in%c("LSC","LESC","LE","CE","Cj","CSSC","CF")]

comb_mat2 <- as.data.frame(comb_mat2)
comb_mat2$Motif <- rownames(comb_mat2)

matz <- as.matrix(lakorna)
m2f<- m2f[m2f$Factor%in%rownames(matz),]

comb_mat2 <- merge(x = comb_mat2, y = m2f, by = "Motif", all = TRUE)

comb_mat2 <- na.omit(comb_mat2)
comb_mat2$Factor <- toupper(comb_mat2$Factor)
comb_mat2 <- comb_mat2[match(unique(comb_mat2$Factor), comb_mat2$Factor),]
comb_mat22 <- comb_mat2
#comb_mat22$Motif <- NULL

matz <- matz[rownames(matz)%in%comb_mat22$Factor,]

matz <- matz[,c("LSC","LESC","LE","CE","Cj","CSSC","CF")]

comb_mat22 <- comb_mat22[comb_mat22$Factor%in%rownames(matz),]

rownames(comb_mat22)<-comb_mat22$Factor
comb_mat22$Factor<-NULL
comb_mat22$Motif <- NULL

comb_mat22 <- comb_mat22[rownames(matz),]
# Correlation score of rows
sel_vec <- sapply(1:nrow(as.matrix(comb_mat22)), function(i) cor(as.matrix(comb_mat22)[i,], as.matrix(matz)[i,]))
#sel_vec <- sel_vec > 0.6
comb_mat22$cor <- sel_vec
comb_mat22$Factor <- rownames(comb_mat22)

comb_mat22 <- comb_mat22[,c("cor","Factor")]
comb_mat_cor <- merge(x = comb_mat2, y = comb_mat22, by = "Factor", all = TRUE)
# mark the most correlating TFs

library(dplyr)
comb_mat_cor <- comb_mat_cor %>% arrange(factor(Motif, as.character(unique(Motif))), desc(cor))
comb_mat_cor <- Reduce(rbind,                                 # Top N highest values by group
                    by(comb_mat_cor,
                       comb_mat_cor["Motif"],
                       head,
                       n = 1))

rownames(comb_mat_cor) <- comb_mat_cor$Motif
comb_mat_cor$Motif <- NULL

# label the rows according to the corresponding TFs
comb_mat2 <- comb_mat

# Select the top 10 fators in scATAC motif prediciton
vec2 <- NULL
for (i in colnames(comb_mat2)[1:7]) {
        sel <- comb_mat2[order(comb_mat2[,i],decreasing = T),]
        vec <- rownames(sel)[1:10]
        print(vec)
        vec2 <- c(vec2,vec)
}
vec2 <- unique(vec2)

comb_mat2 <- comb_mat2[rownames(comb_mat2) %in% vec2,]

# split combined mat in two for the scales
comb_mat3 <- comb_mat2[,colnames(comb_mat2) %in%c("LSC_scepia","LESC_scepia","LE_scepia","CE_scepia","Cj_scepia","CSSC_scepia","CF_scepia")]

comb_mat2 <- comb_mat2[,colnames(comb_mat2) %in%c("LSC","LESC","LE","CE","Cj","CSSC","CF")]

ht1 = Heatmap(comb_mat2, col = f2, cluster_columns = F,cluster_rows = T, name = "Relative motif accesibility")

roword <- row_order(ht1)
comb_mat2 <- comb_mat2[roword,]

comb_mat3 <- comb_mat3[rownames(comb_mat2),]

ht1 = Heatmap(comb_mat2, col = f2, cluster_columns = F,cluster_rows = F, name = "Relative motif accesibility")

f4 = colorRamp2(c(-2e-3, 0, 2e-3), c("grey39", "#EEEEEE", "orange"), space = "RGB")
ht2 = Heatmap(comb_mat3, col = f4, cluster_columns = F,cluster_rows = F, name = "Relative motif accesibility Scepia")


comb_mat_cor <- comb_mat_cor[rownames(comb_mat2),]
vec2 <- na.omit(comb_mat_cor$Factor)
vec1 <- which(!is.na(comb_mat_cor$Factor),arr.ind=TRUE)

# # Adding the motif names on the right
# v <- rownames(comb_mat2)
# 
# # Specify motif type
# library(stringr)
# v2 <- v
# v2 <- v2 %>% str_replace("GM.5.0.", "")
# v2 <- gsub("\\..*","",v2)
# 
# f3 <-viridis(length(unique(v2)),option = "H",alpha = 0.8)

pdf(paste(resultsdir,'/complexheatmap_motifs_marked_motifs_scepia_facts.pdf',sep="/") ,width=8,height=10,paper='special')
htlist <-ht1 + ht2 + Heatmap(v2, col = f3,name = "Motif type", width = unit(0.5, "cm"), heatmap_legend_param = list(
        at = unique(v2),
        labels = unique(v2),
        title = "Motif type",
        col = f3,
        legend_height = unit(4, "cm"))) + rowAnnotation(link = anno_mark(at =  vec1,labels = vec2))
draw(htlist, column_title = "inferred motifs comparison")
dev.off()

#df2 <- df2[df2$rownames.matz.%in%data_new1$rownames.matz.,]



# label the rows according to the corresponding TFs
comb_mat4<-comb_mat2[vec1,]
comb_mat5<-comb_mat3[vec1,]

mat_RNA <- lakorna[rownames(lakorna) %in% vec2,]
mat_RNA <- mat_RNA[vec2,]
# # Adding the motif names on the right
# v <- rownames(comb_mat2)
# 
# # Specify motif type
# library(stringr)
# v2 <- v
# v2 <- v2 %>% str_replace("GM.5.0.", "")
# v2 <- gsub("\\..*","",v2)
# 
ht1 = Heatmap(comb_mat4, col = f2, cluster_columns = F,cluster_rows = F, name = "Relative motif accesibility")

ht2 = Heatmap(comb_mat5, col = f4, cluster_columns = F,cluster_rows = F, name = "Relative motif accesibility Scepia")

f5 = colorRamp2(c(-2, 0, 2), c("blue", "#EEEEEE", "red"), space = "RGB")

mat_RNA <- mat_RNA[,c("LSC","LESC","LE","CE","Cj","CSSC","CF")]

ht3 = Heatmap(mat_RNA, col = f5, cluster_columns = F,cluster_rows = F, name = "Relative expression")

# Adding the motif names on the right
v <- rownames(comb_mat4)

# Specify motif type
library(stringr)
v2 <- v
v2 <- v2 %>% str_replace("GM.5.0.", "")
v2 <- gsub("\\..*","",v2)

f3 <-viridis(length(unique(v2)),option = "H",alpha = 0.8)

comb_mat_cor <- comb_mat_cor[rownames(comb_mat4),]
vec2 <- na.omit(comb_mat_cor$Factor)
vec1 <- which(!is.na(comb_mat_cor$Factor),arr.ind=TRUE)

pdf(paste(resultsdir,'/complexheatmap_motifs_marked_motifs_scepia_facts_RNA.pdf',sep="/") ,width=7.5,height=7.5,paper='special')
htlist <-ht1 + ht2 + ht3 + Heatmap(v2, col = f3,name = "Motif type", width = unit(0.5, "cm"), heatmap_legend_param = list(
        at = unique(v2),
        labels = unique(v2),
        title = "Motif type",
        col = f3,
        legend_height = unit(4, "cm"))) + rowAnnotation(link = anno_mark(at =  vec1,labels = vec2))
draw(htlist, column_title = "inferred motifs comparison")
dev.off()

# annotate the methods at the bottom
ha_met <- columnAnnotation(foo = anno_block(gp = gpar(fill = "lightgrey"),
                               labels = "scATAC-seq", 
                               labels_gp = gpar(col = "black", fontsize = 12)))

ha_met2 <- columnAnnotation(foo = anno_block(gp = gpar(fill = "lightgrey"),
                                            labels = "Scepia", 
                                            labels_gp = gpar(col = "black", fontsize = 12)))
ha_met3 <- columnAnnotation(foo = anno_block(gp = gpar(fill = "lightgrey"),
                                            labels = "scRNA-seq", 
                                            labels_gp = gpar(col = "black", fontsize = 12)))

ht1 = Heatmap(comb_mat4, col = f2, cluster_columns = F,cluster_rows = F, name = "Relative motif accesibility",top_annotation = ha_met)

ht2 = Heatmap(comb_mat5, col = f4, cluster_columns = F,cluster_rows = F, name = "Relative motif accesibility Scepia",top_annotation = ha_met2)

f5 = colorRamp2(c(-2, 0, 2), c("blue", "#EEEEEE", "red"), space = "RGB")

mat_RNA <- mat_RNA[,c("LSC","LESC","LE","CE","Cj","CSSC","CF")]

ht3 = Heatmap(mat_RNA, col = f5, cluster_columns = F,cluster_rows = F, name = "Relative expression",top_annotation = ha_met3)

# Adding the motif names on the right
v <- rownames(comb_mat4)

# Specify motif type
library(stringr)
v2 <- v
v2 <- v2 %>% str_replace("GM.5.0.", "")
v2 <- gsub("\\..*","",v2)

f3 <-viridis(length(unique(v2)),option = "H",alpha = 0.8)

comb_mat_cor <- comb_mat_cor[rownames(comb_mat4),]
vec2 <- na.omit(comb_mat_cor$Factor)
vec1 <- which(!is.na(comb_mat_cor$Factor),arr.ind=TRUE)

pdf(paste(resultsdir,'/complexheatmap_motifs_marked_motifs_scepia_facts_RNA_mtehods.pdf',sep="/") ,width=7.5,height=8,paper='special')
htlist <-ht1 + ht2 + ht3 + Heatmap(v2, col = f3,name = "Motif type", width = unit(0.5, "cm"), heatmap_legend_param = list(
  at = unique(v2),
  labels = unique(v2),
  title = "Motif type",
  col = f3,
  legend_height = unit(4, "cm"))) + rowAnnotation(link = anno_mark(at =  vec1,labels = vec2))
draw(htlist, column_title = "inferred motifs comparison")
dev.off()
#################################
# SCEPIA_full
#############################################################
# Show the top TF for each motif
SCEPIA_full<- read.table(file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/pseudobulk_02052022_meta.tsv",sep = '\t', header = TRUE,row.names = 1)
SCEPIA_full <- as.data.frame(t(SCEPIA_full))
SCEPIA_full$Motif <- rownames(SCEPIA_full)
SCEPIA_full$Motif <- NULL
comb_mat <- as.matrix(SCEPIA_full)

comb_mat2 <- comb_mat

# Select the top 10 fators in SCEPIA prediciton
vec2 <- NULL
for (i in colnames(comb_mat2)) {
  sel <- comb_mat2[order(comb_mat2[,i],decreasing = T),]
  vec <- rownames(sel)[1:10]
  print(vec)
  vec2 <- c(vec2,vec)
}
vec2 <- unique(vec2)

comb_mat2 <- comb_mat2[rownames(comb_mat2) %in% vec2,]

  
colord <- c("LSC_scepia","LESC_scepia","LE_scepia","CE_scepia","Cj_scepia","CSSC_scepia","CF_scepia","IC_scepia","Mel_scepia","Ves_scepia")
comb_mat2 <- comb_mat2[,colnames(comb_mat2) %in%colord]
comb_mat2 <- comb_mat2[, colord]
comb_mat3 <- comb_mat2

f4 = colorRamp2(c(-2e-3, 0, 2e-3), c("grey39", "#EEEEEE", "orange"), space = "RGB")
ht2 = Heatmap(comb_mat3, col = f4, cluster_columns = F,cluster_rows = F, name = "Relative motif accesibility Scepia")

vec2 <- rownames(comb_mat2)
vec1 <- which(rownames(comb_mat2) %in% vec2, arr.ind = T)

# Adding the motif names on the right
v <- rownames(comb_mat2)

# Specify motif type
library(stringr)
v2 <- v
v2 <- v2 %>% str_replace("GM.5.0.", "")
v2 <- gsub("\\..*","",v2)

f3 <-viridis(length(unique(v2)),option = "H",alpha = 0.8)

pdf(paste(resultsdir,'/complexheatmap_motifs_marked_motifs_scepia_all.pdf',sep="/") ,width=10,height=15,paper='special')
htlist <-ht2 + Heatmap(v2, col = f3,name = "Motif type", width = unit(0.5, "cm"), heatmap_legend_param = list(
  at = unique(v2),
  labels = unique(v2),
  title = "Motif type",
  col = f3,
  legend_height = unit(4, "cm"))) + rowAnnotation(link = anno_mark(at =  vec1,labels = vec2))
draw(htlist, column_title = "inferred motifs comparison")
dev.off()

comb_mat2 <- as.data.frame(comb_mat2)
comb_mat2$Motif <- rownames(comb_mat2)

matz <- as.matrix(lakorna)
m2f<- m2f[m2f$Factor%in%rownames(matz),]

comb_mat2 <- merge(x = comb_mat2, y = m2f, by = "Motif", all = TRUE)

comb_mat2 <- na.omit(comb_mat2)
comb_mat2$Factor <- toupper(comb_mat2$Factor)
comb_mat2 <- comb_mat2[match(unique(comb_mat2$Factor), comb_mat2$Factor),]
comb_mat22 <- comb_mat2

# reload matz
## load in the z-score normalized dataset:
lakorna <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/Z-score_datasets20220225/Zscoretable_markers_split2.tsv", sep="\t",row.names=1)

lakorna <- na.omit(lakorna)

matz <- as.matrix(lakorna)

# for filtering
matz <- matz[rownames(matz)%in%comb_mat22$Factor,]

comb_mat22 <- comb_mat22[comb_mat22$Factor%in%rownames(matz),]

rownames(comb_mat22)<-comb_mat22$Factor
comb_mat22$Factor<-NULL
comb_mat22$Motif <- NULL

comb_mat22 <- comb_mat22[rownames(matz),]

# Correlation score of rows
sel_vec <- sapply(1:nrow(as.matrix(comb_mat22)), function(i) cor(as.matrix(comb_mat22)[i,], as.matrix(matz)[i,]))
#sel_vec <- sel_vec > 0.6
comb_mat22$cor <- sel_vec
comb_mat22$Factor <- rownames(comb_mat22)

comb_mat22 <- comb_mat22[,c("cor","Factor")]
comb_mat_cor <- merge(x = comb_mat2, y = comb_mat22, by = "Factor", all = TRUE)
# mark the most correlating TFs

library(dplyr)
comb_mat_cor <- comb_mat_cor %>% arrange(factor(Motif, as.character(unique(Motif))), desc(cor))
comb_mat_cor <- Reduce(rbind,                                 # Top N highest values by group
                       by(comb_mat_cor,
                          comb_mat_cor["Motif"],
                          head,
                          n = 1))

rownames(comb_mat_cor) <- comb_mat_cor$Motif
comb_mat_cor$Motif <- NULL

# label the rows according to the corresponding TFs
comb_mat2 <- comb_mat
comb_mat2 <- comb_mat2[,colnames(comb_mat2) %in%colord]

# Select the top 10 fators in scATAC motif prediciton
vec2 <- NULL
for (i in colnames(comb_mat2)[1:7]) {
  sel <- comb_mat2[order(comb_mat2[,i],decreasing = T),]
  vec <- rownames(sel)[1:10]
  print(vec)
  vec2 <- c(vec2,vec)
}
vec2 <- unique(vec2)

comb_mat2 <- comb_mat2[rownames(comb_mat2) %in% vec2,]

comb_mat3 <- comb_mat2[rownames(comb_mat2),]

f4 = colorRamp2(c(-2e-3, 0, 2e-3), c("grey39", "#EEEEEE", "orange"), space = "RGB")
ht2 = Heatmap(comb_mat3, col = f4, cluster_columns = F,cluster_rows = F, name = "Relative motif accesibility Scepia")

comb_mat_cor <- comb_mat_cor[rownames(comb_mat2),]
vec2 <- na.omit(comb_mat_cor$Factor)
vec1 <- which(!is.na(comb_mat_cor$Factor),arr.ind=TRUE)

# # Adding the motif names on the right
v <- rownames(comb_mat2)
# 
# # Specify motif type
library(stringr)
v2 <- v
v2 <- v2 %>% str_replace("GM.5.0.", "")
v2 <- gsub("\\..*","",v2)
# 
f3 <-viridis(length(unique(v2)),option = "H",alpha = 0.8)

pdf(paste(resultsdir,'/complexheatmap_motifs_marked_motifs_scepia_facts_all.pdf',sep="/") ,width=8,height=10,paper='special')
htlist <-ht2 + Heatmap(v2, col = f3,name = "Motif type", width = unit(0.5, "cm"), heatmap_legend_param = list(
  at = unique(v2),
  labels = unique(v2),
  title = "Motif type",
  col = f3,
  legend_height = unit(4, "cm"))) + rowAnnotation(link = anno_mark(at =  vec1,labels = vec2))
draw(htlist, column_title = "inferred motifs comparison")
dev.off()

# label the rows according to the corresponding TFs
comb_mat5<-comb_mat3[vec1,]

mat_RNA <- lakorna[rownames(lakorna) %in% vec2,]
mat_RNA <- mat_RNA[vec2,]
mat_RNA <- as.matrix(mat_RNA)

colord2 <- c("LSC","LESC","LE","CE","Cj","CSSC","CF","IC","Mel","Ves")
mat_RNA <- mat_RNA[,colord2]
comb_mat5<- comb_mat5[,colord]

ht2 = Heatmap(comb_mat5, col = f4, cluster_columns = F,cluster_rows = F, name = "Relative motif accesibility Scepia")

f5 = colorRamp2(c(-2, 0, 2), c("blue", "#EEEEEE", "red"), space = "RGB")

ht3 = Heatmap(mat_RNA, col = f5, cluster_columns = F,cluster_rows = F, name = "Relative expression")

# Adding the motif names on the right
v <- rownames(comb_mat5)

# Specify motif type
library(stringr)
v2 <- v
v2 <- v2 %>% str_replace("GM.5.0.", "")
v2 <- gsub("\\..*","",v2)

f3 <-viridis(length(unique(v2)),option = "H",alpha = 0.8)

comb_mat_cor <- comb_mat_cor[rownames(comb_mat5),]
vec2 <- na.omit(comb_mat_cor$Factor)
vec1 <- which(!is.na(comb_mat_cor$Factor),arr.ind=TRUE)

pdf(paste(resultsdir,'/complexheatmap_motifs_marked_motifs_scepia_facts_RNA_all.pdf',sep="/") ,width=8,height=10,paper='special')
htlist <-ht2 + ht3 + Heatmap(v2, col = f3,name = "Motif type", width = unit(0.5, "cm"), heatmap_legend_param = list(
  at = unique(v2),
  labels = unique(v2),
  title = "Motif type",
  col = f3,
  legend_height = unit(4, "cm"))) + columnAnnotation(link = anno_mark(at =  vec1,labels = vec2))
draw(htlist, column_title = "inferred motifs comparison")
dev.off()

comb_mat5<- t(comb_mat5)

ht2 = Heatmap(comb_mat5, col = f4, cluster_columns = F,cluster_rows = F, name = "Relative motif accesibility Scepia")

mat_RNA<- t(mat_RNA)
ht3 = Heatmap(mat_RNA, col = f5, cluster_columns = F,cluster_rows = F, name = "Relative expression")
v3 <- t(v2)
anno <- anno_mark(at =  vec1,labels = vec2,which="column")

pdf(paste(resultsdir,'/complexheatmap_motifs_marked_motifs_scepia_facts_RNA_all_rot.pdf',sep="/") ,width=12,height=5,paper='special')
htlist <-Heatmap(v3, col = f3,name = "Motif type", heatmap_legend_param = list(
  at = unique(v3),
  labels = unique(v3),
  title = "Motif type",
  col = f3,
  legend_height = unit(4, "cm")))%v% ht2 %v% ht3
draw(htlist, column_title = "inferred motifs comparison")
dev.off()
#################################
# Scenic motif enrichment analysis normalized enrichment score (NES) test
Scenic_IC <- read.table(file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/scenic/26042022_IC/reg.csv",sep = ",", header = TRUE)

Scenic_IC <- Scenic_IC[-1:-3,1:6]
names(Scenic_IC) <-c("TF","Motif","AUC","NES","MSQ","OI")
Scenic_IC <- Scenic_IC[Scenic_IC$OI>=1,]
Scenic_IC <- Scenic_IC %>% group_by(Motif) %>% arrange(desc(NES),.by_group = T) %>% slice_head(n=1)

Scenic_Mel <- read.table(file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/scenic/26042022_Mel/reg.csv",sep = ',', header = TRUE)
Scenic_Mel <- Scenic_Mel[-1:-3,1:6]
names(Scenic_Mel) <-c("TF","Motif","AUC","NES","MSQ","OI")
Scenic_Mel <- Scenic_Mel[Scenic_Mel$OI>=1,]
Scenic_Mel <- Scenic_Mel %>% group_by(Motif) %>% arrange(desc(NES),.by_group = T) %>% slice_head(n=1)

Scenic_LSC <- read.table(file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/scenic/26042022_LSC/reg.csv",sep = ',', header = TRUE)
Scenic_LSC <- Scenic_LSC[-1:-3,1:6]
names(Scenic_LSC) <-c("TF","Motif","AUC","NES","MSQ","OI")
Scenic_LSC <- Scenic_LSC[Scenic_LSC$OI>=1,]
Scenic_LSC <- Scenic_LSC %>% group_by(Motif) %>% arrange(desc(NES),.by_group = T) %>% slice_head(n=1)


Scenic_IC <- Scenic_IC[,c("Motif","NES")]
Scenic_Mel <- Scenic_Mel[,c("Motif","NES")]
Scenic_LSC <- Scenic_LSC[,c("Motif","NES")]

Scenic_full <- merge(Scenic_IC,Scenic_Mel, by = 'Motif', all = TRUE)
Scenic_full <-merge(Scenic_full,Scenic_LSC, by = 'Motif', all = TRUE)

rownames(Scenic_full)<- Scenic_full$Motif
Scenic_full$Motif<- NULL
names(Scenic_full)<- c("Scenic_IC","Scenic_Mel","Scenic_LSC")


Scenic_full[is.na(Scenic_full)] <- 0

# Select the top 10 fators in scenic motif prediciton
vec2 <- NULL
for (i in colnames(Scenic_full)) {
  sel <- Scenic_full[order(Scenic_full[,i],decreasing = T),]
  vec <- rownames(sel)[1:10]
  print(vec)
  vec2 <- c(vec2,vec)
}
vec2 <- unique(vec2)

Scen_mat <- as.matrix(sapply(Scenic_full, as.numeric))  
rownames(Scen_mat) <- rownames(Scenic_full)
Scen_mat2 <-Scen_mat[rownames(Scen_mat) %in% vec2,]

Z <- t(scale(t(Scen_mat2)))

# Print the matrix of the Z-score NES
f6 = colorRamp2(c(-1, 0, 1), c("grey39", "#EEEEEE", "orange"), space = "RGB")

Z <- Z[,c("Scenic_LSC","Scenic_Mel","Scenic_IC")]
ht4 = Heatmap(Z, col = f6, cluster_columns = F,cluster_rows = T, name = "Z-score NES")

pdf(paste(resultsdir,'/complexheatmap_motifs_scenic.pdf',sep="/") ,width=6,height=6,paper='special')
htlist <-ht4
draw(htlist, column_title = "Z-score of SCENIC NES")
dev.off()

# Adding the TFs to the motifs
lakorna <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/Z-score_datasets20220225/Zscoretable_markers_split2.tsv", sep="\t",row.names=1)

lakorna <- na.omit(lakorna)

matz <- as.matrix(lakorna)

# Scenic motif enrichment analysis normalized enrichment score (NES) test
Scenic_IC <- read.table(file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/scenic/26042022_IC/reg.csv",sep = ",", header = TRUE)

Scenic_IC <- Scenic_IC[-1:-3,1:6]
names(Scenic_IC) <-c("TF","Motif","AUC","NES","MSQ","OI")
Scenic_IC <- Scenic_IC[Scenic_IC$OI>=1,]
Scenic_IC <- Scenic_IC %>% group_by(Motif) %>% arrange(desc(NES),.by_group = T) %>% slice_head(n=1)

Scenic_Mel <- read.table(file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/scenic/26042022_Mel/reg.csv",sep = ',', header = TRUE)
Scenic_Mel <- Scenic_Mel[-1:-3,1:6]
names(Scenic_Mel) <-c("TF","Motif","AUC","NES","MSQ","OI")
Scenic_Mel <- Scenic_Mel[Scenic_Mel$OI>=1,]
Scenic_Mel <- Scenic_Mel %>% group_by(Motif) %>% arrange(desc(NES),.by_group = T) %>% slice_head(n=1)

Scenic_LSC <- read.table(file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/scenic/26042022_LSC/reg.csv",sep = ',', header = TRUE)
Scenic_LSC <- Scenic_LSC[-1:-3,1:6]
names(Scenic_LSC) <-c("TF","Motif","AUC","NES","MSQ","OI")
Scenic_LSC <- Scenic_LSC[Scenic_LSC$OI>=1,]
Scenic_LSC <- Scenic_LSC %>% group_by(Motif) %>% arrange(desc(NES),.by_group = T) %>% slice_head(n=1)

Scenic_IC <- Scenic_IC[,c("TF","Motif")]
Scenic_Mel <- Scenic_Mel[,c("TF","Motif")]
Scenic_LSC <- Scenic_LSC[,c("TF","Motif")]


total <- rbind(Scenic_IC, Scenic_Mel,Scenic_LSC)
total <- unique.data.frame(total) 

total<- total[total$TF%in%rownames(matz),]

Scenic_full2 <- Scenic_full
Scenic_full2$Motif <- rownames(Scenic_full2)
comb_mat2 <- merge(x = Scenic_full2, y = total, by = "Motif", all = TRUE)

comb_mat2 <- na.omit(comb_mat2)
comb_mat2$TF <- toupper(comb_mat2$TF)
comb_mat2 <- comb_mat2[match(unique(comb_mat2$TF), comb_mat2$TF),]
comb_mat22 <- comb_mat2
#comb_mat22$Motif <- NULL

matz <- matz[rownames(matz)%in%comb_mat22$TF,]

matz <- matz[,c("LSC","Mel","IC")]

comb_mat22 <- comb_mat22[comb_mat22$TF%in%rownames(matz),]

rownames(comb_mat22)<-comb_mat22$TF
comb_mat22$TF<-NULL
comb_mat22$Motif <- NULL



comb_mat22 <- comb_mat22[rownames(matz),c("Scenic_LSC","Scenic_Mel","Scenic_IC")]
rn <- rownames(comb_mat22)
comb_mat22 <- as.matrix(sapply(comb_mat22, as.numeric))  
rownames(comb_mat22) <- rn


# Correlation score of rows
comb_mat22 <- t(scale(t(comb_mat22)))
sel_vec <- sapply(1:nrow(as.matrix(comb_mat22)), function(i) cor(as.matrix(comb_mat22)[i,], as.matrix(matz)[i,]))
#sel_vec <- sel_vec > 0.6
comb_mat22 <- as.data.frame(comb_mat22)
comb_mat22$cor <- sel_vec
comb_mat22$TF <- rownames(comb_mat22)

comb_mat22 <- comb_mat22[,c("cor","TF")]
comb_mat2 <- t(scale(t(comb_mat2)))
comb_mat_cor <- merge(x = comb_mat2, y = comb_mat22, by = "TF", all = TRUE)
# mark the most correlating TFs

library(dplyr)
comb_mat_cor <- comb_mat_cor %>% arrange(factor(Motif, as.character(unique(Motif))), desc(cor))
comb_mat_cor <- Reduce(rbind,                                 # Top N highest values by group
                       by(comb_mat_cor,
                          comb_mat_cor["Motif"],
                          head,
                          n = 1))

rownames(comb_mat_cor) <- comb_mat_cor$Motif
comb_mat_cor$Motif <- NULL

comb_mat_cor2<- comb_mat2[comb_mat2$Motif %in% rownames(Z),]
rownames(comb_mat_cor2)<-comb_mat_cor2$Motif
comb_mat_cor2$Motif<-NULL
Z<-unique.data.frame(Z)


comb_mat_cor2 <- comb_mat_cor[rownames(comb_mat_cor) %in% rownames(Scen_mat2),]

Z2 <- comb_mat_cor[rownames(comb_mat_cor) %in% rownames(Z),]

comb_mat_cor2 <- comb_mat_cor2[rownames(Z),]

vec1 <- Z2$TF


order_TFs <- rownames(Z2)
vec2 <- order(match(rownames(Z),order_TFs))[1:length(order_TFs)]

# Print the heatmap with factor associated with the motifs
ht4 = Heatmap(Z, col = f6, cluster_columns = F,cluster_rows = T, name = "Z-score NES")

pdf(paste(resultsdir,'/complexheatmap_motifs_scenic_motifs.pdf',sep="/") ,width=6,height=6,paper='special')
htlist <-ht4+ rowAnnotation(link = anno_mark(at =  vec2,labels = vec1))
draw(htlist, column_title = "Z-score of SCENIC NES")
dev.off()


# Scenic motif enrichment analysis normalized enrichment score (NES) test
Scenic_IC <- read.table(file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/scenic/26042022_IC/reg.csv",sep = ",", header = TRUE)

Scenic_IC <- Scenic_IC[-1:-3,1:6]
names(Scenic_IC) <-c("TF","Motif","AUC","NES","MSQ","OI")
Scenic_IC <- Scenic_IC[Scenic_IC$OI>=1,]
Scenic_IC <- Scenic_IC %>% group_by(TF) %>% arrange(desc(NES),.by_group = T) %>% slice_head(n=1)

Scenic_Mel <- read.table(file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/scenic/26042022_Mel/reg.csv",sep = ',', header = TRUE)
Scenic_Mel <- Scenic_Mel[-1:-3,1:6]
names(Scenic_Mel) <-c("TF","Motif","AUC","NES","MSQ","OI")
Scenic_Mel <- Scenic_Mel[Scenic_Mel$OI>=1,]
Scenic_Mel <- Scenic_Mel %>% group_by(TF) %>% arrange(desc(NES),.by_group = T) %>% slice_head(n=1)

Scenic_LSC <- read.table(file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/scenic/26042022_LSC/reg.csv",sep = ',', header = TRUE)
Scenic_LSC <- Scenic_LSC[-1:-3,1:6]
names(Scenic_LSC) <-c("TF","Motif","AUC","NES","MSQ","OI")
Scenic_LSC <- Scenic_LSC[Scenic_LSC$OI>=1,]
Scenic_LSC <- Scenic_LSC %>% group_by(TF) %>% arrange(desc(NES),.by_group = T) %>% slice_head(n=1)


Scenic_IC <- Scenic_IC[,c("TF","NES")]
Scenic_Mel <- Scenic_Mel[,c("TF","NES")]
Scenic_LSC <- Scenic_LSC[,c("TF","NES")]

Scenic_full <- merge(Scenic_IC,Scenic_Mel, by = 'TF', all = TRUE)
Scenic_full <-merge(Scenic_full,Scenic_LSC, by = 'TF', all = TRUE)

rownames(Scenic_full)<- Scenic_full$TF
Scenic_full$TF<- NULL
names(Scenic_full)<- c("Scenic_IC","Scenic_Mel","Scenic_LSC")


Scenic_full[is.na(Scenic_full)] <- 0

# Select the top 10 fators in scenic TF prediciton
vec2 <- NULL
for (i in colnames(Scenic_full)) {
  sel <- Scenic_full[order(Scenic_full[,i],decreasing = T),]
  vec <- rownames(sel)[1:10]
  print(vec)
  vec2 <- c(vec2,vec)
}
vec2 <- unique(vec2)

Scen_mat <- as.matrix(sapply(Scenic_full, as.numeric))  
rownames(Scen_mat) <- rownames(Scenic_full)
Scen_mat2 <-Scen_mat[rownames(Scen_mat) %in% vec2,]

Z <- t(scale(t(Scen_mat2)))

# Print the matrix of the Z-score NES
f6 = colorRamp2(c(-1, 0, 1), c("grey39", "#EEEEEE", "orange"), space = "RGB")

Z <- Z[,c("Scenic_LSC","Scenic_Mel","Scenic_IC")]
ht4 = Heatmap(Z, col = f6, cluster_columns = F,cluster_rows = T, name = "Z-score NES")

pdf(paste(resultsdir,'/complexheatmap_tfs_scenic.pdf',sep="/") ,width=6,height=6,paper='special')
htlist <-ht4
draw(htlist, column_title = "Z-score of SCENIC NES")
dev.off()

# Adding the TF expression
lakorna <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/Z-score_datasets20220225/Zscoretable_markers_split2.tsv", sep="\t",row.names=1)

lakorna <- na.omit(lakorna)

matz <- as.matrix(lakorna)

matz <- matz[rownames(matz)%in%rownames(Z),]

matz <- matz[,c("LSC","Mel","IC")]

matz<-matz[rownames(Z),]

ht4 = Heatmap(Z, col = f6, cluster_columns = F,cluster_rows = T, name = "Z-score NES")
Z <- Z[row_order(ht4),]
matz<-matz[rownames(Z),]

ht4 = Heatmap(Z, col = f6, cluster_columns = F,cluster_rows = F, name = "Z-score NES")
ht5 = Heatmap(matz, col = f1, cluster_columns = F,cluster_rows = F, name = "Relative expression")

pdf(paste(resultsdir,'/complexheatmap_tfs_scenic_expression.pdf',sep="/") ,width=6,height=6,paper='special')
htlist <-ht4 +ht5
draw(htlist, column_title = "Z-score of SCENIC NES")
dev.off()


#### analysis of epigenomics factors
epig_factors <- read.table(file = "genes.txt",sep = '\t',header = T)

epig_factors <- epig_factors %>% group_by(Function) %>% arrange(desc(Function),.by_group = T)

epig_factors2 <- epig_factors[epig_factors$Protein.complex!="-",]

mat_RNA <- lakorna[rownames(lakorna) %in% epig_factors2$HGNC.approved.symbol,]

mat_RNA <- mat_RNA[order(match(rownames(mat_RNA),epig_factors2$HGNC.approved.symbol )), , drop = FALSE]

mat_RNA <- mat_RNA[,c("LSC","LESC","LE","CE","Cj","CSSC","CF")]
#mat_RNA <- mat_RNA[,]
mat_RNA <- as.matrix(mat_RNA)
ht3 = Heatmap(mat_RNA, col = f5, cluster_columns = F,cluster_rows = T, name = "Relative expression")
ht3
