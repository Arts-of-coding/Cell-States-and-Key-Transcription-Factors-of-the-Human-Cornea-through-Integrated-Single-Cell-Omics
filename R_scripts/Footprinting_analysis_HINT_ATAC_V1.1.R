# Deseq2 significance calculation peaks
##based upon http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
library(dplyr)
library(stringr)
library(gridExtra)

# loading in the count matrix and coldata files

# For the deseq2 matrix vs ESC
workdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/hint_atac/"

## setting up results directory
dateoftoday <- gsub("-", "", as.character(Sys.Date()))
resultsdir <- paste0(workdir, dateoftoday)
system(paste("mkdir -p ", resultsdir))

motif_p53 <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/hint_atac/DiffFootprinting_3/Lineplots/GM.5.0.p53.0001.txt', sep = '\t', header = TRUE)
motif_p53 <- motif_p53[51:150,]
rownames(motif_p53)<--49:50

# Convert df to long format
# install.packages("reshape")
library(reshape)

df <- data.frame(x = seq_along(motif_p53[, 1]),
                 motif_p53)

# Long format
df <- melt(df, id.vars = "x")
df$x <- rep(-49:50,length(unique(df$variable)))

motif <- "GM.5.0.p53.0001"
calc_motif <- function(i) {
  number <- read.table(file = paste0("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/hint_atac/MotifMatching/",i,"_mpbs.bed"))
  sum(number$V4==motif)
}

cells <- c("LSC","LESC","LE","Cj","CE","CSSC","CF")
df2 <- sapply(cells, calc_motif)

library(ggplot2)

colsordered<-c("#83FF52CC","#9AFE42CC" ,"#26EDA5CC","#3E3891CC","#455CD0CC","#F3C83ACC","#FE992CCC")#,"#FEAA33CC","#FE992CCC","#36AAF9CC","#7A0403CC","#FABA39CC","#1AE4B6CC","#A41301CC","#900A01CC",NA)

lp <- ggplot(df, aes(x = x, y = value, color = variable)) +
  geom_line()+
  scale_color_manual(values = colsordered,name = "variable")

i <- "LSC"

calc_ends <- function(i) {
  df[df$x==45&df$variable==i,]$value
}

cells <- c("LSC","LESC","LE","Cj","CE","CSSC","CF")
df3 <- sapply(cells, calc_ends)

annotation <- data.frame(
  x = rep(53,7),
  y = df3,
  label = df2
)

# Determine the motif length to determine the bar where the motif is located
number <- read.table(file = paste0("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/hint_atac/MotifMatching/",i,"_mpbs.bed"))
df_length<-number[number$V4==motif,c("V2","V3")]
df_length <-df_length[1,]
move_from_0 <-(df_length$V3-df_length$V2)/2

# Inspect plot below to see where the line ends 
# 0.265 for LE: LE > LESC
annotation$y <- c(0.3,0.255,0.265,0.24,0.23,0.20,0.18)

pdf(paste0(resultsdir,'/hint_atac_p53.pdf') ,width=6,height=5,paper='special')

print(lp + geom_text(data=annotation, aes( x=x+2, y=y, label=label),                 
            color="black", 
            size=4) +
  geom_vline(xintercept=c(0-move_from_0,0+move_from_0), linetype="longdash")+
    #xlim(-50, 50)+
    #scale_x_continuous(limits=c(-51, 51)) +
#annotate("rect", xmin = 0-move_from_0, xmax = 0+move_from_0, ymin = min(df$value), ymax = max(df$value)+0.05,alpha = .2) +
  labs(y= "Number of reads", x = "Distance from motif center")+ theme_minimal(base_size = 16))

print(lp +
        geom_vline(xintercept=c(0-move_from_0,0+move_from_0), linetype="longdash")+
        #xlim(-50, 50)+
        #scale_x_continuous(limits=c(-51, 51)) +
        #annotate("rect", xmin = 0-move_from_0, xmax = 0+move_from_0, ymin = min(df$value), ymax = max(df$value)+0.05,alpha = .2) +
        labs(y= "Number of reads", x = "Distance from motif center")+ theme_minimal(base_size = 16))

annotation_plot <- annotation["label"]
names(annotation_plot) <- "#"

tt2 <- ttheme_minimal(base_size = 6,core=list(fg_params=list(hjust=1, x=0.9)),
                      rowhead=list(fg_params=list(hjust=1, x=0.95)))
g2 <-tableGrob(t(annotation_plot), theme=tt2)
g2$widths <- unit(rep(0.6, ncol(g2)), "npc")

print(lp +
        geom_vline(xintercept=c(0-move_from_0,0+move_from_0), linetype="longdash")+
        #xlim(-50, 50)+
        #scale_x_continuous(limits=c(-51, 51)) +
        #annotate("rect", xmin = 0-move_from_0, xmax = 0+move_from_0, ymin = min(df$value), ymax = max(df$value)+0.05,alpha = .2) +
        labs(y= "Number of reads", x = "Distance from motif center")+ theme_minimal(base_size = 16)
      +
        #annotate(geom = "table", x = -50, y = 0.20, label = list(annotation_plot), 
                # vjust = 1, hjust = 0)+ 
        annotation_custom(g2, xmin=25, xmax=35, ymin=0.075, ymax=0.125))


dev.off()

###################################
# PAX6
motif <- "GM.5.0.Homeodomain_Paired_box.0002"
motif_pax6 <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/hint_atac/DiffFootprinting_3/Lineplots/GM.5.0.Homeodomain_Paired_box.0002.txt', sep = '\t', header = TRUE)
motif_pax6 <- motif_pax6[51:150,]
rownames(motif_pax6)<--49:50

# Convert df to long format
# install.packages("reshape")
library(reshape)

df <- NULL
df2 <- NULL
df3 <- NULL
df <- data.frame(x = seq_along(motif_pax6[, 1]),
                 motif_pax6)

# Long format
df <- melt(df, id.vars = "x")
df$x <- rep(-49:50,length(unique(df$variable)))

df2 <- sapply(cells, calc_motif)

calc_ends <- function(i) {
  df[df$x==45&df$variable==i,]$value
}

cells <- c("LSC","LESC","LE","Cj","CE","CSSC","CF")
df3 <- sapply(cells, calc_ends)

annotation <- data.frame(
  x = rep(53,7),
  y = df3,
  label = df2
)

# Determine the motif length to determine the bar where the motif is located
number <- read.table(file = paste0("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/hint_atac/MotifMatching/",i,"_mpbs.bed"))
df_length<-number[number$V4==motif,c("V2","V3")]
df_length <-df_length[1,]
move_from_0 <-(df_length$V3-df_length$V2)/2

# Inspect plot below to see where the line ends 
annotation$y <- c(0.22,0.21,0.24,0.23,0.26,0.25,0.20)

lp <- ggplot(df, aes(x = x, y = value, color = variable)) +
  geom_line()+
  scale_color_manual(values = colsordered,name = "variable")

pdf(paste0(resultsdir,'/hint_atac_pax6.pdf') ,width=6,height=5,paper='special')
print(lp + geom_text(data=annotation, aes( x=x+2, y=y, label=label),                 
                     color="black", 
                     size=4) +
        geom_vline(xintercept=c(0-move_from_0,0+move_from_0), linetype="longdash")+
        labs(y= "Number of reads", x = "Distance from motif center")+ theme_minimal(base_size = 16))
print(lp +
        geom_vline(xintercept=c(0-move_from_0,0+move_from_0), linetype="longdash")+
        #xlim(-50, 50)+
        #scale_x_continuous(limits=c(-51, 51)) +
        #annotate("rect", xmin = 0-move_from_0, xmax = 0+move_from_0, ymin = min(df$value), ymax = max(df$value)+0.05,alpha = .2) +
        labs(y= "Number of reads", x = "Distance from motif center")+ theme_minimal(base_size = 16))

annotation_plot <- annotation["label"]
names(annotation_plot) <- "#"

tt2 <- ttheme_minimal(base_size = 6,core=list(fg_params=list(hjust=1, x=0.9)),
                      rowhead=list(fg_params=list(hjust=1, x=0.95)))
g2 <-tableGrob(t(annotation_plot), theme=tt2)
g2$widths <- unit(rep(0.6, ncol(g2)), "npc")

print(lp +
        geom_vline(xintercept=c(0-move_from_0,0+move_from_0), linetype="longdash")+
        #xlim(-50, 50)+
        #scale_x_continuous(limits=c(-51, 51)) +
        #annotate("rect", xmin = 0-move_from_0, xmax = 0+move_from_0, ymin = min(df$value), ymax = max(df$value)+0.05,alpha = .2) +
        labs(y= "Number of reads", x = "Distance from motif center")+ theme_minimal(base_size = 16)
      +
        #annotate(geom = "table", x = -50, y = 0.20, label = list(annotation_plot), 
        # vjust = 1, hjust = 0)+ 
        annotation_custom(g2, xmin=25, xmax=35, ymin=0.1, ymax=0.15))

dev.off()

###################################
# smad3
motif <- "GM.5.0.SMAD.0007"
motif_tab <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/hint_atac/DiffFootprinting/Lineplots/GM.5.0.SMAD.0007.txt', sep = '\t', header = TRUE)
motif_tab <- motif_tab[51:150,]
rownames(motif_tab)<--49:50

# Convert df to long format
# install.packages("reshape")
library(reshape)

df <- NULL
df2 <- NULL
df3 <- NULL
df <- data.frame(x = seq_along(motif_tab[, 1]),
                 motif_tab)

# Long format
df <- melt(df, id.vars = "x")
df$x <- rep(-49:50,length(unique(df$variable)))

df2 <- sapply(cells, calc_motif)

calc_ends <- function(i) {
  df[df$x==45&df$variable==i,]$value
}

cells <- c("LSC","LESC","LE","Cj","CE","CSSC","CF")
df3 <- sapply(cells, calc_ends)

annotation <- data.frame(
  x = rep(53,7),
  y = df3,
  label = df2
)

# Determine the motif length to determine the bar where the motif is located
number <- read.table(file = paste0("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/hint_atac/MotifMatching/",i,"_mpbs.bed"))
df_length<-number[number$V4==motif,c("V2","V3")]
df_length <-df_length[1,]
move_from_0 <-(df_length$V3-df_length$V2)/2

# Inspect plot below to see where the line ends 
annotation$y <- c(0.22,0.21,0.24,0.23,0.26,0.25,0.20)

lp <- ggplot(df, aes(x = x, y = value, color = variable)) +
  geom_line()+
  scale_color_manual(values = colsordered,name = "variable")

pdf(paste0(resultsdir,'/hint_atac_smad3.pdf') ,width=6,height=5,paper='special')
print(lp + geom_text(data=annotation, aes( x=x+2, y=y, label=label),                 
                     color="black", 
                     size=4) +
        geom_vline(xintercept=c(0-move_from_0,0+move_from_0), linetype="longdash")+
        labs(y= "Number of reads", x = "Distance from motif center")+ theme_minimal(base_size = 16))
print(lp +
        geom_vline(xintercept=c(0-move_from_0,0+move_from_0), linetype="longdash")+
        #xlim(-50, 50)+
        #scale_x_continuous(limits=c(-51, 51)) +
        #annotate("rect", xmin = 0-move_from_0, xmax = 0+move_from_0, ymin = min(df$value), ymax = max(df$value)+0.05,alpha = .2) +
        labs(y= "Number of reads", x = "Distance from motif center")+ theme_minimal(base_size = 16))

annotation_plot <- annotation["label"]
names(annotation_plot) <- "#"

tt2 <- ttheme_minimal(base_size = 6,core=list(fg_params=list(hjust=1, x=0.9)),
                      rowhead=list(fg_params=list(hjust=1, x=0.95)))
g2 <-tableGrob(t(annotation_plot), theme=tt2)
g2$widths <- unit(rep(0.6, ncol(g2)), "npc")

print(lp +
        geom_vline(xintercept=c(0-move_from_0,0+move_from_0), linetype="longdash")+
        #xlim(-50, 50)+
        #scale_x_continuous(limits=c(-51, 51)) +
        #annotate("rect", xmin = 0-move_from_0, xmax = 0+move_from_0, ymin = min(df$value), ymax = max(df$value)+0.05,alpha = .2) +
        labs(y= "Number of reads", x = "Distance from motif center")+ theme_minimal(base_size = 16)
      +
        #annotate(geom = "table", x = -50, y = 0.20, label = list(annotation_plot), 
        # vjust = 1, hjust = 0)+ 
        annotation_custom(g2, xmin=25, xmax=35, ymin=0.075, ymax=0.125))
dev.off()

###################################
# twist1
motif <- "GM.5.0.bHLH.0065"
motif_tab <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/hint_atac/DiffFootprinting/Lineplots/GM.5.0.bHLH.0065.txt', sep = '\t', header = TRUE)
motif_tab <- motif_tab[51:150,]
rownames(motif_tab)<--49:50

# Convert df to long format
# install.packages("reshape")
library(reshape)

df <- NULL
df2 <- NULL
df3 <- NULL
df <- data.frame(x = seq_along(motif_tab[, 1]),
                 motif_tab)

# Long format
df <- melt(df, id.vars = "x")
df$x <- rep(-49:50,length(unique(df$variable)))

df2 <- sapply(cells, calc_motif)

calc_ends <- function(i) {
  df[df$x==45&df$variable==i,]$value
}

cells <- c("LSC","LESC","LE","Cj","CE","CSSC","CF")
df3 <- sapply(cells, calc_ends)

annotation <- data.frame(
  x = rep(53,7),
  y = df3,
  label = df2
)

# Determine the motif length to determine the bar where the motif is located
number <- read.table(file = paste0("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/hint_atac/MotifMatching/",i,"_mpbs.bed"))
df_length<-number[number$V4==motif,c("V2","V3")]
df_length <-df_length[1,]
move_from_0 <-(df_length$V3-df_length$V2)/2

# Inspect plot below to see where the line ends 
annotation$y <- c(0.10,0.11,0.14,0.12,0.13,0.21,0.22)

lp <- ggplot(df, aes(x = x, y = value, color = variable)) +
  geom_line()+
  scale_color_manual(values = colsordered,name = "variable")

pdf(paste0(resultsdir,'/hint_atac_twist1.pdf') ,width=6,height=5,paper='special')
print(lp + geom_text(data=annotation, aes( x=x+2, y=y, label=label),                 
                     color="black", 
                     size=4) +
        geom_vline(xintercept=c(0-move_from_0,0+move_from_0), linetype="longdash")+
        labs(y= "Number of reads", x = "Distance from motif center")+ theme_minimal(base_size = 16))
print(lp +
        geom_vline(xintercept=c(0-move_from_0,0+move_from_0), linetype="longdash")+
        #xlim(-50, 50)+
        #scale_x_continuous(limits=c(-51, 51)) +
        #annotate("rect", xmin = 0-move_from_0, xmax = 0+move_from_0, ymin = min(df$value), ymax = max(df$value)+0.05,alpha = .2) +
        labs(y= "Number of reads", x = "Distance from motif center")+ theme_minimal(base_size = 16))

annotation_plot <- annotation["label"]
names(annotation_plot) <- "#"

tt2 <- ttheme_minimal(base_size = 6,core=list(fg_params=list(hjust=1, x=0.9)),
                      rowhead=list(fg_params=list(hjust=1, x=0.95)))
g2 <-tableGrob(t(annotation_plot), theme=tt2)
g2$widths <- unit(rep(0.6, ncol(g2)), "npc")

print(lp +
        geom_vline(xintercept=c(0-move_from_0,0+move_from_0), linetype="longdash")+
        #xlim(-50, 50)+
        #scale_x_continuous(limits=c(-51, 51)) +
        #annotate("rect", xmin = 0-move_from_0, xmax = 0+move_from_0, ymin = min(df$value), ymax = max(df$value)+0.05,alpha = .2) +
        labs(y= "Number of reads", x = "Distance from motif center")+ theme_minimal(base_size = 16)
      +
        #annotate(geom = "table", x = -50, y = 0.20, label = list(annotation_plot), 
        # vjust = 1, hjust = 0)+ 
        annotation_custom(g2, xmin=25, xmax=35, ymin=0.065, ymax=0.115))
dev.off()

###################################
# grhl1
motif <- "GM.5.0.Grainyhead.0002"
motif_tab <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/hint_atac/DiffFootprinting/Lineplots/GM.5.0.Grainyhead.0002.txt', sep = '\t', header = TRUE)
motif_tab <- motif_tab[51:150,]
rownames(motif_tab)<--49:50

# Convert df to long format
# install.packages("reshape")
library(reshape)

df <- NULL
df2 <- NULL
df3 <- NULL
df <- data.frame(x = seq_along(motif_tab[, 1]),
                 motif_tab)

# Long format
df <- melt(df, id.vars = "x")
df$x <- rep(-49:50,length(unique(df$variable)))

df2 <- sapply(cells, calc_motif)

calc_ends <- function(i) {
  df[df$x==45&df$variable==i,]$value
}

cells <- c("LSC","LESC","LE","Cj","CE","CSSC","CF")
df3 <- sapply(cells, calc_ends)

annotation <- data.frame(
  x = rep(53,7),
  y = df3,
  label = df2
)

# Determine the motif length to determine the bar where the motif is located
number <- read.table(file = paste0("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/hint_atac/MotifMatching/",i,"_mpbs.bed"))
df_length<-number[number$V4==motif,c("V2","V3")]
df_length <-df_length[1,]
move_from_0 <-(df_length$V3-df_length$V2)/2

# Inspect plot below to see where the line ends 
annotation$y <- c(0.22,0.28,0.25,0.30,0.33,0.15,0.13)

lp <- ggplot(df, aes(x = x, y = value, color = variable)) +
  geom_line()+
  scale_color_manual(values = colsordered,name = "variable")

pdf(paste0(resultsdir,'/hint_atac_grhl1.pdf') ,width=6,height=5,paper='special')
print(lp + geom_text(data=annotation, aes( x=x+2, y=y, label=label),                 
                     color="black", 
                     size=4) +
        geom_vline(xintercept=c(0-move_from_0,0+move_from_0), linetype="longdash")+
        labs(y= "Number of reads", x = "Distance from motif center")+ theme_minimal(base_size = 16))
print(lp +
        geom_vline(xintercept=c(0-move_from_0,0+move_from_0), linetype="longdash")+
        #xlim(-50, 50)+
        #scale_x_continuous(limits=c(-51, 51)) +
        #annotate("rect", xmin = 0-move_from_0, xmax = 0+move_from_0, ymin = min(df$value), ymax = max(df$value)+0.05,alpha = .2) +
        labs(y= "Number of reads", x = "Distance from motif center")+ theme_minimal(base_size = 16))

annotation_plot <- annotation["label"]
names(annotation_plot) <- "#"

tt2 <- ttheme_minimal(base_size = 6,core=list(fg_params=list(hjust=1, x=0.9)),
                      rowhead=list(fg_params=list(hjust=1, x=0.95)))
g2 <-tableGrob(t(annotation_plot), theme=tt2)
g2$widths <- unit(rep(0.6, ncol(g2)), "npc")

print(lp +
        geom_vline(xintercept=c(0-move_from_0,0+move_from_0), linetype="longdash")+
        #xlim(-50, 50)+
        #scale_x_continuous(limits=c(-51, 51)) +
        #annotate("rect", xmin = 0-move_from_0, xmax = 0+move_from_0, ymin = min(df$value), ymax = max(df$value)+0.05,alpha = .2) +
        labs(y= "Number of reads", x = "Distance from motif center")+ theme_minimal(base_size = 16)
      +
        #annotate(geom = "table", x = -50, y = 0.20, label = list(annotation_plot), 
        # vjust = 1, hjust = 0)+ 
        annotation_custom(g2, xmin=25, xmax=35, ymin=0.075, ymax=0.125))
dev.off()

###################################
# NFIC
motif <- "GM.5.0.SMAD.0002"
motif_tab <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/hint_atac/DiffFootprinting_0707/Lineplots/GM.5.0.SMAD.0002.txt', sep = '\t', header = TRUE)
motif_tab <- motif_tab[51:150,]
rownames(motif_tab)<--49:50

# Convert df to long format
# install.packages("reshape")
library(reshape)

df <- NULL
df2 <- NULL
df3 <- NULL
df <- data.frame(x = seq_along(motif_tab[, 1]),
                 motif_tab)

# Long format
df <- melt(df, id.vars = "x")
df$x <- rep(-49:50,length(unique(df$variable)))

calc_motif <- function(i) {
  number <- read.table(file = paste0("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/hint_atac/MotifMatching0707/",i,"_mpbs.bed"))
  sum(number$V4==motif)
}

df2 <- sapply(cells, calc_motif)

calc_ends <- function(i) {
  df[df$x==45&df$variable==i,]$value
}

cells <- c("LSC","LESC","LE","Cj","CE","CSSC","CF")
df3 <- sapply(cells, calc_ends)

annotation <- data.frame(
  x = rep(53,7),
  y = df3,
  label = df2
)

# Determine the motif length to determine the bar where the motif is located
number <- read.table(file = paste0("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/hint_atac/MotifMatching0707/",i,"_mpbs.bed"))
df_length<-number[number$V4==motif,c("V2","V3")]
df_length <-df_length[1,]
move_from_0 <-(df_length$V3-df_length$V2)/2

# Inspect plot below to see where the line ends 
annotation$y <- c(0.20,0.14,0.15,0.18,0.17,0.16,0.19)

lp <- ggplot(df, aes(x = x, y = value, color = variable)) +
  geom_line()+
  scale_color_manual(values = colsordered,name = "variable")

pdf(paste0(resultsdir,'/hint_atac_NFIC.pdf') ,width=6.5,height=5,paper='special')
print(lp + geom_text(data=annotation, aes( x=x+2, y=y, label=label),                 
                     color="black", 
                     size=4) +
        geom_vline(xintercept=c(0-move_from_0,0+move_from_0), linetype="longdash")+
        labs(y= "Number of reads", x = "Distance from motif center")+ theme_minimal(base_size = 16))
print(lp +
        geom_vline(xintercept=c(0-move_from_0,0+move_from_0), linetype="longdash")+
        #xlim(-50, 50)+
        #scale_x_continuous(limits=c(-51, 51)) +
        #annotate("rect", xmin = 0-move_from_0, xmax = 0+move_from_0, ymin = min(df$value), ymax = max(df$value)+0.05,alpha = .2) +
        labs(y= "Number of reads", x = "Distance from motif center")+ theme_minimal(base_size = 16))

annotation_plot <- annotation["label"]
names(annotation_plot) <- "#"

tt2 <- ttheme_minimal(base_size = 6,core=list(fg_params=list(hjust=1, x=0.9)),
                      rowhead=list(fg_params=list(hjust=1, x=0.95)))
g2 <-tableGrob(t(annotation_plot), theme=tt2)
g2$widths <- unit(rep(0.6, ncol(g2)), "npc")

print(lp +
        geom_vline(xintercept=c(0-move_from_0,0+move_from_0), linetype="longdash")+
        #xlim(-50, 50)+
        #scale_x_continuous(limits=c(-51, 51)) +
        #annotate("rect", xmin = 0-move_from_0, xmax = 0+move_from_0, ymin = min(df$value), ymax = max(df$value)+0.05,alpha = .2) +
        labs(y= "Number of reads", x = "Distance from motif center")+ theme_minimal(base_size = 16)
      +
        #annotate(geom = "table", x = -50, y = 0.20, label = list(annotation_plot), 
        # vjust = 1, hjust = 0)+ 
        annotation_custom(g2, xmin=25, xmax=35, ymin=0.075, ymax=0.125))
dev.off()

###################################
# OTX1
motif <- "GM.5.0.Homeodomain.0097"
motif_tab <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/hint_atac/DiffFootprinting_0707/Lineplots/GM.5.0.Homeodomain.0097.txt', sep = '\t', header = TRUE)
motif_tab <- motif_tab[51:150,]
rownames(motif_tab)<--49:50

# Convert df to long format
# install.packages("reshape")
library(reshape)

df <- NULL
df2 <- NULL
df3 <- NULL
df <- data.frame(x = seq_along(motif_tab[, 1]),
                 motif_tab)

# Long format
df <- melt(df, id.vars = "x")
df$x <- rep(-49:50,length(unique(df$variable)))

calc_motif <- function(i) {
  number <- read.table(file = paste0("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/hint_atac/MotifMatching0707/",i,"_mpbs.bed"))
  sum(number$V4==motif)
}

df2 <- sapply(cells, calc_motif)

calc_ends <- function(i) {
  df[df$x==45&df$variable==i,]$value
}

cells <- c("LSC","LESC","LE","Cj","CE","CSSC","CF")
df3 <- sapply(cells, calc_ends)

annotation <- data.frame(
  x = rep(53,7),
  y = df3,
  label = df2
)

# Determine the motif length to determine the bar where the motif is located
number <- read.table(file = paste0("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/hint_atac/MotifMatching0707/",i,"_mpbs.bed"))
df_length<-number[number$V4==motif,c("V2","V3")]
df_length <-df_length[1,]
move_from_0 <-(df_length$V3-df_length$V2)/2

# Inspect plot below to see where the line ends 
annotation$y <- c(0.18,0.15,0.16,0.17,0.19,0.13,0.14)

lp <- ggplot(df, aes(x = x, y = value, color = variable)) +
  geom_line()+
  scale_color_manual(values = colsordered,name = "variable")

pdf(paste0(resultsdir,'/hint_atac_OTX1.pdf') ,width=6,height=5,paper='special')
print(lp + geom_text(data=annotation, aes( x=x+2, y=y, label=label),                 
                     color="black", 
                     size=4) +
        geom_vline(xintercept=c(0-move_from_0,0+move_from_0), linetype="longdash")+
        labs(y= "Number of reads", x = "Distance from motif center")+ theme_minimal(base_size = 16))
print(lp +
        geom_vline(xintercept=c(0-move_from_0,0+move_from_0), linetype="longdash")+
        #xlim(-50, 50)+
        #scale_x_continuous(limits=c(-51, 51)) +
        #annotate("rect", xmin = 0-move_from_0, xmax = 0+move_from_0, ymin = min(df$value), ymax = max(df$value)+0.05,alpha = .2) +
        labs(y= "Number of reads", x = "Distance from motif center")+ theme_minimal(base_size = 16))

annotation_plot <- annotation["label"]
names(annotation_plot) <- "#"

tt2 <- ttheme_minimal(base_size = 6,core=list(fg_params=list(hjust=1, x=0.9)),
                      rowhead=list(fg_params=list(hjust=1, x=0.95)))
g2 <-tableGrob(t(annotation_plot), theme=tt2)
g2$widths <- unit(rep(0.6, ncol(g2)), "npc")

print(lp +
        geom_vline(xintercept=c(0-move_from_0,0+move_from_0), linetype="longdash")+
        #xlim(-50, 50)+
        #scale_x_continuous(limits=c(-51, 51)) +
        #annotate("rect", xmin = 0-move_from_0, xmax = 0+move_from_0, ymin = min(df$value), ymax = max(df$value)+0.05,alpha = .2) +
        labs(y= "Number of reads", x = "Distance from motif center")+ theme_minimal(base_size = 16)
      +
        #annotate(geom = "table", x = -50, y = 0.20, label = list(annotation_plot), 
        # vjust = 1, hjust = 0)+ 
        annotation_custom(g2, xmin=25, xmax=35, ymin=0.05, ymax=0.1))
dev.off()