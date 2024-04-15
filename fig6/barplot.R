

# load libraries

library("Seurat")
library("dplyr")
library("ggplot2")
library("ComplexHeatmap")
library("circlize")
library("scales")
library("RColorBrewer")

DTC <- c(0,1)
DTC_renamed <- c("non-DTC","DTC")

# select palette

cols <- c(brewer.pal(4,"Greys")[2],brewer.pal(4,"Oranges")[4])
names(cols) <- DTC_renamed


### EXP1B

# directory and file names

DIR <- "../analysis/gene_expression"

EXP <- "1B_GEX"
EXP_NAME <- "EXP1B"

SAMPLES <- c("T3","T9","T13","T17")
SAMPLES_RENAMED <- c("d3","d9","d13","d17")

OBJECT <- file.path(DIR, EXP, "treated_filt/hvg_pca_clust/object.Rds")

BARPLOT <- paste0("barplot_CB_DTC_", EXP_NAME)

# load data and compute stats

object <- readRDS(OBJECT)

object@meta.data <- object@meta.data %>% mutate(sample.name = SAMPLES_RENAMED[match(object@meta.data$sample.name,SAMPLES)])
object@meta.data$sample.name <- factor(object@meta.data$sample.name, levels = SAMPLES_RENAMED)
object@meta.data <- object@meta.data %>% mutate(DRC_sc = DTC_renamed[match(object@meta.data$DRC_sc,DTC)])
object@meta.data$DRC_sc <- factor(object@meta.data$DRC_sc, levels = DTC_renamed)

n <- length(SAMPLES)
m <- length(DTC)
df <- as.data.frame(matrix(0, nrow = n*m, ncol = 3))
colnames(df) <- c("sample", "class", "count")
for (i in 1:length(SAMPLES)) {
    for (j in 1:length(DTC)) {
        c <- sum(object@meta.data$sample.name == SAMPLES_RENAMED[i] & object@meta.data$DRC_sc == DTC_renamed[j])
        df[m*(i-1)+j,] <- c(SAMPLES_RENAMED[i], DTC_renamed[j], c)
    }
}
df$class <- factor(df$class, levels = DTC_renamed)
df$sample <- factor(df$sample, levels = SAMPLES_RENAMED)
df$count <- as.numeric(df$count)

# figure with barplot colored by DTC

g <- ggplot(data = df, aes(x = sample, y = count, fill = class)) + theme_classic() + geom_bar(stat = "identity")
g <- g + xlab("") + ylab("Number of cells") + scale_fill_manual(values = cols)

BARPLOT_PLOT <- paste0(BARPLOT, ".pdf")
pdf(BARPLOT_PLOT, width = 3.5, height = 2.5)
g
dev.off()


### EXP1D

# directory and file names

DIR <- "../analysis/gene_expression"

EXP <- "1D_GEX"
EXP_NAME <- "EXP1D"

SAMPLES <- c("T5_2_Multiseq","T7_2_Multiseq","T9_2_Multiseq","T11_2_Multiseq","T13_2_Multiseq","T15_2_Multiseq")
SAMPLES_RENAMED <- c("d5","d7","d9","d11","d13","d15")

OBJECT <- file.path(DIR, EXP, "treated_filt/hvg_pca_clust/object.Rds")

BARPLOT <- paste0("barplot_CB_DTC_", EXP_NAME)

# load data and compute stats

object <- readRDS(OBJECT)

object@meta.data <- object@meta.data %>% mutate(sample.name = SAMPLES_RENAMED[match(object@meta.data$sample.name,SAMPLES)])
object@meta.data$sample.name <- factor(object@meta.data$sample.name, levels = SAMPLES_RENAMED)
object@meta.data <- object@meta.data %>% mutate(DRC_sc = DTC_renamed[match(object@meta.data$DRC_sc,DTC)])
object@meta.data$DRC_sc <- factor(object@meta.data$DRC_sc, levels = DTC_renamed)

n <- length(SAMPLES)
m <- length(DTC)
df <- as.data.frame(matrix(0, nrow = n*m, ncol = 3))
colnames(df) <- c("sample", "class", "count")
for (i in 1:length(SAMPLES)) {
    for (j in 1:length(DTC)) {
        c <- sum(object@meta.data$sample.name == SAMPLES_RENAMED[i] & object@meta.data$DRC_sc == DTC_renamed[j])
        df[m*(i-1)+j,] <- c(SAMPLES_RENAMED[i], DTC_renamed[j], c)
    }
}
df$class <- factor(df$class, levels = DTC_renamed)
df$sample <- factor(df$sample, levels = SAMPLES_RENAMED)
df$count <- as.numeric(df$count)

# figure with barplot colored by DTC

g <- ggplot(data = df, aes(x = sample, y = count, fill = class)) + theme_classic() + geom_bar(stat = "identity")
g <- g + xlab("") + ylab("Number of cells") + scale_fill_manual(values = cols)

BARPLOT_PLOT <- paste0(BARPLOT, ".pdf")
pdf(BARPLOT_PLOT, width = 4, height = 2.5)
g
dev.off()


### DTC vs clusters

INFO <- "../tables/DATA_FRAME_CLONES_treated.tsv"

SAMPLES_1B <- c("T13","T17")
SAMPLES_1B_RENAMED <- paste0("d",c(13,17),"_2")

SAMPLES_1D <- c("T13_2_Multiseq","T15_2_Multiseq")
SAMPLES_1D_RENAMED <- paste0("d",c(13,15),"_1")

SAMPLES <- c(SAMPLES_1D, SAMPLES_1B)
SAMPLES_RENAMED <- c(SAMPLES_1D_RENAMED, SAMPLES_1B_RENAMED)

CL_1B <- "cluster_1B_P0_filt_renamed"
CL_1D <- "cluster_1D_P0_filt_renamed"

PLOT_1B <- "barplot_lineage_multi_cl_1B.pdf"
PLOT_1D <- "barplot_lineage_multi_cl_1D.pdf"

PLOT_1B_AUG <- "barplot_lineage_augmented_multi_cl_1B.pdf"
PLOT_1D_AUG <- "barplot_lineage_augmented_multi_cl_1D.pdf"

# collect the normalized abundance of each GBC type in samples

df_info <- read.table(INFO, sep = "\t", header = TRUE)
idx <- match(paste0("cpt.",SAMPLES), colnames(df_info))
colnames(df_info)[idx] <- SAMPLES_RENAMED

df <- as.data.frame(matrix(NA, nrow = 0, ncol = 6))
colnames(df) <- c("cl_1B","cl_1D","lineage","cpt","sample","lineage_augmented")
i <- 1
for (s in SAMPLES_RENAMED) {
    df_tmp <- df_info[,c(CL_1B,CL_1D,"lineage",s)]
    v <- rep(s, nrow(df_tmp))
    w <- df_info[,paste0("lineage_augmented.", SAMPLES[i])]
    df_tmp <- cbind(df_tmp, v, w)
    colnames(df_tmp) <- colnames(df)
    df <- rbind(df, df_tmp)
    i <- i+1
}

idx <- which(df$cl_1B == -1)
df$cl_1B[idx] <- rep("unmapped", length(idx))
idx <- which(df$cl_1D == -1)
df$cl_1D[idx] <- rep("unmapped", length(idx))

cl_1B <- sort(unique(df$cl_1B))
cl_1D <- sort(unique(df$cl_1D))
ncl_1B <- length(cl_1B)
ncl_1D <- length(cl_1D)

samples <- sort(unique(df$sample))
df$sample <- factor(df$sample, levels = SAMPLES_RENAMED)

lin <- sort(unique(df$lineage[!is.na(df$lineage)]))

df <- df[!is.na(df$lineage_augmented),]
df$lineage_augmented <- factor(df$lineage_augmented, levels = lin)

df_lin <- df[!is.na(df$lineage),]
df_lin$lineage <- factor(df_lin$lineage, levels = lin)

df$cl_1B <- factor(df$cl_1B, levels = rev(cl_1B))
df$cl_1D <- factor(df$cl_1D, levels = rev(cl_1D))

df_lin$cl_1B <- factor(df_lin$cl_1B, levels = rev(cl_1B))
df_lin$cl_1D <- factor(df_lin$cl_1D, levels = rev(cl_1D))

col_1B <- rev(c(brewer.pal(n = ncl_1B-1, name = "Set1"),"gray90"))
names(col_1B) <- rev(cl_1B)
col_1D <- rev(c(brewer.pal(n = ncl_1D-1, name = "Set2"),"gray90"))
names(col_1D) <- rev(cl_1D)

g <- ggplot(data = df_lin, aes(x = lineage, y = cpt, fill = cl_1B)) 
g <- g + theme_classic() + ylab("CPT fraction")
g <- g + geom_bar(stat = "identity", position = "fill") + scale_fill_manual(values = col_1B, drop = FALSE)
g <- g + facet_wrap(~ sample, ncol = 4)

pdf(PLOT_1B, width = 8, height = 2.5)
print(g)
dev.off()

g <- ggplot(data = df_lin, aes(x = lineage, y = cpt, fill = cl_1D)) 
g <- g + theme_classic() + ylab("CPT fraction")
g <- g + geom_bar(stat = "identity", position = "fill") + scale_fill_manual(values = col_1D, drop = FALSE)
g <- g + facet_wrap(~ sample, ncol = 4)

pdf(PLOT_1D, width = 8, height = 2.5)
print(g)
dev.off()

g <- ggplot(data = df, aes(x = lineage_augmented, y = cpt, fill = cl_1B)) 
g <- g + theme_classic() + ylab("CPT fraction")
g <- g + geom_bar(stat = "identity", position = "fill") + scale_fill_manual(values = col_1B, drop = FALSE)
g <- g + facet_wrap(~ sample, ncol = 4)

pdf(PLOT_1B_AUG, width = 8, height = 2.5)
print(g)
dev.off()

g <- ggplot(data = df, aes(x = lineage_augmented, y = cpt, fill = cl_1D)) 
g <- g + theme_classic() + ylab("CPT fraction")
g <- g + geom_bar(stat = "identity", position = "fill") + scale_fill_manual(values = col_1D, drop = FALSE)
g <- g + facet_wrap(~ sample, ncol = 4)

pdf(PLOT_1D_AUG, width = 8, height = 2.5)
print(g)
dev.off()


q()

