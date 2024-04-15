
# load libraries

library("Seurat")
library("dplyr")
library("ggplot2")
library("ComplexHeatmap")
library("circlize")
library("scales")
library("RColorBrewer")


### EXP1B

# directory and file names

EXP_NAME <- "EXP1B"

INFO <- paste0("../tables/DATA_FRAME_",EXP_NAME,"_treated_filt.tsv")

UMAP <- paste0("umap_DTC_treated_",EXP_NAME,".pdf")

# load data

df <- read.table(INFO, sep = "\t", header = TRUE)

# select palette

col <- c(brewer.pal(4,"Greys")[2],brewer.pal(4,"Oranges")[4])
names(col) <- c(0,1)

# figure with umap colored by cluster

df$sample.name <- factor(df$sample.name, levels = unique(df$sample.name))
df$DRC_bulk_fisher <- factor(df$DRC_bulk_fisher, levels = c(0,1))

g <- ggplot(df, aes(x = -UMAP_1, y = UMAP_2, color = DRC_bulk_fisher)) + theme_classic() + geom_point(size = 0.2)
g <- g + scale_color_manual(values = col) + facet_wrap( ~ sample.name, nrow = 1)

pdf(UMAP, width = 10, height = 3)
g
dev.off()


UMAP <- paste0("umap_lineage_treated_",EXP_NAME,".pdf")

# load data

df <- read.table(INFO, sep = "\t", header = TRUE)
df$lineage[is.na(df$lineage)] <- rep(0, sum(is.na(df$lineage)))

# select palette

col <- c(brewer.pal(4,"Greys")[2],hue_pal()(2))

# figure with umap colored by lineage

df$sample.name <- factor(df$sample.name, levels = unique(df$sample.name))
df$lineage <- factor(df$lineage, levels = c(0,1,2))

g <- ggplot(df, aes(x = -UMAP_1, y = UMAP_2, color = lineage)) + theme_classic() + geom_point(size = 0.2)
g <- g + scale_color_manual(values = col) + facet_wrap( ~ sample.name, nrow = 1)

pdf(UMAP, width = 10, height = 3)
g
dev.off()



### EXP1D

# directory and file names

EXP_NAME <- "EXP1D"

INFO <- paste0("../tables/DATA_FRAME_",EXP_NAME,"_treated_filt.tsv")

UMAP <- paste0("umap_DTC_treated_",EXP_NAME,".pdf")

# load data

df <- read.table(INFO, sep = "\t", header = TRUE)

# select palette

col <- c(brewer.pal(4,"Greys")[2],brewer.pal(4,"Oranges")[4])
names(col) <- c(0,1)

# figure with umap colored by cluster

df$sample.name <- factor(df$sample.name, levels = unique(df$sample.name))
df$DRC_bulk_fisher <- factor(df$DRC_bulk_fisher, levels = c(0,1))

g <- ggplot(df, aes(x = -UMAP_1, y = UMAP_2, color = DRC_bulk_fisher)) + theme_classic() + geom_point(size = 0.2)
g <- g + scale_color_manual(values = col) + facet_wrap( ~ sample.name, nrow = 1)

pdf(UMAP, width = 14, height = 3)
g
dev.off()


UMAP <- paste0("umap_lineage_treated_",EXP_NAME,".pdf")

# load data

df <- read.table(INFO, sep = "\t", header = TRUE)
df$lineage[is.na(df$lineage)] <- rep(0, sum(is.na(df$lineage)))

# select palette

col <- c(brewer.pal(4,"Greys")[2],hue_pal()(2))

# figure with umap colored by lineage

df$sample.name <- factor(df$sample.name, levels = unique(df$sample.name))
df$lineage <- factor(df$lineage, levels = c(0,1,2))

g <- ggplot(df, aes(x = -UMAP_1, y = UMAP_2, color = lineage)) + theme_classic() + geom_point(size = 0.2)
g <- g + scale_color_manual(values = col) + facet_wrap( ~ sample.name, nrow = 1)

pdf(UMAP, width = 14, height = 3)
g
dev.off()


### EXP1B time-course treated

# directory and file names

EXP_NAME <- "EXP1B"

INFO <- paste0("../tables/DATA_FRAME_",EXP_NAME,"_treated_filt.tsv")

UMAP <- paste0("umap_cl_treated_",EXP_NAME,".pdf")

# load data

df <- read.table(INFO, sep = "\t", header = TRUE)
cl_levels <- sort(unique(df$cl_renamed))

# figure with umap colored by cluster

df$cl_renamed <- factor(df$cl_renamed, levels = cl_levels)

g <- ggplot(df, aes(x = -UMAP_1, y = UMAP_2, color = cl_renamed)) + theme_classic() + geom_point(size = 0.2) + ggtitle("exp2")

pdf(UMAP, width = 5, height = 4.5)
g
dev.off()


### EXP1D time-course treated

# directory and file names

EXP_NAME <- "EXP1D"

INFO <- paste0("../tables/DATA_FRAME_",EXP_NAME,"_treated_filt.tsv")

UMAP <- paste0("umap_cl_treated_",EXP_NAME,".pdf")


# load data

df <- read.table(INFO, sep = "\t", header = TRUE)
cl_levels <- sort(unique(df$cl_renamed))

# figure with umap colored by cluster

df$cl_renamed <- factor(df$cl_renamed, levels = cl_levels)

g <- ggplot(df, aes(x = -UMAP_1, y = UMAP_2, color = cl_renamed)) + theme_classic() + geom_point(size = 0.2) + ggtitle("exp1")

pdf(UMAP, width = 5, height = 4.5)
g
dev.off()



q()

