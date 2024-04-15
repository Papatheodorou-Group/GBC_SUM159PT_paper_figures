
# load libraries

library("ggplot2")
library("circlize")
library("scales")
library("RColorBrewer")


### EXP1B P0 untreated

# directory and file names

EXP_NAME <- "EXP1B"

INFO <- paste0("../tables/DATA_FRAME_",EXP_NAME,"_P0_filt.tsv")

UMAP <- paste0("umap_cl_P0_",EXP_NAME,".pdf")

# load data

df <- read.table(INFO, sep = "\t", header = TRUE)
cl_levels <- sort(unique(df$cl_renamed))

# select palette

col_cl <- brewer.pal(n = length(cl_levels), name = "Set1")

# figure with umap colored by cluster

df$cl_renamed <- factor(df$cl_renamed, levels = cl_levels)

g <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = cl_renamed)) + theme_classic() + geom_point(size = 0.5)
g <- g + scale_color_manual(values = col_cl) + ggtitle("T0")

pdf(UMAP, width = 5, height = 4.5)
g
dev.off()


### EXP1D P0 untreated

# directory and file names

EXP_NAME <- "EXP1D"

INFO <- paste0("../tables/DATA_FRAME_",EXP_NAME,"_P0_filt.tsv")

UMAP <- paste0("umap_cl_P0_",EXP_NAME,".pdf")

# load data

df <- read.table(INFO, sep = "\t", header = TRUE)
cl_levels <- sort(unique(df$cl_renamed))

# select palette

col_cl <- brewer.pal(n = length(cl_levels), name = "Set2")

# figure with umap colored by cluster

df$cl_renamed <- factor(df$cl_renamed, levels = cl_levels)

g <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = cl_renamed)) + theme_classic() + geom_point(size = 0.5)
g <- g + scale_color_manual(values = col_cl) + ggtitle("T1")

pdf(UMAP, width = 5, height = 4.5)
g
dev.off()



q()


