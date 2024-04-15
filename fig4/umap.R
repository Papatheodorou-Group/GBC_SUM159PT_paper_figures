
# load libraries

library("ggplot2")
library("circlize")
library("scales")
library("RColorBrewer")


### EXP1C

# directory and file names

EXP_NAME <- "EXP1C"

INFO <- paste0("../tables/DATA_FRAME_",EXP_NAME,".tsv")

UMAP_RNA <- paste0("umap_cl_",EXP_NAME,"_rna.pdf")
UMAP_ATAC <- paste0("umap_cl_",EXP_NAME,"_atac.pdf")

# load data

df <- read.table(INFO, sep = "\t", header = TRUE)
cl_levels <- sort(unique(df$cl_renamed))

# select palette

col_cl <- brewer.pal(n = length(cl_levels), name = "Dark2")

# figure with umap colored by cluster

df$cl_renamed <- factor(df$cl_renamed, levels = cl_levels)

g <- ggplot(df, aes(x = UMAP_1_rna, y = UMAP_2_rna, color = cl_renamed)) + theme_classic() + geom_point(size = 0.5)
g <- g + scale_color_manual(values = col_cl) + ggtitle("T0")

pdf(UMAP_RNA, width = 5.2, height = 4.5)
g
dev.off()

g <- ggplot(df, aes(x = UMAP_1_atac, y = UMAP_2_atac, color = cl_renamed)) + theme_classic() + geom_point(size = 0.5)
g <- g + scale_color_manual(values = col_cl) + ggtitle("T0")

pdf(UMAP_ATAC, width = 5.2, height = 4.5)
g
dev.off()


### EXP1E

# directory and file names

EXP_NAME <- "EXP1E"

INFO <- paste0("../tables/DATA_FRAME_",EXP_NAME,".tsv")

UMAP_RNA <- paste0("umap_cl_",EXP_NAME,"_rna.pdf")
UMAP_ATAC <- paste0("umap_cl_",EXP_NAME,"_atac.pdf")

# load data

df <- read.table(INFO, sep = "\t", header = TRUE)
cl_levels <- sort(unique(df$cl_renamed))

# select palette

col_cl <- brewer.pal(n = length(cl_levels), name = "Dark2")

# figure with umap colored by cluster

df$cl_renamed <- factor(df$cl_renamed, levels = cl_levels)

g <- ggplot(df, aes(x = UMAP_1_rna, y = UMAP_2_rna, color = cl_renamed)) + theme_classic() + geom_point(size = 0.5)
g <- g + scale_color_manual(values = col_cl) + ggtitle("T0")

pdf(UMAP_RNA, width = 5.2, height = 4.5)
g
dev.off()

g <- ggplot(df, aes(x = UMAP_1_atac, y = UMAP_2_atac, color = cl_renamed)) + theme_classic() + geom_point(size = 0.5)
g <- g + scale_color_manual(values = col_cl) + ggtitle("T0")

pdf(UMAP_ATAC, width = 5.2, height = 4.5)
g
dev.off()



q()


