
# load libraries

library("ggplot2")
library("circlize")
library("scales")
library("RColorBrewer")

TIC_BIN <- c(0,1)


# select palettes

col_unsel <- brewer.pal(4,"Greys")[2]
col_tic <- brewer.pal(4,"Blues")[4]
col_sel <- "black"
cols_tic <- c(col_unsel, col_tic)
names(cols_tic) <- TIC_BIN


### EXP1B

# directory and file names

EXP_NAME <- "EXP1B"

INFO <- paste0("../tables/DATA_FRAME_",EXP_NAME,"_P0_filt.tsv")

UMAP <- paste0("umap_TIC_P0_",EXP_NAME,".pdf")

# load data

df <- read.table(INFO, sep = "\t", header = TRUE)

# figure with umap colored by TIC

df$TIC_fisher <- factor(df$TIC_fisher, levels = TIC_BIN)
df <- df[order(df$TIC_fisher),]

g <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = TIC_fisher)) + theme_classic() + geom_point()
g <- g + scale_color_manual(values = cols_tic) + ggtitle("TIC")

pdf(UMAP, width = 5, height = 4.5)
g
dev.off()


### EXP1D

# directory and file names

EXP_NAME <- "EXP1D"

INFO <- paste0("../tables/DATA_FRAME_",EXP_NAME,"_P0_filt.tsv")

UMAP <- paste0("umap_TIC_P0_",EXP_NAME,".pdf")

# load data

df <- read.table(INFO, sep = "\t", header = TRUE)

# figure with umap colored by TIC

df$TIC_fisher <- factor(df$TIC_fisher, levels = TIC_BIN)
df <- df[order(df$TIC_fisher),]

g <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = TIC_fisher)) + theme_classic() + geom_point()
g <- g + scale_color_manual(values = cols_tic) + ggtitle("TIC")

pdf(UMAP, width = 5, height = 4.5)
g
dev.off()



### EXP1B 

# directory and file names

EXP_NAME <- "EXP1B"

INFO <- paste0("../tables/DATA_FRAME_",EXP_NAME,"_P0_filt.tsv")

UMAP <- paste0("umap_TIC-treatedTIC-both_P0_",EXP_NAME,".pdf")

# load data

df <- read.table(INFO, sep = "\t", header = TRUE)

# select color for highlighting

cols_sel <- c(col_unsel, col_sel)
names(cols_sel) <- TIC_BIN
pt_size_sel <- rep(0.5,nrow(df))
pt_size_sel[df$sel == 1] <- rep(1,sum(df$sel == 1))

# figure with umap colored by selected cells

v <- as.numeric(df$TIC_fisher + df$treated_TIC_fisher > 0)
w <- rep("none", nrow(df))
w[df$TIC_fisher == 1] <- rep("TIC-only", sum(df$TIC_fisher == 1))
w[df$TIC_fisher == 1 & df$treated_TIC_fisher == 1] <- rep("both", sum(df$TIC_fisher == 1 & df$treated_TIC_fisher == 1))
w[df$TIC_fisher == 0 & df$treated_TIC_fisher == 1] <- rep("treatedTIC-only", sum(df$TIC_fisher == 0 & df$treated_TIC_fisher == 1))
df <- cbind(df, sel = v, class = w)

df <- df[order(df$sel),]
df <- df[df$class != "none",]

df$sel <- factor(df$sel, levels = TIC_BIN)
df$class <- factor(df$class, levels = c("TIC-only","treatedTIC-only","both"))

g <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = sel)) + theme_classic() + geom_point()
g <- g + scale_color_manual(values = cols_sel) + facet_wrap(~ class, ncol = 3)

pdf(UMAP, width = 8, height = 3)
g
dev.off()


### EXP1D

# directory and file names

EXP_NAME <- "EXP1D"

INFO <- paste0("../tables/DATA_FRAME_",EXP_NAME,"_P0_filt.tsv")

UMAP <- paste0("umap_TIC-treatedTIC-both_P0_",EXP_NAME,".pdf")

# load data

df <- read.table(INFO, sep = "\t", header = TRUE)

# select color for highlighting

cols_sel <- c(col_unsel, col_sel)
names(cols_sel) <- TIC_BIN
pt_size_sel <- rep(0.5,nrow(df))
pt_size_sel[df$sel == 1] <- rep(1,sum(df$sel == 1))

# figure with umap colored by selected cells

v <- as.numeric(df$TIC_fisher + df$treated_TIC_fisher > 0)
w <- rep("none", nrow(df))
w[df$TIC_fisher == 1] <- rep("TIC-only", sum(df$TIC_fisher == 1))
w[df$TIC_fisher == 1 & df$treated_TIC_fisher == 1] <- rep("both", sum(df$TIC_fisher == 1 & df$treated_TIC_fisher == 1))
w[df$TIC_fisher == 0 & df$treated_TIC_fisher == 1] <- rep("treatedTIC-only", sum(df$TIC_fisher == 0 & df$treated_TIC_fisher == 1))
df <- cbind(df, sel = v, class = w)

df <- df[order(df$sel),]
df <- df[df$class != "none",]

df$sel <- factor(df$sel, levels = TIC_BIN)
df$class <- factor(df$class, levels = c("TIC-only","treatedTIC-only","both"))

g <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = sel)) + theme_classic() + geom_point()
g <- g + scale_color_manual(values = cols_sel) + facet_wrap(~ class, ncol = 3)

pdf(UMAP, width = 8, height = 3)
g
dev.off()



q()


