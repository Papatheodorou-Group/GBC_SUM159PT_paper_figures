
# load libraries

library("ggplot2")
library("circlize")
library("scales")
library("RColorBrewer")

TIC_BIN <- c(0,1)


### DTC

# select palettes

col_unsel <- brewer.pal(4,"Greys")[2]
col_tic <- brewer.pal(4,"Oranges")[4]
cols_tic <- c(col_unsel, col_tic)
names(cols_tic) <- TIC_BIN
col_sel <- "black"
cols_sel <- c(col_unsel, col_sel)
names(cols_sel) <- TIC_BIN

# directory and file names

EXP_NAME <- "EXP1B"

INFO <- paste0("../tables/DATA_FRAME_",EXP_NAME,"_P0_filt.tsv")

UMAP <- paste0("umap_DTC_P0_",EXP_NAME,".pdf")

# load data

df <- read.table(INFO, sep = "\t", header = TRUE)

# select color for highlighting

pt_size_tic <- rep(0.5,nrow(df))
pt_size_tic[df$tic == 1] <- rep(1,sum(df$tic == 1))

# figure with umap colored by TIC

df$DRC_sc <- factor(df$DRC_sc, levels = TIC_BIN)
df <- df[order(df$DRC_sc),]

g <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = DRC_sc)) + theme_classic() + geom_point()
g <- g + scale_color_manual(values = cols_tic) + ggtitle("DTC")

pdf(UMAP, width = 5, height = 4.5)
g
dev.off()


# directory and file names

EXP_NAME <- "EXP1D"

INFO <- paste0("../tables/DATA_FRAME_",EXP_NAME,"_P0_filt.tsv")

UMAP <- paste0("umap_DTC_P0_",EXP_NAME,".pdf")

# load data

df <- read.table(INFO, sep = "\t", header = TRUE)

# select color for highlighting

pt_size_tic <- rep(0.5,nrow(df))
pt_size_tic[df$tic == 1] <- rep(1,sum(df$tic == 1))

# figure with umap colored by TIC

df$DRC_sc <- factor(df$DRC_sc, levels = TIC_BIN)
df <- df[order(df$DRC_sc),]

g <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = DRC_sc)) + theme_classic() + geom_point()
g <- g + scale_color_manual(values = cols_tic) + ggtitle("DTC")

pdf(UMAP, width = 5, height = 4.5)
g
dev.off()


### TIC-DTC comparison

# directory and file names

EXP_NAME <- "EXP1B"

INFO <- paste0("../tables/DATA_FRAME_",EXP_NAME,"_P0_filt.tsv")

UMAP <- paste0("umap_TIC-DTC-both_P0_",EXP_NAME,".pdf")

# load data

df <- read.table(INFO, sep = "\t", header = TRUE)

# select color for highlighting

pt_size_tic <- rep(0.5,nrow(df))
pt_size_tic[df$tic == 1] <- rep(1,sum(df$tic == 1))

# figure with umap colored by selected cells

v <- as.numeric(df$TIC_fisher + df$DRC_sc > 0)
w <- rep("none", nrow(df))
w[df$TIC_fisher == 1] <- rep("TIC-only", sum(df$TIC_fisher == 1))
w[df$DRC_sc == 1 & df$TIC_fisher == 1] <- rep("both", sum(df$DRC_sc == 1 & df$TIC_fisher == 1))
w[df$DRC_sc == 1 & df$TIC_fisher == 0] <- rep("DTC-only", sum(df$DRC_sc == 1 & df$TIC_fisher == 0))
df <- cbind(df, sel = v, class = w)

df <- df[order(df$sel),]
df <- df[df$class != "none",]

df$sel <- factor(df$sel, levels = TIC_BIN)
df$class <- factor(df$class, levels = c("TIC-only","DTC-only","both"))

g <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = sel)) + theme_classic() + geom_point()
g <- g + scale_color_manual(values = cols_sel) + facet_wrap(~ class, ncol = 3)

pdf(UMAP, width = 8, height = 3)
g
dev.off()


# directory and file names

EXP_NAME <- "EXP1D"

INFO <- paste0("../tables/DATA_FRAME_",EXP_NAME,"_P0_filt.tsv")

UMAP <- paste0("umap_TIC-DTC-both_P0_",EXP_NAME,".pdf")

# load data

df <- read.table(INFO, sep = "\t", header = TRUE)

pt_size_tic <- rep(0.5,nrow(df))
pt_size_tic[df$tic == 1] <- rep(1,sum(df$tic == 1))

# figure with umap colored by selected cells

v <- as.numeric(df$TIC_fisher + df$DRC_sc > 0)
w <- rep("none", nrow(df))
w[df$TIC_fisher == 1] <- rep("TIC-only", sum(df$TIC_fisher == 1))
w[df$DRC_sc == 1 & df$TIC_fisher == 1] <- rep("both", sum(df$DRC_sc == 1 & df$TIC_fisher == 1))
w[df$DRC_sc == 1 & df$TIC_fisher == 0] <- rep("DTC-only", sum(df$DRC_sc == 1 & df$TIC_fisher == 0))
df <- cbind(df, sel = v, class = w)

df <- df[order(df$sel),]
df <- df[df$class != "none",]

df$sel <- factor(df$sel, levels = TIC_BIN)
df$class <- factor(df$class, levels = c("TIC-only","DTC-only","both"))

g <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = sel)) + theme_classic() + geom_point()
g <- g + scale_color_manual(values = cols_sel) + facet_wrap(~ class, ncol = 3)

pdf(UMAP, width = 8, height = 3)
g
dev.off()


# directory and file names

EXP_NAME <- "EXP1B"

INFO <- paste0("../tables/DATA_FRAME_",EXP_NAME,"_P0_filt.tsv")

UMAP <- paste0("umap_TIC-DTCbulk-both_P0_",EXP_NAME,".pdf")

# load data

df <- read.table(INFO, sep = "\t", header = TRUE)

pt_size_tic <- rep(0.5,nrow(df))
pt_size_tic[df$tic == 1] <- rep(1,sum(df$tic == 1))

# figure with umap colored by selected cells

v <- as.numeric(df$TIC_fisher + df$DRC_bulk_fisher > 0)
w <- rep("none", nrow(df))
w[df$TIC_fisher == 1] <- rep("TIC-only", sum(df$TIC_fisher == 1))
w[df$DRC_bulk_fisher == 1 & df$TIC_fisher == 1] <- rep("both", sum(df$DRC_bulk_fisher == 1 & df$TIC_fisher == 1))
w[df$DRC_bulk_fisher == 1 & df$TIC_fisher == 0] <- rep("DTC-only", sum(df$DRC_bulk_fisher == 1 & df$TIC_fisher == 0))
df <- cbind(df, sel = v, class = w)

df <- df[order(df$sel),]
df <- df[df$class != "none",]

df$sel <- factor(df$sel, levels = TIC_BIN)
df$class <- factor(df$class, levels = c("TIC-only","DTC-only","both"))

g <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = sel)) + theme_classic() + geom_point()
g <- g + scale_color_manual(values = cols_sel) + facet_wrap(~ class, ncol = 3)

pdf(UMAP, width = 8, height = 3)
g
dev.off()


# directory and file names

EXP_NAME <- "EXP1D"

INFO <- paste0("../tables/DATA_FRAME_",EXP_NAME,"_P0_filt.tsv")

UMAP <- paste0("umap_TIC-DTCbulk-both_P0_",EXP_NAME,".pdf")

# load data

df <- read.table(INFO, sep = "\t", header = TRUE)

pt_size_tic <- rep(0.5,nrow(df))
pt_size_tic[df$tic == 1] <- rep(1,sum(df$tic == 1))

# figure with umap colored by selected cells

v <- as.numeric(df$TIC_fisher + df$DRC_bulk_fisher > 0)
w <- rep("none", nrow(df))
w[df$TIC_fisher == 1] <- rep("TIC-only", sum(df$TIC_fisher == 1))
w[df$DRC_bulk_fisher == 1 & df$TIC_fisher == 1] <- rep("both", sum(df$DRC_bulk_fisher == 1 & df$TIC_fisher == 1))
w[df$DRC_bulk_fisher == 1 & df$TIC_fisher == 0] <- rep("DTC-only", sum(df$DRC_bulk_fisher == 1 & df$TIC_fisher == 0))
df <- cbind(df, sel = v, class = w)

df <- df[order(df$sel),]
df <- df[df$class != "none",]

df$sel <- factor(df$sel, levels = TIC_BIN)
df$class <- factor(df$class, levels = c("TIC-only","DTC-only","both"))

g <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = sel)) + theme_classic() + geom_point()
g <- g + scale_color_manual(values = cols_sel) + facet_wrap(~ class, ncol = 3)

pdf(UMAP, width = 8, height = 3)
g
dev.off()


# directory and file names

EXP_NAME <- "EXP1B"

INFO <- paste0("../tables/DATA_FRAME_",EXP_NAME,"_P0_filt.tsv")

UMAP <- paste0("umap_treatedTIC-DTC-both_P0_",EXP_NAME,".pdf")

# load data

df <- read.table(INFO, sep = "\t", header = TRUE)

pt_size_tic <- rep(0.5,nrow(df))
pt_size_tic[df$tic == 1] <- rep(1,sum(df$tic == 1))

# figure with umap colored by selected cells

v <- as.numeric(df$treated_TIC_fisher + df$DRC_sc > 0)
w <- rep("none", nrow(df))
w[df$treated_TIC_fisher == 1] <- rep("TIC-only", sum(df$treated_TIC_fisher == 1))
w[df$DRC_sc == 1 & df$treated_TIC_fisher == 1] <- rep("both", sum(df$DRC_sc == 1 & df$treated_TIC_fisher == 1))
w[df$DRC_sc == 1 & df$treated_TIC_fisher == 0] <- rep("DTC-only", sum(df$DRC_sc == 1 & df$treated_TIC_fisher == 0))
df <- cbind(df, sel = v, class = w)

df <- df[order(df$sel),]
df <- df[df$class != "none",]

df$sel <- factor(df$sel, levels = TIC_BIN)
df$class <- factor(df$class, levels = c("TIC-only","DTC-only","both"))

g <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = sel)) + theme_classic() + geom_point()
g <- g + scale_color_manual(values = cols_sel) + facet_wrap(~ class, ncol = 3)

pdf(UMAP, width = 8, height = 3)
g
dev.off()


# directory and file names

EXP_NAME <- "EXP1B"

INFO <- paste0("../tables/DATA_FRAME_",EXP_NAME,"_P0_filt.tsv")

UMAP <- paste0("umap_treatedTIC-DTCbulk-both_P0_",EXP_NAME,".pdf")

# load data

df <- read.table(INFO, sep = "\t", header = TRUE)

pt_size_tic <- rep(0.5,nrow(df))
pt_size_tic[df$tic == 1] <- rep(1,sum(df$tic == 1))

# figure with umap colored by selected cells

v <- as.numeric(df$treated_TIC_fisher + df$DRC_bulk_fisher > 0)
w <- rep("none", nrow(df))
w[df$treated_TIC_fisher == 1] <- rep("TIC-only", sum(df$treated_TIC_fisher == 1))
w[df$DRC_bulk_fisher == 1 & df$treated_TIC_fisher == 1] <- rep("both", sum(df$DRC_bulk_fisher == 1 & df$treated_TIC_fisher == 1))
w[df$DRC_bulk_fisher == 1 & df$treated_TIC_fisher == 0] <- rep("DTC-only", sum(df$DRC_bulk_fisher == 1 & df$treated_TIC_fisher == 0))
df <- cbind(df, sel = v, class = w)

df <- df[order(df$sel),]
df <- df[df$class != "none",]

df$sel <- factor(df$sel, levels = TIC_BIN)
df$class <- factor(df$class, levels = c("TIC-only","DTC-only","both"))

g <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = sel)) + theme_classic() + geom_point()
g <- g + scale_color_manual(values = cols_sel) + facet_wrap(~ class, ncol = 3)

pdf(UMAP, width = 8, height = 3)
g
dev.off()


# directory and file names

EXP_NAME <- "EXP1D"

INFO <- paste0("../tables/DATA_FRAME_",EXP_NAME,"_P0_filt.tsv")

UMAP <- paste0("umap_treatedTIC-DTCbulk-both_P0_",EXP_NAME,".pdf")

# load data

df <- read.table(INFO, sep = "\t", header = TRUE)

pt_size_tic <- rep(0.5,nrow(df))
pt_size_tic[df$tic == 1] <- rep(1,sum(df$tic == 1))

# figure with umap colored by selected cells

v <- as.numeric(df$treated_TIC_fisher + df$DRC_bulk_fisher > 0)
w <- rep("none", nrow(df))
w[df$treated_TIC_fisher == 1] <- rep("TIC-only", sum(df$treated_TIC_fisher == 1))
w[df$DRC_bulk_fisher == 1 & df$treated_TIC_fisher == 1] <- rep("both", sum(df$DRC_bulk_fisher == 1 & df$treated_TIC_fisher == 1))
w[df$DRC_bulk_fisher == 1 & df$treated_TIC_fisher == 0] <- rep("DTC-only", sum(df$DRC_bulk_fisher == 1 & df$treated_TIC_fisher == 0))
df <- cbind(df, sel = v, class = w)

df <- df[order(df$sel),]
df <- df[df$class != "none",]

df$sel <- factor(df$sel, levels = TIC_BIN)
df$class <- factor(df$class, levels = c("TIC-only","DTC-only","both"))

g <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = sel)) + theme_classic() + geom_point()
g <- g + scale_color_manual(values = cols_sel) + facet_wrap(~ class, ncol = 3)

pdf(UMAP, width = 8, height = 3)
g
dev.off()


q()


