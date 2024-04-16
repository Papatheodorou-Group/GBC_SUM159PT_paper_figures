

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

EXP_NAME <- "EXP1B"

SAMPLES <- c("T3","T9","T13","T17")
SAMPLES_RENAMED <- c("d3","d9","d13","d17")

OBJECT <- "../seurat_objects/timecourse_EXP1B.Rds"

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

EXP_NAME <- "EXP1D"

SAMPLES <- c("T5_2_Multiseq","T7_2_Multiseq","T9_2_Multiseq","T11_2_Multiseq","T13_2_Multiseq","T15_2_Multiseq")
SAMPLES_RENAMED <- c("d5","d7","d9","d11","d13","d15")

OBJECT <- "../seurat_objects/timecourse_EXP1D.Rds"

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



q()

