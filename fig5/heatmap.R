
# load libraries

library("ComplexHeatmap")
library("circlize")
library("scales")
library("RColorBrewer")
library("ggrepel")
    

### EXP1C

SC_EXP_NAME <- "1C_ARC"
TABLE <- "../tables/DATA_FRAME_EXP1C.tsv"

OUTDIR <- paste0("../wes_matrices/out_", SC_EXP_NAME, "-WES_intersect")

NORM_COUNT_TABLE <- file.path(OUTDIR, "matrix.tsv")
COL_ANNOTATION <- file.path(OUTDIR, "cl_annotation.txt")
ROW_ANNOTATION <- file.path(OUTDIR, "logfc_annotation.tsv")
    
NORM_COUNT_HEATMAP <- paste0("heatmap_WES_", SC_EXP_NAME, ".pdf")

# load data
    
M <- as.matrix(read.table(NORM_COUNT_TABLE, sep = "\t", header = TRUE))
row_annotation <- as.matrix(read.table(ROW_ANNOTATION))
colnames(row_annotation) <- paste0("WES",1:ncol(row_annotation))
col_annotation <- read.table(COL_ANNOTATION)[,1]

data <- read.table(TABLE, sep = "\t", header = TRUE)
idx <- match(gsub("\\.","-",colnames(M)),data$cell_ID)
data <- data[idx,]

ordering <- order(data$cl_renamed)
col_annotation <- sort(data$cl_renamed)
idx <- order(apply(row_annotation,1,mean),decreasing = TRUE)
M <- M[idx,ordering]
row_annotation <- row_annotation[idx,]

M <- t(scale(t(M)))
 
# generate norm count heatmap
  
cl <- unique(col_annotation)
ncl <- length(cl)
   
col_fun_cl <- brewer.pal(n = length(cl), name = "Dark2")
names(col_fun_cl) <- cl
col_fun <- colorRamp2(c(-2,0,2), c("white","gray80","black"))
    
# NB: do not put DTC annotation because we don't see anything
    
row_split <- rep("depleted p.t.", nrow(row_annotation))
idx <- which(apply(row_annotation,1,mean) > 0)
row_split[idx] <- rep("amplified p.t.", length(idx))
   
labelfont <- rep("plain",nrow(M))
labelfont[1] <- "bold"

col_fun_ann <- colorRamp2(c(min(row_annotation),0,max(row_annotation)), c("blue","gray90","red"))

ha <- HeatmapAnnotation(log2fc = row_annotation, which = "row", col = list(log2fc = col_fun_ann), 
                       show_legend = TRUE, show_annotation_name = TRUE) 

hb <- HeatmapAnnotation(cluster = col_annotation, which = "column", col = list(cluster = col_fun_cl),
                        show_legend = TRUE, show_annotation_name = FALSE)
                        
h <- Heatmap(M, name = "Z-score", col = col_fun, left_annotation = ha, top_annotation = hb,
             column_split = col_annotation, row_split = row_split, row_names_gp = gpar(fontface = labelfont),
             cluster_rows = FALSE, show_column_names = FALSE, column_title = "DNA accessibility at day 0",
             row_title = "reproducible CNV regions")
    
pdf(NORM_COUNT_HEATMAP, width = 7.5, height = 5)
draw(h)
dev.off()


TOP <- 10

OUTDIR <- paste0("../wes_matrices/out_double_treatment_", SC_EXP_NAME, "-WES")

NORM_COUNT_TABLE <- file.path(OUTDIR, "matrix.tsv")
COL_ANNOTATION <- file.path(OUTDIR, "cl_annotation.txt")
ROW_ANNOTATION <- file.path(OUTDIR, "logfc_annotation.tsv")
    
NORM_COUNT_HEATMAP <- paste0("heatmap_WES_double_treatment_", SC_EXP_NAME, ".pdf")

# load data
    
M <- as.matrix(read.table(NORM_COUNT_TABLE, sep = "\t", header = TRUE))
row_annotation <- as.matrix(read.table(ROW_ANNOTATION))
col_annotation <- read.table(COL_ANNOTATION)[,1]

data <- read.table(TABLE, sep = "\t", header = TRUE)
idx <- match(gsub("\\.","-",colnames(M)),data$cell_ID)
data <- data[idx,]

ordering <- order(data$cl_renamed)
col_annotation <- sort(data$cl_renamed)
idx <- order(apply(row_annotation,1,mean),decreasing = TRUE)
M <- M[idx,ordering]
row_annotation <- row_annotation[idx]

M <- t(scale(t(M)))
 
# generate norm count heatmap
  
cl <- unique(col_annotation)
ncl <- length(cl)
   
col_fun_cl <- brewer.pal(n = length(cl), name = "Dark2")
names(col_fun_cl) <- cl
col_fun <- colorRamp2(c(-2,0,2), c("white","gray80","black"))
    
# NB: do not put DTC annotation because we don't see anything
    
row_split <- rep("depleted p.t.", length(row_annotation))
idx <- which(row_annotation > 0)
row_split[idx] <- rep("amplified p.t.", length(idx))
   
labelfont <- rep("plain",nrow(M))
labelfont[1] <- "bold"

col_fun_ann <- colorRamp2(c(-0.5,0,0.5), c("blue","gray90","red"))

val <- row_annotation[order(abs(row_annotation),decreasing=TRUE)[TOP]] # min abs value for the row names to be shown

ha <- HeatmapAnnotation(log2fc = row_annotation, which = "row", col = list(log2fc = col_fun_ann), 
                        show_legend = TRUE, show_annotation_name = TRUE) 

hb <- HeatmapAnnotation(cluster = col_annotation, which = "column", col = list(cluster = col_fun_cl),
                        show_legend = TRUE, show_annotation_name = FALSE)
                        
idx_val <- which(abs(row_annotation) > val)
gene_labels_row <- rowAnnotation(highlight = anno_mark(at = idx_val, labels = rownames(M)[idx_val], 
                                 side = 'right', labels_gp = gpar(fontsize = 9)))
                        
h <- Heatmap(M, name = "Z-score", col = col_fun, top_annotation = hb, left_annotation = ha, 
             column_split = col_annotation, row_split = row_split, row_names_gp = gpar(fontface = labelfont),
             cluster_rows = FALSE, show_column_names = FALSE, column_title = "DNA accessibility at day 0",
             show_row_names = FALSE, cluster_column_slices = FALSE, right_annotation = gene_labels_row,
             row_title = "CNV regions")
    
pdf(NORM_COUNT_HEATMAP, width = 6, height = 5)
draw(h)
dev.off()
    

### EXP1E

SC_EXP_NAME <- "1E_ARC"
TABLE <- "../tables/DATA_FRAME_EXP1E.tsv"

OUTDIR <- paste0("../wes_matrices/out_", SC_EXP_NAME, "-WES_intersect")

NORM_COUNT_TABLE <- file.path(OUTDIR, "matrix.tsv")
COL_ANNOTATION <- file.path(OUTDIR, "cl_annotation.txt")
ROW_ANNOTATION <- file.path(OUTDIR, "logfc_annotation.tsv")
    
NORM_COUNT_HEATMAP <- paste0("heatmap_WES_", SC_EXP_NAME, ".pdf")
    
# load data
    
M <- as.matrix(read.table(NORM_COUNT_TABLE, sep = "\t", header = TRUE))
row_annotation <- as.matrix(read.table(ROW_ANNOTATION))
colnames(row_annotation) <- paste0("WES",1:ncol(row_annotation))
col_annotation <- read.table(COL_ANNOTATION)[,1]

data <- read.table(TABLE, sep = "\t", header = TRUE)
idx <- match(gsub("\\.","-",colnames(M)),data$cell_ID)
data <- data[idx,]

ordering <- order(data$cl_renamed)
col_annotation <- sort(data$cl_renamed)
idx <- order(apply(row_annotation,1,mean),decreasing = TRUE)
M <- M[idx,ordering]
row_annotation <- row_annotation[idx,]

M <- t(scale(t(M)))
 
# generate norm count heatmap
  
cl <- unique(col_annotation)
ncl <- length(cl)
   
col_fun_cl <- brewer.pal(n = length(cl), name = "Dark2")
names(col_fun_cl) <- cl
col_fun <- colorRamp2(c(-2,0,2), c("white","gray80","black"))
    
# NB: do not put DTC annotation because we don't see anything
    
row_split <- rep("depleted p.t.", nrow(row_annotation))
idx <- which(apply(row_annotation,1,mean) > 0)
row_split[idx] <- rep("amplified p.t.", length(idx))
   
labelfont <- rep("plain",nrow(M))
labelfont[1] <- "bold"

col_fun_ann <- colorRamp2(c(min(row_annotation),0,max(row_annotation)), c("blue","gray90","red"))

ha <- HeatmapAnnotation(log2fc = row_annotation, which = "row", col = list(log2fc = col_fun_ann), 
                       show_legend = TRUE, show_annotation_name = TRUE) 

hb <- HeatmapAnnotation(cluster = col_annotation, which = "column", col = list(cluster = col_fun_cl),
                        show_legend = TRUE, show_annotation_name = FALSE)
                        
h <- Heatmap(M, name = "Z-score", col = col_fun, left_annotation = ha, top_annotation = hb,
             column_split = col_annotation, row_split = row_split, row_names_gp = gpar(fontface = labelfont),
             cluster_rows = FALSE, show_column_names = FALSE, column_title = "DNA accessibility at day 0",
             row_title = "reproducible CNV regions")
    
pdf(NORM_COUNT_HEATMAP, width = 7.5, height = 5)
draw(h)
dev.off()
    


TOP <- 10

OUTDIR <- paste0("../wes_matrices/out_double_treatment_", SC_EXP_NAME, "-WES")

NORM_COUNT_TABLE <- file.path(OUTDIR, "matrix.tsv")
COL_ANNOTATION <- file.path(OUTDIR, "cl_annotation.txt")
ROW_ANNOTATION <- file.path(OUTDIR, "logfc_annotation.tsv")
    
NORM_COUNT_HEATMAP <- paste0("heatmap_WES_double_treatment_", SC_EXP_NAME, ".pdf")

# load data
    
M <- as.matrix(read.table(NORM_COUNT_TABLE, sep = "\t", header = TRUE))
row_annotation <- as.matrix(read.table(ROW_ANNOTATION))
col_annotation <- read.table(COL_ANNOTATION)[,1]

data <- read.table(TABLE, sep = "\t", header = TRUE)
idx <- match(gsub("\\.","-",colnames(M)),data$cell_ID)
data <- data[idx,]

ordering <- order(data$cl_renamed)
col_annotation <- sort(data$cl_renamed)
idx <- order(apply(row_annotation,1,mean),decreasing = TRUE)
M <- M[idx,ordering]
row_annotation <- row_annotation[idx]

M <- t(scale(t(M)))
 
# generate norm count heatmap
  
cl <- unique(col_annotation)
ncl <- length(cl)
   
col_fun_cl <- brewer.pal(n = length(cl), name = "Dark2")
names(col_fun_cl) <- cl
col_fun <- colorRamp2(c(-2,0,2), c("white","gray80","black"))
    
# NB: do not put DTC annotation because we don't see anything
    
row_split <- rep("depleted p.t.", length(row_annotation))
idx <- which(row_annotation > 0)
row_split[idx] <- rep("amplified p.t.", length(idx))
   
labelfont <- rep("plain",nrow(M))
labelfont[1] <- "bold"

col_fun_ann <- colorRamp2(c(-0.5,0,0.5), c("blue","gray90","red"))

val <- row_annotation[order(abs(row_annotation),decreasing=TRUE)[TOP]] # min abs value for the row names to be shown

ha <- HeatmapAnnotation(log2fc = row_annotation, which = "row", col = list(log2fc = col_fun_ann), 
                        show_legend = TRUE, show_annotation_name = TRUE) 

hb <- HeatmapAnnotation(cluster = col_annotation, which = "column", col = list(cluster = col_fun_cl),
                        show_legend = TRUE, show_annotation_name = FALSE)
                        
idx_val <- which(abs(row_annotation) > val)
gene_labels_row <- rowAnnotation(highlight = anno_mark(at = idx_val, labels = rownames(M)[idx_val], 
                                 side = 'right', labels_gp = gpar(fontsize = 9)))
                        
h <- Heatmap(M, name = "Z-score", col = col_fun, top_annotation = hb, left_annotation = ha, 
             column_split = col_annotation, row_split = row_split, row_names_gp = gpar(fontface = labelfont),
             cluster_rows = FALSE, show_column_names = FALSE, column_title = "DNA accessibility at day 0",
             show_row_names = FALSE, cluster_column_slices = FALSE, right_annotation = gene_labels_row,
             row_title = "CNV regions")
    
pdf(NORM_COUNT_HEATMAP, width = 6, height = 5)
draw(h)
dev.off()
    


q()

