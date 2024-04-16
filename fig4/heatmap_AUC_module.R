
# load libraries

library("pROC")
library("RColorBrewer")
library("circlize")
library("ComplexHeatmap")

PHENOTYPE <- c("DRC_bulk_fisher","TIC_fisher","treated_TIC_fisher")
PHENOTYPE_RENAMED <- c("DTC_invitro","TIC","DTC_invivo")
TOP <- 240


### EXP1C

SC_EXP_NAME <- "1C_ARC"
TABLE <- "../tables/DATA_FRAME_EXP1C.tsv"

TOPIC_PAIRS <- "../multiome/IDR_1C_ARC-1E_ARC.csv"
IN_TABLE <- "../multiome/module_cell_AUC_1C.tsv"

OUT_TABLE <- "AUC_module_cell_AUC_1C.tsv"
HEATMAP_AUC <- "heatmap_AUC_module_cell_AUC_1C.pdf"

auc_all <- read.table(IN_TABLE, sep = "\t", header = TRUE)

TOP <- min(TOP, nrow(auc_all))
auc_all <- auc_all[1:TOP,]

topic_pairs <- read.csv2(TOPIC_PAIRS, sep = ";", header = TRUE)
topic_pairs <- topic_pairs[order(topic_pairs$repr.score, decreasing = TRUE),]
repr_score <- topic_pairs$repr.score[1:TOP]
module_size <- topic_pairs$IDR5[1:TOP]

data <- read.table(TABLE, sep = "\t", header = TRUE)
idx <- match(gsub("\\.","-",colnames(auc_all)),data$cell_ID)
data <- data[idx,]

col_annotation <- data$cl_renamed
cl <- sort(unique(col_annotation))
ncl <- length(cl)

nphe <- length(PHENOTYPE)

TOP <- nrow(auc_all)

auc_cl <- matrix(NA, nrow = TOP, ncol = ncl+nphe)
for (i in 1:TOP) {
    v <- NULL
    for (c in cl) {
        is.cl <- as.numeric(col_annotation == c)
        cell.auc <- as.numeric(auc_all[i,])
        x <- auc(response = is.cl, predictor = cell.auc, direction = "<")
        v <- c(v,x)
    }
    for (j in PHENOTYPE) {
        is.phe <- data[[j]]
        cell.auc <- as.numeric(auc_all[i,])
        x <- auc(response = is.phe, predictor = cell.auc, direction = "<")
        v <- c(v,x)   
    }
    auc_cl[i,] <- v
}
rownames(auc_cl) <- rownames(auc_all)
colnames(auc_cl) <- c(cl,PHENOTYPE_RENAMED)

write.table(auc_cl, file = OUT_TABLE, sep = "\t", quote = FALSE)

col_fun_auc <- colorRamp2(c(0.5,1), c("gray95","black"))
col_fun_cl <- brewer.pal(n = length(cl), name = "Dark2")
col_fun_phe <- brewer.pal(n = length(PHENOTYPE), name = "Dark2")
col_fun <- c(col_fun_cl, col_fun_phe)
names(col_fun) <- c(cl, PHENOTYPE_RENAMED)
idr_col_fun <- circlize::colorRamp2(breaks = c(0, 0.2, 0.5, 0.8, 1), colors = c("gray95","yellow","orange", "red", "darkred"))

ha <- HeatmapAnnotation(repr.score = repr_score, which = "row", col = list(repr.score = idr_col_fun),
                        show_legend = TRUE, show_annotation_name = FALSE)

hb <- HeatmapAnnotation(cluster = c(cl,PHENOTYPE_RENAMED), which = "column", col = list(cluster = col_fun),
                        show_legend = TRUE, show_annotation_name = FALSE)
                        
hc <- HeatmapAnnotation(size = anno_barplot(module_size, bar_width = 1, border = FALSE, gp = gpar(fill = "white")), which = "row", 
                        show_annotation_name = TRUE)
                        
h <- Heatmap(auc_cl, name = "AUC", col = col_fun_auc, 
             top_annotation = hb, left_annotation = ha, right_annotation = hc, row_names_side = "left",
             cluster_rows = FALSE, show_column_names = TRUE, show_row_names = TRUE, column_names_side = "top",
             cluster_columns = FALSE, row_title = paste("top", TOP, "modules"), row_names_gp = gpar(fontsize = 10))

pdf(HEATMAP_AUC, width = 4.3, height = 4.6)
draw(h)
dev.off()


### EXP1E

SC_EXP_NAME <- "1E_ARC"
TABLE <- "../tables/DATA_FRAME_EXP1E.tsv"

IN_TABLE <- "../multiome/module_cell_AUC_1E.tsv"

OUT_TABLE <- "AUC_module_cell_AUC_1E.tsv"
HEATMAP_AUC <- "heatmap_AUC_module_cell_AUC_1E.pdf"

auc_all <- read.table(IN_TABLE, sep = "\t", header = TRUE)

TOP <- min(TOP, nrow(auc_all))
auc_all <- auc_all[1:TOP,]

topic_pairs <- read.csv2(TOPIC_PAIRS, sep = ";", header = TRUE)
topic_pairs <- topic_pairs[order(topic_pairs$repr.score, decreasing = TRUE),]
repr_score <- topic_pairs$repr.score[1:TOP]
module_size <- topic_pairs$IDR5[1:TOP]

data <- read.table(TABLE, sep = "\t", header = TRUE)
idx <- match(gsub("\\.","-",colnames(auc_all)),data$cell_ID)
data <- data[idx,]

col_annotation <- data$cl_renamed
cl <- sort(unique(col_annotation))
ncl <- length(cl)

nphe <- length(PHENOTYPE)

TOP <- nrow(auc_all)

auc_cl <- matrix(NA, nrow = TOP, ncol = ncl+nphe)
for (i in 1:TOP) {
    v <- NULL
    for (c in cl) {
        is.cl <- as.numeric(col_annotation == c)
        cell.auc <- as.numeric(auc_all[i,])
        x <- auc(response = is.cl, predictor = cell.auc, direction = "<")
        v <- c(v,x)
    }
    for (j in PHENOTYPE) {
        is.phe <- data[[j]]
        cell.auc <- as.numeric(auc_all[i,])
        x <- auc(response = is.phe, predictor = cell.auc, direction = "<")
        v <- c(v,x)   
    }
    auc_cl[i,] <- v
}
rownames(auc_cl) <- rownames(auc_all)
colnames(auc_cl) <- c(cl,PHENOTYPE_RENAMED)

write.table(auc_cl, file = OUT_TABLE, sep = "\t", quote = FALSE)

col_fun_auc <- colorRamp2(c(0.5,1), c("gray95","black"))
col_fun_cl <- brewer.pal(n = length(cl), name = "Dark2")
col_fun_phe <- brewer.pal(n = length(PHENOTYPE), name = "Dark2")
col_fun <- c(col_fun_cl, col_fun_phe)
names(col_fun) <- c(cl, PHENOTYPE_RENAMED)
idr_col_fun <- circlize::colorRamp2(breaks = c(0, 0.2, 0.5, 0.8, 1), colors = c("gray95","yellow","orange", "red", "darkred"))

ha <- HeatmapAnnotation(repr.score = repr_score, which = "row", col = list(repr.score = idr_col_fun),
                        show_legend = TRUE, show_annotation_name = FALSE)

hb <- HeatmapAnnotation(cluster = c(cl,PHENOTYPE_RENAMED), which = "column", col = list(cluster = col_fun),
                        show_legend = TRUE, show_annotation_name = FALSE)
                        
hc <- HeatmapAnnotation(size = anno_barplot(module_size, bar_width = 1, border = FALSE, gp = gpar(fill = "white")), which = "row", 
                        show_annotation_name = TRUE)
                        
h <- Heatmap(auc_cl, name = "AUC", col = col_fun_auc, 
             top_annotation = hb, left_annotation = ha, right_annotation = hc, row_names_side = "left",
             cluster_rows = FALSE, show_column_names = TRUE, show_row_names = TRUE, column_names_side = "top",
             cluster_columns = FALSE, row_title = paste("top", TOP, "modules"), row_names_gp = gpar(fontsize = 10))

pdf(HEATMAP_AUC, width = 4.3, height = 4.6)
draw(h)
dev.off()


q()

