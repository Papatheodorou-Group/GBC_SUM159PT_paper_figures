
# load libraries

library("pROC")
library("RColorBrewer")
library("circlize")
library("ComplexHeatmap")


### EXP1C

SC_EXP_NAME <- "1C_ARC"
TABLE <- "../tables/DATA_FRAME_EXP1C.tsv"

TOPIC_PAIRS <- "../multiome/IDR_1C_ARC-1E_ARC.csv"
IN_TABLE <- "../multiome/module_cell_AUC_1C.tsv"
IN_AUC_TABLE <- "AUC_module_cell_AUC_1C.tsv"

HEATMAP_SCALED <- "heatmap_module_cell_AUC_1C_scaled_selected.pdf"

MIN_SIZE <- 20 # minimum number of regions per module to keep it

auc_all <- read.table(IN_TABLE, sep = "\t", header = TRUE)

topic_pairs <- read.csv2(TOPIC_PAIRS, sep = ";", header = TRUE)
topic_pairs <- topic_pairs[order(topic_pairs$repr.score, decreasing = TRUE),]

module_size <- topic_pairs$IDR5
idx <- which(module_size >= MIN_SIZE)
topic_pairs <- topic_pairs[idx,]
auc_all <- auc_all[idx,]

TOP <- min(20,nrow(auc_all))
auc_all <- auc_all[1:TOP,]
repr_score <- topic_pairs$repr.score[1:TOP]
module_size <- topic_pairs$IDR5[1:TOP]

auc_cl <- read.table(IN_AUC_TABLE, sep = "\t", header = TRUE)
auc_cl <- cbind(auc_cl[idx,1:3], auc_cl[idx,c("TIC","DTC_invitro")])
auc_cl <- auc_cl[1:TOP,]
colnames(auc_cl) <- c(paste0("S",1:3), "TIC", "DTC in vitro")

data <- read.table(TABLE, sep = "\t", header = TRUE)
idx <- match(gsub("\\.","-",colnames(auc_all)),data$cell_ID)
data <- data[idx,]

ordering <- order(data$cl_renamed)
col_annotation <- sort(data$cl_renamed)
auc_all <- auc_all[,ordering]

cl <- sort(unique(col_annotation))
ncl <- length(cl)

col_fun <- colorRamp2(c(0.5,0.7), c("white","black"))
col_fun_auc <- colorRamp2(c(0,0.2,0.5,0.8,1), c("darkblue","blue","gray95","red","darkred"))
col_fun_scaled <- colorRamp2(c(0,1), c("white","black"))
col_fun_cl <- brewer.pal(n = length(cl), name = "Dark2")
names(col_fun_cl) <- cl
idr_col_fun <- circlize::colorRamp2(breaks = c(0, 0.2, 0.5, 0.8, 1), colors = c("gray95","yellow","orange", "red", "darkred"))

# heatmaps with all clusters

auc_all_scaled <- t(scale(t(auc_all)))

ha <- HeatmapAnnotation(repr.score = repr_score, which = "row", col = list(repr.score = idr_col_fun),
                        show_legend = TRUE, show_annotation_name = FALSE)

hb <- HeatmapAnnotation(cluster = col_annotation, which = "column", col = list(cluster = col_fun_cl),
                        show_legend = TRUE, show_annotation_name = FALSE)
                        
hc <- HeatmapAnnotation(size = anno_barplot(module_size, bar_width = 1, border = FALSE, gp = gpar(fill = "white")), which = "row", 
                        show_annotation_name = TRUE)

h1 <- Heatmap(auc_all_scaled, name = "Zscore(AUC)", col = col_fun_scaled, column_split = col_annotation, 
             top_annotation = hb, left_annotation = ha, row_names_side = "left",
             cluster_rows = FALSE, cluster_columns = TRUE, show_column_names = FALSE, show_row_names = TRUE, 
             cluster_column_slices = FALSE, show_row_dend = TRUE,
             row_names_gp = gpar(fontsize = 10), width = unit(9, "cm")) 
             
h2 <- Heatmap(auc_cl, name = "AUC.cl", col = col_fun_auc, right_annotation = hc, 
             column_title = "cluster\nAUC", column_title_gp = gpar(fontsize = 10),
             cluster_rows = FALSE, show_column_names = TRUE, show_row_names = FALSE, column_names_side = "top",
             cluster_columns = FALSE, width = unit(2, "cm"))        

pdf(HEATMAP_SCALED, width = 7, height = 5.2)
draw(h1+h2, ht_gap = unit(0,"cm"))
dev.off()  


### EXP1E

SC_EXP_NAME <- "1E_ARC"
TABLE <- "../tables/DATA_FRAME_EXP1E.tsv"

IN_TABLE <- "../multiome/module_cell_AUC_1E.tsv"
IN_AUC_TABLE <- "AUC_module_cell_AUC_1E.tsv"

HEATMAP_SCALED <- "heatmap_module_cell_AUC_1E_scaled_selected.pdf"

MIN_SIZE <- 20 # minimum number of regions per module to keep it

auc_all <- read.table(IN_TABLE, sep = "\t", header = TRUE)

topic_pairs <- read.csv2(TOPIC_PAIRS, sep = ";", header = TRUE)
topic_pairs <- topic_pairs[order(topic_pairs$repr.score, decreasing = TRUE),]

module_size <- topic_pairs$IDR5
idx <- which(module_size >= MIN_SIZE)
topic_pairs <- topic_pairs[idx,]
auc_all <- auc_all[idx,]

TOP <- min(20,nrow(auc_all))
auc_all <- auc_all[1:TOP,]
repr_score <- topic_pairs$repr.score[1:TOP]
module_size <- topic_pairs$IDR5[1:TOP]

auc_cl <- read.table(IN_AUC_TABLE, sep = "\t", header = TRUE)
auc_cl <- cbind(auc_cl[idx,1:3], auc_cl[idx,c("TIC","DTC_invitro")])
auc_cl <- auc_cl[1:TOP,]
colnames(auc_cl) <- c(paste0("S",1:3), "TIC", "DTC in vitro")

data <- read.table(TABLE, sep = "\t", header = TRUE)
idx <- match(gsub("\\.","-",colnames(auc_all)),data$cell_ID)
data <- data[idx,]

ordering <- order(data$cl_renamed)
col_annotation <- sort(data$cl_renamed)
auc_all <- auc_all[,ordering]

cl <- sort(unique(col_annotation))
ncl <- length(cl)

col_fun <- colorRamp2(c(0.5,0.7), c("white","black"))
col_fun_auc <- colorRamp2(c(0,0.2,0.5,0.8,1), c("darkblue","blue","gray95","red","darkred"))
col_fun_scaled <- colorRamp2(c(0,1), c("white","black"))
col_fun_cl <- brewer.pal(n = length(cl), name = "Dark2")
names(col_fun_cl) <- cl
idr_col_fun <- circlize::colorRamp2(breaks = c(0, 0.2, 0.5, 0.8, 1), colors = c("gray95","yellow","orange", "red", "darkred"))

# heatmaps with all clusters

auc_all_scaled <- t(scale(t(auc_all)))

ha <- HeatmapAnnotation(repr.score = repr_score, which = "row", col = list(repr.score = idr_col_fun),
                        show_legend = TRUE, show_annotation_name = FALSE)

hb <- HeatmapAnnotation(cluster = col_annotation, which = "column", col = list(cluster = col_fun_cl),
                        show_legend = TRUE, show_annotation_name = FALSE)
                        
hc <- HeatmapAnnotation(size = anno_barplot(module_size, bar_width = 1, border = FALSE, gp = gpar(fill = "white")), which = "row", 
                        show_annotation_name = TRUE)

h1 <- Heatmap(auc_all_scaled, name = "Zscore(AUC)", col = col_fun_scaled, column_split = col_annotation, 
             top_annotation = hb, left_annotation = ha, row_names_side = "left",
             cluster_rows = FALSE, cluster_columns = TRUE, show_column_names = FALSE, show_row_names = TRUE, 
             cluster_column_slices = FALSE, show_row_dend = TRUE,
             row_names_gp = gpar(fontsize = 10), width = unit(9, "cm")) 
             
h2 <- Heatmap(auc_cl, name = "AUC.cl", col = col_fun_auc, right_annotation = hc, 
             column_title = "cluster\nAUC", column_title_gp = gpar(fontsize = 10),
             cluster_rows = FALSE, show_column_names = TRUE, show_row_names = FALSE, column_names_side = "top",
             cluster_columns = FALSE, width = unit(2, "cm"))        

pdf(HEATMAP_SCALED, width = 7, height = 5.2)
draw(h1+h2, ht_gap = unit(0,"cm"))
dev.off()  


q()

