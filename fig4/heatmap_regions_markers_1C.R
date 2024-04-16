
TOPIC_PREFIX <- "../multiome/IDR_1C_ARC-1E_ARC_"

TOPIC_1C <- c("2","31","27","15")
TOPIC_1E <- c("16","15","14","19")

CLUSTER_ENRICH <- list(c("S1","S3"),c("S2"),c("S1"),c("S3"))

OBJECT <- "../seurat_objects/P0_multiome_EXP1C.Rds"

SC_EXP_NAME <- "1C_ARC"
TABLE <- "../tables/DATA_FRAME_EXP1C.tsv"

# TOP <- 50
TOP_PERC <- 5
TOP_MIN <- 20

library("Seurat")  
library("Signac")
library("GenomicRanges")
library("RColorBrewer")
library("circlize")
library("ComplexHeatmap")

object <- readRDS(OBJECT)
object <- RunTFIDF(object, assay = "ATAC")

# TODO: only select the markers associated to the populations a models is enriched in -> DISABLED

data_sel <- as.data.frame(matrix(NA, nrow = 0, ncol = 15))
idx_cl_enrich <- c()
TOP <- rep(0,length(TOPIC_1C))
for (i in 1:length(TOPIC_1C)) {
    
    filename <- paste0(TOPIC_PREFIX, TOPIC_1C[i], "-", TOPIC_1E[i], ".csv")
    data <- read.csv2(filename, sep = ";", header = TRUE)

    # sort by markers first, so that the region with the marker is reported on top if probs are equal
    idx <- which(unlist(lapply(data$cluster_50kb, function(x) sum(strsplit(x, split="/") %in% CLUSTER_ENRICH[[i]]))) > 0)
    num_markers <- rep(0, nrow(data))
#    num_markers[idx] <- apply(data[idx,c("marker_500kb","marker_50kb","marker_5kb")], 1, function(x) sum(!is.na(x)))
    num_markers[idx] <- apply(data[idx,c("marker_50kb","marker_5kb")], 1, function(x) sum(!is.na(x)))
#    num_markers <- apply(data[,c("marker_500kb","marker_50kb","marker_5kb")], 1, function(x) sum(!is.na(x)))
    data <- data[order(num_markers, decreasing = TRUE),]
    data <- data[order(data$gIDR, decreasing = TRUE),]
#    data <- data[order(data$IDR, decreasing = TRUE),]
    
    TOP[i] <- max(TOP_MIN, floor(nrow(data)*TOP_PERC/100))
    if (i == 1) {
        idx_cl_enrich <- c(idx_cl_enrich, idx[idx <= TOP[i]])
    } else {
        idx_cl_enrich <- c(idx_cl_enrich, sum(TOP[1:i-1]) + idx[idx <= TOP[i]])
    }
    
    data_sel <- rbind(data_sel, data[1:TOP[i],])
    
}

# granges for regions with top IDR
gr <- GRanges(seqnames = Rle(data_sel$chr), 
              ranges = IRanges(data_sel$start, end = data_sel$end, names = 1:nrow(data_sel)), 
              strand = rep("*",nrow(data_sel)))
                  
# granges for all regions
gr_atac <- object@assays$ATAC@ranges
# ov_loc <- subsetByOverlaps(gr_atac, gr, select = "first")
ov_loc <- findOverlaps(gr_atac, gr)

# TFIDF matrix 
# object <- subset(object, features = rownames(object)[ov_loc@from])
M <- object@assays$ATAC@data
M <- M[ov_loc@from,]
# M <- t(scale(t(M)))

idx <- order(ov_loc@to) 
M <- M[idx,] # rows ordered according to the IDR data frame
which_info <- ov_loc@to[idx] # row index in the IDR data frame

# idx_lab <- which(!is.na(data_sel$marker_50kb)) # indexes in the IDR data frame that contain the markers
idx_lab <- idx_cl_enrich
idx_lab_M <- which(which_info %in% idx_lab)
which_ann <- idx_lab[match(which_info[idx_lab_M], idx_lab)] # indexes in the IDR data frame that contain the markers matched with M rows
markers <- data_sel$marker_50kb[which_ann]
regions_with_marker <- paste(data_sel$chr, data_sel$start, data_sel$end, sep = "-")[which_ann]

at_markers <- match(which_ann, which_info) # M row indexes that contain IDR markers

n1 <- sum(which_info <= TOP[1])
n2 <- sum(which_info <= sum(TOP[1:2])) - n1
n3 <- sum(which_info <= sum(TOP[1:3])) - n1 - n2
n4 <- sum(which_info <= sum(TOP)) - n1 - n2 - n3
module <- c(rep("Module 1",n1),rep("Module 2",n2),rep("Module 3",n3),rep("Module 4",n4))

data <- read.table(TABLE, sep = "\t", header = TRUE)
idx <- match(gsub("\\.","-",colnames(M)),data$cell_ID)
data <- data[idx,]

ordering <- order(data$cl_renamed)
col_annotation <- sort(data$cl_renamed)
M <- M[,ordering]

cl <- unique(col_annotation)
ncl <- length(cl)
   
col_fun <- colorRamp2(c(0,1), c("gray90","black"))
col_fun_cl <- brewer.pal(n = length(cl), name = "Dark2")
names(col_fun_cl) <- cl

gene_labels <- rowAnnotation(highlight = anno_mark(at = at_markers, labels = markers, 
                         which = "row", side = 'right', labels_gp = gpar(fontsize = 9)))
region_labels <- rowAnnotation(highlight = anno_mark(at = at_markers, labels = regions_with_marker, 
                         which = "row", side = 'left', labels_gp = gpar(fontsize = 9)))

hb <- HeatmapAnnotation(cluster = col_annotation, which = "column", col = list(cluster = col_fun_cl),
                        show_legend = TRUE, show_annotation_name = FALSE)
                        
h <- Heatmap(as.matrix(M), name = "TFIDF", col = col_fun, row_split = module, column_split = col_annotation, top_annotation = hb,
                     cluster_rows = FALSE, show_column_names = FALSE, show_row_names = FALSE,
                     # right_annotation = gene_labels, 
                     left_annotation = region_labels,
                     cluster_columns = TRUE, cluster_column_slices = FALSE, cluster_row_slices = FALSE,
                     show_row_dend = FALSE, show_column_dend = FALSE, row_title = "top 5% reproducible regions")

pdf("heatmap_region_markers_1C.pdf", width = 7, height = 4.5)
draw(h)
dev.off()



