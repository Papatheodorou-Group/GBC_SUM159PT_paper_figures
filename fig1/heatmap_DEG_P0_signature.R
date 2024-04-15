
# directory and file names

DIR1 <- "../analysis/gene_expression/1B_GEX/P0_ccRegress_filt_2ndRound/DEA"
DIR2 <- "../analysis/gene_expression/1D_GEX/P0_ccRegress_filt/DEA"
SIGN_DIR <- "../analysis/gene_expression/signatures/DEG_P0_ccRegress"

SIGN_1 <- file.path(SIGN_DIR,"signature1.txt")
SIGN_2 <- file.path(SIGN_DIR,"signature2.txt")
SIGN_3 <- file.path(SIGN_DIR,"signature3.txt")

METHOD <- "MAST"

CLMODE1 <- "clusters_pca_vst_top1000_k40_res0.5"
CLMODE2 <- "clusters_pca_vst_top1000_k40_res0.5"

SIGNATURE <- "DEG_P0_signature.tsv"
HEATMAP_DATA <- "heatmap_DEG_P0_signature.tsv"
HEATMAP <- "heatmap_DEG_P0_signature.pdf"
HEATMAP_HORIZ <- "heatmap_DEG_P0_signature_horizontal.pdf"

CLS1 <- c("6","2","5")
CLS2 <- c("6","2","4")
CLS1_NAMES <- paste0("T0",c(""," ", "   "))
CLS2_NAMES <- paste0("T1",c(""," ", "   "))

# LOG2FC <- 0.5
LOG2FC <- 0
PVAL <- 0.05

# load libraries

library("ggplot2")
library("ggrepel")
library("ComplexHeatmap")
library("circlize")
library("scales")
library("RColorBrewer")

# load data and compile matrix

sign1 <- read.table(SIGN_1)[,1]
sign2 <- read.table(SIGN_2)[,1]
sign3 <- read.table(SIGN_3)[,1]

genes <- unique(c(sign1, sign2, sign3))

# extract the genes and the log2FC, for each cluster
# set 0 if not DE
M <- matrix(0, ncol = length(CLS1_NAMES)+length(CLS2_NAMES), nrow = length(genes))
rownames(M) <- genes
colnames(M) <- unlist(lapply(1:3, function(x) c(CLS1_NAMES[x],CLS2_NAMES[x])))
i <- 1
for (cl in CLS1) {
    DEG_TABLE <- file.path(DIR1, CLMODE1, METHOD, paste0("DEG_MAST_cl",cl,"-all.tsv"))
    deg_table <- read.table(DEG_TABLE, sep = "\t", header = TRUE)
    genes_cl <- deg_table$geneID
    common <- genes_cl[genes_cl %in% genes]
    M[match(common,genes),i] <- deg_table$avg_log2FC[match(common,genes_cl)]
    i <- i + 2
}
i <- 2
for (cl in CLS2) {
    DEG_TABLE <- file.path(DIR2, CLMODE2, METHOD, paste0("DEG_MAST_cl",cl,"-all.tsv"))
    deg_table <- read.table(DEG_TABLE, sep = "\t", header = TRUE)
    genes_cl <- deg_table$geneID
    common <- genes_cl[genes_cl %in% genes]
    M[match(common,genes),i] <- deg_table$avg_log2FC[match(common,genes_cl)]
    i <- i + 2
}

write.table(M, file = HEATMAP_DATA, sep = "\t", quote = FALSE)

# load the signature and map it to the matrix

# select only the genes in the signatures
gene_labels <- c( head(read.table(SIGN_1)[,1],20), 
            head(read.table(SIGN_2)[,1],5), 
            head(read.table(SIGN_3)[,1],5) )

gene_labels_id <- match(gene_labels, rownames(M))

# select palette for clusters

col1 <- brewer.pal(n = length(CLS1_NAMES), name = "Set1")
col2 <- brewer.pal(n = length(CLS2_NAMES), name = "Set2")
cols <- c(col1,col2)
names(cols) <- c(CLS1_NAMES,CLS2_NAMES)
cols <- cols[match(names(cols),colnames(M))]

# select palette for matrix entries (accuracy)

col_fun <- colorRamp2(c(0,max(M)/3,max(M)), c("gray90","red","darkred"))

# figure with F1

split <- c()
genes_part <- genes
i <- 1
for (s in list(sign1,sign2,sign3)) {
    split <- c(split, rep(i, sum(s %in% genes_part)))
    genes_part <- genes_part[!(genes_part %in% s)]
    i <- i+1
}

# annotate only the genes in common between experiments
anno <- HeatmapAnnotation(highlight = anno_mark(at = gene_labels_id, labels = gene_labels, side = 'right', labels_gp = gpar(fontsize = 9)), 
           which = "row")
ha <- HeatmapAnnotation(cluster = colnames(M), which = "column", col = list(cluster = cols), show_legend = FALSE)
h <- Heatmap(M, name = "log2(FC)", col = col_fun, top_annotation = ha, right_annotation = anno,
	cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, show_column_names = TRUE, 
	column_names_side = "top", # column_names_rot = 45, 
	column_split = c(1,1,2,2,3,3), 
	row_split = split, show_parent_dend_line = FALSE,
	column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 9))

pdf(HEATMAP, width = 3, height = 5)
draw(h)
dev.off()

CLS1_NAMES <- paste0(c(""," ", "   "),"T0")
CLS2_NAMES <- paste0(c(""," ", "   "),"T1")

colnames(M) <- unlist(lapply(1:3, function(x) c(CLS1_NAMES[x],CLS2_NAMES[x])))

# select palette for clusters

col1 <- brewer.pal(n = length(CLS1_NAMES), name = "Set1")
col2 <- brewer.pal(n = length(CLS2_NAMES), name = "Set2")
cols <- c(col1,col2)
names(cols) <- c(CLS1_NAMES,CLS2_NAMES)
cols <- cols[match(names(cols),colnames(M))]

# horizontal layout
anno <- HeatmapAnnotation(highlight = anno_mark(at = gene_labels_id, labels = gene_labels, side = 'top', 
           labels_rot = 60, labels_gp = gpar(fontsize = 9)), which = "column")
ha <- HeatmapAnnotation(cluster = colnames(M), which = "row", col = list(cluster = cols), show_legend = FALSE)
h <- Heatmap(t(M), name = "log2(FC)", col = col_fun, left_annotation = ha, top_annotation = anno,
	cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = TRUE, show_column_names = FALSE, 
	row_names_side = "left", 
	row_split = c(1,1,2,2,3,3), 
	column_split = split, 
	column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 9))

pdf(HEATMAP_HORIZ, width = 6, height = 3)
draw(h)
dev.off()


q()
