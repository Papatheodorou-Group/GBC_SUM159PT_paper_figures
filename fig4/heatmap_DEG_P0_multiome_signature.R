
# directory and file names

DIR <- "../analysis/multiome/signatures"
PREFIX1 <- "1C_ARC_seurat.Rds_pca_rna_0.4_clusters_wilcox_marker_genes_cluster"
PREFIX2 <- "1E_ARC_seurat.Rds_pca_rna_0.6_clusters_wilcox_marker_genes_cluster"

SIGN_1 <- file.path(DIR,"1C_1E_signatures_S1.txt")
SIGN_2 <- file.path(DIR,"1C_1E_signatures_S2.txt")
SIGN_3 <- file.path(DIR,"1C_1E_signatures_S3.txt")

SIGN_1_COMM <- file.path(DIR,"Multiome_scRNA_common_signatures_S1.txt")
SIGN_2_COMM <- file.path(DIR,"Multiome_scRNA_common_signatures_S2.txt")
SIGN_3_COMM <- file.path(DIR,"Multiome_scRNA_common_signatures_S3.txt")

HEATMAP_DATA <- "heatmap_DEG_P0_multiome_signature.tsv"
HEATMAP <- "heatmap_DEG_P0_multiome_signature.pdf"
HEATMAP_HORIZ <- "heatmap_DEG_P0_multiome_signature_horizontal.pdf"

CLS1 <- c("6","3","4")
CLS2 <- c("5","3","4")
CLS1_NAMES <- paste0("rep1",c(""," ", "   "))
CLS2_NAMES <- paste0("rep2",c(""," ", "   "))

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

sign1 <- as.character(read.table(SIGN_1)[,1])
sign2 <- as.character(read.table(SIGN_2)[,1])
sign3 <- as.character(read.table(SIGN_3)[,1])

sign1_comm <- as.character(read.table(SIGN_1_COMM)[,1])
sign2_comm <- as.character(read.table(SIGN_2_COMM)[,1])
sign3_comm <- as.character(read.table(SIGN_3_COMM)[,1])

genes <- unique(c(sign1, sign2, sign3))
sign_comm <- c(sign1_comm, sign2_comm, sign3_comm)

# extract the genes and the log2FC, for each cluster
# set 0 if not DE
M <- matrix(0, ncol = length(CLS1_NAMES)+length(CLS2_NAMES), nrow = length(genes))
rownames(M) <- genes
colnames(M) <- unlist(lapply(1:3, function(x) c(CLS1_NAMES[x],CLS2_NAMES[x])))
i <- 1
for (cl in CLS1) {
    DEG_TABLE <- file.path(DIR, paste0(PREFIX1,cl,".tsv"))
    deg_table <- read.table(DEG_TABLE, sep = "\t", header = TRUE)
    genes_cl <- deg_table$geneID
    common <- genes_cl[genes_cl %in% genes]
    M[match(common,genes),i] <- deg_table$avg_log2FC[match(common,genes_cl)]
    i <- i + 2
}
i <- 2
for (cl in CLS2) {
    DEG_TABLE <- file.path(DIR, paste0(PREFIX2,cl,".tsv"))
    deg_table <- read.table(DEG_TABLE, sep = "\t", header = TRUE)
    genes_cl <- deg_table$geneID
    common <- genes_cl[genes_cl %in% genes]
    M[match(common,genes),i] <- deg_table$avg_log2FC[match(common,genes_cl)]
    i <- i + 2
}

write.table(M, file = HEATMAP_DATA, sep = "\t", quote = FALSE)

# load the signature and map it to the matrix

# select only the genes in the signatures
gene_labels <- c( head(sign1,20), 
            head(sign2,5), 
            head(sign3,5) )

gene_labels_id <- match(gene_labels, rownames(M))

# highlight the markers found in scRNA-Seq too

gene_labels_highlight <- which(gene_labels %in% sign_comm)

# select palette for matrix entries (accuracy)

col_fun <- colorRamp2(c(0,max(M)/3,max(M)), c("gray90","red","darkred"))

# split by signature

split <- c()
genes_part <- genes
i <- 1
for (s in list(sign1,sign2,sign3)) {
    split <- c(split, rep(i, sum(s %in% genes_part)))
    genes_part <- genes_part[!(genes_part %in% s)]
    i <- i+1
}

fontface <- rep("plain",length(gene_labels))
fontface[gene_labels_highlight] <- rep("bold",length(gene_labels_highlight))

# annotate only the genes in common between experiments
anno <- HeatmapAnnotation(highlight = anno_mark(at = gene_labels_id, labels = gene_labels, side = 'right', 
           labels_gp = gpar(fontsize = 9, fontface = fontface)), 
           which = "row")
h <- Heatmap(M, name = "log2(FC)", col = col_fun, right_annotation = anno,
	cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, show_column_names = TRUE, 
	column_names_side = "top", # column_names_rot = 45, 
	column_split = c(1,1,2,2,3,3), 
	row_split = split, show_parent_dend_line = FALSE, column_names_rot = 45,
	column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 9))

pdf(HEATMAP, width = 3, height = 5)
draw(h)
dev.off()

CLS1_NAMES <- paste0(c(""," ", "   "),"rep1")
CLS2_NAMES <- paste0(c(""," ", "   "),"rep2")

colnames(M) <- unlist(lapply(1:3, function(x) c(CLS1_NAMES[x],CLS2_NAMES[x])))

# horizontal layout
anno <- HeatmapAnnotation(highlight = anno_mark(at = gene_labels_id, labels = gene_labels, side = 'top', 
           labels_rot = 60, labels_gp = gpar(fontsize = 9, fontface = fontface)), which = "column")
h <- Heatmap(t(M), name = "log2(FC)", col = col_fun, 
	cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = TRUE, show_column_names = FALSE, 
	row_names_side = "left", top_annotation = anno,
	row_split = c(1,1,2,2,3,3), 
	column_split = split, 
	column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 9))

pdf(HEATMAP_HORIZ, width = 6, height = 2.8)
draw(h)
dev.off()


q()
