
# directory and file names

DIR <- "../analysis/gene_expression/1D_GEX/treated_filt/DEA"

METHOD <- "MAST"

CLMODE <- "clusters_pca_vst_top5000_k30_res0.1"

HEATMAP_DATA <- "heatmap_DEG_EXP1D_treated_cl.tsv"
HEATMAP <- "heatmap_DEG_EXP1D_treated_cl.pdf"

CLS <- c(2,1,0)
CLS_NAMES <- c(1,2,3)

# LOG2FC <- 0.5
LOG2FC <- 0
PVAL <- 0.05
TOP <- 10

# load libraries

library("ggplot2")
library("ggrepel")
library("ComplexHeatmap")
library("circlize")
library("scales")
library("RColorBrewer")

# collect all DEGs first
genes <- NULL
top_genes <- NULL
first_top_gene <- NULL
i <- 1
for (cl in CLS) {
    first_top_gene <- c(first_top_gene, length(genes)+1)
    DEG_TABLE <- file.path(DIR, CLMODE, METHOD, paste0("DEG_MAST_cl",cl,"-all.tsv"))
    deg_table <- read.table(DEG_TABLE, sep = "\t", header = TRUE)
    idx <- which(deg_table$avg_log2FC > LOG2FC & deg_table$p_val_adj < PVAL)
    v <- 1:length(idx)
    rank1 <- v[order(deg_table$avg_log2FC[idx], decreasing = TRUE)]
    rank2 <- v[order(deg_table$p_val_adj[idx])]
    w <- apply(matrix(c(rank1,rank2), nrow = length(idx), ncol = 2), 1, mean)
    avg_rank <- v[order(w)]
    ordered_genes <- deg_table$geneID[idx[avg_rank]]
    common_genes <- genes[genes %in% ordered_genes]
    if (i > 1) {
        for (j in 2:length(first_top_gene)) {
            first_top_gene[j] <- first_top_gene[j] - sum(genes[first_top_gene[j-1]:(first_top_gene[j]-1)] %in% common_genes) + 1
        }
    }
    genes <- c(genes[!(genes %in% common_genes)], ordered_genes[!(ordered_genes %in% common_genes)])
    i <- i + 1
}
first_top_gene <- c(first_top_gene, i)
genes <- unique(genes)
top_genes <- NULL
for (i in first_top_gene) {
    top_genes <- c(top_genes, genes[i:(i+TOP-1)])
}

# extract the genes and the log2FC, for each cluster
# set 0 if not DE
M <- matrix(0, ncol = length(CLS_NAMES), nrow = length(genes))
colnames(M) <- CLS_NAMES
rownames(M) <- genes
i <- 1
for (cl in CLS) {
    DEG_TABLE <- file.path(DIR, CLMODE, METHOD, paste0("DEG_MAST_cl",cl,"-all.tsv"))
    deg_table <- read.table(DEG_TABLE, sep = "\t", header = TRUE)
    idx <- which(deg_table$geneID %in% genes)
    genes_cl <- deg_table$geneID[idx]
    log2FC <- deg_table$avg_log2FC[idx]
    idx2 <- match(genes_cl, genes)
    M[idx2,i] <- log2FC
    i <- i + 1
}

write.table(M, file = HEATMAP_DATA, sep = "\t", quote = FALSE)

# load the signature and map it to the matrix

# select only the genes in the signatures
gene_labels <- top_genes
gene_labels_id <- match(gene_labels, rownames(M))

# select palette for clusters

cols <- hue_pal()(length(CLS))
names(cols) <- CLS_NAMES
cols <- cols[match(names(cols),colnames(M))]

col_fun <- colorRamp2(c(min(M),0,max(M)/3,max(M)), c("blue","gray90","red","darkred"))

# annotate only the genes in common between experiments
anno <- HeatmapAnnotation(highlight = anno_mark(at = gene_labels_id, labels = gene_labels, side = 'right', labels_gp = gpar(fontsize = 9)), 
           which = "row")
ha <- HeatmapAnnotation(cluster = colnames(M), which = "column", col = list(cluster = cols), show_legend = FALSE)
h <- Heatmap(M, name = "log2(FC)", col = col_fun, top_annotation = ha, right_annotation = anno,
	cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, show_column_names = TRUE, 
	column_names_side = "top", # column_names_rot = 45, 
	column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 9))

pdf(HEATMAP, width = 2.5, height = 5)
draw(h)
dev.off()


q()
