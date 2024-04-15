
# load libraries

library("circlize")
library("scales")
library("RColorBrewer")

CL_FIELD <- "cl_renamed"
ANN_FIELD <- c("DRC_bulk","DRC_bulk_wilcox","DRC_bulk_fisher","DRC_bulk_gam","DRC_sc",
               "TIC","TIC_wilcox","TIC_fisher","TIC_gam",
               "treated_TIC","treated_TIC_wilcox","treated_TIC_fisher","treated_TIC_gam")

### EXP1B

# files

TABLE <- "../tables/DATA_FRAME_EXP1B_P0_filt.tsv"
             
DOTPLOT <- "odds_ratio_TIC_cl_EXP1B.pdf"
OUT_TABLE <- "odds_ratio_TIC_cl_EXP1B.tsv"

df <- read.table(TABLE, sep = "\t", header = TRUE)

clusters <- sort(unique(df[,CL_FIELD]))
ncl <- length(clusters)
col <- brewer.pal(n = ncl, name = "Set1")

# compute the odds ratio for each annotation and for each cluster

df_odds_ratio <- as.data.frame(matrix(NA, nrow = 0, ncol = length(clusters)))

ncol <- 5
nrow <- ceiling(length(ANN_FIELD)/ncol)
pdf(DOTPLOT, width = 9, height = 2*nrow)
par(mfrow = c(nrow,ncol))
for (ann in ANN_FIELD) {

    odds_ratio <- c()
    for (cl in clusters) {
    
        a <- sum(df[,CL_FIELD] == cl & df[,ann] == 1)
        b <- sum(df[,CL_FIELD] != cl & df[,ann] == 1)
        c <- sum(df[,CL_FIELD] == cl & df[,ann] == 0)
        d <- sum(df[,CL_FIELD] != cl & df[,ann] == 0)
        
        # o <- (a*d)/(b*c)
        o <- log2((a*d)/(b*c))
        odds_ratio <- c(odds_ratio, o)
    }
    
    df_odds_ratio <- rbind(df_odds_ratio, odds_ratio)
    
    plot(x = clusters, y = odds_ratio, col = col, xlab = "clusters", ylab = "log-odds ratio", main = ann, pch = 16, cex = 1.6)
    abline(h = 0, lty = 2, col = "gray50")
}
dev.off()

rownames(df_odds_ratio) <- ANN_FIELD
colnames(df_odds_ratio) <- clusters

write.table(df_odds_ratio, file = OUT_TABLE, sep = "\t", quote = FALSE)


### EXP1D

TABLE <- "../tables/DATA_FRAME_EXP1D_P0_filt.tsv"
             
DOTPLOT <- "odds_ratio_TIC_cl_EXP1D.pdf"
OUT_TABLE <- "odds_ratio_TIC_cl_EXP1D.tsv"

df <- read.table(TABLE, sep = "\t", header = TRUE)

clusters <- sort(unique(df[,CL_FIELD]))
ncl <- length(clusters)
col <- brewer.pal(n = ncl, name = "Set2")

# compute the odds ratio for each annotation and for each cluster

df_odds_ratio <- as.data.frame(matrix(NA, nrow = 0, ncol = length(clusters)))

ncol <- 5
nrow <- ceiling(length(ANN_FIELD)/ncol)
pdf(DOTPLOT, width = 9, height = 2*nrow)
par(mfrow = c(nrow,ncol))
for (ann in ANN_FIELD) {

    odds_ratio <- c()
    for (cl in clusters) {
    
        a <- sum(df[,CL_FIELD] == cl & df[,ann] == 1)
        b <- sum(df[,CL_FIELD] != cl & df[,ann] == 1)
        c <- sum(df[,CL_FIELD] == cl & df[,ann] == 0)
        d <- sum(df[,CL_FIELD] != cl & df[,ann] == 0)
        
        # o <- (a*d)/(b*c)
        o <- log2((a*d)/(b*c))
        odds_ratio <- c(odds_ratio, o)
    }
    
    df_odds_ratio <- rbind(df_odds_ratio, odds_ratio)
    
    plot(x = clusters, y = odds_ratio, col = col, xlab = "clusters", ylab = "log-odds ratio", main = ann, pch = 16, cex = 1.6)
    abline(h = 0, lty = 2, col = "gray50")
}
dev.off()

rownames(df_odds_ratio) <- ANN_FIELD
colnames(df_odds_ratio) <- clusters

write.table(df_odds_ratio, file = OUT_TABLE, sep = "\t", quote = FALSE)

q()


