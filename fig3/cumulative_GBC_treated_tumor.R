
# directory and file names

TABLE <- "../tables/DATA_FRAME_GBC_IN_VIVO.tsv"

SAMPLES <- c("P0_R4,P0_PERT2,P0_PERT3","A2PT_R4,A1PT_PERT1,A5V2PT_PERT1,B1PT_PERT3,B2PT_PTX_Per3,B3PT_PTX_Per3")

TYPE <- c("Parental","Treated")

CUMPLOT <- "cumulative_GBC_treated_tumor.pdf"

# load libraries

library("ggplot2")
library("scales")
library("RColorBrewer")

# load data

df <- read.table(TABLE, sep = "\t", header = TRUE)

# compute range

all_samples <- unlist(strsplit(paste(SAMPLES, collapse = ","), split = ","))
max_gbc_count <- max(colSums(df[,paste0("cpm.",all_samples)] > 0))

df_count <- as.data.frame(matrix(NA, nrow = 0, ncol = 6))

for (i in 1:length(SAMPLES)) {

    samples <- unlist(strsplit(SAMPLES[i], split = ","))

    # compute the cumulative distribution
    
    for (sample in samples) {
        sample_id <- paste0("cpm.",sample)
        count <- df[[sample_id]]
        count <- sort(count[count > 0], decreasing = TRUE)
        s <- rep(sample, length(count))
        c <- 1:length(count)
        cum <- cumsum(count)/sum(count)
        t <- rep(TYPE[i], length(cum))
        df1 <- data.frame(sample = s, clone.id = c, clone.id.frac = c/length(count), GBC.count = cum, count = count, type = t)
        df_count <- rbind(df_count, df1)
    }
    
}
colnames(df_count) <- c("sample", "clone.id", "clone.id.frac", "GBC.count", "count", "type")

# select palette

colors <- c("gray80","black")

# generate figure

df_count$type <- factor(df_count$type, levels = TYPE)

g <- ggplot(data = df_count, aes(x = clone.id, y = GBC.count, color = type)) + theme_classic() + geom_point(alpha = 0.5, size=1) 
g <- g + xlab("GBC sorted by CPM") + ylab("Cumulative CPM fraction")
g <- g + scale_x_continuous(trans = log10_trans(), limits = c(1,max_gbc_count))
g <- g + theme(legend.title=element_blank()) + scale_color_manual(values = colors) + theme(legend.position="top")

pdf(CUMPLOT, width = 2.5, height = 2.8)
g
dev.off()



q()

