
# directory and file names

TABLE <- "../tables/DATA_FRAME_GBC_IN_VIVO.tsv"

SAMPLES <- c("P0_R4","P0_PERT2","P0_PERT3")

CUMPLOT <- "cumulative_GBC_tumor_P0.pdf"

# load libraries

library("ggplot2")
library("scales")
library("RColorBrewer")

# load data

df <- read.table(TABLE, sep = "\t", header = TRUE)
df_count <- as.data.frame(matrix(NA, nrow = 0, ncol = length(SAMPLES)))

for (sample in SAMPLES) {
    sample_id <- paste0("cpm.",sample)
    count <- df[[sample_id]]
    count <- sort(count[count > 0], decreasing = TRUE)
    s <- rep(sample, length(count))
    c <- 1:length(count)
    cum <- cumsum(count)/sum(count)
    df1 <- data.frame(sample = s, clone.id = c, clone.id.frac = c/length(count), GBC.count = cum, count = count)
    df_count <- rbind(df_count, df1)
}
colnames(df_count) <- c("sample", "clone.id", "clone.id.frac", "GBC.count", "count")

# generate figure

g <- ggplot(data = df_count, aes(x = clone.id, y = GBC.count)) + theme_classic() + geom_point(alpha = 0.5, size=1) 
g <- g + xlab("GBC sorted by CPM") + ylab("Cumulative CPM fraction") + ggtitle("Parental")
g <- g + theme(legend.title=element_blank()) + scale_x_continuous(trans = log10_trans())

pdf(CUMPLOT, width = 2.7, height = 2.8)
g
dev.off()


q()




q()
