
# directory and file names

TABLE <- "../tables/DATA_FRAME_GBC_IN_VIVO.tsv"

SAMPLES <- c("P0_R4","P0_PERT2","P0_PERT3",
             "A3PT_PERT3","C3M1_PERT2","C3PT_PERT1","D1PT_PERT3","D4PT_PERT3","D5PT_PER3","C1PT_PTX","C2PT_PTX","D3PT_NT")
CLASSES <- c(rep("Parental",3),rep("Tumour",9))
CLASSES_UNIQ <- c("Parental","Tumour")

BOXPLOT <- "boxplot_GBC_count_tumor_P0.pdf"

# load libraries

library("ggplot2")
library("scales")
library("RColorBrewer")

# load data

df <- read.table(TABLE, sep = "\t", header = TRUE)
idx <- match(paste0("cpm.",SAMPLES), colnames(df))
count <- colSums(df[,idx] > 0)
df_count <- data.frame(sample = SAMPLES, class = CLASSES, count = count)

df_count$class <- factor(df_count$class, levels = CLASSES_UNIQ)

# generate figure (abs)

g <- ggplot(data = df_count, aes(x = class, y = count)) + theme_classic() + geom_boxplot() 
g <- g + xlab("") + ylab("GBC count") 
g <- g + theme(axis.text.x = element_text(angle = 45, hjust=1))
g <- g + theme(legend.title=element_blank()) 

pdf(BOXPLOT, width = 1.5, height = 2.8)
g
dev.off()




q()
