
# load libraries

library("ggplot2")
library("scales")
library("RColorBrewer")

### EXP1

# directory and file names

EXP <- "EXP1"
TABLE <- paste0("../tables/DATA_FRAME_GBC_IN_VITRO_", EXP, ".tsv")

SAMPLES <- c("D0","S1_P0_1","S2_P0_2","P_end_A","P_end_B","P_end_C","D13","D17","D24","D15_A","D15_B","D15_C")
CLASSES <- c(rep("Parental",3),rep("Untreated",3),rep("Treated",6))
CLASSES_UNIQ <- c("Parental","Untreated","Treated")

BOXPLOT <- paste0("boxplot_GBC_count_invitro_",EXP,".pdf")



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

pdf(BOXPLOT, width = 1.8, height = 2.8)
g
dev.off()


### EXP4 

# directory and file names

EXP <- "EXP4"
TABLE <- paste0("../tables/DATA_FRAME_GBC_IN_VITRO_", EXP, ".tsv")

SAMPLES <- c("P0","P_end_1","P_end_2","P_end_3","D15_A","D15_B","D15_C","D15_D","D15_E")
CLASSES <- c(rep("Parental",1),rep("Untreated",3),rep("Treated",5))
CLASSES_UNIQ <- c("Parental","Untreated","Treated")

BOXPLOT <- paste0("boxplot_GBC_count_invitro_",EXP,".pdf")

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

pdf(BOXPLOT, width = 1.8, height = 2.8)
g
dev.off()


### EXP5

# directory and file names

EXP <- "EXP5"
TABLE <- paste0("../tables/DATA_FRAME_GBC_IN_VITRO_", EXP, ".tsv")

SAMPLES <- c("P5OR_GBC_P0_1","P5OR_GBC_P0_2","P5OR_GBC_PACLI_D17","P5OR_GBC_PACLI_D20")
CLASSES <- c(rep("Parental",2),rep("Treated",2))
CLASSES_UNIQ <- c("Parental","Treated")

BOXPLOT <- paste0("boxplot_GBC_count_invitro_",EXP,".pdf")

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
