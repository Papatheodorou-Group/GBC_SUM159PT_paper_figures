             
# load libraries

library("ggplot2")
library("scales")
library("RColorBrewer")

# select palette

colors <- c("gray80","gray50","black")


### EXP1

# directory and file names

TABLE <- "../tables/DATA_FRAME_GBC_IN_VITRO_EXP1.tsv"
EXP <- "EXP1"

SAMPLES <- c("D0,S1_P0_1,S2_P0_2","P_end_A,P_end_B,P_end_C","D13,D17,D24,D15_A,D15_B,D15_C")

TYPE <- c("Parental","Untreated","Treated")

CUMPLOT <- c(paste0("cumulative_GBC_untreated_invitro_",EXP,".pdf"),
             paste0("cumulative_GBC_treated_invitro_",EXP,".pdf"))

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
    colnames(df_count) <- c("sample", "clone.id", "clone.id.frac", "GBC.count", "count", "type")
    
}

# generate figures

df_count$type <- factor(df_count$type, levels = TYPE)

idx <- which(df_count$type %in% c("Parental","Untreated"))
g <- ggplot(data = df_count[idx,], aes(x = clone.id, y = GBC.count, color = type)) + theme_classic() + geom_point(alpha = 0.5, size=1) 
g <- g + xlab("GBC sorted by CPM") + ylab("Cumulative CPM fraction")
g <- g + scale_x_continuous(trans = log10_trans(), limits = c(1,max_gbc_count))
g <- g + theme(legend.title=element_blank()) + scale_color_manual(values = colors) + theme(legend.position="top")

pdf(CUMPLOT[1], width = 2.5, height = 2.8)
print(g)
dev.off()

df_count$type <- factor(df_count$type, levels = TYPE)

idx <- which(df_count$type %in% c("Parental","Treated"))
g <- ggplot(data = df_count[idx,], aes(x = clone.id, y = GBC.count, color = type)) + theme_classic() + geom_point(alpha = 0.5, size=1) 
g <- g + xlab("GBC sorted by CPM") + ylab("Cumulative CPM fraction")
g <- g + scale_x_continuous(trans = log10_trans(), limits = c(1,max_gbc_count))
g <- g + theme(legend.title=element_blank()) + scale_color_manual(values = colors) + theme(legend.position="top")

pdf(CUMPLOT[2], width = 2.5, height = 2.8)
print(g)
dev.off()

q()


# select palette

colors <- c("gray80","black")

### EXP4

# directory and file names

TABLE <- "../tables/DATA_FRAME_GBC_IN_VITRO_EXP4.tsv"
EXP <- "EXP4"

SAMPLES <- c("P0","P_end_1,P_end_2,P_end_3","D15_A,D15_B,D15_C,D15_D,D15_E")

TYPE <- c("Parental","Untreated","Treated")

CUMPLOT <- c(paste0("cumulative_GBC_untreated_invitro_",EXP,".pdf"),
             paste0("cumulative_GBC_treated_invitro_",EXP,".pdf"))

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
    colnames(df_count) <- c("sample", "clone.id", "clone.id.frac", "GBC.count", "count", "type")
    
}

# generate figures

df_count$type <- factor(df_count$type, levels = TYPE)

idx <- which(df_count$type %in% c("Parental","Untreated"))
g <- ggplot(data = df_count[idx,], aes(x = clone.id, y = GBC.count, color = type)) + theme_classic() + geom_point(alpha = 0.5, size=1) 
g <- g + xlab("GBC sorted by CPM") + ylab("Cumulative CPM fraction")
g <- g + scale_x_continuous(trans = log10_trans(), limits = c(1,max_gbc_count))
g <- g + theme(legend.title=element_blank()) + scale_color_manual(values = colors) + theme(legend.position="top")

pdf(CUMPLOT[1], width = 2.5, height = 2.8)
print(g)
dev.off()

df_count$type <- factor(df_count$type, levels = TYPE)

idx <- which(df_count$type %in% c("Parental","Treated"))
g <- ggplot(data = df_count[idx,], aes(x = clone.id, y = GBC.count, color = type)) + theme_classic() + geom_point(alpha = 0.5, size=1) 
g <- g + xlab("GBC sorted by CPM") + ylab("Cumulative CPM fraction")
g <- g + scale_x_continuous(trans = log10_trans(), limits = c(1,max_gbc_count))
g <- g + theme(legend.title=element_blank()) + scale_color_manual(values = colors) + theme(legend.position="top")

pdf(CUMPLOT[2], width = 2.5, height = 2.8)
print(g)
dev.off()

q()


### EXP5

# directory and file names

TABLE <- "../tables/DATA_FRAME_GBC_IN_VITRO_EXP5.tsv"
EXP <- "EXP5"

SAMPLES <- c("P5OR_GBC_P0_1,P5OR_GBC_P0_2","P5OR_GBC_PACLI_D17,P5OR_GBC_PACLI_D20")

TYPE <- c("Parental","Treated")

CUMPLOT <- paste0("cumulative_GBC_invitro_",EXP,".pdf")

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
    colnames(df_count) <- c("sample", "clone.id", "clone.id.frac", "GBC.count", "count", "type")
    
}

# generate figures

df_count$type <- factor(df_count$type, levels = TYPE)

g <- ggplot(data = df_count, aes(x = clone.id, y = GBC.count, color = type)) + theme_classic() + geom_point(alpha = 0.5, size=1) 
g <- g + xlab("GBC sorted by CPM") + ylab("Cumulative CPM fraction")
g <- g + scale_x_continuous(trans = log10_trans(), limits = c(1,max_gbc_count))
g <- g + theme(legend.title=element_blank()) + scale_color_manual(values = colors) + theme(legend.position="top")

pdf(CUMPLOT, width = 2.5, height = 2.8)
print(g)
dev.off()

q()



