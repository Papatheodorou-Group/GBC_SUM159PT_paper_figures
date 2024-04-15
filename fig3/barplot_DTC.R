
# load libraries

library("ggplot2")
library("dplyr")
library("RColorBrewer")


col_bar <- c(brewer.pal(4,"Greys")[2],brewer.pal(4,"Oranges")[4])

### DTC

INFO <- "../tables/DATA_FRAME_GBC_IN_VITRO_EXP1.tsv"

SAMPLES_P0 <- c("D0","S1_P0_1","S2_P0_2")
SAMPLES_P0_RENAMED <- paste0("p",1:3)

SAMPLES_UNTREATED <- c("P_end_A","P_end_B","P_end_C")         
SAMPLES_UNTREATED_RENAMED <- paste0("ut",1:3)

SAMPLES_TREATED <- c("D13","D17","D24","D15_A","D15_B","D15_C")
SAMPLES_TREATED_RENAMED <- paste0("t",1:6)

SAMPLES <- c(SAMPLES_P0, SAMPLES_UNTREATED, SAMPLES_TREATED)
SAMPLES_RENAMED <- c(SAMPLES_P0_RENAMED, SAMPLES_UNTREATED_RENAMED, SAMPLES_TREATED_RENAMED)

TYPE <- c("DRC_bulk","DRC_bulk_fisher","DRC_bulk_gam","DRC_bulk_wilcox")

PLOT <- "barplot_DTC_multi.pdf"

PLOT_BOOL <- "barplot_DTC_GBC_multi.pdf"

# collect the abundance of each GBC type in samples

df_info <- read.table(INFO, sep = "\t", header = TRUE)
idx <- match(paste0("cpm.",SAMPLES), colnames(df_info))
colnames(df_info)[idx] <- SAMPLES_RENAMED

df <- as.data.frame(matrix(0, nrow = 0, ncol = 5))
colnames(df) <- c("type", "sample", "class", "count", "GBC.count")

# p0 + untreated + treated
for (t in TYPE) {

    for (s in SAMPLES_RENAMED) {
    
        sel <- sum(df_info[df_info[[t]] == 1,s]) 
        non_sel <- sum(df_info[df_info[[t]] == 0,s])
        
        sel_bool <- sum(df_info[df_info[[t]] == 1,s] > 0)
        non_sel_bool <- sum(df_info[df_info[[t]] == 0,s] > 0)
        
        df_tmp <- data.frame(type = rep(t,2), sample = rep(s,2), class = c(1,0), 
                             count = c(sel,non_sel), GBC.count = c(sel_bool,non_sel_bool))
        df <- rbind(df, df_tmp) 
    }   
}

df$sample <- factor(df$sample, levels = SAMPLES_RENAMED)
df$class <- factor(df$class, levels = c(0,1))

g <- ggplot(data = df, aes(x = sample, y = count, fill = class)) + theme_classic() + ylab("CPM fraction")
g <- g + geom_bar(stat = "identity", position = "fill") + scale_fill_manual(values = col_bar)
g <- g + facet_wrap(~ type, ncol = 4)

pdf(PLOT, width = 12, height = 2)
print(g)
dev.off()

g <- ggplot(data = df, aes(x = sample, y = GBC.count, fill = class)) + theme_classic() + ylab("GBC fraction")
g <- g + geom_bar(stat = "identity", position = "fill") + scale_fill_manual(values = col_bar)
g <- g + facet_wrap(~ type, ncol = 4)

pdf(PLOT_BOOL, width = 12, height = 2)
print(g)
dev.off()


### DTC time-course

INFO <- "../tables/DATA_FRAME_EXP1D_treated_filt.tsv"

SAMPLES <- c("T5_2_Multiseq","T7_2_Multiseq","T9_2_Multiseq","T11_2_Multiseq","T13_2_Multiseq","T15_2_Multiseq")
SAMPLES_RENAMED <- paste0("d",c(5,7,9,11,13,15))

PLOT <- "barplot_DTC_time_course.pdf"

# collect the abundance of each GBC type in samples

df <- read.table(INFO, sep = "\t", header = TRUE)
df$sample.name <- SAMPLES_RENAMED[match(df$sample.name,SAMPLES)]

df$DRC_bulk_fisher[is.na(df$DRC_bulk_fisher)] <- rep(0,sum(is.na(df$DRC_bulk_fisher)))

df$sample.name <- factor(df$sample.name, levels = SAMPLES_RENAMED)
df$DRC_bulk_fisher <- factor(df$DRC_bulk_fisher, levels = c(0,1))

g <- ggplot(data = df, aes(x = sample.name, fill = DRC_bulk_fisher)) + theme_classic() + ylab("cell fraction")
g <- g + geom_bar(stat = "count", position = "fill") + scale_fill_manual(values = col_bar)

pdf(PLOT, width = 3.5, height = 2)
print(g)
dev.off()

q()

