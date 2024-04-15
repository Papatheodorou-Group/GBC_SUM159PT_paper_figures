
# load libraries

library("ggplot2")
library("dplyr")
library("RColorBrewer")

### tumour by TICs

INFO <- "../tables/DATA_FRAME_GBC_IN_VIVO.tsv"

SAMPLES_P0 <- c("P0_R4","P0_PERT2","P0_PERT3")
SAMPLES_P0_RENAMED <- paste0("p",1:3)

SAMPLES_TUMORS <- c("A3PT_PERT3","C3M1_PERT2","C3PT_PERT1","D1PT_PERT3","D4PT_PERT3","D5PT_PER3","C1PT_PTX","C2PT_PTX","D3PT_NT")         
SAMPLES_TUMORS_RENAMED <- paste0("t",1:9)

SAMPLES_TUMORS_TREATED <- c("A2PT_R4","A1PT_PERT1","A5V2PT_PERT1","B1PT_PERT3","B2PT_PTX_Per3","B3PT_PTX_Per3")
SAMPLES_TUMORS_TREATED_RENAMED <- paste0("tt",1:6)

SAMPLES <- c(SAMPLES_P0, SAMPLES_TUMORS, SAMPLES_TUMORS_TREATED)
SAMPLES_RENAMED <- c(SAMPLES_P0_RENAMED, SAMPLES_TUMORS_RENAMED, SAMPLES_TUMORS_TREATED_RENAMED)

TYPE_UNTR <- c("TIC","TIC_fisher","TIC_gam","TIC_wilcox")
TYPE_TR <- c("treated_TIC","treated_TIC_fisher","treated_TIC_gam","treated_TIC_wilcox")
TYPE <- c(TYPE_UNTR,TYPE_TR)

CL_1B <- "cluster_1B_P0_filt_renamed"
CL_1D <- "cluster_1D_P0_filt_renamed"

PLOT_UNTR <- "barplot_TIC_untreated_multi.pdf"
PLOT_TR <- "barplot_TIC_treated_multi.pdf"

PLOT_UNTR_BOOL <- "barplot_TIC_untreated_GBC_multi.pdf"
PLOT_TR_BOOL <- "barplot_TIC_treated_GBC_multi.pdf"

col_bar <- c(brewer.pal(4,"Greys")[2],brewer.pal(4,"Blues")[4])

# collect the abundance of each GBC type in samples

df_info <- read.table(INFO, sep = "\t", header = TRUE)
idx <- match(paste0("cpm.",SAMPLES), colnames(df_info))
colnames(df_info)[idx] <- SAMPLES_RENAMED

df <- as.data.frame(matrix(0, nrow = 0, ncol = 5))
colnames(df) <- c("type", "sample", "class", "count", "GBC.count")

# p0 + untreated tumors
for (t in TYPE_UNTR) {

    for (s in c(SAMPLES_P0_RENAMED, SAMPLES_TUMORS_RENAMED)) {
    
        sel <- sum(df_info[df_info[[t]] == 1,s]) 
        non_sel <- sum(df_info[df_info[[t]] == 0,s])
        
        sel_bool <- sum(df_info[df_info[[t]] == 1,s] > 0)
        non_sel_bool <- sum(df_info[df_info[[t]] == 0,s] > 0)
        
        df_tmp <- data.frame(type = rep(t,2), sample = rep(s,2), class = c(1,0), 
                             count = c(sel,non_sel), GBC.count = c(sel_bool,non_sel_bool))
        df <- rbind(df, df_tmp) 
    }   
}

# p0 + treated tumor
for (t in TYPE_TR) {

    for (s in c(SAMPLES_P0_RENAMED, SAMPLES_TUMORS_TREATED_RENAMED)) {
    
        sel <- sum(df_info[df_info[[t]] == 1,s]) 
        non_sel <- sum(df_info[df_info[[t]] == 0,s])
        
        sel_bool <- sum(df_info[df_info[[t]] == 1,s] > 0)
        non_sel_bool <- sum(df_info[df_info[[t]] == 0,s] > 0)
        
        df_tmp <- data.frame(type = rep(t,2), sample = rep(s,2), class = c(1,0), 
                             count = c(sel,non_sel), GBC.count = c(sel_bool,non_sel_bool))
        df <- rbind(df, df_tmp) 
    }   
}

df$class <- factor(df$class, levels = c(0,1))

g <- ggplot(data = df[df$type %in% TYPE_UNTR,], aes(x = sample, y = count, fill = class)) + theme_classic() + ylab("CPM fraction")
g <- g + geom_bar(stat = "identity", position = "fill") + scale_fill_manual(values = col_bar)
g <- g + facet_wrap(~ type, ncol = 4)

pdf(PLOT_UNTR, width = 10, height = 2)
print(g)
dev.off()

g <- ggplot(data = df[df$type %in% TYPE_TR,], aes(x = sample, y = count, fill = class)) + theme_classic() + ylab("CPM fraction")
g <- g + geom_bar(stat = "identity", position = "fill") + scale_fill_manual(values = col_bar)
g <- g + facet_wrap(~ type, ncol = 4)

pdf(PLOT_TR, width = 10, height = 2)
print(g)
dev.off()

g <- ggplot(data = df[df$type %in% TYPE_UNTR,], aes(x = sample, y = GBC.count, fill = class)) + theme_classic() + ylab("GBC fraction")
g <- g + geom_bar(stat = "identity", position = "fill") + scale_fill_manual(values = col_bar)
g <- g + facet_wrap(~ type, ncol = 4)

pdf(PLOT_UNTR_BOOL, width = 10, height = 2)
print(g)
dev.off()

g <- ggplot(data = df[df$type %in% TYPE_TR,], aes(x = sample, y = GBC.count, fill = class)) + theme_classic() + ylab("GBC fraction")
g <- g + geom_bar(stat = "identity", position = "fill") + scale_fill_manual(values = col_bar)
g <- g + facet_wrap(~ type, ncol = 4)

pdf(PLOT_TR_BOOL, width = 10, height = 2)
print(g)
dev.off()


### tumour by clusters

# files

INFO <- "../tables/DATA_FRAME_GBC_IN_VIVO.tsv"

SAMPLES_P0 <- c("P0_R4","P0_PERT2","P0_PERT3")
SAMPLES_P0_RENAMED <- paste0("p",1:3)

SAMPLES_TUMORS <- c("A3PT_PERT3","C3M1_PERT2","C3PT_PERT1","D1PT_PERT3","D4PT_PERT3","D5PT_PER3","C1PT_PTX","C2PT_PTX","D3PT_NT")         
SAMPLES_TUMORS_RENAMED <- paste0("t",1:9)

SAMPLES_TUMORS_TREATED <- c("A2PT_R4","A1PT_PERT1","A5V2PT_PERT1","B1PT_PERT3","B2PT_PTX_Per3","B3PT_PTX_Per3")
SAMPLES_TUMORS_TREATED_RENAMED <- paste0("tt",1:6)

SAMPLES_UNTREATED <- c(SAMPLES_P0, SAMPLES_TUMORS)
SAMPLES_UNTREATED_RENAMED <- c(SAMPLES_P0_RENAMED, SAMPLES_TUMORS_RENAMED)

SAMPLES_TREATED <- c(SAMPLES_P0, SAMPLES_TUMORS_TREATED)
SAMPLES_TREATED_RENAMED <- c(SAMPLES_P0_RENAMED, SAMPLES_TUMORS_TREATED_RENAMED)

SAMPLES <- c(SAMPLES_P0, SAMPLES_TUMORS, SAMPLES_TUMORS_TREATED)
SAMPLES_RENAMED <- c(SAMPLES_P0_RENAMED, SAMPLES_TUMORS_RENAMED, SAMPLES_TUMORS_TREATED_RENAMED)

TYPE_UNTR <- c("TIC","TIC_fisher","TIC_gam","TIC_wilcox")
TYPE_TR <- c("treated_TIC","treated_TIC_fisher","treated_TIC_gam","treated_TIC_wilcox")
TYPE <- c(TYPE_UNTR,TYPE_TR)

CL_1B <- "cluster_1B_P0_filt_renamed"
CL_1D <- "cluster_1D_P0_filt_renamed"

DF <- "barplot_TIC_multi_cl_NEW.tsv"

PLOT_UNTR_1B <- "barplot_TIC_untreated_multi_cl_1B_NEW.pdf"
PLOT_UNTR_1D <- "barplot_TIC_untreated_multi_cl_1D_NEW.pdf"
PLOT_TR_1B <- "barplot_TIC_treated_multi_cl_1B_NEW.pdf"
PLOT_TR_1D <- "barplot_TIC_treated_multi_cl_1D_NEW.pdf"

# collect the normalized abundance of each GBC type in samples

df_info <- read.table(INFO, sep = "\t", header = TRUE)
idx <- match(paste0("cpm.",SAMPLES), colnames(df_info))
colnames(df_info)[idx] <- SAMPLES_RENAMED

# collect the abundance of each cluster in TIC in samples

idx <- which(df_info[[CL_1B]]==-1)
df_info[idx,CL_1B] <- rep("undef",length(idx))
idx <- which(df_info[[CL_1D]]==-1)
df_info[idx,CL_1D] <- rep("undef",length(idx))

cl_1B <- sort(unique(df_info[[CL_1B]]))
cl_1D <- sort(unique(df_info[[CL_1D]]))
ncl_1B <- length(cl_1B)
ncl_1D <- length(cl_1D)

df <- as.data.frame(matrix(0, nrow = 0, ncol = 5))
colnames(df) <- c("type", "sample", "class", "sc_sample", "count")

# untreated tumors
for (t in TYPE_UNTR) {

    for (s in SAMPLES_UNTREATED_RENAMED) {
    
        sel_1B <- sapply(cl_1B, function(x) sum(df_info[df_info[[t]] == 1 & df_info[[CL_1B]] == x,s]))
        sel_1D <- sapply(cl_1D, function(x) sum(df_info[df_info[[t]] == 1 & df_info[[CL_1D]] == x,s]))
        
        df_tmp <- data.frame(type = rep(t,ncl_1B+ncl_1D), sample = rep(s,ncl_1B+ncl_1D), class = c(cl_1B,cl_1D), 
                             sc_sample = c(rep("1B",ncl_1B),rep("1D",ncl_1D)), count = c(sel_1B,sel_1D))
        df <- rbind(df, df_tmp) 
    }   
}
# NEW: collect all clones, without TIC selection
for (s in SAMPLES_UNTREATED_RENAMED) {
    
    sel_1B <- sapply(cl_1B, function(x) sum(df_info[df_info[[CL_1B]] == x,s]))
    sel_1D <- sapply(cl_1D, function(x) sum(df_info[df_info[[CL_1D]] == x,s]))

    df_tmp <- data.frame(type = rep("all_untreated",ncl_1B+ncl_1D), sample = rep(s,ncl_1B+ncl_1D), class = c(cl_1B,cl_1D), 
                     sc_sample = c(rep("1B",ncl_1B),rep("1D",ncl_1D)), count = c(sel_1B,sel_1D))
    df <- rbind(df, df_tmp) 
} 

# treated tumors
for (t in TYPE_TR) {

    for (s in SAMPLES_TREATED_RENAMED) {
    
        sel_1B <- sapply(cl_1B, function(x) sum(df_info[df_info[[t]] == 1 & df_info[[CL_1B]] == x,s])) 
        sel_1D <- sapply(cl_1D, function(x) sum(df_info[df_info[[t]] == 1 & df_info[[CL_1D]] == x,s]))
        
        df_tmp <- data.frame(type = rep(t,ncl_1B+ncl_1D), sample = rep(s,ncl_1B+ncl_1D), class = c(cl_1B,cl_1D), 
                             sc_sample = c(rep("1B",ncl_1B),rep("1D",ncl_1D)), count = c(sel_1B,sel_1D))
        df <- rbind(df, df_tmp) 
    }   
}
for (s in SAMPLES_TREATED_RENAMED) {

    sel_1B <- sapply(cl_1B, function(x) sum(df_info[df_info[[CL_1B]] == x,s])) 
    sel_1D <- sapply(cl_1D, function(x) sum(df_info[df_info[[CL_1D]] == x,s]))

    df_tmp <- data.frame(type = rep("all_treated",ncl_1B+ncl_1D), sample = rep(s,ncl_1B+ncl_1D), class = c(cl_1B,cl_1D), 
                     sc_sample = c(rep("1B",ncl_1B),rep("1D",ncl_1D)), count = c(sel_1B,sel_1D))
    df <- rbind(df, df_tmp) 
}

classes <- unique(c(cl_1B,cl_1D))
df$class <- factor(df$class, levels = rev(classes))

# sort samples by decreasing number of detected GBCs
idx <- order(colSums(df_info[,SAMPLES_RENAMED] > 0), decreasing = TRUE)
ordered_samples <- SAMPLES_RENAMED[idx]
df$sample <- factor(df$sample, levels = ordered_samples)

write.table(df, file = DF, quote = FALSE, sep = "\t", row.names = FALSE)

col_1B <- rev(c(brewer.pal(n = ncl_1B-1, name = "Set1"),"gray90"))
names(col_1B) <- rev(cl_1B)
col_1D <- rev(c(brewer.pal(n = ncl_1D-1, name = "Set2"),"gray90"))
names(col_1D) <- rev(cl_1D)

g <- ggplot(data = df[df$type %in% c(TYPE_UNTR,"all_untreated") & df$sc_sample == "1B",], aes(x = sample, y = count, fill = class)) 
g <- g + theme_classic() + ylab("CPM fraction")
g <- g + geom_bar(stat = "identity", position = "fill") + scale_fill_manual(values = col_1B)
g <- g + facet_wrap(~ type, ncol = 5)

pdf(PLOT_UNTR_1B, width = 10, height = 2.5)
print(g)
dev.off()

g <- ggplot(data = df[df$type %in% c(TYPE_UNTR,"all_untreated") & df$sc_sample == "1D",], aes(x = sample, y = count, fill = class)) 
g <- g + theme_classic() + ylab("CPM fraction")
g <- g + geom_bar(stat = "identity", position = "fill") + scale_fill_manual(values = col_1D)
g <- g + facet_wrap(~ type, ncol = 5)

pdf(PLOT_UNTR_1D, width = 10, height = 2.5)
print(g)
dev.off()

g <- ggplot(data = df[df$type %in% c(TYPE_TR,"all_treated") & df$sc_sample == "1B",], aes(x = sample, y = count, fill = class)) 
g <- g + theme_classic() + ylab("CPM fraction")
g <- g + geom_bar(stat = "identity", position = "fill") + scale_fill_manual(values = col_1B)
g <- g + facet_wrap(~ type, ncol = 5)

pdf(PLOT_TR_1B, width = 8, height = 2.5)
print(g)
dev.off()

g <- ggplot(data = df[df$type %in% c(TYPE_TR,"all_treated") & df$sc_sample == "1D",], aes(x = sample, y = count, fill = class)) 
g <- g + theme_classic() + ylab("CPM fraction")
g <- g + geom_bar(stat = "identity", position = "fill") + scale_fill_manual(values = col_1D)
g <- g + facet_wrap(~ type, ncol = 5)

pdf(PLOT_TR_1D, width = 8, height = 2.5)
print(g)
dev.off()


q()

