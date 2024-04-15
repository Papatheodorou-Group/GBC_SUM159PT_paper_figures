
# load libraries

library("Seurat")
library("Signac")
library("pROC")

DIR <- "../analysis/multiome"



### EXP1C

EXP_NAME <- "EXP1C"

DF <- file.path("../tables", paste0("DATA_FRAME_",EXP_NAME,".tsv"))

AUC <- "module_cell_AUC_1C.tsv"

META <- c("expr.GBC.num","expr.GBC.list","DRC_bulk","DRC_bulk_wilcox","DRC_bulk_fisher","DRC_bulk_gam",
          "TIC","TIC_wilcox","TIC_fisher","TIC_gam",
          "treated_TIC","treated_TIC_wilcox","treated_TIC_fisher","treated_TIC_gam","DRC_sc") 
          
CLUSTER_S1 <- "1"
CLUSTER_S3 <- "3"

ROC_MODULE20 <- "ROC_EXP1C_Module20_DRC_S3.pdf"
ROC_MODULE1 <- "ROC_EXP1C_Module1_TIC_S1_S3.pdf"

df_info <- read.table(DF, sep = "\t", header = TRUE)
auc_info <- read.table(AUC, sep = "\t", header = TRUE)
module20 <- as.numeric(auc_info["Module 20",])
module1 <- as.numeric(auc_info["Module 1",])

auc_DRC <- auc(response = df_info$DRC_bulk_fisher, predictor = module20, direction = "<")
auc_S3 <- auc(response = (df_info$cl_renamed == "3"), predictor = module20, direction = "<")

roc_DRC <- roc(response = df_info$DRC_bulk_fisher, predictor = module20, direction = "<")
roc_S3 <- roc(response = (df_info$cl_renamed == "3"), predictor = module20, direction = "<")

pdf(ROC_MODULE20, width = 5, height = 5)
plot(roc_DRC, col = "gray50", main = "Module 20")
lines(roc_S3, col = "black")
text(paste("DTC",round(auc_DRC,2)), x = 0.15, y = 0.2, cex = 1.5, col = "gray50")
text(paste("S3",round(auc_S3,2)), x = 0.15, y = 0.1, cex = 1.5, col = "black")
dev.off()

auc_TIC <- auc(response = df_info$TIC_fisher, predictor = module1, direction = "<")
auc_S1 <- auc(response = (df_info$cl_renamed == "1"), predictor = module1, direction = "<")
auc_S3 <- auc(response = (df_info$cl_renamed == "3"), predictor = module1, direction = "<")

roc_TIC <- roc(response = df_info$TIC_fisher, predictor = module1, direction = "<")
roc_S1 <- roc(response = (df_info$cl_renamed == "1"), predictor = module1, direction = "<")
roc_S3 <- roc(response = (df_info$cl_renamed == "3"), predictor = module1, direction = "<")

pdf(ROC_MODULE1, width = 5, height = 5)
plot(roc_TIC, col = "gray50", main = "Module 1")
lines(roc_S1, col = "red")
lines(roc_S3, col = "black")
text(paste("TIC",round(auc_TIC,2)), x = 0.15, y = 0.3, cex = 1.5, col = "gray50")
text(paste("S1",round(auc_S1,2)), x = 0.15, y = 0.2, cex = 1.5, col = "red")
text(paste("S3",round(auc_S3,2)), x = 0.15, y = 0.1, cex = 1.5, col = "black")
dev.off()



EXP <- "1C_ARC"

OBJECT <- file.path(DIR, EXP, paste0(EXP,"_seurat.Rds"))
DF <- file.path("../tables", paste0("DATA_FRAME_",EXP_NAME,".tsv"))
          
TOPIC <- "Topic_15"
CLUSTER <- "3"

ROC <- "ROC_EXP1C_DRC_S3.pdf"

object <- readRDS(OBJECT)
df_info <- read.table(DF, sep = "\t", header = TRUE)

idx <- match(df_info$cell_ID, colnames(object))
topic <- object@reductions$cisTopic@cell.embeddings[,TOPIC]
topic <- topic[idx]

auc_DRC <- auc(response = df_info$DRC_bulk_fisher, predictor = topic, direction = "<")
auc_S3 <- auc(response = (df_info$cl_renamed == CLUSTER), predictor = topic, direction = "<")

roc_DRC <- roc(response = df_info$DRC_bulk_fisher, predictor = topic, direction = "<")
roc_S3 <- roc(response = (df_info$cl_renamed == CLUSTER), predictor = topic, direction = "<")

pdf(ROC, width = 5, height = 5)
plot(roc_DRC)
lines(roc_S3, col = "red")
text(round(auc_DRC,2), x = 0.2, y = 0.3, cex = 1.5)
text(round(auc_S3,2), x = 0.2, y = 0.2, cex = 1.5)
dev.off()



### EXP1E

EXP_NAME <- "EXP1E"

DF <- file.path("../tables", paste0("DATA_FRAME_",EXP_NAME,".tsv"))

AUC <- "module_cell_AUC_1E.tsv"

META <- c("expr.GBC.num","expr.GBC.list","DRC_bulk","DRC_bulk_wilcox","DRC_bulk_fisher","DRC_bulk_gam",
          "TIC","TIC_wilcox","TIC_fisher","TIC_gam",
          "treated_TIC","treated_TIC_wilcox","treated_TIC_fisher","treated_TIC_gam","DRC_sc") 

CLUSTER_S1 <- "1"
CLUSTER_S3 <- "3"

ROC_MODULE20 <- "ROC_EXP1E_Module20_DRC_S3.pdf"
ROC_MODULE1 <- "ROC_EXP1E_Module1_TIC_S1_S3.pdf"

df_info <- read.table(DF, sep = "\t", header = TRUE)
auc_info <- read.table(AUC, sep = "\t", header = TRUE)
module20 <- as.numeric(auc_info["Module 20",])
module1 <- as.numeric(auc_info["Module 1",])

auc_DRC <- auc(response = df_info$DRC_bulk_fisher, predictor = module20, direction = "<")
auc_S3 <- auc(response = (df_info$cl_renamed == "3"), predictor = module20, direction = "<")

roc_DRC <- roc(response = df_info$DRC_bulk_fisher, predictor = module20, direction = "<")
roc_S3 <- roc(response = (df_info$cl_renamed == "3"), predictor = module20, direction = "<")

pdf(ROC_MODULE20, width = 5, height = 5)
plot(roc_DRC, col = "gray50", main = "Module 20")
lines(roc_S3, col = "black")
text(paste("DTC",round(auc_DRC,2)), x = 0.15, y = 0.2, cex = 1.5, col = "gray50")
text(paste("S3",round(auc_S3,2)), x = 0.15, y = 0.1, cex = 1.5, col = "black")
dev.off()

auc_TIC <- auc(response = df_info$TIC_fisher, predictor = module1, direction = "<")
auc_S1 <- auc(response = (df_info$cl_renamed == "1"), predictor = module1, direction = "<")
auc_S3 <- auc(response = (df_info$cl_renamed == "3"), predictor = module1, direction = "<")

roc_TIC <- roc(response = df_info$TIC_fisher, predictor = module1, direction = "<")
roc_S1 <- roc(response = (df_info$cl_renamed == "1"), predictor = module1, direction = "<")
roc_S3 <- roc(response = (df_info$cl_renamed == "3"), predictor = module1, direction = "<")

pdf(ROC_MODULE1, width = 5, height = 5)
plot(roc_TIC, col = "gray50", main = "Module 1")
lines(roc_S1, col = "red")
lines(roc_S3, col = "black")
text(paste("TIC",round(auc_TIC,2)), x = 0.15, y = 0.3, cex = 1.5, col = "gray50")
text(paste("S1",round(auc_S1,2)), x = 0.15, y = 0.2, cex = 1.5, col = "red")
text(paste("S3",round(auc_S3,2)), x = 0.15, y = 0.1, cex = 1.5, col = "black")
dev.off()



EXP <- "1E_ARC"

OBJECT <- file.path(DIR, EXP, paste0(EXP,"_seurat.Rds"))
DF <- file.path("../tables", paste0("DATA_FRAME_",EXP_NAME,".tsv"))
          
TOPIC <- "Topic_19"
CLUSTER <- "2"

ROC <- "ROC_EXP1E_DRC_S3.pdf"

library("Seurat")
library("Signac")
library("pROC")

object <- readRDS(OBJECT)
df_info <- read.table(DF, sep = "\t", header = TRUE)

idx <- match(df_info$cell_ID, colnames(object))
topic <- object@reductions$cisTopic@cell.embeddings[,TOPIC]
topic <- topic[idx]

auc_DRC <- auc(response = df_info$DRC_bulk_fisher, predictor = topic, direction = "<")
auc_S3 <- auc(response = (df_info$cl_renamed == CLUSTER), predictor = topic, direction = "<")

roc_DRC <- roc(response = df_info$DRC_bulk_fisher, predictor = topic, direction = "<")
roc_S3 <- roc(response = (df_info$cl_renamed == CLUSTER), predictor = topic, direction = "<")

pdf(ROC, width = 5, height = 5)
plot(roc_DRC)
lines(roc_S3, col = "red")
text(round(auc_DRC,2), x = 0.2, y = 0.3, cex = 1.5)
text(round(auc_S3,2), x = 0.2, y = 0.2, cex = 1.5)
dev.off()

q()


