## 16_cellchat.R - ligand-receptor: brain microenvironment -> tumor cells
suppressPackageStartupMessages({library(Seurat); library(CellChat); library(dplyr)})

seu <- readRDS("D:/BCBM_Project/data/raw/GSE324453_seurat.rds")
seu <- NormalizeData(seu, verbose=FALSE)

# downsample for memory
set.seed(1)
n <- ncol(seu)
keep <- sample(colnames(seu), min(25000, n))
seu <- subset(seu, cells=keep)
cat("working cells:", ncol(seu), "\n")

# use existing cell type labels as identities
Idents(seu) <- factor(seu$common_cell_clusters)
levels(Idents(seu))

# Create CellChat
cellchat <- createCellChat(object=seu, group.by="common_cell_clusters", assay="RNA")
cellchat <- addMeta(cellchat, meta=seu@meta.data)
cellchat <- setIdent(cellchat, ident.use="common_cell_clusters")
cat("cell groups:", levels(cellchat@idents), "\n")

CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat, do.fast=FALSE)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells=5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# signaling TO tumor cells
df <- subsetCommunication(cellchat)
to_tumor <- df[grepl("Tumor", df$target), c("source","target","pathway_name","interaction_name","prob","pval")]
to_tumor <- to_tumor[order(-to_tumor$prob), ]
cat("\n=== top signaling TO tumor cells ===\n")
print(head(to_tumor, 30), row.names=FALSE)
write.csv(to_tumor, "D:/BCBM_Project/results/tables/cellchat_to_tumor.csv", row.names=FALSE)

# specifically astrocyte/microglia(myeloid)/neuron -> tumor
focus <- to_tumor[to_tumor$source %in% c("Astrocytes","Myeloid cells","Neurons","ODG","Endothelial cells"), ]
cat("\n=== microenvironment -> tumor (top 25) ===\n")
print(head(focus,25), row.names=FALSE)
write.csv(focus, "D:/BCBM_Project/results/tables/cellchat_microenv_to_tumor.csv", row.names=FALSE)
cat("\nsaved.\n")
