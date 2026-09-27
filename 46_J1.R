## 46_J1_similarity.R
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
seu <- readRDS("D:/BCBM_Project/data/raw/GSE324453_seurat.rds")
sc1 <- LayerData(seu, layer="counts"); ct1 <- as.character(seu$common_cell_clusters)
map1 <- c("Astrocytes"="Astrocyte","Neurons"="Neuron","ODG"="Oligodendrocyte",
  "Tumor cells"="Tumor","Tumor cells (prol.)"="Tumor","T/NK cells"="Tcell",
  "B cells"="Bcell","Plasma cells"="Plasma","Myeloid cells"="Myeloid",
  "Endothelial cells"="Endothelial","Pericytes"="Pericyte","Stromal cells"="Stroma")
lab1 <- unname(map1[ct1])
base <- "D:/BCBM_Project/data/raw/Wu_etal_2021_BRCA_scRNASeq/"
genes <- read.table(paste0(base,"count_matrix_genes.tsv"),stringsAsFactors=FALSE)$V1
bar <- read.table(paste0(base,"count_matrix_barcodes.tsv"),stringsAsFactors=FALSE)$V1
sc2 <- readMM(paste0(base,"count_matrix_sparse.mtx")); rownames(sc2)<-genes; colnames(sc2)<-bar
meta <- read.csv(paste0(base,"metadata.csv"), row.names=1)
map2 <- c("Cancer Epithelial"="Tumor","Normal Epithelial"="NormalEpi","CAFs"="CAF","PVL"="PVL",
  "Endothelial"="Endothelial","T-cells"="Tcell","B-cells"="Bcell","Plasmablasts"="Plasma","Myeloid"="Myeloid")
lab2 <- unname(map2[meta$celltype_major])
common <- intersect(rownames(sc1), rownames(sc2))
ref <- cbind(sc1[common,], sc2[common,]); labels <- c(lab1, lab2)
# subsample 500 per type for speed
set.seed(1)
keep <- unlist(lapply(split(seq_along(labels), labels), function(ii) sample(ii, min(length(ii),500))))
ref <- ref[,keep]; labels <- labels[keep]
cat("subsampled cells:",length(labels),"\n")
# average log expression per type
cpm <- function(M){t(t(M)/colSums(M))*1e6}
avg <- list()
for(typ in unique(labels)){
  idx <- which(labels==typ)
  avg[[typ]] <- rowMeans(log2(cpm(ref[,idx,drop=FALSE])+1))
}
A <- do.call(cbind, avg)
cat("cell types:",colnames(A),"\n")
cat("cell counts:\n"); print(table(labels))
# spearman correlation
cm <- cor(A, method="spearman")
print(round(cm,3))
# astrocyte similarity ranking
ast <- cm["Astrocyte",]
ast <- sort(ast, decreasing=TRUE)
cat("\n=== Astrocyte similarity to other types ===\n")
print(round(ast,3))
write.csv(cm, "D:/BCBM_Project/results/tables/J1_similarity_matrix.csv")
cat("saved\n")
