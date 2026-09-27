## 18_tf_trajectory.R
## (a) TF activity: score key neural vs EMT/EMT TFs in tumor cells
## (b) Slingshot trajectory on tumor cells, ordered by adaptation program.
suppressPackageStartupMessages({library(Seurat); library(dplyr)})

seu <- readRDS("D:/BCBM_Project/data/raw/GSE324453_seurat.rds")
seu <- NormalizeData(seu, verbose=FALSE)

adapt <- readLines("D:/BCBM_Project/results/tables/geneset_adaptation.txt")
comp  <- readLines("D:/BCBM_Project/results/tables/geneset_competence.txt")
M <- LayerData(seu, layer="data")
colmean <- function(genes){ g <- intersect(genes, rownames(M)); setNames(colMeans(M[g,,drop=FALSE]), colnames(M)) }
seu$adapt <- colmean(adapt); seu$comp <- colmean(comp)

# ---- TF activity (target-based regulon scoring; curated high-confidence targets) ----
regulons <- list(
  PAX6   = c("PAX6","TUBB3","NEUROD1","SOX2","STMN2","GAP43","NCAM1","NRCAM"),
  ZIC2   = c("ZIC2","SOX2","PAX6","NES","SOX1","HESX1"),
  SOX2   = c("SOX2","NES","NES","PAX6","ZIC2","SOX1","SOX21","FOXG1"),
  ASCL1  = c("ASCL1","NEUROG1","NEUROD1","DELTA1","HES6","STMN2"),
  NEUROD1= c("NEUROD1","TUBB3","SYT1","SNAP25","GAP43","CHGA","SCG2"),
  MYC    = c("MYC","CCND1","NCL","NPM1","RPL23A","ODC1","MCM4"),
  E2F1   = c("E2F1","MKI67","CCNB1","CCNA2","MCM5","PCNA","BIRC5"),
  ZEB1   = c("ZEB1","VIM","CDH2","FN1","SNAI2","ZEB2","ITGB1"),
  SNAI2  = c("SNAI2","VIM","CDH2","FN1","MMP2","ZEB1"),
  TWIST1 = c("TWIST1","VIM","CDH2","AKT2","TNIK")
)
for (nm in names(regulons)) {
  seu[[paste0("TF_",nm)]] <- colmean(regulons[[nm]])
}

# tumor cells only
tum <- subset(seu, subset = grepl("Tumor", common_cell_clusters))
cat("tumor cells:", ncol(tum), "\n")
agg <- data.frame(
  adapt = tum$adapt,
  comp = tum$comp,
  t(scale(t(as.matrix(tum@meta.data[, paste0("TF_",names(regulons))]))))
)
# correlation of TF scores with adaptation
cors <- sapply(paste0("TF_",names(regulons)), function(v) cor(tum$adapt, tum@meta.data[[v]], use="complete"))
cat("\n=== correlation of TF activity with adaptation score (tumor cells) ===\n")
print(sort(cors, decreasing=TRUE))

# ---- Slingshot trajectory on a downsampled tumor UMAP ----
set.seed(1)
ts <- subset(tum, cells = sample(colnames(tum), min(8000, ncol(tum))))
ts <- FindVariableFeatures(ts, verbose=FALSE)
ts <- ScaleData(ts, verbose=FALSE)
ts <- RunPCA(ts, npcs=15, verbose=FALSE)
ts <- RunUMAP(ts, dims=1:15, verbose=FALSE)
# clusters
ts <- FindNeighbors(ts, dims=1:15, verbose=FALSE)
ts <- FindClusters(ts, resolution=0.5, verbose=FALSE)

suppressPackageStartupMessages(library(slingshot))
Idents(ts) <- factor(ts$seurat_clusters)
sce <- as.SingleCellExperiment(ts)
colLabels(sce) <- Idents(ts)
sce <- slingshot(sce, clusterLabels = colLabels(sce), reducedDim = "UMAP",
                 start.clus = NULL)
ts$pseudotime <- slingPseudotime(sce)[,1]
cat("\npseudotime range:", range(ts$pseudotime, na.rm=TRUE), "\n")

# correlation of pseudotime with program scores
cat("cor(pseudotime, adapt):", round(cor(ts$pseudotime, ts$adapt, use="complete"),3), "\n")
cat("cor(pseudotime, comp):", round(cor(ts$pseudotime, ts$comp, use="complete"),3), "\n")

# TF activity along pseudotime bins
ts$pt_bin <- cut(ts$pseudotime, breaks=quantile(ts$pseudotime, probs=0:5/5, na.rm=TRUE), include.lowest=TRUE)
pt_tf <- ts@meta.data %>% group_by(pt_bin) %>%
  summarise(across(starts_with("TF_"), mean, na.rm=TRUE))
cat("\n=== TF activity along pseudotime ===\n"); print(as.data.frame(pt_tf), digits=2)

write.csv(data.frame(barcode=colnames(tum), adapt=tum$adapt, comp=tum$comp),
          "D:/BCBM_Project/results/tables/tumor_adapt_scores.csv", row.names=FALSE)
write.csv(sort(cors, decreasing=TRUE), "D:/BCBM_Project/results/tables/tf_cor_adaptation.csv")
write.csv(pt_tf, "D:/BCBM_Project/results/tables/tf_along_pseudotime.csv", row.names=FALSE)
cat("\nsaved.\n")
