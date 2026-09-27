## 37_background_calibration.R - BayesPrism on normal-breast pseudobulk to estimate floor
suppressPackageStartupMessages({library(Seurat); library(BayesPrism); library(Matrix); library(openxlsx)})
set.seed(1)
# Build pseudobulk from GSE176078 normal epithelial + CAF + PVL + endothelial (no brain, no tumor)
base <- "D:/BCBM_Project/data/raw/Wu_etal_2021_BRCA_scRNASeq/"
genes <- read.table(paste0(base,"count_matrix_genes.tsv"),stringsAsFactors=FALSE)$V1
bar <- read.table(paste0(base,"count_matrix_barcodes.tsv"),stringsAsFactors=FALSE)$V1
sc <- readMM(paste0(base,"count_matrix_sparse.mtx")); rownames(sc)<-genes; colnames(sc)<-bar
meta <- read.csv(paste0(base,"metadata.csv"), row.names=1)
normal_types <- c("Normal Epithelial","CAFs","PVL","Endothelial")
normal_idx <- which(meta$celltype_major %in% normal_types)
cat("normal cells:",length(normal_idx),"\n")
# 45 pseudobulk samples, 3000 cells each
N<-3000; nsamp<-45
pseudo <- matrix(0,nrow=nrow(sc),ncol=nsamp); rownames(pseudo)<-genes
for(k in 1:nsamp){
  idx <- sample(normal_idx, N, replace=TRUE)
  pseudo[,k] <- rowSums(sc[,idx,drop=FALSE])
}
colnames(pseudo)<-paste0("NORM",1:nsamp)
cat("pseudobulk built\n")

# Reference: same combined as module C
seu <- readRDS("D:/BCBM_Project/data/raw/GSE324453_seurat.rds")
sc1 <- LayerData(seu, layer="counts"); ct1 <- as.character(seu$common_cell_clusters)
map1 <- c("Astrocytes"="Astrocyte","Neurons"="Neuron","ODG"="Oligodendrocyte",
  "Tumor cells"="Tumor","Tumor cells (prol.)"="Tumor","T/NK cells"="Tcell",
  "B cells"="Bcell","Plasma cells"="Plasma","Myeloid cells"="Myeloid",
  "Endothelial cells"="Endothelial","Pericytes"="Pericyte","Stromal cells"="Stroma")
lab1 <- unname(map1[ct1])
map2 <- c("Cancer Epithelial"="Tumor","Normal Epithelial"="NormalEpi","CAFs"="CAF","PVL"="PVL",
  "Endothelial"="Endothelial","T-cells"="Tcell","B-cells"="Bcell","Plasmablasts"="Plasma","Myeloid"="Myeloid")
lab2 <- unname(map2[meta$celltype_major])
common <- intersect(rownames(sc1), rownames(sc))
ref <- cbind(sc1[common,], sc[common,]); labels <- c(lab1, lab2)
keep <- unlist(lapply(split(seq_along(labels), labels), function(ii) sample(ii, min(length(ii),250))))
ref <- ref[,keep]; labels <- labels[keep]

# map pseudobulk genes
pm <- read.csv("D:/BCBM_Project/results/tables/program_decomposition.csv")
# pseudobulk already gene-symbol rownames
common2 <- intersect(rownames(ref), rownames(pseudo))
ref2 <- t(as.matrix(ref[common2,]))
mix2 <- t(as.matrix(pseudo[common2,]))
cat("ref:",dim(ref2)," mix:",dim(mix2),"\n")
prism <- new.prism(reference=ref2, mixture=mix2, input.type="counts",
  cell.type.labels=labels, cell.state.labels=labels, key="Tumor")
prism <- run.prism(prism, n.cores=1, gibbs.control=list(chain.length=150, burn.in=50, thinning=2))
frac <- get.fraction(prism, which.theta="final", state.or.type="type")
write.csv(frac, "D:/BCBM_Project/results/tables/bp_background_floor.csv")
brain <- frac[,"Neuron"]+frac[,"Astrocyte"]+frac[,"Oligodendrocyte"]
cat("\n=== BACKGROUND FLOOR (normal breast pseudobulk, n=45) ===\n")
cat("brain fraction median:",median(brain)," IQR:",IQR(brain),
    " q95:",quantile(brain,0.95),"\n")
cat("tumor fraction median:",median(frac[,"Tumor"]),"\n")
