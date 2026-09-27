## 33_bp_combined.R
suppressPackageStartupMessages({library(Seurat); library(BayesPrism); library(Matrix); library(openxlsx)})
set.seed(1)
seu <- readRDS("D:/BCBM_Project/data/raw/GSE324453_seurat.rds")
sc1 <- LayerData(seu, layer="counts")
ct1 <- as.character(seu$common_cell_clusters)
map1 <- c("Astrocytes"="Astrocyte","Neurons"="Neuron","ODG"="Oligodendrocyte",
  "Tumor cells"="Tumor","Tumor cells (prol.)"="Tumor","T/NK cells"="Tcell",
  "B cells"="Bcell","Plasma cells"="Plasma","Myeloid cells"="Myeloid",
  "Endothelial cells"="Endothelial","Pericytes"="Pericyte","Stromal cells"="Stroma")
lab1 <- unname(map1[ct1])
base <- "D:/BCBM_Project/data/raw/Wu_etal_2021_BRCA_scRNASeq/"
genes <- read.table(paste0(base,"count_matrix_genes.tsv"),stringsAsFactors=FALSE)$V1
bar <- read.table(paste0(base,"count_matrix_barcodes.tsv"),stringsAsFactors=FALSE)$V1
sc2 <- readMM(paste0(base,"count_matrix_sparse.mtx"))
rownames(sc2) <- genes; colnames(sc2) <- bar
meta <- read.csv(paste0(base,"metadata.csv"), row.names=1)
map2 <- c("Cancer Epithelial"="Tumor","Normal Epithelial"="NormalEpi","CAFs"="CAF","PVL"="PVL",
  "Endothelial"="Endothelial","T-cells"="Tcell","B-cells"="Bcell","Plasmablasts"="Plasma","Myeloid"="Myeloid")
lab2 <- unname(map2[meta$celltype_major])
common <- intersect(rownames(sc1), rownames(sc2))
ref <- cbind(sc1[common,], sc2[common,]); labels <- c(lab1, lab2)
keep <- unlist(lapply(split(seq_along(labels), labels), function(ii) sample(ii, min(length(ii),250))))
ref <- ref[,keep]; labels <- labels[keep]
cat("ref cells:",length(labels),"\n"); print(table(labels))
bulk <- openxlsx::read.xlsx("D:/BCBM_Project/data/raw/GSE184869_expression.xlsx", sheet="log2TMMCPM")
rownames(bulk) <- bulk$ensembl_gene_id; bulk$ensembl_gene_id <- NULL
pm <- read.csv("D:/BCBM_Project/results/tables/program_decomposition.csv")
emap <- setNames(pm$symbol, pm$ensembl)
bulk$symbol <- emap[rownames(bulk)]; bulk <- bulk[!is.na(bulk$symbol),]
g <- aggregate(.~symbol, data=bulk, FUN=mean)
rownames(g) <- g$symbol; g$symbol <- NULL
mixture <- t(as.matrix(g))
common2 <- intersect(rownames(ref), colnames(mixture))
ref2 <- t(as.matrix(ref[common2,])); mix2 <- mixture[, common2]
cat("ref:",dim(ref2)," mix:",dim(mix2),"\n")
prism <- new.prism(reference=ref2, mixture=mix2, input.type="counts",
  cell.type.labels=labels, cell.state.labels=labels, key="Tumor")
prism <- run.prism(prism, n.cores=1, gibbs.control=list(chain.length=150, burn.in=50, thinning=2))
frac <- get.fraction(prism, which.theta="final", state.or.type="type")
write.csv(frac, "D:/BCBM_Project/results/tables/bp_combined_fractions.csv")
cat("\n=== mean fractions ===\n"); print(round(colMeans(frac),3))
