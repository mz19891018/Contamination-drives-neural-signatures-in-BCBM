## 26_bayesprism.R - deconvolve GSE184869 bulk using GSE324453 scRNA reference
suppressPackageStartupMessages({library(Seurat); library(BayesPrism); library(Matrix)})
set.seed(1)

seu <- readRDS("D:/BCBM_Project/data/raw/GSE324453_seurat.rds")
sc <- LayerData(seu, layer="counts")   # genes x cells
ct <- as.character(seu$common_cell_clusters)
# subsample reference cells to ~6000 total for speed
set.seed(1)
keep <- unlist(lapply(split(seq_along(ct), ct), function(ii) sample(ii, min(length(ii), 250))))
sc <- sc[, keep]; ct <- ct[keep]
cat("reference cells after subsample:", length(ct), "\n")

bulk <- openxlsx::read.xlsx("D:/BCBM_Project/data/raw/GSE184869_expression.xlsx", sheet="log2TMMCPM")
rownames(bulk) <- bulk$ensembl_gene_id; bulk$ensembl_gene_id <- NULL
pm <- read.csv("D:/BCBM_Project/results/tables/program_decomposition.csv")
emap <- setNames(pm$symbol, pm$ensembl)
bulk$symbol <- emap[rownames(bulk)]
bulk <- bulk[!is.na(bulk$symbol),]
g <- aggregate(.~symbol, data=bulk, FUN=mean)
rownames(g) <- g$symbol; g$symbol <- NULL
mixture <- t(as.matrix(g))   # samples x genes
cat("bulk:", dim(mixture), " sc:", dim(sc), "\n")

common <- intersect(rownames(sc), colnames(mixture))
sc2 <- t(as.matrix(sc[common,]))   # cells x genes
mix2 <- mixture[, common]
cat("common genes:", length(common), " ref:", dim(sc2), "\n")

prism <- new.prism(
  reference = sc2,
  mixture = mix2,
  input.type = "counts",
  cell.type.labels = ct,
  cell.state.labels = ct,
  key = NULL
)
prism <- run.prism(prism, n.cores = 1,
                   gibbs.control = list(chain.length = 150, burn.in = 50, thinning = 2))
theta <- get.fraction(prism, which.theta="final", state.or.type="type")
write.csv(theta, "D:/BCBM_Project/results/tables/bp_cell_fractions.csv")
cat("\n=== mean cell-type fraction across 90 samples ===\n")
print(round(colMeans(theta),3))
cat("\nsaved.\n")
