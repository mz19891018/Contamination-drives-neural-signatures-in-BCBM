## 10b_singlecell.R - manual module scoring (Seurat v5 layer-safe)
suppressPackageStartupMessages({library(Seurat); library(dplyr)})

seu <- readRDS("D:/BCBM_Project/data/raw/GSE324453_seurat.rds")
cat("dim:", dim(seu), "\n")
seu <- NormalizeData(seu, verbose=FALSE)

adapt <- readLines("D:/BCBM_Project/results/tables/geneset_adaptation.txt")
comp  <- readLines("D:/BCBM_Project/results/tables/geneset_competence.txt")
ag <- intersect(adapt, rownames(seu)); cg <- intersect(comp, rownames(seu))
cat("present adapt:", length(ag), " comp:", length(cg), "\n")

expr <- t(as.matrix(LayerData(seu, layer="data")))  # cells x genes

# z-score each gene across all cells, then mean per set (module score)
z <- scale(expr)
scorem <- function(genes) {
  g <- intersect(genes, colnames(z))
  rowMeans(z[, g, drop=FALSE], na.rm=TRUE)
}
seu$prog_adapt <- scorem(ag)
seu$prog_comp  <- scorem(cg)

meta <- seu@meta.data
agg <- meta %>% group_by(common_cell_clusters) %>%
  summarise(n=n(), adapt=mean(prog_adapt,na.rm=TRUE),
            comp=mean(prog_comp,na.rm=TRUE)) %>% arrange(desc(adapt))
cat("\n=== mean program score by cell type ===\n")
print(as.data.frame(agg), digits=3)

# focus tumor cells
tum <- meta[grep("Tumor", meta$common_cell_clusters), ]
cat("\n=== Tumor cells only ===\n")
cat("n tumor cells:", nrow(tum), "\n")
cat("adapt score range:", round(range(tum$prog_adapt),3), "\n")
cat("comp  score range:", round(range(tum$prog_comp),3), "\n")

write.csv(meta[, c("common_cell_clusters","sampleDemuxIDstrict","prog_adapt","prog_comp")],
          "D:/BCBM_Project/results/tables/singlecell_program_scores.csv", row.names=FALSE)

# save a compact summary
write.csv(agg, "D:/BCBM_Project/results/tables/singlecell_by_celltype.csv", row.names=FALSE)
cat("\nsaved.\n")
