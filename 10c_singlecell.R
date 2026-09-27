## 10c_singlecell.R - memory-efficient sparse module scoring
suppressPackageStartupMessages({library(Seurat); library(dplyr); library(Matrix)})

seu <- readRDS("D:/BCBM_Project/data/raw/GSE324453_seurat.rds")
cat("dim:", dim(seu), "\n")
seu <- NormalizeData(seu, verbose=FALSE)

adapt <- readLines("D:/BCBM_Project/results/tables/geneset_adaptation.txt")
comp  <- readLines("D:/BCBM_Project/results/tables/geneset_competence.txt")
ag <- intersect(adapt, rownames(seu)); cg <- intersect(comp, rownames(seu))
cat("present adapt:", length(ag), " comp:", length(cg), "\n")

# sparse normalized matrix: genes x cells
M <- LayerData(seu, layer="data")
cat("class:", class(M), "\n")

# module score = mean log-normalized expression of the gene set per cell (cells x 1)
score <- function(genes) {
  g <- intersect(genes, rownames(M))
  colMeans(M[g, , drop=FALSE])
}
seu$prog_adapt <- score(ag)
seu$prog_comp  <- score(cg)

meta <- seu@meta.data
agg <- meta %>% group_by(common_cell_clusters) %>%
  summarise(n=n(), adapt=mean(prog_adapt,na.rm=TRUE),
            comp=mean(prog_comp,na.rm=TRUE)) %>% arrange(desc(adapt))
cat("\n=== mean program score by cell type ===\n")
print(as.data.frame(agg), digits=3)

tum <- meta[grep("Tumor", meta$common_cell_clusters), ]
cat("\n=== Tumor cells only (n=", nrow(tum), ") ===\n", sep="")
cat("adapt mean:", round(mean(tum$prog_adapt),3), " range:", round(range(tum$prog_adapt),3), "\n")
cat("comp  mean:", round(mean(tum$prog_comp),3), " range:", round(range(tum$prog_comp),3), "\n")

# tumor subclusters: split tumor cells by median adaptation
tum$adapt_hi <- tum$prog_adapt > median(tum$prog_adapt)
cat("\nadapt-high tumor cells vs adapt-low:\n")
cat("  comp score in adapt-high:", round(mean(tum$prog_comp[tum$adapt_hi]),3),
    " vs adapt-low:", round(mean(tum$prog_comp[!tum$adapt_hi]),3), "\n")

write.csv(meta[, c("common_cell_clusters","sampleDemuxIDstrict","prog_adapt","prog_comp")],
          "D:/BCBM_Project/results/tables/singlecell_program_scores.csv", row.names=FALSE)
write.csv(agg, "D:/BCBM_Project/results/tables/singlecell_by_celltype.csv", row.names=FALSE)
cat("\nsaved.\n")
