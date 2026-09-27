## 10_singlecell.R
## Load GSE324453 Seurat object, inspect annotations, score programs in malignant cells.
suppressPackageStartupMessages({library(Seurat); library(dplyr)})

rds <- "D:/BCBM_Project/data/raw/GSE324453_seurat.rds"
seu <- readRDS(rds)
cat("dim:", dim(seu), "\n")
cat("meta colnames:\n"); print(colnames(seu@meta.data))
cat("\n=== first rows of meta ===\n"); print(head(seu@meta.data, 3))

# guess cell-type / annotation columns
ann_cols <- grep("cell|type|cluster|ident|orig|sample|patient", colnames(seu@meta.data),
                 value=TRUE, ignore.case=TRUE)
cat("\n=== candidate annotation columns ===\n"); print(ann_cols)
for (cc in ann_cols) {
  cat("\n---", cc, "---\n"); print(table(seu@meta.data[[cc]], useNA="ifany"))
}

# program gene sets (from Python output)
adapt <- readLines("D:/BCBM_Project/results/tables/geneset_adaptation.txt")
comp  <- readLines("D:/BCBM_Project/results/tables/geneset_competence.txt")
cat("\nadapt genes:", length(adapt), " comp genes:", length(comp), "\n")

# which genes present in the object
ag <- intersect(adapt, rownames(seu)); cg <- intersect(comp, rownames(seu))
cat("present adapt:", length(ag), " present comp:", length(cg), "\n")

# Seurat AddModuleScore
seu <- AddModuleScore(seu, features=list(adapt=ag, comp=cg), name="prog_")
cat("\n=== mean program score by annotation ===\n")
for (cc in ann_cols) {
  cat("\n--- by", cc, "---\n")
  agg <- seu@meta.data %>% group_by(.data[[cc]]) %>%
    summarise(n=n(), adapt=mean(prog_adapt, na.rm=TRUE),
              comp=mean(prog_comp, na.rm=TRUE)) %>% arrange(desc(adapt))
  print(as.data.frame(agg))
}
# save scores
out <- data.frame(barcode=colnames(seu),
                  prog_adapt=seu$prog_adapt, prog_comp=seu$prog_comp,
                  seu@meta.data[, ann_cols, drop=FALSE])
write.csv(out, "D:/BCBM_Project/results/tables/singlecell_program_scores.csv", row.names=FALSE)
cat("\nsaved scores.\n")
