## 20_nichenet_full.R - official NicheNet ligand activity
suppressPackageStartupMessages({library(Seurat); library(nichenetr); library(dplyr); library(Matrix)})

seu <- readRDS("D:/BCBM_Project/data/raw/GSE324453_seurat.rds")
seu <- NormalizeData(seu, verbose=FALSE)
Idents(seu) <- factor(seu$common_cell_clusters)

adapt <- readLines("D:/BCBM_Project/results/tables/geneset_adaptation.txt")

ligand_target_matrix <- readRDS("D:/BCBM_Project/data/raw/ligand_target_matrix.rds")
lr_network <- readRDS("D:/BCBM_Project/data/raw/lr_network_human.rds")
cat("ltm dim:", dim(ligand_target_matrix), " lr_network:", nrow(lr_network), "\n")

# mean expression per cell type, genes x groups
M <- LayerData(seu, layer="data")
ct <- seu$common_cell_clusters
groups <- c("Astrocytes","Myeloid cells","Neurons","ODG","Endothelial cells","Tumor cells","Tumor cells (prol.)")
expr_mean <- sapply(groups, function(g) rowMeans(M[, ct==g, drop=FALSE]))
rownames(expr_mean) <- rownames(M)

# expressed ligands/receivers threshold: > 0.5 log-normalized
expressed_genes <- rownames(expr_mean)[rowMeans(expr_mean)>0.5]

# Proper NicheNet: receiver genes-of-interest = genes in TUMOR cells correlated with adaptation score
tcells <- which(ct %in% c("Tumor cells","Tumor cells (prol.)"))
ag <- intersect(adapt, rownames(M))
adapt_score <- setNames(colMeans(M[ag, tcells, drop=FALSE]), colnames(M)[tcells])
Mt <- M[, tcells, drop=FALSE]
# correlation of each expressed gene with adaptation score (subset for speed)
cand <- intersect(expressed_genes, rownames(M))
# use genes in the NicheNet target matrix columns
mat_targets <- colnames(ligand_target_matrix)
goi_candidates <- intersect(cand, mat_targets)
cat("genes tested for GOI:", length(goi_candidates), "\n")
Mt_sub <- as.matrix(Mt[goi_candidates,])
gcors <- apply(Mt_sub, 1, function(g) cor(g, adapt_score, method="spearman"))
gcors <- sort(gcors, decreasing=TRUE)
goi <- head(names(gcors), min(100, length(gcors)))
cat("GOI (top100 adapted genes in tumor):", length(goi), "\n")
cat("top GOI:", paste(head(goi,10), collapse=", "), "\n")

targets <- goi
background <- intersect(expressed_genes, mat_targets)
cat("target genes:", length(targets), " background:", length(background), "\n")

# candidate ligands: expressed by senders
senders <- c("Astrocytes","Myeloid cells","Neurons","ODG","Endothelial cells")
ligands_senders <- rownames(expr_mean)[rowMeans(expr_mean[,senders,drop=FALSE])>0.5]
candidate_ligands <- intersect(intersect(lr_network$from, ligands_senders), rownames(ligand_target_matrix))
cat("candidate ligands from senders:", length(candidate_ligands), "\n")

# activity = sum of NicheNet predicted weights of a ligand on the target set
submat <- ligand_target_matrix[candidate_ligands, targets, drop=FALSE]
lig_act <- sort(setNames(rowSums(submat), rownames(submat)), decreasing=TRUE)
lig_act_df <- data.frame(ligand=names(lig_act), activity=as.numeric(lig_act))

# which senders express the top ligands?
top <- head(lig_act_df, 25)
top$best_sender <- sapply(top$ligand, function(l){
  vals <- expr_mean[l, senders]; names(which.max(vals))
})
cat("\n=== NicheNet-ranked ligands inducing the ADAPTATION program ===\n")
print(top, row.names=FALSE)
write.csv(lig_act_df, "D:/BCBM_Project/results/tables/nichenet_full_ligand_activity.csv", row.names=FALSE)

# top ligand-target links for top ligands
top_lig <- head(lig_act_df$ligand,10)
links <- lapply(top_lig, function(l){
  o <- order(ligand_target_matrix[l, targets], decreasing=TRUE)[1:8]
  data.frame(ligand=l, target=targets[o], weight=ligand_target_matrix[l,targets[o]])
})
links <- do.call(rbind, links)
print(head(links,30), row.names=FALSE)
write.csv(links, "D:/BCBM_Project/results/tables/nichenet_full_top_targets.csv", row.names=FALSE)
cat("\nsaved.\n")
