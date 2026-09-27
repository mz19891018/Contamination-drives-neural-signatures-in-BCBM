## 19_niche.R - NicheNet-style causal ligand ranking
## Goal: which sender ligands/receptors best predict the tumor adaptation program?
## Uses CellChatDB human LR pairs; ranks ligands by (sender expression * receptor-adaptation coupling).
suppressPackageStartupMessages({library(Seurat); library(CellChat); library(dplyr); library(Matrix)})

seu <- readRDS("D:/BCBM_Project/data/raw/GSE324453_seurat.rds")
seu <- NormalizeData(seu, verbose=FALSE)
M <- LayerData(seu, layer="data")   # genes x cells, sparse

adapt <- readLines("D:/BCBM_Project/results/tables/geneset_adaptation.txt")
ag <- intersect(adapt, rownames(M))
adapt_score <- setNames(colMeans(M[ag,,drop=FALSE]), colnames(M))

ct <- seu$common_cell_clusters
senders <- c("Astrocytes","Myeloid cells","Neurons","ODG","Endothelial cells")
targets <- c("Tumor cells","Tumor cells (prol.)")

# CellChatDB LR pairs
db <- CellChatDB.human$interaction
lr <- db[, c("ligand","receptor","pathway_name")]
lr <- lr[lr$ligand %in% rownames(M) & lr$receptor %in% rownames(M), ]

# mean ligand expression per sender type
lig_mat <- sapply(senders, function(s){
  cells <- which(ct==s)
  rowMeans(M[lr$ligand, cells, drop=FALSE])
})
rownames(lig_mat) <- lr$ligand

# receptor expression & adaptation coupling in tumor cells
tcells <- which(ct %in% targets)
adapt_t <- adapt_score[tcells]
# for each receptor, correlation with adaptation score
rec_expr <- as.matrix(M[lr$receptor, tcells, drop=FALSE])
rec_cor <- sapply(seq_len(nrow(rec_expr)), function(i){
  x <- rec_expr[i,]
  if(sd(x)==0) return(0)
  cor(x, adapt_t, method="spearman")
})
names(rec_cor) <- lr$receptor

# score each LR pair: mean sender ligand level * |receptor-adaptation correlation|
lr$sender <- apply(lig_mat[match(lr$ligand, rownames(lig_mat)),,drop=FALSE],1,max)
lr$rec_cor <- rec_cor[match(lr$receptor, names(rec_cor))]
lr$score <- lr$sender * abs(lr$rec_cor)
lr$sign <- sign(lr$rec_cor)

# focus on pairs where receptor positively tracks adaptation (induction)
res <- lr[order(-lr$score), ]
res <- res[!duplicated(res$interaction_name <- paste(res$ligand,res$receptor,sep="_")),]
cat("\n=== TOP ligand->receptor pairs coupled to tumor adaptation ===\n")
print(head(res[,c("ligand","receptor","pathway_name","sender","rec_cor","score")], 25), row.names=FALSE)

write.csv(res, "D:/BCBM_Project/results/tables/nichenet_style_ligand_ranking.csv", row.names=FALSE)

# aggregate by sender ligand
agg <- res %>% group_by(ligand) %>% summarise(top_pathway=first(pathway_name),
       max_sender=max(sender), best_cor=first(rec_cor), score=max(score)) %>% arrange(desc(score))
cat("\n=== TOP ligands by sender expression x adaptation coupling ===\n")
print(head(agg, 20), row.names=FALSE)
write.csv(agg, "D:/BCBM_Project/results/tables/nichenet_style_ligand_agg.csv", row.names=FALSE)
cat("\nsaved.\n")
