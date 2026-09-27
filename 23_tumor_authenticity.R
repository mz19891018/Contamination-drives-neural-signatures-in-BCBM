## 23_tumor_authenticity.R - is adaptation real in tumor cells, or ambient/background?
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
set.seed(42)
seu <- readRDS("D:/BCBM_Project/data/raw/GSE324453_seurat.rds")
seu <- NormalizeData(seu, verbose=FALSE)
M <- LayerData(seu, layer="data")   # genes x cells
ct <- seu$common_cell_clusters

adapt <- readLines("D:/BCBM_Project/results/tables/geneset_adaptation.txt")
comp  <- readLines("D:/BCBM_Project/results/tables/geneset_competence.txt")

# negative-control tissue markers (NOT brain)
lung  <- c("SFTPC","SFTPB","SFTPA1","AGER","FOXJ1","MUC5AC","SCGB1A1","NKX2-1")
liver <- c("ALB","AFP","TTR","APOB","HNF4A","CYP3A4","SERPINA1","FGB")
prolif<- c("MKI67","TOP2A","PCNA","CCNB1","CDK1","MCM2")

score <- function(genes, cells){
  g <- intersect(genes, rownames(M))
  if(length(g)<3) return(NA)
  rowMeans(M[g, cells, drop=FALSE])
}

tum <- which(ct %in% c("Tumor cells","Tumor cells (prol.)"))
neur<- which(ct=="Neurons")

res <- data.frame(
  set = c("Adaptation","Competence","Lung(neg)","Liver(neg)","Proliferation"),
  tumor_mean = c(mean(score(adapt,tum)), mean(score(comp,tum)), mean(score(lung,tum)),
                 mean(score(liver,tum)), mean(score(prolif,tum))),
  neuron_mean = c(mean(score(adapt,neur)), mean(score(comp,neur)), mean(score(lung,neur)),
                  mean(score(liver,neur)), mean(score(prolif,neur)))
)
print(res, row.names=FALSE)

# random gene-set null: match adaptation size, same overall expression level
allg <- rownames(M)
base <- rowMeans(M)
n <- length(intersect(adapt, allg))
rand_scores <- replicate(500, {
  pick <- sample(allg, n)
  mean(rowMeans(M[pick, tum, drop=FALSE]))
})
obs <- mean(score(adapt,tum))
p_rand <- mean(rand_scores >= obs)
cat(sprintf("\nAdaptation in tumor: %.4f ; random-gene null mean=%.4f sd=%.4f ; empirical p=%.3f\n",
            obs, mean(rand_scores), sd(rand_scores), p_rand))

comp_obs <- mean(score(comp,tum))
n2 <- length(intersect(comp, allg))
rand2 <- replicate(500, mean(rowMeans(M[sample(allg,n2), tum, drop=FALSE])))
cat(sprintf("Competence in tumor: %.4f ; random null mean=%.4f ; empirical p=%.3f\n",
            comp_obs, mean(rand2), mean(rand2>=comp_obs)))

write.csv(res, "D:/BCBM_Project/results/tables/tumor_authenticity_scores.csv", row.names=FALSE)
saveRDS(list(adapt_rand=rand_scores, comp_rand=rand2), "D:/BCBM_Project/results/tables/tumor_authenticity_null.rds")
cat("saved.\n")
