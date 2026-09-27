## 29_mixing_simulation.R - in silico mixing: how much normal brain recreates observed logFC?
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
set.seed(20260926)

seu <- readRDS("D:/BCBM_Project/data/raw/GSE324453_seurat.rds")
cnt <- LayerData(seu, layer="counts")  # genes x cells
ct  <- as.character(seu$common_cell_clusters)
don <- as.character(seu$sampleDemuxIDstrict)

# tumor donors with >=500 tumor cells
tum_idx <- which(ct %in% c("Tumor cells","Tumor cells (prol.)"))
tab <- table(don[tum_idx])
good_don <- names(tab)[tab>=500]
cat("good tumor donors (>=500):", length(good_don), "\n")

# contaminating brain cells
brain_idx <- which(ct %in% c("Neurons","Astrocytes","ODG"))
# relative abundance for mixing
bt <- table(ct[brain_idx]); print(bt)
w_neuron <- bt["Neurons"]/sum(bt); w_odg <- bt["ODG"]/sum(bt); w_ast <- bt["Astrocytes"]/sum(bt)
cat("weights neuron/odg/ast:", round(c(w_neuron,w_odg,w_ast),3), "\n")

neuron_idx <- which(ct=="Neurons"); odg_idx <- which(ct=="ODG"); ast_idx <- which(ct=="Astrocytes")

# observed logFC from real paired DE
obs <- read.csv("D:/BCBM_Project/results/tables/program_decomposition.csv")
obs_lfc <- setNames(obs$mean_logFC_BM_vs_BP, obs$symbol)
obs_up <- rownames(obs[obs$program!="Intermediate" & obs$mean_logFC_BM_vs_BP>0.3,])
obs_up_sym <- obs$symbol[obs$program!="Intermediate" & obs$mean_logFC_BM_vs_BP>0.3]
cat("observed up genes to match:", length(obs_up_sym), "\n")

N <- 3000
fs <- c(0,0.02,0.05,0.075,0.10,0.125,0.15,0.20,0.25,0.30,0.40,0.50)
nrep <- 6
f_list <- c()
rho_list <- c()
dose_C1QL1 <- c(); dose_ATP1A2 <- c(); dose_TUBB4A <- c()

for (f in fs){
  rhos <- c()
  d1<-c(); d2<-c(); d3<-c()
  for(rep in 1:nrep){
    # build 45 pairs, same tumor donor within pair
    lfc <- rep(0, nrow(cnt)); names(lfc) <- rownames(cnt)
    P_mat <- matrix(0,nrow=nrow(cnt),ncol=45)
    B_mat <- P_mat
    for(k in 1:45){
      d <- sample(good_don,1)
      t_idx <- which(don==d & ct %in% c("Tumor cells","Tumor cells (prol.)"))
      # sample 3000 tumor cells WITH REPLACEMENT
      tumor_draw <- sample(t_idx, N, replace=TRUE)
      n_brain <- round(N*f)
      # brain draw
      nb_n <- round(n_brain*w_neuron); nb_o <- round(n_brain*w_odg); nb_a <- n_brain-nb_n-nb_o
      b_draw <- c(sample(neuron_idx,nb_n,replace=TRUE),
                  sample(odg_idx,nb_o,replace=TRUE),
                  sample(ast_idx,nb_a,replace=TRUE))
      P <- rowSums(cnt[,tumor_draw,drop=FALSE])
      B <- P + rowSums(cnt[,b_draw,drop=FALSE])
      P_mat[,k] <- P; B_mat[,k] <- B
    }
    # CPM log2
    cpm <- function(M){ M <- t(t(M)/colSums(M))*1e6; log2(M+1) }
    Pc <- cpm(P_mat); Bc <- cpm(B_mat)
    sim_lfc <- rowMeans(Bc)-rowMeans(Pc)
    names(sim_lfc) <- rownames(cnt)
    # spearman on observed up genes
    common <- intersect(obs_up_sym, names(sim_lfc))
    o <- obs_lfc[common]; s <- sim_lfc[common]
    rhos <- c(rhos, suppressWarnings(cor(o,s,method="spearman")))
    d1<-c(d1, sim_lfc["C1QL1"]); d2<-c(d2, sim_lfc["ATP1A2"]); d3<-c(d3, sim_lfc["TUBB4A"])
  }
  f_list<-c(f_list,f); rho_list<-c(rho_list,mean(rhos,na.rm=TRUE))
  dose_C1QL1<-c(dose_C1QL1,mean(d1,na.rm=TRUE))
  dose_ATP1A2<-c(dose_ATP1A2,mean(d2,na.rm=TRUE))
  dose_TUBB4A<-c(dose_TUBB4A,mean(d3,na.rm=TRUE))
  cat(sprintf("f=%.3f  spearman=%.3f  C1QL1=%.2f ATP1A2=%.2f TUBB4A=%.2f\n",
              f, mean(rhos,na.rm=TRUE), mean(d1,na.rm=TRUE),mean(d2,na.rm=TRUE),mean(d3,na.rm=TRUE)),flush=TRUE)
}

res <- data.frame(f=fs, spearman=rho_list, C1QL1=dose_C1QL1, ATP1A2=dose_ATP1A2, TUBB4A=dose_TUBB4A)
write.csv(res, "D:/BCBM_Project/results/tables/mixing_simulation.csv", row.names=FALSE)
fstar <- fs[which.max(rho_list)]
cat("\n>>> f* (best Spearman) =", fstar, " rho=",max(rho_list),"\n")
cat(">>> BayesPrism estimated brain fraction in BM = 0.16\n")
