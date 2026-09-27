## 36_sim_full_logFC.R - save full simulated logFC per f (1 rep)
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
set.seed(20260926)
seu <- readRDS("D:/BCBM_Project/data/raw/GSE324453_seurat.rds")
cnt <- LayerData(seu, layer="counts")
ct <- as.character(seu$common_cell_clusters); don <- as.character(seu$sampleDemuxIDstrict)
tum_idx <- which(ct %in% c("Tumor cells","Tumor cells (prol.)"))
tab <- table(don[tum_idx]); good <- names(tab)[tab>=500]
neuron_idx <- which(ct=="Neurons"); odg_idx <- which(ct=="ODG"); ast_idx <- which(ct=="Astrocytes")
bt <- table(ct[ct %in% c("Neurons","Astrocytes","ODG")])
wn <- bt["Neurons"]/sum(bt); wo <- bt["ODG"]/sum(bt); wa <- bt["Astrocytes"]/sum(bt)
N<-3000
fs <- c(0,0.02,0.05,0.075,0.10,0.125,0.15,0.20,0.25,0.30,0.40,0.50)
out <- matrix(0,nrow=nrow(cnt),ncol=length(fs)); rownames(out)<-rownames(cnt); colnames(out)<-as.character(fs)
for(i in seq_along(fs)){
  f<-fs[i]
  Pm<-matrix(0,nrow=nrow(cnt),ncol=45); Bm<-Pm
  for(k in 1:45){
    d<-sample(good,1); ti<-which(don==d & ct %in% c("Tumor cells","Tumor cells (prol.)"))
    td<-sample(ti,N,replace=TRUE)
    nb<-round(N*f); nn<-round(nb*wn); no<-round(nb*wo); na<-nb-nn-no
    bd<-c(sample(neuron_idx,nn,replace=TRUE),sample(odg_idx,no,replace=TRUE),sample(ast_idx,na,replace=TRUE))
    Pm[,k]<-rowSums(cnt[,td,drop=FALSE]); Bm[,k]<-Pm[,k]+rowSums(cnt[,bd,drop=FALSE])
  }
  cpmf<-function(M){M<-t(t(M)/colSums(M))*1e6; log2(M+1)}
  out[,i]<-rowMeans(cpmf(Bm))-rowMeans(cpmf(Pm))
  cat("f=",f,"done\n",flush=TRUE)
}
saveRDS(out,"D:/BCBM_Project/results/tables/sim_full_logFC.rds")
cat("saved\n")
