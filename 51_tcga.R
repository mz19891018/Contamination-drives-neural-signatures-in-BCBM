## 51_tcga_background.R - TCGA-BRCA independent background
suppressPackageStartupMessages({library(recount3); library(Seurat); library(BayesPrism); library(Matrix); library(SummarizedExperiment)})
cat("Fetching TCGA-BRCA...\n")
rse <- create_rse_manual(organism="human", project="BRCA", project_home="data_sources/tcga",
  type="gene", annotation="gencode_v29")
cat("TCGA dims:", dim(rse), "\n")
# primary solid tumor samples
cold <- as.data.frame(colData(rse))
print(table(cold$tcga.gdc_cases.samples.sample_type, useNA="ifany")[1:5])
prim <- cold$tcga.gdc_cases.samples.sample_type == "Primary Tumor"
cat("primary tumors:", sum(prim), "\n")
# random 100
set.seed(1)
idx <- sample(which(prim), min(50, sum(prim)))
rse2 <- rse[,idx]
counts <- assay(rse2, "raw_counts")
# gene symbols from gencode
gn <- rowData(rse2)$gene_name
rownames(counts) <- gn
counts <- counts[!is.na(gn) & gn!="",]
counts <- rowsum(counts, group=rownames(counts))
cat("counts matrix:", dim(counts), "\n")

# reference (same combined as main)
seu <- readRDS("D:/BCBM_Project/data/raw/GSE324453_seurat.rds")
sc1 <- LayerData(seu, layer="counts"); ct1 <- as.character(seu$common_cell_clusters)
map1 <- c("Astrocytes"="Astrocyte","Neurons"="Neuron","ODG"="Oligodendrocyte",
  "Tumor cells"="Tumor","Tumor cells (prol.)"="Tumor","T/NK cells"="Tcell",
  "B cells"="Bcell","Plasma cells"="Plasma","Myeloid cells"="Myeloid",
  "Endothelial cells"="Endothelial","Pericytes"="Pericyte","Stromal cells"="Stroma")
lab1 <- unname(map1[ct1])
base <- "D:/BCBM_Project/data/raw/Wu_etal_2021_BRCA_scRNASeq/"
genes2 <- read.table(paste0(base,"count_matrix_genes.tsv"),stringsAsFactors=FALSE)$V1
bar2 <- read.table(paste0(base,"count_matrix_barcodes.tsv"),stringsAsFactors=FALSE)$V1
sc2 <- readMM(paste0(base,"count_matrix_sparse.mtx")); rownames(sc2)<-genes2; colnames(sc2)<-bar2
meta <- read.csv(paste0(base,"metadata.csv"), row.names=1)
map2 <- c("Cancer Epithelial"="Tumor","Normal Epithelial"="NormalEpi","CAFs"="CAF","PVL"="PVL",
  "Endothelial"="Endothelial","T-cells"="Tcell","B-cells"="Bcell","Plasmablasts"="Plasma","Myeloid"="Myeloid")
lab2 <- unname(map2[meta$celltype_major])
common <- intersect(rownames(sc1), rownames(sc2))
ref <- cbind(sc1[common,], sc2[common,]); labels <- c(lab1, lab2)
keep <- unlist(lapply(split(seq_along(labels), labels), function(ii) sample(ii, min(length(ii),250))))
ref <- ref[,keep]; labels <- labels[keep]
common2 <- intersect(rownames(ref), rownames(counts))
ref2 <- t(as.matrix(ref[common2,]))
mix2 <- t(as.matrix(counts[common2,]))
cat("ref:",dim(ref2)," mix:",dim(mix2),"\n")
prism <- new.prism(reference=ref2, mixture=mix2, input.type="counts",
  cell.type.labels=labels, cell.state.labels=labels, key="Tumor")
prism <- run.prism(prism, n.cores=1, gibbs.control=list(chain.length=150, burn.in=50, thinning=2), update.gibbs=FALSE)
frac <- get.fraction(prism, which.theta="first", state.or.type="type")
write.csv(frac, "D:/BCBM_Project/results/tables/tcga_background_fractions.csv")
brain <- frac[,"Neuron"]+frac[,"Astrocyte"]+frac[,"Oligodendrocyte"]
cat("\n=== TCGA-BRCA BACKGROUND (n=",nrow(frac),") ===\n")
cat("brain total: median=",median(brain)," IQR=",IQR(brain)," q95=",quantile(brain,0.95),"\n")
cat("by type: Neuron median=",median(frac[,"Neuron"])," Astrocyte=",median(frac[,"Astrocyte"]),
    " ODG=",median(frac[,"Oligodendrocyte"]),"\n")
cat("Tumor median=",median(frac[,"Tumor"]),"\n")
