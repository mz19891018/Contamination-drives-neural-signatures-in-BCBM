cat("R:", as.character(getRversion()), "\n")
cat("BiocManager:", requireNamespace("BiocManager", quietly=TRUE), "\n")
for (p in c("slingshot","TrajectoryUtils","SingleCellExperiment","decoupleR","dorothea","UCell","fgsea","pheatmap")) {
  cat(p, ":", requireNamespace(p, quietly=TRUE), "\n")
}
