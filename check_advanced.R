## check advanced R packages
for (p in c("CellChat","decoupleR","dorothea","SingleCellExperiment","UCell","fgsea")) {
  cat(p, ":", requireNamespace(p, quietly=TRUE), "\n")
}
