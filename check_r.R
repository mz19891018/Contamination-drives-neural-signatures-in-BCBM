cat("R", as.character(getRversion()), "\n")
for (pkg in c("Seurat","dplyr","ggplot2","hdf5r")) {
  cat(pkg, ":", requireNamespace(pkg, quietly=TRUE), "\n")
}
