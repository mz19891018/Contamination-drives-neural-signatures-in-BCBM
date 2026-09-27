for (p in c("nichenetr","tidyverse","matrixStats")) {
  cat(p, ":", requireNamespace(p, quietly=TRUE), "\n")
}
