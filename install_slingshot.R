options(repos = c(CRAN = "https://cloud.r-project.org"))
suppressMessages(BiocManager::install(c("slingshot","SingleCellExperiment"), ask=FALSE, update=FALSE, quiet=FALSE))
for (p in c("slingshot","SingleCellExperiment")) cat(p, ":", requireNamespace(p, quietly=TRUE), "\n")
