options(timeout=600)
library(GEOquery)
g <- getGEO("GSE17607", destdir="D:/BCBM_Project/data/raw", getGPL=FALSE)
cat("done\n")
