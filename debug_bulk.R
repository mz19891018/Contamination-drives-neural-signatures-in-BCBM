suppressPackageStartupMessages({library(readxl)})
bulk <- read_excel("D:/BCBM_Project/data/raw/GSE184869_expression.xlsx", sheet="log2TMMCPM")
bulk <- as.data.frame(bulk); rownames(bulk) <- bulk[[1]]; bulk[[1]]<-NULL
ens2sym <- read.csv("D:/BCBM_Project/results/tables/program_decomposition.csv")
m <- setNames(ens2sym$symbol, ens2sym$ensembl)
sym <- m[rownames(bulk)]
cat("N total:",length(sym)," NA:",sum(is.na(sym))," empty:",sum(sym=="",na.rm=TRUE),"\n")
cat("unique non-empty:",length(unique(sym[!is.na(sym)&sym!=""])),"\n")
print(head(table(sym,useNA="ifany")))
