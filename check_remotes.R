for (p in c("remotes","devtools")) cat(p, ":", requireNamespace(p, quietly=TRUE), "\n")
