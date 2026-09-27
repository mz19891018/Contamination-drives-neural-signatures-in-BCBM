suppressPackageStartupMessages(library(BayesPrism))
print(args(BayesPrism:::valid.gibbs.control))
cat("\n--- default gibbs control: ---\n")
print(BayesPrism:::valid.gibbs.control(list()))
