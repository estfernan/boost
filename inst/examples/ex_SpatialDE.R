library(boost)

## load sample dataset
data(mob)

## scale gene expression levels and extract one gene
mob.t <- normalize.st(mob, scaling.method = "TSS")
g     <- mob.t[, "Apoe"]

## run the statistical test
res <- SpatialDE(g, mob.spots, gene.name = "Apoe")
print(res)
