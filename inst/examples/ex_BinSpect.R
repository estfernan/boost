library(boost)

## load sample dataset
data(mob)

## extract a sample gene, dichotomise expression levels, and get spatial network
g <- binarize.st(mob, "Apoe", cluster.method = "GMC")
A <- get.neighbors(mob.spots, 4, method = "distance")

## run the statistical test
res <- binSpect(g, A, do.fisher.test = FALSE, gene.name = "Apoe")
print(res)
