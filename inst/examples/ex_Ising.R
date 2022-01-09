library(boost)

## load sample dataset
data(mob)

## extract a sample gene, dichotomise expression levels, and get spatial network
g <- binarize.st(mob, "Apoe", cluster.method = "GMC")
A <- get.neighbors(mob.spots, 4, method = "distance")

## fit the model
res <- BOOST.Ising(g, A, gene.name = "Apoe", n.iter = 500)
print(res)
