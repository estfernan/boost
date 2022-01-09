library(boost)

## load sample dataset
data(mob)

## estimate the size factor and extract one sample gene
s <- get.size.factor(mob, estimation.method = "TSS")
g <- mob[, "Apoe"]

## run the statistical test
res <- SPARK(g, mob.spots, s, gene.name = "Apoe")
print(res)
