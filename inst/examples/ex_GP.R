\dontrun{
library(boost)

## load sample dataset
data(mob)

## extract a sample gene and get size factor
g <- mob[, "Apoe"]
s <- get.size.factor(mob, estimation.method = "TSS")

## fit the model
res <- BOOST.GP(g, mob.spots, size.factor = s, gene.name = "Apoe")
print(res)
}
