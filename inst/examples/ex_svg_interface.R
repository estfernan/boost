## Methods available in boost and the SVG criteria used in STAr
svg.methods()

## Gene-filtering criteria used before each method in STAr
svg.filter.criteria()

\dontrun{
data(mob)
flt <- filter.st(mob, mob.spots, method = "SPARK")
sf  <- get.size.factor(flt$count, estimation.method = "TSS")

## single-gene results converted to the standard format
fit <- SPARK(flt$count[, "Penk"], flt$spots, size.factor = sf, gene.name = "Penk")
as.svg.result(fit)
}
