##
## Source code is from the edgeR package (https://bioinf.wehi.edu.au/edgeR/) and
##  downloaded from https://rdrr.io/bioc/edgeR/src/R/calcNormFactors.R
##  on April 12, 2018.
##
## Citations from authors:
##
##   Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor package
##   for differential expression analysis of digital gene expression data.
##   Bioinformatics 26, 139-140
##
##   McCarthy DJ, Chen Y and Smyth GK (2012). Differential expression analysis of
##   multifactor RNA-Seq experiments with respect to biological variation.
##   Nucleic Acids Research 40, 4288-4297
##
##   Chen Y, Lun ATL, Smyth GK (2016). From reads to genes to pathways:
##   differential expression analysis of RNA-Seq experiments using Rsubread and
##   the edgeR quasi-likelihood pipeline. F1000Research 5, 1438
##

##' Scale normalization of RNA-Seq data
##'
##' @noRd
##'
calcNormFactors_copy <- function(
  object, lib.size = NULL, method = c("TMM", "RLE", "upperquartile", "none"),
  refColumn = NULL,
  logratioTrim = 0.3, sumTrim = 0.05,
  doWeighting = TRUE, Acutoff = -1e10, p = 0.75,
  ...
)
{
  x <- as.matrix(object)

  if (any(is.na(x)))
  {
    stop("NA counts not permitted")
  }

  if (is.null(lib.size))
  {
    lib.size <- colSums(x)
  }

  if (any(is.na(lib.size)))
  {
    stop("NA lib.sizes not permitted")
  }

  method <- match.arg(method)
  allzero <- .rowSums(x > 0, nrow(x), ncol(x)) == 0

  if (any(allzero))
  {
    x <- x[!allzero, , drop = FALSE]
  }

  if (nrow(x) == 0 || ncol(x) == 1)
  {
    method = "none"
  }

  f <- switch(method,
    TMM = {
      f75 <- .calcFactorQuantile_copy(data = x, lib.size = lib.size, p = 0.75)

      if (is.null(refColumn))
      {
        refColumn <- which.min(abs(f75 - mean(f75)))
      }

      if (length(refColumn) == 0 | refColumn < 1 | refColumn > ncol(x))
      {
        refColumn <- 1
      }

      f <- rep(NA, ncol(x))

      for (i in 1:ncol(x))
      {
        f[i] <- .calcFactorWeighted_copy(
          obs = x[, i], ref = x[, refColumn],
          libsize.obs = lib.size[i],
          libsize.ref = lib.size[refColumn],
          logratioTrim = logratioTrim, sumTrim = sumTrim,
          doWeighting = doWeighting, Acutoff = Acutoff
        )
      }

      # f
      f * lib.size
    },
    RLE = .calcFactorRLE_copy(x),
    upperquartile = .calcFactorQuantile_copy(x, lib.size, p = p),
    none = rep(1, ncol(x))
  )

  return(f)
}

.calcFactorRLE_copy <- function(data)
{
  gm <- exp(rowMeans(log(data)))
  apply(data, 2, function(u) median((u/gm)[gm > 0]))
}

.calcFactorQuantile_copy <- function(data, lib.size, p = 0.75)
{
  y <- t(t(data)/lib.size)
  apply(y, 2, function(x) quantile(x, p = p))
}

.calcFactorWeighted_copy <- function(
  obs, ref, libsize.obs = NULL, libsize.ref = NULL,
  logratioTrim = 0.3, sumTrim = 0.05,
  doWeighting = TRUE, Acutoff = -1e10
)
{
  obs <- as.numeric(obs)
  ref <- as.numeric(ref)
  if (is.null(libsize.obs)) {
    nO <- sum(obs)
  } else {
    nO <- libsize.obs
  }
  if (is.null(libsize.ref)) {
    nR <- sum(ref)
  } else {
    nR <- libsize.ref
  }
  logR <- log2((obs/nO)/(ref/nR))
  absE <- (log2(obs/nO) + log2(ref/nR))/2
  v <- (nO - obs)/nO/obs + (nR - ref)/nR/ref
  fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
  logR <- logR[fin]
  absE <- absE[fin]
  v <- v[fin]

  if (max(abs(logR)) < 1e-6)
  {
    return(1)
  }

  n <- length(logR)
  loL <- floor(n*logratioTrim) + 1
  hiL <- n + 1 - loL
  loS <- floor(n*sumTrim) + 1
  hiS <- n + 1 - loS
  keep <- (rank(logR) >= loL & rank(logR) <= hiL) & (rank(absE) >= loS & rank(absE) <= hiS)

  if (doWeighting)
  {
    f <- sum(logR[keep]/v[keep], na.rm = TRUE)/sum(1/v[keep], na.rm = TRUE)
  }
  else
  {
    f <- mean(logR[keep], na.rm = TRUE)
  }

  if (is.na(f))
  {
    f <- 0
  }

  return(2^f)
}
