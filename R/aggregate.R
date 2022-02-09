##
## R package boost by Esteban Fernández, Xi Jiang, Suhana Bedi, and Qiwei Li
## Copyright (C) 2021
##
## This file is part of the R package boost.
##
## The R package boost is free software: You can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by the
## Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package boost is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY without even the implied warranty of MERCHANTABILITY
## or FITNESS FOR A PARTICULAR PURPOSE.
##

##' Aggregate Rankings
##'
##' Generate an aggregated ranking from multiple base rankers.
##'
##' See Schimek, M. G. et al. (2015) for more information on the procedure.
##'
##' @param data An \eqn{n}-by-\eqn{p} integer matrix \eqn{Y} that denotes the
##'   absolute gene expression count table. Each entry is the read count for
##'   gene \eqn{j} collected at spot \eqn{i}.
##' @param data An \eqn{n}-by-\eqn{(p+1)} numeric matrix that represents a
##'   table of \eqn{p}-values with \eqn{n} observations and
##'   \eqn{p+1} variables. Note that the first column should be the gene name
##'   and the columns following after should be the methods
##'   (i.e., base rankers).
##' @param K A numeric value that indicates how many top-\eqn{K} genes to sort
##'   out in each base ranker.
##' @param method A character string that specifies the rank aggregation method.
##' @param ties.method A character string that specifies how ties should
##'   be treated (see "Details" in [base::rank()]).
##'
##' @return A numeric matrix that represents a table with genes and their rank.
##'
##' @references
##'
##' Li, X., Wang, X., & Xiao, G. (2019). A comparative study of rank
##' aggregation methods for partial and top ranked lists in genomic
##' applications. _Briefings in bioinformatics_, _20_(1), 178–189.
##' <https://doi.org/10.1093/bib/bbx101>.
##'
##' Schimek M., Budinska E., Kugler K., Svendova V., Ding J., Lin S.
##' (2015). “TopKLists: a comprehensive R package for statistical inference,
##' stochastic aggregation, and visualization of multiple omics ranked lists.”
##' _Stat Appl Genet Mol Biol_, 311-6.
##' <http://www.degruyter.com/doi/10.1515/sagmb-2014-0093>.
##'
##' @export
##' @importFrom TopKLists Borda MC
##'
rank.aggregation <- function(
  data, K,
  method = c('GEO', 'MC2'),
  ties.method = c("average", "first", "last", "random", "max", "min")
)
{
  p <- dim(data)[2] - 1
  n <- dim(data)[1]

  ## Order data by gene name
  data <- data[order(data$gene), ]

  ## data Pre-processing for rank aggregation ================================

  ## Sort out rows with missing values for all methods
  mm   <- which(rowSums(is.na(data[, 2:(p + 1)])) == p)
  temp <- data[mm, ]      # Save these data into particular file
  data <- data[-mm, ]     # remove these rows


  ## Sort out genes in the top K in any of the lists
  gene_top <- unique(c(sapply(2:(p + 1), function(x) data$gene[order(data[, x])][1:K])))

  ## Genes are not in the top K in any of the five lists
  gene_tt <- data$gene[!(data$gene %in% gene_top)] # %nin%
  tt <- data[!(data$gene %in% gene_top), ] # %nin%
  temp <- rbind(temp, tt)

  ## Generate final dataset used for rank aggregation
  data <- data[data$gene %in% gene_top, ]
  data <- data[order(data$gene), ]

  ## Generate the list containing individual ranked lists.
  input <- lapply(2:(p + 1), function(x) data$gene[order(data[, x], na.last = NA)])

  ## Rank aggregation using geometric mean (GEO) =============================
  if (method == 'GEO')
  {
    borda = Borda(input)

    ## GEO rank result
    score <- borda$Scores[, 3][order(borda$TopK[, 3])]
    data$rank <- rank(score, ties.method = ties.method)

    ## Genes with all missing results are assigned with the missing values
    temp$rank <- NA

    ## Gene not with all missing results but not ranked are assigned with maximum rank plus one.
    temp$rank[temp$gene %in% gene_tt] <- length(gene_top) + 1
  }

  ## Rank aggregation using MC2 ==============================================
  else if (method == 'MC2')
  {
    ## MC2 rank result
    MCO <- MC(input = input)
    score <- MCO$MC2.Prob[order(MCO$MC2.TopK)]
    data$rank <- rank(-score, ties.method = ties.method)  # Better rank (small values) with large probability

    ## Genes with all missing results are assigned with the missing values
    temp$rank <- NA

    ## Gene not with all missing results but not ranked are assigned with maximum rank plus one.
    temp$rank[temp$gene %in% gene_tt] <- length(gene_top) + 1
  }

  else {
    stop("value passed to 'method' is not a valid option")
  }

  ## Combine data
  output <- rbind(data, temp)
  output <- output[order(output$gene), c(1, (p + 2))]

  ## Export the results ======================================================

  return(output)
}
