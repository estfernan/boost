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

##' Performance of Gene Identification
##'
##' Calculate the performance of spatially variable (SV) gene identification
##' on simulated data.
##'
##' @param predictor A numeric vector of length \eqn{n} that denotes the
##'   \eqn{p}-values or Bayes factors (BFs).
##' @param truth A logical vector of length \eqn{n} that represents the ground
##'   truth corresponding to the predictor.
##' @param predictor.type A character string that specifies whether
##'   \eqn{p}-values of Bayes factors (BFs) were provided.
##' @param threshold A numeric value that specifies the cutoff for
##'   defining SV genes.
##'
##' @return A list object that contains six performance metrics
##'   (Sensitivity, Specificity, F1_score, FDR, AUC, and MCC).
##'
##' @export
##' @keywords metrics
##'
##' @importFrom pROC auc roc
##'
compute.metrics <- function(
  predictor, truth,
  predictor.type = c("BF", "p-value"),
  threshold = NULL
)
{
  if (is.vector(predictor) == FALSE)
  {
    stop("value passed to 'predictor' is not a valid option")
  }

  if (is.vector(truth) == FALSE)
  {
    stop("value passed to 'truth' is not a valid option")
  }

  if (length(predictor) != length(truth))
  {
    stop("values passed to 'predictor' and 'truth' have different lengths")
  }

  if (!(predictor.type %in% c("BF", "p-value")))
  {
    stop("value passed to 'predictor.type' is not a valid option")
  }

  if (is.null(threshold))
  {
    if (predictor.type == "BF")
    {
      threshold <- 150
    }
    else if (predictor.type == "p-value")
    {
      threshold <- 0.05
    }
  }

  AUC <- auc(roc(truth, predictor))

  if (predictor.type == "BF")
  {
    predictor_binary <- predictor >= threshold
  }
  else if (predictor.type == "p-value")
  {
    predictor_binary <- predictor <= threshold
  }

  confusion_matrix <- table(predictor_binary,truth)

  if (nrow(confusion_matrix) == 1 & rownames(confusion_matrix)[1] == "TRUE")
  {
    confusion_matrix <- rbind(rep(0, 2), confusion_matrix)
    row.names(confusion_matrix) <- c("FALSE","TRUE")
  }

  if (nrow(confusion_matrix) == 1 & rownames(confusion_matrix)[1] == "FALSE")
  {
    confusion_matrix <- rbind(confusion_matrix, rep(0, 2))
    row.names(confusion_matrix) <- c("FALSE","TRUE")
  }

  TN <- confusion_matrix[1]
  FP <- confusion_matrix[2]
  FN <- confusion_matrix[3]
  TP <- confusion_matrix[4]

  Sensitivity <- ifelse(TP == 0, 0, TP/(TP + FN))
  Specificity <- ifelse(TN == 0, 0, TN/(TN + FP))

  F1_score <- ifelse(TP == 0, 0, 2*TP/(2*TP + FP + FN))
  FDR      <- ifelse(FP == 0, 0, FP/(FP + TP))

  MCC <- ifelse((TP + FP) == 0 | (TP + FN) == 0 | (TN + FP) == 0 | (TN + FN) == 0,
                0,
                (TP*TN - FP*FN) / (sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))))

  obj <- list(Sensitivity = Sensitivity,
              Specificity = Specificity,
              F1_score    = F1_score,
              FDR         = FDR,
              AUC         = AUC,
              MCC         = MCC)

  return(obj)
}
