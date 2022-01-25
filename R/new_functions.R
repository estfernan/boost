rank.aggregation <- function(data, K, method = c('GEO', 'MC2'), ties.method = ties.method = c("average", "first", "last", "random", "max", "min"){
  
  library(Hmisc)
  library(TopKLists)
  
  p <- dim(data)[2] - 1
  n <- dim(data)[1]
  
  ## Order data by gene name
  data <- data[order(data$gene), ]
  
  
  ## data Pre-processing for rank aggregation ==================================
  ## Sort out rows with missing values for all methods
  mm <- which(rowSums(is.na(data[, 2:(p + 1)])) == p)
  temp <- data[mm, ]     # Save these data into particular file
  data <- data[-mm, ]     # remove these rows
  
  
  ## Sort out genes in the top K in any of the lists
  gene_top <- unique(c(sapply(2:(p + 1), function(x) data$gene[order(data[, x])][1:K])))
  
  ## Genes are not in the top K in any of the five lists
  gene_tt <- data$gene[data$gene %nin% gene_top]
  tt <- data[data$gene %nin% gene_top, ]
  temp <- rbind(temp, tt)
  
  ## Generate final dataset used for rank aggregation
  data <- data[data$gene %in% gene_top, ]
  data <- data[order(data$gene), ]    
  
  ## Generate the list containing individual ranked lists.
  input <- lapply(2:(p + 1), function(x) data$gene[order(data[, x], na.last = NA)])
  
  
  if (method == 'GEO'){
  
  ## Rank aggregation using geometric mean (GEO) ===============================
  borda = Borda(input)
  
  ## GEO rank result
  score <- borda$Scores[, 3][order(borda$TopK[, 3])]
  data$rank <- rank(score, ties.method = ties.method)
  
  ## Genes with all missing results are assigned with the missing values
  temp$rank <- NA

  ## Gene not with all missing results but not ranked are assigned with maximum rank plus one.
  temp$rank[temp$gene %in% gene_tt] <- length(gene_top) + 1

}

else if (method == 'MC2'){
 
  ## Rank aggregation using MC2 ================================================
  ## MC2 rank result
  MCO <- MC(input = input)
  score <- MCO$MC2.Prob[order(MCO$MC2.TopK)]
  data$rank <- rank(-score, ties.method = ties.method)  # Better rank (small values) with large probability
  
  ## Genes with all missing results are assigned with the missing values
  temp$rank <- NA
  
  ## Gene not with all missing results but not ranked are assigned with maximum rank plus one.
  temp$rank[temp$gene %in% gene_tt] <- length(gene_top) + 1
  

}
  
else {stop("value passed to 'method' is not a valid option")}
  
  ## Combine data
  output <- rbind(data, temp)
  output <- output[order(output$gene), c(1, (p + 2))]    
  
  
  ## Export the results ========================================================
  return(output)
}  


compute.metrics <- function(predictor, truth, predictor.type = c('BF', 'p-value'), threshold = NULL){
  if (is.vector(predictor) == FALSE){
    stop("value passed to 'predictor' is not a valid option")
  }
  if (is.vector(truth) == FALSE){
    stop("value passed to 'truth' is not a valid option")
  }
  if (length(predictor) != length(truth)){
    stop("values passed to 'predictor' and 'truth' have different length")
  }
  if (!(predictor.type %in% c('BF', 'p-value'))){
    stop("value passed to 'predictor.type' is not a valid option")
  }
  if (is.null(threshold)){
    if (predictor.type == 'BF'){
      threshold <- 150
    }
    else if (predictor.type == 'p-value'){
      threshold <- 0.05
    }
  }
  
  library(pROC)
  
  AUC <- auc(roc(truth, predictor))
  
  if (predictor.type == 'BF'){
    predictor_binary <- predictor >= threshold
  }
  else if (predictor.type == 'p-value'){
    predictor_binary <- predictor <= threshold
  }
  confusion_matrix <- table(predictor_binary,truth)

  if(nrow(confusion_matrix)==1 &rownames(confusion_matrix)[1]=='TRUE'){
    confusion_matrix = rbind( c(0,0),confusion_matrix )
    row.names(confusion_matrix) = c("FALSE","TRUE")
  }
  if(nrow(confusion_matrix)==1& rownames(confusion_matrix)[1]=='FALSE'){
    confusion_matrix = rbind( confusion_matrix,c(0,0) )
    row.names(confusion_matrix) = c("FALSE","TRUE")
  }
  
  TN <- confusion_matrix[1]
  FP <- confusion_matrix[2]
  FN <- confusion_matrix[3]
  TP <- confusion_matrix[4]
  
  Sensitivity <- ifelse(TP==0,0,TP/(TP+FN))
  Specificity <- ifelse(TN==0,0,TN/(TN+FP))
  F1_score <- ifelse(TP==0,0,2*TP/(2*TP+FP+FN))
  FDR <- ifelse(FP==0,0,FP/(FP+TP))
  MCC <- ifelse((TP+FP)==0|(TP+FN)==0|(TN+FP)==0|(TN+FN)==0,0,(TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  return(list(Sensitivity =Sensitivity,Specificity = Specificity,F1_score = F1_score,FDR = FDR, AUC = AUC,MCC = MCC ))
}

  
