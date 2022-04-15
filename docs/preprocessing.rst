Data Pre-processing
=============================

Before SV gene detection, there are several steps to prepare the data. 

Data Filtration
------------------------------

Spots with fewer than ten total counts across all genes are excluded and genes with more than 90% zero read counts across all spots are dropped. 

::
    # filter data
    filter.st <- function(count, sample_info, min_total = 10,min_percentage = 0.1){
      gene_num <- ncol(count)
      sample_num <- nrow(count)
      if(sum(rowSums(count) < min_total) == 0){
        if (sum(colSums(count == 0) > (1-min_percentage)*sample_num) > 0){
          sample_f <- sample_info
          count_f <- count[,-(which(colSums(count == 0) > (1-min_percentage)*sample_num))]
        }
        else{
          sample_f <- sample_info
          count_f <- count
        }}
      else{
        if (sum(colSums(count[-which(rowSums(count)<min_total),] == 0) > (1-min_percentage)*sample_num) > 0){
          sample_f <- sample_info[-which(rowSums(count)<min_total),]
          count_f <- count[-which(rowSums(count)<min_total),-(which(colSums(count[-which(rowSums(count)<min_total),] == 0) > (1-min_percentage)*sample_num))]
        }
        else{
          sample_f <- sample_info[-which(rowSums(count)<min_total),]
          count_f <- count[-which(rowSums(count)<min_total),]
        }
      }
      return(list(sample_f,count_f))
    }


    filter_result <- filter.st(count, sample_info)
    count <- filter_result[[2]]
    sample_info <- filter_result[[1]]

    print(paste0('Dimension of count matrix: ', dim(count)[1],', ', dim(count)[2]))
    ## [1] "Dimension of count matrix: 260, 11360"
    
    print(paste0('Dimension of location information matrix: ', dim(sample_info)[1],', ', dim(sample_info)[2]))
    ## [1] "Dimension of location information matrix: 260, 2"

After filtration, the Mouse Olfactory Bulb data has 260 sample points and 11360 genes.


Size Factor Estimation
-----------------------------
Size factor is the input for some SV gene detection approaches. In boost package, we use get.size.factor function to estimate it. The inputs of this function are count matrix Y and size factor estimation method ("TSS", "Q75", "RLE" and "TMM"). Output is a vector with length $n$, representing the estimated size factor for spots. In this example, we choose TSS (total sum scale) as the estimation method. 

::
	size_factor <- get.size.factor(count, estimation.method = "TSS")
	size_factor <- size_factor/mean(size_factor)


Expression Counts Normalization
------------------------------------

Some methods need normalized gene expression levels as input. normalize.st is the function for count data to adjust for the library size, stablize the variance, and do the log-transformation. It provides seven normalization methods ("TSS", "Q75", "RLE", "TMM", "A-VST", "N-VST" and "log-VST"). The output is the normalized expression level matrix, which has the same shape as the input count matrix $Y$.

::
    normalized_count <- normalize.st(count, scaling.method = "TSS")
    print(normalized_count[1:10, 1:10])



