
Simulation Study - Example
======================================

SMP-Gym website (https://lce.biohpc.swmed.edu/smp_gym) provides a comprehensive comparison of different SV gene identification methods on simulated data. Here we provide an example to demonstrate how to implement these methods on the simulated data. 

Load Data
------------------------
In this example, we choose the first replicate under the setting of linear pattern and no zero-inflation.
::
    # download data from SMP-Gym website
    url <- 'https://lce.biohpc.swmed.edu/smp_gym/download/simulated_data/linear_pattern/data/linear_pattern_zero_0_replicate_1.zip'
    download.file(url, destfile= "simulation_example.zip", mode = "wb")
    unzip("simulation_example.zip")

    # read data
    count <- read.csv('linear_pattern_zero_0_replicate_1.count.csv')
    sample_info <- read.csv('linear_pattern_zero_0_replicate_1.loc.csv')
    ground_truth <- read.csv('linear_pattern_zero_0_replicate_1.gamma.csv')$x

    # select first 10 genes to run
    count <- count[, 1:10]
    ground_truth <- ground_truth[1:10]


Data Preparation
------------------------

Same as the data pre-processing step in real data example, we need to estimate the size factor for BOOST-GP and SPARK and get normalized expression counts for SpatialDE, BOOST-Ising and BinSpect. 
::
    # size factor estimation
    size_factor <- get.size.factor(count, estimation.method = "TSS")
    size_factor <- size_factor/mean(size_factor)

    # normalize count
    normalized_count_tss <- normalize.st(count, scaling.method = "TSS")
    normalized_count_log_vst <- normalize.st(count, scaling.method = "log-VST")

    # get neighbor information
    neighbor_info <- get.neighbors(sample_info, n.neighbors = 4, method = "distance")


Run SV Gene Detection Methods
--------------------------------------
All methods in this package require the expression level of one gene as input. To implement the SV gene identification for the whole dataset, we need to run these functions gene by gene. We can simply use 'for' loop to implement it. 
::
    n_gene <- dim(count)[2]

    # create a matrix to save results
    result <- matrix(0, nrow = n_gene, ncol = 6)
    colnames(result) <- c('BOOST-GP', 'BOOST-Ising', 'BinSpect-rank', 'BinSpect-kmeans','SPARK', 'SpatialDE')

    # run all methods and save SV gene measurements for all methods
    for (gene_index in 1:n_gene){

      abs_expr <- count[, gene_index]
      norm_expr_log_vst <- normalized_count_log_vst[, gene_index]
      norm_expr_tss <- normalized_count_tss[, gene_index]
      norm_expr_tss_log <- log(norm_expr_tss + 1)

      # run BOOST-GP
      result_boost_gp_temp <- BOOST.GP(abs_expr, sample_info, size.factor = size_factor, gene.name=colnames(count)[gene_index], n.iter = 1000)
      result[gene_index, 'BOOST-GP'] <- result_boost_gp_temp$measures$BF

      # run BOOST-Ising
      binary_expr_gmc <- binarize.st(normalized_count_tss, colnames(count)[gene_index], cluster.method =  "GMC")
      result_boost_ising_temp <- BOOST.Ising(binary_expr_gmc, neighbor_info, gene.name = colnames(count)[gene_index])
      result[gene_index, 'BOOST-Ising'] <- result_boost_ising_temp$measures$BF.neg

      # run BinSpect-rank
      binary_expr_rank <- binarize.st(normalized_count_tss, colnames(count)[gene_index], cluster.method =  "rank")
      result_binspect_rank_temp <- binSpect(binary_expr_rank, neighbor_info, gene.name = colnames(count)[gene_index])
      result[gene_index, 'BinSpect-rank'] <- result_binspect_rank_temp$measures$p.val

      # run BinSpect-kmeans
      binary_expr_kmeans <- binarize.st(log(normalized_count_tss + 1), colnames(count)[gene_index], cluster.method =  "k-means")
      result_binspect_kmeans_temp <- binSpect(binary_expr_kmeans, neighbor_info, gene.name = colnames(count)[gene_index])
      result[gene_index, 'BinSpect-kmeans'] <- result_binspect_kmeans_temp$measures$p.val

      # run SPARK
      result_spark_temp <- SPARK(abs_expr, sample_info, size.factor = size_factor, gene.name = colnames(count)[gene_index])
      result[gene_index, 'SPARK'] <- result_spark_temp$measures$p.value

      # run SpatialDE
      result_spatialde_temp <- SpatialDE(norm_expr_log_vst, sample_info, gene.name = colnames(count)[gene_index])
      result[gene_index, 'SpatialDE'] <- result_spatialde_temp$measures$p.value
    }


Rank Aggregation
---------------------------------

SMP-Gym provides the results for two rank aggregation methods: GEO and MC2. In boost package, we can conduct the rank aggregation via the function 'rank.aggregation'. This function aggregates rankings from :math:`m` base rankers to generate an aggregated ranking using GEO or MC2 method. Inputs are 1) data with the first column 'gene' records the gene names, 2) K - Sort out top-K genes in each base ranker; 3) method: 'GEO' or 'MC2'; 4) ties.method - a character string specifying how ties are treated.
::
    # create data frame for rank aggregation
    result_df <- data.frame(gene = colnames(count), BOOST_GP =rank(-result[, 'BOOST-GP'], ties.method = "random"))
    result_df$BOOST_Ising <- rank(-result[, 'BOOST-Ising'], ties.method = "random")
    result_df$BinSpect_rank <- rank(result[, 'BinSpect-rank'], ties.method = "random")
    result_df$BinSpect_kmeans <- rank(result[, 'BinSpect-kmeans'], ties.method = "random")
    result_df$SPARK <- rank(result[, 'SPARK'], ties.method = "random")
    result_df$SpatialDE <- rank(result[, 'SpatialDE'], ties.method = "random")

    # rank aggregation
    rank_result <- rank.aggregation(result_df, n_gene, method = 'GEO', ties.method = "random")
    
    print(rank_result)
    ##    gene rank
    ## 6    V6    1
    ## 10  V10    2
    ## 5    V5    3
    ## 3    V3    4
    ## 2    V2    5
    ## 4    V4    6
    ## 9    V9    7
    ## 1    V1    8
    ## 8    V8    9
    ## 7    V7   10

Output is a table with genes and their rank. Gene 'V6' ranks first, which is a SV gene we generate in this simulated data. 

Compute Performace Metrics
--------------------------------

SMP-Gym applies six metrics to comprehensively quantify the performance of SV gene identification for each method. 

* Sensitivity: measure the proportion of correctly identified SV genes across all SV genes in the studied data replicate. Sensitivity ranges from 0 to 1, large sensitivity value corresponds to better classifier model performance. 
* Specificity: measure the proportion of correctly identified non-SV genes across all non-SV genes in the studied data replicate. Specificity ranges from 0 to 1, large specificity value denotes high ability of model to correctly classify non-SV genes. 
* F1-score: harmonic mean of the precision and recall, which simultaneously evaluates the ability of the model to detect true SV genes across all SV genes identified and all true SV genes in the dataset. F1-score ranges from 0 to 1 and higher value indicates better performance. 
* False discovery rate (FDR): calculates the ratio of the number of SV genes detected incorrectly to the total number of SV genes detected, which ranges from 0 to 1. Lower FDR indicates better performance. 
* AUC: Area under the receiver operating characteristic curve (ROC curve), which measures the model performance under all possible thresholds. AUC has range from 0 to 1. When AUC gets closer to 1, model has a better performance. 
* Matthews Correlation Coefficient (MCC): a robust measure to evaluate model performance under imbalance issue, which incorporates all elements in the confusion matrix. MCC ranges from -1 to 1, and when MCC approaches to 1, model has perfect classification ability.

In boost package, we can compute these metrics via the function 'compute.metrics'. Outputs include six performance measurements:
::
    # Compute performance metrics for six methods
    compute.metrics(result[, 'BOOST-GP'], ground_truth, predictor.type = 'BF', threshold = 150)
    ## $Sensitivity
    ## [1] 0
    ## 
    ## $Specificity
    ## [1] 1
    ## 
    ## $F1_score
    ## [1] 0
    ## 
    ## $FDR
    ## [1] 0
    ## 
    ## $AUC
    ## Area under the curve: 1
    ## 
    ## $MCC
    ## [1] 0
    
    compute.metrics(result[, 'BOOST-Ising'], ground_truth, predictor.type = 'BF', threshold = 150)
    ## $Sensitivity
    ## [1] 1
    ## 
    ## $Specificity
    ## [1] 1
    ## 
    ## $F1_score
    ## [1] 1
    ## 
    ## $FDR
    ## [1] 0
    ## 
    ## $AUC
    ## Area under the curve: 1
    ## 
    ## $MCC
    ## [1] 1
    
    compute.metrics(p.adjust(result[, 'BinSpect-rank'], 'BH'), ground_truth, predictor.type = 'p-value', threshold = 0.05)
    ## $Sensitivity
    ## [1] 1
    ## 
    ## $Specificity
    ## [1] 0.8888889
    ## 
    ## $F1_score
    ## [1] 0.6666667
    ## 
    ## $FDR
    ## [1] 0.5
    ## 
    ## $AUC
    ## Area under the curve: 1
    ## 
    ## $MCC
    ## [1] 0.6666667
    
    compute.metrics(p.adjust(result[, 'BinSpect-kmeans'], 'BH'), ground_truth, predictor.type = 'p-value', threshold = 0.05)
    ## $Sensitivity
    ## [1] 1
    ## 
    ## $Specificity
    ## [1] 1
    ## 
    ## $F1_score
    ## [1] 1
    ## 
    ## $FDR
    ## [1] 0
    ## 
    ## $AUC
    ## Area under the curve: 1
    ## 
    ## $MCC
    ## [1] 1
    
    compute.metrics(p.adjust(result[, 'SPARK'], 'BH'), ground_truth, predictor.type = 'p-value', threshold = 0.05)
    ## $Sensitivity
    ## [1] 1
    ## 
    ## $Specificity
    ## [1] 0.8888889
    ## 
    ## $F1_score
    ## [1] 0.6666667
    ## 
    ## $FDR
    ## [1] 0.5
    ## 
    ## $AUC
    ## Area under the curve: 1
    ## 
    ## $MCC
    ## [1] 0.6666667
    
    compute.metrics(p.adjust(result[, 'SpatialDE'], 'BH'), ground_truth, predictor.type = 'p-value', threshold = 0.05)
    ## $Sensitivity
    ## [1] 1
    ## 
    ## $Specificity
    ## [1] 1
    ## 
    ## $F1_score
    ## [1] 1
    ## 
    ## $FDR
    ## [1] 0
    ## 
    ## $AUC
    ## Area under the curve: 1
    ## 
    ## $MCC
    ## [1] 1
