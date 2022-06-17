Load Data
==========================

The current version of boost requires two input data:

* The gene expression count matrix :math:`Y: n \times p` (:math:`n` - number of spots; :math:`p` - number of genes)
* The location information matrix :math:`T: n \times 2`. It includes the :math:`x` and :math:`y` coordinate for each sample point.

Both data should be stored in R matrix format. For gene expression count matrix :math:`Y`, column names should be gene names.

In the following, we conduct the analysis on Mouse Olfactory Bulb data (replicate 11) as an example to show the functions in boost package.

Files are originally retrieved from the `STAr website <https://lce.biohpc.swmed.edu/star/>`_ . 

First, we load the packages needed for loading data.
::
        library(utils)
        library(stats)
        library(grDevices)

We directly retrieve the data from the website. 
::
        # Download data from STAr website
        url <- 'https://lce.biohpc.swmed.edu/star/download/st/mouse_olfactory_bulb_st/processed_data/mouse_olfactory_bulb_replicate_11.zip'
        download.file(url, destfile= "mob.zip", mode = "wb")
        unzip("mob.zip")

        # read data
        count <- read.csv('mouse_olfactory_bulb_replicate_11.count.csv')
        count <- as.matrix(count)
        sample_info <- read.csv('mouse_olfactory_bulb_replicate_11.loc.csv')
        sample_info <- as.matrix(sample_info[, 1:2])
        gene_name <- read.csv('mouse_olfactory_bulb_replicate_11.gene_name.csv')
        
        # set column names of count as gene names 
        colnames(count) <- gene_name$x

        print(dim(count))
        ## [1]   262 16218
        print(dim(sample_info))
        ## [1] 262   2
        
The Mouse Olfactory Bulb data (replicate 11) has 262 sample points and 16218 genes. The 'sample_info' table records the x and y coordinates of all sample points. 
::
        print(count[1:10, 1:10])
        ##       Nrf1 Zbtb5 Ccnl1 Lrrfip1 Bbs1 Lix1 Whrn Ate1 Ubac1 Rab34
        ##  [1,]    1     1     1       2    1    2    1    1     2     1
        ##  [2,]    0     0     3       2    2    7    0    2     3     2
        ##  [3,]    0     1     1       0    0    0    1    0     0     0
        ##  [4,]    1     0     1       0    4    6    1    4     3     1
        ##  [5,]    0     0     0       3    0    2    0    2     4     0
        ##  [6,]    0     0     0       5    0    1    1    0     2     0
        ##  [7,]    0     0     2       1    3    4    0    3     2     0
        ##  [8,]    1     1     3       0    1    3    0    4     0     0
        ##  [9,]    0     0     2       0    0    4    0    1     0     2
        ## [10,]    0     0     2       3    0    2    0    1     0     0
       
       print(head(sample_info))
       ##           X      Y
       ## [1,] 16.920  9.015
       ## [2,] 16.945 11.075
       ## [3,] 16.970 10.118
       ## [4,] 16.939 12.132
       ## [5,] 16.949 13.055
       ## [6,] 16.942 15.088


