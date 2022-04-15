Load Data
==========================

The current version of boost requires two input data:

* The gene expression count matrix :math:`Y: n \times p` (:math:`n` - number of spots; :math:`p` - number of genes)
* The location information matrix :math:`T: n \times 2`. It includes the :math:`x` and :math:`y` coordinate for each sample point.

Both data should be stored in R matrix format. For gene expression count matrix :math:`Y`, column names should be gene names.

In the following, we conduct the analysis on Mouse Olfactory Bulb data as an example to show the functions of boost package.

Files originally retrieved from the `SMP-Gym website <https://lce.biohpc.swmed.edu/smp_gym/explorer.php>`_ .:


        library(utils)
        library(stats)
        library(grDevices)

        # Download data from SMP-Gym website
        url <- 'https://lce.biohpc.swmed.edu/smp_gym/download/st/mouse_olfactory_bulb_st/processed_data/mouse_olfactory_bulb_replicate_11.zip'
        download.file(url, destfile= "mob.zip", mode = "wb")
        unzip("mob.zip")

        # read data
        count <- read.csv('mouse_olfactory_bulb_replicate_11.count.csv')
        count <- as.matrix(count)
        sample_info <- read.csv('mouse_olfactory_bulb_replicate_11.loc.csv')
        sample_info <- as.matrix(sample_info[, 1:2])
        gene_name <- read.csv('mouse_olfactory_bulb_replicate_11.gene_name.csv')

        colnames(count) <- gene_name$x

        print(dim(count))
        print(count[1:10, 1:10])
        print(dim(sample_info))
        print(head(sample_info))
        
