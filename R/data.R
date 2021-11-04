#' Simulated count matrix from one batch
#'
#' This matrix was simulated from the Negative Binomial(NB) mixture distribution. 
#' More specifically, we simulated 3,000 genes and five group of cells with 400 cells each. 
#' Different cell groups have different mean and dispersion parameters in the NB distribution. 
#' One group of cells was treated as the baseline group with mean parameters
#' generated from truncated Normal with mean and standard deviation 0.2, 0.1 and 
#' dispersion generated from truncated with mean and standard deviation 0.5, 0.1. 
#' For each of the other four groups, 25 genes had mean parameters shifted 
#' from the baseline group with log2 fold change 2.4.
#' Dispersion parameters of feature genes were shifted by 
#' a normalized factor with mean 0 and standard deviation 0.05.
#' The true cell type labels are in 
#' 'one_batch_matrix_celltype", and the cell clustering results generated from 
#' two clustering methods are in "one_batch_matrix_label_df".
#'
#' @format one_batch_matrix is a matrix with 3,000 genes (rows) and 2,000 cells (columns). 
#' The rows are named with g1, g2, ..., g3000 representing gene 1 to gene 3000; 
#' the columns are named with rc1, rc2, ..., rc2000 representing raw count 1 to 
#' raw count 2000. 
#' 
#' 
"one_batch_matrix"

#' Cell type labels of simulated count matrix from one batch
#'
#' This is the true label set corresponds to the cells in the one_batch_matrix. 
#'
#' @format A vector of 2,000 characters indicating the cell types of the 
#' reference cells. For each cell type, there are 400 cells with the name 
#' "type1-5". 
#' 
#' 
"one_batch_matrix_celltype"

#' Clustering labels for simulated one-batch single-cell count matrix
#'
#' This is the label sets correspond to the cells in the one_batch_matrix. 
#' Cells were clustered by K-Means on first 200 PCs and Seurat v3.1.5 
#' with 2000 (default) feature genes. The number of clusters are set to 
#' be {2,3,..., 7}. 
#'
#' @format A data frame of 2,000 rows and 12 columns. Each row represents a cell. 
#' Each column represents a clustering method and the corresponding number of clusters. 
"one_batch_matrix_label_df"




#' Simulated count matrix from two batches
#'
#' 
#' This matrix was simulated from the Negative Binomial(NB) mixture distribution. 
#' More specifically, we first simulated one batch of cells with 1,000 cells in total. 
#' Within the batch, we simulated 3,000 genes and five group of cells with 400 cells each. 
#' One group of cells was treated as the baseline group. The mean parameters
#' are generated from truncated Normal with mean and standard deviation 0.2, 0.1 and 
#' dispersion generated from truncated with mean and standard deviation 0.5, 0.1. 
#' For each of the other four groups, 25 genes had mean parameters shifted 
#' from the baseline group with log2 fold change 2.6.
#' Dispersion parameters of feature genes were shifted by 
#' a normalized factor with mean 0 and standard deviation 0.05.
#' To generate the second batch of cells, we shifted the mean parameters of 
#' randomly selected 20% genes across all cell types by 0.1. 
#'
#' @format A matrix with 3,000 genes (rows) and 2,000 cells (columns). 
#' The cells are from 5 cell type, and each cell type contains 400 cells. 
#' Within each cell type, 200 of cells are from batch 1, and the other cells are 
#' from batch 2. See Note for details of simulation setting. 
#' 
#' 
"two_batch_matrix"



#' Batch labels of simulated count matrix from two batches
#'
#' This is the batch label set corresponds to the cells in the two_batch_matrix. 
#'
#' @format A vector of 2,000 characters indicating the batch labels of cells 
#' in two_batch_matrix. 
#' For each batch, there are 1,000 cells. 
"two_batch_matrix_batch"

#' Cell type labels of simulated count matrix from two batches
#'
#' This is the true label set corresponds to the cells in the two_batch_matrix. 
#'
#' @format A vector of 2,000 characters indicating the cell types of cells 
#' in two_batch_matrix. 
#' For each cell type, there are 400 cells with the name "type1-5". 
"two_batch_matrix_celltype"


#' Clustering labels for simulated two-batch single-cell count matrix
#'
#' This is the label sets correspond to the cells in the two_batch_matrix. 
#' Cells in two batches of the two_batch_matrix was first integrated by Seurat v4. 
#' The dataset after integration was then clustered by K-Means on first 200 PCs 
#' and Seurat with 2000 (default) feature genes.
#' The number of clusters are set to be {2,3,..., 10}. 
#'
#' @format A data frame of 2,000 rows and 18 columns. Each row represents a cell. 
#' Each column represents a clustering method and the corresponding 
#' number of clusters. 
"two_batch_matrix_label_df"

