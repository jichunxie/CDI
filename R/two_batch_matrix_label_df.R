#' Clustering labels for simulated two-batch single-cell count matrix
#'
#' This is the label sets correspond to the cells in the two_batch_matrix. 
#' Cells in two batches of the two_batch_matrix was first integrated by Seurat v4. 
#' The dataset after integration was then clustered by K-Means on first 200 PCs 
#' and Seurat with 2000 (default) feature genes.
#' The number of clusters are set to be {2,3,..., 10}. 
#'
#' @format A data frame of 2,000 rows and 18 columns. Each row represents a cell. 
#' Each column represents a clustering method and the corresponding number of clusters. 
"two_batch_matrix_label_df"
