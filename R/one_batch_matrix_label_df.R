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
