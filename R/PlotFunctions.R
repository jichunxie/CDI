#' Visualize CDI values via a lineplot
#' 
#' This function visualize CDI outputs with a lineplot. 
#' The x-axis is for the number of clusters. 
#' The y-axis is for the CDI values. 
#' Different colors represent different clustering methods.
#' The red triangle represents the optimal clustering result corresponding to the smallest CDI value. 
#' The purple star represents the CDI value for the benchmark (main) cell type label set. 
#' he brown star represents the CDI value for the benchmark sub-type label set. 
#'
#' @param cdi_dataframe A data frame of CDI values. Each row represents one clustering 
#' method and the number of clusters combination. The columns include "Cluster_method", 
#' "CDI_AIC", "CDI_BIC", and "N_cluster". 
#' @param cdi_type A string indication the type of CDI. It can be either "CDI_AIC" or "CDI_BIC".
#' @param benchmark_celltype_cdi A list of the output from the CalculateCDI function 
#' on the benchmark cell type label set. Default is null. 
#' @param benchmark_celltype_ncluster A number indicating the number of cell types 
#' in the benchmark cell type label set. Default is null. 
#' @param benchmark_maintype_cdi A list of the output from the CalculateCDI function 
#' on the benchmark main type label set. Default is null. 
#' @param benchmark_maintype_ncluster A number indicating the number of cell types 
#' in the benchmark main type label set. Default is null. 
#' @param clustering_method A vector of characters indicating the corresponding clustering 
#' method for each label set. If this is provided, the 
#' lineplot will be colored differently by different clustering methods. 
#' This color scheme can also be obtained if the column names in 
#' candidate_label_dataframe of CalculateCDI are provided 
#' with the form "[ClusteringMethod]_k[NumberOfClusters]" such as "KMeans_k5".  
#' @param show_axis_names A bool value indicating whether the axis names should be shown or not in the plot. 
#' @param show_method_legend A bool value indicating whether the legend of methods should be shown or not in the plot. 
#' @importFrom grDevices rgb
#' @importFrom ggsci scale_color_d3
#' @importFrom ggplot2 ggplot geom_point geom_line theme theme_bw aes_string element_blank
#' @return A ggplot object.
#' @examples
#' ## This data frame is generated here for demonstration only:
#' CDI_return <- data.frame(Cluster_method = c("HC", "HC",  "HC", "KMeans", "KMeans", "KMeans"),
#' N_cluster = c(2,4,6,2,4,6),
#' CDI_AIC = c(150, 200, 250, 220, 160, 180),
#' CDI_BIC = c(170, 210, 280, 250, 180, 200))
#' CDI_lineplot(CDI_return, cdi_type = "CDI_AIC")
#' benchmark_cdi_return <- list(CDI_AIC = 150, CDI_BIC = 170)
#' benchmark_ncelltype <- 3
#' CDI_lineplot(CDI_return, 
#'             cdi_type = "CDI_AIC", 
#'             benchmark_celltype_cdi = benchmark_cdi_return, 
#'             benchmark_celltype_ncluster = benchmark_ncelltype)
#' @export

# provide a vector of methods
CDI_lineplot <- function(cdi_dataframe, 
	cdi_type, 
	benchmark_celltype_cdi = NULL, 
	benchmark_celltype_ncluster = NULL, 
	benchmark_maintype_cdi = NULL, 
	benchmark_maintype_ncluster = NULL,
	clustering_method = NULL,
	show_axis_names = TRUE, 
	show_method_legend = TRUE){
  min_indx <- which.min(cdi_dataframe[, cdi_type])
  if(!is.null(clustering_method)){
    cdi_dataframe["Cluster_method"] <- clustering_method
  } else if(sum(is.na(cdi_dataframe$Cluster_method))){
    cdi_dataframe[is.na(cdi_dataframe$Cluster_method), "Cluster_method"] <- "Unknown_method"
  } 
  p1 <- ggplot(cdi_dataframe, aes_string(x = "N_cluster", y = cdi_type, group = "Cluster_method")) + 
    geom_point(aes_string(color = "Cluster_method"), size = 4, alpha = 0.7) + 
    geom_line(aes_string(color = "Cluster_method"), size = 2, alpha = 0.7) + 
    geom_point(data = cdi_dataframe[min_indx, ], 
    	colour = grDevices::rgb(210,85,62, max = 255), 
    	shape = 17, size = 4) + labs(x = "Number of clusters") + 
  	scale_color_d3() + theme_bw()
  if(show_axis_names == FALSE){
    p1 <- p1 + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  }
  if(show_method_legend == FALSE){
    p1 <- p1 + theme(legend.position = "none")
  }
  if((missing(benchmark_celltype_cdi) | missing(benchmark_celltype_ncluster)) & (missing(benchmark_maintype_cdi) | missing(benchmark_maintype_ncluster))){
    return(p1)
  } 
  if(cdi_type == "CDI_AIC"){
    if(!is.null(benchmark_celltype_cdi)){
      p1 <- p1 + geom_point(data = data.frame(
      	N_cluster = benchmark_celltype_ncluster, 
      	CDI_AIC = benchmark_celltype_cdi$CDI_AIC, 
      	Cluster_method = "Benchmark"), 
      	colour = "purple3",
      	shape = "*", 
      	size = 15) 
    }
    if(!is.null(benchmark_maintype_cdi)){
      p1 <- p1 + geom_point(data = data.frame(
      	N_cluster = benchmark_maintype_ncluster, 
      	CDI_AIC = benchmark_celltype_cdi$CDI_AIC, 
      	Cluster_method = "Benchmark"), 
      	colour="#7E6148FF", 
      	shape = "*", 
      	size = 15)
    }
    
  } 
  if(cdi_type == "CDI_BIC"){
    if(!is.null(benchmark_celltype_cdi)){
      p1 <- p1 + geom_point(data = data.frame(
      	N_cluster = benchmark_celltype_ncluster, 
      	CDI_BIC = benchmark_celltype_cdi$CDI_BIC,
      	Cluster_method = "Benchmark"),
      	colour="purple3", 
      	shape = "*", 
      	size = 15)
    }
    if(!is.null(benchmark_maintype_cdi)){
      p1 <- p1 + geom_point(data = data.frame(
      	N_cluster = benchmark_maintype_ncluster,
      	CDI_BIC = benchmark_maintype_cdi$CDI_BIC,
      	Cluster_method = "Benchmark"), 
      	colour="#7E6148FF", 
      	shape = "*", 
      	size = 15)
    }
    
  }
  
   return(p1)
}




#' Heatmap of contingency table 
#' 
#' 
#' If the benchmark cell type label set is available, we can also compare one 
#' candidate label set 
#' (e.g. the optimal clustering label set selected by CDI) with the benchmark 
#' cell type labels. 
#' Here we provide the heatmap of contingency table for comparison.
#' Each row represents a cell type in the benchmark label set, 
#' and each column represents a cluster in the clustering label set. 
#' Each rectangle is color-scaled by the proportion of the cells in the given 
#' cluster coming from the benchmark types. 
#' Each column sums to 1. 

#'
#' @param candidate_label A vector of characters indicating the candidate 
#' clustering label set of cells.
#' @param benchmark_label A vector of characters indicating the benchmark 
#' cell type label set of cells.
#' @param proportion_size A number indicating the label size of proportion 
#' values inside each rectangle. 
#' The label will not be shown if this parameter is set to be FALSE.  
#' @param show_axis_names A bool value indicating whether the axis names 
#' should be shown or not in the plot. 
#' @param show_legend A bool value indicating whether the legend of methods 
#' should be shown or not in the plot.
#' @param rename_candidate_clusters A bool value indicating whether the  
#' candidate cluster names will be changed or not. If true, the candidate
#' cluster names will be changed to the candidate_cluster_names.
#' @param candidate_cluster_names A vector of characters indicating the 
#' candidate cluster names. The order of names will be the same as the order
#' without customizing the candidate cluster names. For example, when 
#' rename_candidate_clusters = FALSE, and the output figure has x axis label
#' (cluster1, cluster0). If rename_candidate_clusters = TRUE, and 
#' candidate_cluster_names is (c1, c2). That means cluster1 -> c1, and cluster0 -> c0.
#' @importFrom ggplot2 ggplot geom_tile geom_text theme theme_classic 
#' aes aes_string element_blank element_line labs scale_fill_gradientn 
#' scale_x_discrete element_text
#' @importFrom grDevices rgb
#' @importFrom reshape2 melt
#' @return A ggplot object.
#' @examples
#' contingency_heatmap(
#' 	benchmark_label = c(rep("type_b", 45), rep("type_a", 145), rep("type_c", 10)),
#' 	candidate_label = paste0('cluster', c(rep(0,10), rep(1, 180), rep(2, 10))))
#' contingency_heatmap(
#' 	benchmark_label = c(rep("type_b", 45), rep("type_a", 145), rep("type_c", 10)),
#' 	candidate_label = paste0('cluster', c(rep(0,10), rep(1, 180), rep(2, 10))),
#' 	rename_candidate_clusters = TRUE, 
#' 	candidate_cluster_names = c('1', '2', '3'))
#' @export

contingency_heatmap <- function(benchmark_label, 
	candidate_label, 
	proportion_size = 4,
	show_axis_names = TRUE,
	show_legend = TRUE,
	rename_candidate_clusters = FALSE,
	candidate_cluster_names = NULL){
  ## create contingency table
  count_mat <- as.matrix(table(data.frame(benchmark_label, candidate_label)))
  ## change from counts to proportions
  col_sum_matrix <- matrix(rep(colSums(count_mat), nrow(count_mat)), nrow = nrow(count_mat), byrow = TRUE)
  ## each row is a candidate class, and each row sums to 1.
  prop_mat <- t(round(count_mat / col_sum_matrix, 2))
  longData <- reshape2::melt(prop_mat)
  ## To avoid the problem of undifined variable
	PropValue = NULL
  colnames(longData)[3] <- "PropValue"
  ## remove zero proportions
  longData <- longData[longData$PropValue != 0, ]
  ## order benchmark cell types according to the composition of candidate cell types
  order_level <- NULL
  for(cell_type in levels(as.factor(longData$benchmark_label))){
    large_prop_indx <- names(which(prop_mat[,cell_type] >= 0.2))
    large_prop_sort <- large_prop_indx[order(prop_mat[large_prop_indx, cell_type], decreasing = TRUE)]
    distinct_indx <- setdiff(large_prop_sort, order_level)
    order_level <- c(order_level, distinct_indx)
  }
  p1 <- ggplot(longData, 
  	aes(x = factor(candidate_label, levels = order_level), y = factor(benchmark_label), fill = PropValue)) + 
    scale_fill_gradientn(limits = c(0,1), 
    	colours = c(rgb(204,204,204, maxColorValue = 255), rgb(25,150,125, maxColorValue = 255))) + 
    labs(x = "Candidate label", y = "Benchmark label", fill = "Proportion") + 
    theme_classic() + geom_tile() + 
    theme(axis.text = element_text(size=10, angle=0, vjust=0.3), axis.line = element_line(size = 1))
  if(rename_candidate_clusters){
  	p1 <- p1 + scale_x_discrete(labels = candidate_cluster_names)
  }
  if(proportion_size){
    p1 <- p1 + geom_text(label = round(longData$PropValue, 2), 
    	size = proportion_size) 
  }
  if(show_axis_names == FALSE){
    p1 <- p1 + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  }
  if(show_legend == FALSE){
    p1 <- p1 + theme(legend.position = "none")
  }
  return(p1)
}



