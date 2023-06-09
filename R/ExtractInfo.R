## -----------------------------------------------------------------------------
##          Internal function -- extract_count
## -----------------------------------------------------------------------------
#' Extract count matrix from input X
#' 
#' This is an internal function extracting count matrix from input X. 
#'  
#' @param X The class of X can be "matrix", "Seurat" object, or "SingleCellExperiment" object. 
#' If X is a matrix, it should be a UMI count matrix where each row represents a gene, and 
#' each column represents a cell. 
#' If X is a Seurat or SingleCellExperiment object, the UMI count matrix should be 
#' stored in count_slot of assays.
#' 
#' @param count_slot A string indicating the location of raw UMI count. 
#' For Seurat object, it is a slot in "RNA" of "assays"; 
#' For SingleCellExperiment object, it is a slot in "assays". 
#' Each row represents a gene, and each column represents a cell. 
#' The genes should be those before feature gene selection.
#' 
#' @return a matrix of extracted UMI counts where each row represents a gene and each column represents a cell. 
#' 
#' @importFrom methods is
#' @importFrom SummarizedExperiment assays
#' @importFrom Seurat GetAssayData
#' @noRd
extract_count <- function(X, count_slot = NULL){
	if(!is(X, "Seurat") & !is(X, "matrix") & !is(X, "SingleCellExperiment")){
		stop("The type of X need to be one of the following: matrix, Seurat object, or SingleCellExperiment object.")
	}
	# (1) Input Seurat object, extract batch and count with selected genes
	if(is(X, "Seurat")){
		if(is.null(count_slot)) stop("count_slot is not specified for Seurat object X.")
		gcmat <- tryCatch(
			GetAssayData(object = X, slot = count_slot), 
			error = function(e) stop("count_slot is not in assays of Seurat object X."))
	 # (2) Input SingleCellExperiment object, extract batch and count with selected genes
	 # Note: unlike Seurat, sce doesn't give an error if is.null(assays(X)[[count_slot]])
	 # That is why try catch is not used here
	} else if(is(X, "SingleCellExperiment")){
		if(is.null(count_slot)) stop("count_slot is not specified for SingleCellExperiment object X.")
		if(!is.null(assays(X)[[count_slot]])){
			gcmat <- assays(X)[[count_slot]]
		} else stop("count_slot is not in assays of SingleCellExperiment object X.")
	# (3) Input matrix, rename X as gcmat
	} else {
		gcmat <- X
	}
	return(as.matrix(gcmat))
}




## -----------------------------------------------------------------------------
##          Internal function -- extract_batch
## -----------------------------------------------------------------------------
#' Extract batch information from input X
#' 
#' This is an internal function extracting batch labels from input X. 
#'  
#' @param X The class of X can be "matrix", "Seurat" object, or "SingleCellExperiment" object. 
#' If X is a matrix, it should be a UMI count matrix where each row represents a gene, and 
#' each column represents a cell. 
#' If X is a Seurat object, the batch labels should be stored in batch_slot of "meta.data". 
#' If X is a SingleCellExperiment object, the batch labels should be stored in batch_slot of "colData". 
#' 
#' @param batch_label A vector of characters indicating the batch labels of the cells. 
#' The length of batch_label should be the same as the number of columns 
#' in the count matrix.The default value is NULL indicating that there is no batch labels 
#' or the batch label should be found in the batch_slot of object X. 
#' 
#' @param batch_slot A string indicating the location of batch labels of cells.
#' For Seurat object, it is a slot in meta.data;
#' For SingleCellExperiment object, it is a slot in "colData". 
#' The default value is NULL indicating that there is no batch information available. 
#' 
#' @return a vector of extracted batch labels "batch_label". 
#' 
#' @importFrom methods is
#' @importFrom SummarizedExperiment colData
#' @importFrom Seurat FetchData
#' @noRd
extract_batch <- function(X, 
	batch_label = NULL, 
	batch_slot = NULL){
	if(!is(X, "Seurat") & !is(X, "matrix") & !is(X, "SingleCellExperiment")){
		stop("The type of X need to be one of the following: matrix, Seurat object, or SingleCellExperiment object.")
	}
	# (1) Input Seurat object, extract batch and count with selected genes
	if(is(X, "Seurat")){
	  # if batch_slot is NULL, there are two cases:
		# Either there is no batch label, or the batch label is specified in batch_label.
		# In both cases we don't need to change batch_label.  
		# So only need to consider when batch_slot is not NULL. 
		if(!is.null(batch_slot)){
			batch_label <- tryCatch(
			FetchData(object = X, vars = batch_slot),
			error = function(e) stop("batch_slot is not in meta.data of Seurat object X."))
		}
	 # (2) Input SingleCellExperiment object, extract batch and count with selected genes
	 # Note: unlike Seurat, sce doesn't give an error if is.null(assays(X)[[count_slot]])
	 # That is why try catch is not used here
	} else if(is(X, "SingleCellExperiment")){
		if(!is.null(batch_slot)){
			if(!is.null(colData(X)[batch_slot])){
				batch_label <- colData(X)[[batch_slot]]
			} else stop("batch_slot is not in SingleCellExperiment object X.")
		}
	}
	# (3) Input matrix: directly return batch_label
	return(as.vector(unlist(batch_label)))
}

