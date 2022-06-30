## -----------------------------------------------------------------------------
##          External function -- extract_sce
## -----------------------------------------------------------------------------
#' Extract information from SingleCellExperiment object
#' 
#' This function extracts count matrix, cell labels, and batch labels from 
#' SingleCellExperiment object. 
#'
#' @param sce_obj SingleCellExperiment object
#' @param info_type A character indicating the type of information to extract. 
#' "count" represents the count matrix; "label" represents cell labels; "batch"
#' represent batch labels. 
#' @param info_slot Slot(s) in assays that stores the count matrix, or slot in 
#' colData that saves cluster / batch labels.
#' @importFrom methods is
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment colData
#' @return A count matrix if info_type = "count"; A data frame if 
#' info_type is "label"; A vector of characters if info_type is "batch".
#'
#' @examples
#' ## Simulate count matrix, batch, and cell type labels
#' ng <- 100; nc <- 100
#' set.seed(1)
#' X <- cbind(
#' 	matrix(
#' 		c(rnbinom(ng*nc/4, size = 1, mu = 0.1),
#' 			rnbinom(ng*nc/4, size = 1, mu = 0.5)),
#' 		nrow = ng, 
#' 		byrow = TRUE),
#' 	matrix(
#' 		c(rnbinom(ng*nc/4, size = 1, mu = 1),
#' 			rnbinom(ng*nc/4, size = 1, mu = 0.5)),
#' 		nrow = ng, 
#' 		byrow = TRUE))
#' colnames(X) <- paste0('c', seq_len(nc))
#' rownames(X) <- paste0('g', seq_len(ng))
#' Batches <- rep(seq_len(2), nc/2)
#' Method1_k2 <- rep(seq_len(2), c(nc/2,nc/2))
#' Method1_k3 <- sample(seq_len(3), nc, replace = TRUE)
#' 
#' library(SingleCellExperiment)
#' sim_sce <- SingleCellExperiment(
#'   list(count = X),
#'   colData = data.frame(
#'     Cell_name = colnames(X),
#' 		 Batch = Batches,
#'   	 Method1_k2 = Method1_k2,
#'   	 Method1_k3 = Method1_k3),
#' 	rowData = data.frame(Gene_name = rownames(X)))
#' 
#' extract_sce(sim_sce, "count", "count")
#' extract_sce(sim_sce, "batch", "Batch")
#' extract_sce(sim_sce, "label", c("Method1_k2", "Method1_k3"))
#' @export

extract_sce <- function(sce_obj, info_type, info_slot){
	if(!is(sce_obj, "SingleCellExperiment")){
		stop("Input must be a SingleCellExperiment object.")
	}
	if(info_type == "count"){
		return(as.matrix(assays(sce_obj)[[info_slot]]))
	} 
	if(info_type == "label"){
		return(as.data.frame(colData(sce_obj)[, info_slot]))
	} 
	if(info_type == "batch"){
			return(as.vector(colData(sce_obj)[, info_slot]))
	}
	return("Please input one of options for info_type: count, label, or batch.")
}


## -----------------------------------------------------------------------------
##          External function -- extract_seurat
## -----------------------------------------------------------------------------


#' Extract information from Seurat object
#' 
#' This function extracts count matrix, cell labels, and batch labels from 
#' Seurat object. 
#'
#' @param seurat_obj Seurat object
#' @param info_type A character indicating the type of information to extract. 
#' "count" represents the count matrix; "label" represents cell labels; "batch"
#' represent batch labels. 
#' @param info_slot Slot(s) in assays$RNA that stores the count matrix, or slot in 
#' meta.data that saves cluster / batch labels.
#' @importFrom Seurat CreateSeuratObject
#' @return A count matrix if info_type is "count"; A data frame if 
#' info_type = "label"; A vector of characters if info_type is "batch".
#'
#' @examples
#' ## Simulate count matrix, batch, and cell type labels
#' ng <- 100; nc <- 100
#' set.seed(1)
#' X <- cbind(
#' 	matrix(
#' 		c(rnbinom(ng*nc/4, size = 1, mu = 0.1),
#' 			rnbinom(ng*nc/4, size = 1, mu = 0.5)),
#' 		nrow = ng, 
#' 		byrow = TRUE),
#' 	matrix(
#' 		c(rnbinom(ng*nc/4, size = 1, mu = 1),
#' 			rnbinom(ng*nc/4, size = 1, mu = 0.5)),
#' 		nrow = ng, 
#' 		byrow = TRUE))
#' colnames(X) <- paste0('c', seq_len(nc))
#' rownames(X) <- paste0('g', seq_len(ng))
#' Batches <- rep(seq_len(2), nc/2)
#' Method1_k2 <- rep(seq_len(2), c(nc/2,nc/2))
#' Method1_k3 <- sample(seq_len(3), nc, replace = TRUE)
#' 
#' library(Seurat)
#' library(SeuratObject)
#' sim_seurat <- CreateSeuratObject(counts = as.data.frame(X))
#' sim_seurat <- AddMetaData(sim_seurat, colnames(X), "Cell_name")
#' sim_seurat <- AddMetaData(sim_seurat, Batches, "Batch")
#' sim_seurat <- AddMetaData(sim_seurat, Method1_k2, "Method1_k2")
#' sim_seurat <- AddMetaData(sim_seurat, Method1_k3, "Method1_k3")
#' 	
#' extract_seurat(sim_seurat, "count", "counts")
#' extract_seurat(sim_seurat, "label", c("Method1_k2", "Method1_k3"))
#' extract_seurat(sim_seurat, "batch", "Batch")
#' @export


extract_seurat <- function(seurat_obj, info_type, info_slot){
	if(!is(seurat_obj, "Seurat")){
		stop("Input must be a Seurat object.")
	}
	if(info_type == "count"){
		return(as.matrix(eval(parse(
			text = paste0("seurat_obj@assays$RNA@", info_slot)))))
	} 
	if(info_type == "label"){
		return(as.data.frame(seurat_obj@meta.data[, info_slot]))
	} 
	if(info_type == "batch"){
			return(as.vector(seurat_obj@meta.data[, info_slot]))
	}
	return("Please input one of options for info_type: count, label, or batch.")
}




