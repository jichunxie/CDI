## -----------------------------------------------------------------------------
##          Internal function -- one_batch_feature_gene_rank
## -----------------------------------------------------------------------------
##  rank genes for one batch
one_batch_feature_gene_rank <- function(
	gcmat = NULL, 
	method = "wds", 
	nfeature = 500,
	zp_threshold = 0.95){
	## count genes and cells
	ng <- nrow(gcmat)
  nc <- ncol(gcmat)
	## Gene selection via wds
	if(method == "wds"){
		zp <- rowMeans(gcmat == 0)
    zp_indx <- c(seq_len(ng))[zp < zp_threshold]
    mu_g <- rowMeans(gcmat); var_g <- rowVars(gcmat)
    phi_est <- (var_g-mu_g)/mu_g^2
    phi_df <- data.frame(indx = seq_len(ng), phi_est = phi_est)
    phi_df[zp_indx, "phi_rank"] <- rank(-phi_est[zp_indx])
    ## Assign genes with NA rank to a large value to remove NA
    phi_df$phi_rank[is.na(phi_df$phi_rank)] <- ng + 1000
    return(phi_df$phi_rank)
	}
	## Gene selection via vst
  if(method == "vst"){
		if(is.null(colnames(gcmat))){
			colnames(gcmat) <- paste0("c", seq_len(nc))
			}
		if(is.null(rownames(gcmat))){
			rownames(gcmat) <- paste0("g", seq_len(ng))
			}
		gene_indx_name <- data.frame(indx = c(seq_len(ng)), grank = rep(ng, ng))
		rownames(gene_indx_name) <- rownames(gcmat)
		Seurat_obj <- CreateSeuratObject(counts = as.data.frame(gcmat))
		Seurat_obj <- AddMetaData(Seurat_obj, colnames(gcmat), "Cell_name")
  	Seurat_obj <- NormalizeData(Seurat_obj, verbose = FALSE)
    Seurat_obj <- FindVariableFeatures(
    	Seurat_obj, 
    	selection.method = "vst", 
    	nfeatures = nfeature, 
    	verbose = FALSE)
    tmp_df <- data.frame(
    	gname = Seurat_obj@assays[["RNA"]]@var.features, 
    	grank = seq_len(nfeature))
    gene_indx_name[as.character(tmp_df$gname), "grank"] <- tmp_df$grank
    return(gene_indx_name$grank)
  }
}


## -----------------------------------------------------------------------------
##        External function -- feature_gene_selection
## -----------------------------------------------------------------------------
#' Select feature genes
#' 
#' This function selects a subset of feature genes that are expected to express 
#' differently across cell types before calculating CDI.
#'
#'
#' @param X The class of X can be "matrix", "Seurat" object, or "SingleCellExperiment" object. 
#' If X is a matrix, it should be a raw UMI count matrix where each row represents a gene, and 
#' each column represents a cell. The genes should only included feature genes 
#' (that are selected by feature_gene_selection function).
#' If X is a Seurat object, users need to specify where the count matrix, feature gene, and
#' batch labels are stored in seurat_count_slot, seurat_feature_slot, and seurat_batch_slot, 
#' respectively. 
#' If X is a SingleCellExperiment object, users need to specify where the count matrix, 
#' feature gene, and batch labels are stored in sce_count_slot, sce_feature_slot, 
#' and sce_batch_slot, respectively.
#' 
#' @param batch_label A vector of characters indicating the batch labels of the cells. 
#' The length of batch_label should be the same as the number of columns 
#' in the count matrix.
#' 
#' @param seurat_count_slot A string indicating the location of raw UMI count. The 
#' default value is "@assays$RNA@counts". Each row represents a gene, and each column
#' represents a cell. The genes should be those before feature gene selection.
#' 
#' @param seurat_batch_slot A string indicating the location of batch labels of 
#' cells. The default value is NULL. 
#' 
#' @param sce_count_slot A string indicating the location of raw UMI count in 
#' sce object. The default value is "count". Each row represents a gene, and each column 
#' represents a cell. The genes should be those before feature gene selection.
#' 
#' @param sce_batch_slot A string indicating the location of batch labels of 
#' cells in sce object. The default value is NULL.
#' 
#' @param method A character indicating the method used to select 
#' feature genes. "wds" (default) represent the working dispersion score 
#' proposed in our study; "vst" is the default method for feature gene selection 
#' in FindVariableFeatures function of Seurat package.
#' 
#' @param nfeature An integer indicating the number of features to select.
#' 
#' @param zp_threshold A number of zero proportion threshold. The range 
#' should be (0,1). Genes with 
#' zero proportion greater than this value will be excluded before feature selection.
#' 
#' @importFrom matrixStats rowVars rowMins
#' @importFrom SeuratObject AddMetaData
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures
#' @importFrom SingleCellExperiment SingleCellExperiment rowData colData
#' @importFrom SummarizedExperiment assays
#' @return A vector of indices of selected feature genes corresponding to the 
#' row indices of input count matrix.
#'
#' @examples
#' ## Simulate a matrix of 100 rows (genes), where the first 50 genes have
#' ## different mean expression levels.
#' ## Apply feature_gene_selection to select genes
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
#' 
#' ## Input: matrix
#' feature_gene_selection(
#' 	X = X,
#' 	batch_label = Batches,
#' 	nfeature = 20,
#' 	zp_threshold = 0.95)
#' 
#' 
#' ## Input: SingleCellExperiment object
#' library(SingleCellExperiment)
#' sim_sce <- SingleCellExperiment(
#'   list(count = X),
#'   colData = data.frame(
#'   	Cell_name = colnames(X),
#' 		batch = Batches),
#' 	rowData = data.frame(gene_name = rownames(X)))
#' selected_genes <- feature_gene_selection(
#'   X = sim_sce,
#'   sce_batch_slot = "batch",
#'   nfeature = 20,
#'   zp_threshold = 0.95)
#' # Record feature genes in sim_sce
#' rowData(sim_sce)["var.features"] <- FALSE
#' rowData(sim_sce)["var.features"][selected_genes,'var.features'] <- TRUE
#' 
#' ## Input: Seurat object
#' library(Seurat)
#' library(SeuratObject)
#' sim_seurat <- CreateSeuratObject(counts = as.data.frame(X))
#' sim_seurat <- AddMetaData(sim_seurat, colnames(X), "Cell_name")
#' sim_seurat <- AddMetaData(sim_seurat, Batches, "batch")
#' selected_genes <- feature_gene_selection(
#' 	X = sim_seurat,
#' 	seurat_batch_slot = "@meta.data$batch",
#' 	nfeature = 20,
#' 	zp_threshold = 0.95)
#' # Record feature genes in sim_seurat
#' sim_seurat@assays$RNA@var.features <- selected_genes
#'                      
#' @references Stuart and Butler et al. Comprehensive Integration of 
#' Single-Cell Data. Cell (2019) [Seurat V3]
#' 
#' @export
#' 
feature_gene_selection <- function(
	X = NULL, 
	batch_label = NULL,
	seurat_count_slot = "@assays$RNA@counts", 
	seurat_batch_slot = NULL,
	sce_count_slot = "count",
	sce_batch_slot = NULL,
	method = "wds",
	nfeature = 500,
	zp_threshold = 0.95){
	X_class = class(X)[1]
	if(!(X_class %in% c("matrix", "Seurat", "SingleCellExperiment"))){
		stop("The type of X need to be one of the following: matrix, Seurat object, or SingleCellExperiment object.")
	}
	if(X_class == "Seurat"){
		  # Seurat generates error if the seurat_count_slot doesn't exist
	    gcmat <- tryCatch(
			as.matrix(eval(parse(text = paste0("X", seurat_count_slot)))), 
			error = function(e) stop("seurat_count_slot is not in Seurat object X."))
	    if(!is.null(seurat_batch_slot)){
	    	batch_label <- tryCatch(
				eval(parse(text = paste0("X", seurat_batch_slot))), 
				error = function(e) stop("seurat_batch_slot is not in Seurat object X."))
	    }
	} else if(X_class == "SingleCellExperiment"){
		if(!is.null(assays(X)[[sce_count_slot]])){
			gcmat <- assays(X)[[sce_count_slot]]
		} else{
			stop("sce_count_slot is not in SingleCellExperiment object X.")
		}
		# if sce_batch_slot is NULL, there are two cases:
		# 1. there is no batch label
		# 2. the batch label is specified in batch_label
		if(!is.null(sce_batch_slot)){
			if(!is.null(colData(X)[sce_batch_slot])){
				batch_label <- colData(X)[[sce_batch_slot]]
			} else{
				stop("sce_batch_slot is not in SingleCellExperiment object X.")
			}
		}
	} else{
		# X_class == "matrix"
		assign("gcmat", X)
	}
	ng <- nrow(gcmat)
	## One batch
	if(is.null(batch_label) | (length(unique(batch_label)) == 1)){
  	feature_bool = (one_batch_feature_gene_rank(
  		gcmat = gcmat, 
  		method = method, 
  		nfeature = nfeature, 
  		zp_threshold = zp_threshold) <= nfeature)
    return(c(seq_len(ng))[feature_bool])
  } else{
  	## Multi-batches
    ## split gcmat columns according to the batch labels
    mat_list <- lapply(
    	split(seq_along(batch_label), batch_label), 
    	function(indx, mat) mat[,indx], 
    	mat = gcmat)
    ## rank genes in each matrix
    gene_rank_list <- lapply(
    	mat_list, 
    	one_batch_feature_gene_rank, 
    	method = method, 
    	nfeature = nfeature, 
    	zp_threshold = zp_threshold)
    ## return top nfeature genes
    gene_rank_mat <- do.call(cbind, gene_rank_list)
    gene_min_rank <- rowMins(gene_rank_mat)
    names(gene_min_rank) <- seq_len(nrow(gcmat))
    as.numeric(names(gene_min_rank)[order(gene_min_rank)[seq_len(nfeature)]])
    return(sort(as.numeric(names(gene_min_rank)[order(gene_min_rank)[seq_len(nfeature)]])))
  }
  
}
