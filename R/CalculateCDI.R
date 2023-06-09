## -----------------------------------------------------------------------------
##          Internal functions
## -----------------------------------------------------------------------------

## Calculate negative sum of NB log-likelihoods
#' @param input_parameter values of (mu, r) in NB distribution
#' @param x_vec a vector of x values as observations from NB distribution
#' @param sc_vec a vector of size factor s_c in NB(s_c*mu, r)
#' @importFrom stats dnbinom
#' @noRd
neg_nb_logsum <- function(input_parameter, x_vec, sc_vec){
  mu <- input_parameter[1]
  r <- input_parameter[2]
  return(-sum(dnbinom(x_vec, size = r, mu = sc_vec * mu, log = TRUE)))
}


## Calculate MLE of NB distribution with size factor 
## per cell-type per gene
nb_size_mle <- function(x_vec, sc_vec){
  if(sum(x_vec) == 0){
    return(list(param = c(0.001, 0.2), 
                value = neg_nb_logsum(c(0.001, 0.2), x_vec, sc_vec)))
  }
  avg <- mean(x_vec); s2 <- var(x_vec)
  init_mu <- avg
  tmp_init_r <- avg^2/(s2 - avg)
  init_r_0 <- ifelse(is.na(tmp_init_r)|(tmp_init_r < 0.1), 0.1, tmp_init_r)
  init_r <- ifelse(init_r_0 > 50, 50, init_r_0)
  nb_est <- nlminb(c(init_mu, init_r), 
                  objective = neg_nb_logsum, 
                  gradient = NULL, 
                  lower = c(1e-6, 1e-6), 
                  upper = c(1e3, 1e6),
                  x_vec = x_vec,
                  sc_vec = sc_vec)
  return(list(param = nb_est$par, value = nb_est$objective))
}



## Calculate the sum of negative likelihood
## all cell-type per gene per batch
single_batch_one_gene_likelihood <- function(
	gvec, 
	candidate_label, 
	ncluster, 
	cell_size_factor){
  gvec_list <- split(gvec, f = factor(candidate_label))
  sc_list <- split(cell_size_factor, f = factor(candidate_label))
  neg_llk_g <- 0
  for(k in seq_len(ncluster)){
    gvec_ct <- gvec_list[[k]]
    sc_ct <- sc_list[[k]]
    neg_llk_g  <- neg_llk_g  + nb_size_mle(gvec_ct, sc_ct)$value
  }
  return(neg_llk_g)
  
}



## Calculate the sum of negative likelihood
## all cell-types per gene all batches
multi_batch_one_gene_likelihood <- function(gvec,
	candidate_label, 
	ncluster,
	size_ct_list,
	batch_ct_list,
	lrt_pval_threshold = 0.01){
  gvec_ct_list <- split(gvec, f = factor(candidate_label))
  lrt_test_pval <- numeric(ncluster)
  total_nllk <- 0
  for(k in seq_len(ncluster)){
    cur_ct_gvec <- gvec_ct_list[[k]]
    cur_ct_size <- size_ct_list[[k]]
    cur_ct_batch <- batch_ct_list[[k]]
    cur_ct_nbatch <- length(unique(cur_ct_batch))
    
    ## fit one NB for all batches within this cluster
    common_nb_nllk <- nb_size_mle(cur_ct_gvec, cur_ct_size)$value
    if(cur_ct_nbatch == 1){
      total_nllk <- total_nllk + common_nb_nllk 
      lrt_test_pval[k] <- 1
    } else{
      sep_nllk <- 0
      gvec_ct_batch_list <- split(cur_ct_gvec, f = cur_ct_batch)
      size_ct_batch_list <- split(cur_ct_size, f = cur_ct_batch)
      
      ## fit separate NB for each batch within this cluster
      for(b in seq_len(cur_ct_nbatch)){
        sep_nllk <- sep_nllk + nb_size_mle(gvec_ct_batch_list[[b]], 
                                        size_ct_batch_list[[b]])$value
      }
      ## likelihood ratio test (LRT) to decide which one to choose
      lrt_test_pval_cur_ct <- round(pchisq(2*(common_nb_nllk - sep_nllk), 
                                          2*(cur_ct_nbatch - 1), 
                                          lower.tail = FALSE), 4)
      lrt_test_pval[k] <- lrt_test_pval_cur_ct
      total_nllk <- total_nllk + ifelse(lrt_test_pval_cur_ct < lrt_pval_threshold,  
                                       sep_nllk, common_nb_nllk)
    }
  }
  return(list(NegLLK = total_nllk, nReject = sum(lrt_test_pval < lrt_pval_threshold)))
  
}


## Calculate the CDI values for one label set
calculate_CDI_oneset <- function(
	sub_gcmat, 
	candidate_label, 
	batch_label = NULL, 
	cell_size_factor, 
	BPPARAM,
	lrt_pval_threshold = 0.01){
  original_ncluster <- length(unique(candidate_label))
  
  #### One-batch scenario
  if(is.null(batch_label) | (length(unique(batch_label)) == 1)){
    ## filter clusters with small number of cells
    min_ncell <- min(table(candidate_label))
    if(min_ncell < 3){
      sub_indx <- c(seq_len(length(candidate_label)))[(candidate_label %in% names(which(table(candidate_label) > 2)))]
      candidate_label <- candidate_label[sub_indx]
      sub_gcmat <- sub_gcmat[,sub_indx]
      cell_size_factor <- cell_size_factor[sub_indx]
    }
    ng <- nrow(sub_gcmat); nc <- ncol(sub_gcmat)
    ## after filtering, it is possible ncluster < original_cluster
    ## for fair comparison, use original_cluster in penalty
    ncluster <- length(unique(candidate_label))
    
    ## calculating log-likelihood
    if(is.null(rownames(sub_gcmat))){
      rownames(sub_gcmat) <- paste0("g", seq_len(ng))
    }
    sub_gclist <- split(sub_gcmat, f = rownames(sub_gcmat))
    neg_llk_list <- bplapply(sub_gclist, 
                            single_batch_one_gene_likelihood,  
                            candidate_label = candidate_label, 
                            ncluster = ncluster, 
                            cell_size_factor = cell_size_factor, 
                            BPPARAM = BPPARAM)
    neg_llk <- sum(unlist(neg_llk_list))
    npara <- ng * original_ncluster * 2
    
    #### Multi-batch scenario
  } else{
    combine_label <- paste0("ct_", candidate_label, "_b_", batch_label)
    min_combine_ncell <- min(table(combine_label))
    
    ## filter clusters with small number of cells
    if(min_combine_ncell < 3){
      sub_indx <- c(seq_len(length(combine_label)))[(combine_label %in% names(which(table(combine_label) > 2)))]
      candidate_label <- candidate_label[sub_indx]
      batch_label <- batch_label[sub_indx]
      sub_gcmat <- sub_gcmat[,sub_indx]
      cell_size_factor <- cell_size_factor[sub_indx]
    }
    ng <- nrow(sub_gcmat); nc <- ncol(sub_gcmat)
    ## after filtering, it is possible ncluster < original_cluster
    ## for fair comparison, use original_cluster in penalty
    ncluster <- length(unique(candidate_label))
    if(is.null(rownames(sub_gcmat))){
      rownames(sub_gcmat) <- paste0("g", seq_len(ng))
    }
    batch_ct_list <- split(batch_label, f = candidate_label)
    size_ct_list <- split(cell_size_factor, f = candidate_label)
    sub_gclist <- split(sub_gcmat, f = rownames(sub_gcmat))
    neg_llk_list <- bplapply(sub_gclist, 
                            multi_batch_one_gene_likelihood,  
                            candidate_label = candidate_label, 
                            ncluster = ncluster, 
                            batch_ct_list = batch_ct_list, 
                            size_ct_list  = size_ct_list,
                            lrt_pval_threshold = lrt_pval_threshold,
                            BPPARAM = BPPARAM)
    neg_llk <- sum(unlist(lapply(neg_llk_list, '[[', 'NegLLK')))
    total_rej <- sum(unlist(lapply(neg_llk_list, '[[', 'nReject')))
    npara <- (ng * original_ncluster + total_rej) * 2
  }
  return(data.frame(
  	CDI_AIC = 2*neg_llk + 2*npara, 
  	CDI_BIC = 2*neg_llk + npara*log(nc), 
  	neg_llk_val =  neg_llk,
  	N_cluster = original_ncluster))
  
}





## -----------------------------------------------------------------------------
##          External function -- size_factor
## -----------------------------------------------------------------------------
#' Size factor of each cell
#' 
#' Different cells have different library sizes. 
#' This function calculates the size factor of each cell in the UMI count matrix 
#' to capture the variation in cell library size. 
#'
#' @param X The class of X can be "matrix", "Seurat" object, or "SingleCellExperiment" object. 
#' If X is a matrix, it should be a raw UMI count matrix where each row represents a gene, and 
#' each column represents a cell. The genes should be those before feature gene selection.
#' If X is a Seurat object or SingleCellExperiment object, users need to specify where the count 
#' matrix is stored in count_slot.
#' 
#' @param count_slot A string indicating the location of raw UMI count. 
#' For Seurat object, it is a slot in "RNA" of "assays"; 
#' For SingleCellExperiment object, it is a slot in "assays". 
#' Each row represents a gene, and each column represents a cell. 
#' The genes should be those before feature gene selection.
#'
#' @importFrom matrixStats colMedians
#' @importFrom SingleCellExperiment SingleCellExperiment rowData colData
#' @importFrom SummarizedExperiment assays
#' @importFrom Seurat GetAssayData 
#' @importFrom methods is 
#' 
#' @return A numeric vector indicating the size factor of the cells. 
#' This should be one of the inputs of the function calculate_CDI.
#' @examples
#' 
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
#' 
#' ## Input: matrix
#' cell_size <- size_factor(X = X)
#' 
#' ## Input: SingleCellExperiment object
#' library(SingleCellExperiment)
#' sim_sce <- SingleCellExperiment(
#'   list(count = X),
#'   colData = data.frame(Cell_name = colnames(X)),
#' 	 rowData = data.frame(Gene_name = rownames(X)))
#' cell_size <- size_factor(X = sim_sce, count_slot = "count")
#' 
#' ## Input: Seurat object
#' library(Seurat)
#' library(SeuratObject)
#' sim_seurat <- CreateSeuratObject(counts = as.data.frame(X))
#' sim_seurat <- AddMetaData(sim_seurat, colnames(X), "Cell_name")
#' cell_size <- size_factor(X = sim_seurat, count_slot = "counts")
#' 
#' @export


size_factor <- function(
	X,
	count_slot = NULL){
	# extract counts
	gcmat = extract_count(X, count_slot)
	# calculate size factors
  gcmat[gcmat == 0] <- 0.5
  nc <- ncol(gcmat)
  log_gcmat <- log(gcmat)
  ref_size <- exp(rowMeans(log_gcmat))
  ratio_to_ref <- sweep(gcmat, 1, ref_size, "/")
  cell_size_factor <- colMedians(ratio_to_ref)
  return(cell_size_factor)
}



## -----------------------------------------------------------------------------
##          External function -- calculate_CDI
## -----------------------------------------------------------------------------

#' Clustering Deviance Index (CDI)
#'
#' This function calculates CDI-AIC and CDI-BIC for each candidate set of cell labels.
#' CDI calculates AIC and BIC of cell-type-specific gene-specific NB model for UMI counts, 
#' where the cell types are based on each candidate label set, 
#' and only the selected subset of genes are considered. 
#' Whether to use CDI-AIC or CDI-BIC depend on the goals. 
#' We suggest using CDI-BIC to select optimal main cell types and using CDI-AIC 
#' to select optimal subtypes, because BIC puts more penalty on the complexity 
#' of models (number of clusters). 
#' 
#' @param X The class of X can be "matrix", "Seurat" object, or "SingleCellExperiment" object. 
#' If X is a matrix, it should be a UMI count matrix where each row represents a gene, and 
#' each column represents a cell. 
#' If X is a Seurat object or SingleCellExperiment object, 
#' users need to specify where the count matrix and 
#' batch labels are stored in count_slot and batch_slot, respectively.
#' If feature_gene_index is NULL, genes in X should only included feature genes 
#' (that are selected by feature_gene_selection function); if feature_gene_index
#' is not NULL, this function will extract a subset of X 
#' with genes indexed by feature_gene_index. 
#' 
#' @param feature_gene_index A vector of unique integers indicating the indices of feature 
#' genes. The default value if NULL, which means all genes in X will be used to 
#' calculate CDI. The integers in feature_gene_index need to be no greater than 
#' the number of genes in X. 
#' 
#' @param cand_lab_df A vector of cluster labels of the cells or 
#' a data frame where each column corresponds to one set of cluster labels of 
#' the cells. This (these) label sets can be clustering results obtained by 
#' any clustering methods. The length (number of rows) of 
#' cand_lab_df should be the same as the number of columns in 
#' the count matrix. 
#' If the column names of label set data frame are provided with the format 
#' "[ClusteringMethod]_k[NumberOfClusters]" such as "KMeans_K5, `calculate_CDI` 
#' will extract the "[ClusteringMethod]" as the Cluster_method. 
#' The clustering method can also be provided in the 
#' argument "clustering_method" for each label set. 
#' 
#' @param cell_size_factor A numeric vector indicating the size factor 
#' of the cells. This should be the output of function size_factor. 
#' The length of cell_size_factor should be the same as the number of columns 
#' in the count matrix. 
#' 
#' @param batch_label A vector of characters indicating the batch labels of the cells. 
#' The length of batch_label should be the same as the number of columns 
#' in the count matrix.
#' 
#' @param count_slot A string indicating the location of raw UMI count. 
#' For Seurat object, it is a slot in "RNA" of "assays"; 
#' For SingleCellExperiment object, it is a slot in "assays". 
#' Each row represents a gene, and each column represents a cell. 
#' The genes should be those before feature gene selection.
#' 
#' @param batch_slot A string indicating the location of batch labels of cells.
#' For Seurat object, it is a slot in meta.data;
#' For SingleCellExperiment object, it is a slot in "colData". 
#' The default value is NULL indicating that there is no batch information available. 
#' 
#' @param lrt_pval_threshold A numeric value within (0, 1) indicating 
#' the p-value threshold for the likelihood ratio test (LRT). If multiple 
#' batches exist, within each cluster and each gene, CDI will test whether 
#' a batch-common NB model or a batch-specific NB model should be fitted 
#' with the LRT. If the p-value is less than this threshold, a batch-specific 
#' NB model will be fitted. Otherwise, a batch-common NB model will be fitted.
#' 
#' @param clustering_method A vector of characters indicating the corresponding clustering 
#' method for each label set. The length of the vector needs to be the same 
#' as the number of columns in cand_lab_df.
#' 
#' @param BPPARAM A \code{\link{BiocParallelParam}} object from the BiocParallel 
#' package. By specifying this argument, users can control over how to perform 
#' the parallel computing. Default is \code{\link{SerialParam}} which uses a 
#' single core.
#'        
#' @importFrom SingleCellExperiment SingleCellExperiment rowData colData
#' @importFrom SummarizedExperiment assays
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom stats nlminb pchisq var
#' @importFrom matrixStats rowMedians
#' @importFrom Seurat GetAssayData FetchData
#' @importFrom methods is
#' 
#' @return calculate_CDI returns a data frame with 5 columns. The columns are 
#' Label_name (name of each label set), Cluster_method (clustering method), CDI-AIC, 
#' CDI-BIC, and N_cluster (number of clusters). Each row corresponds to one set of cell labels.
#' 
#' @examples
## Simulate count matrix, batch, and cell clustering labels
#' ng <- 100; nc <- 100
#' set.seed(1)
#' 
#' # count matrix
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
#' 
#' # batch label
#' Batches <- rep(seq_len(2), nc/2)
#' 
#' # cell clustering labels
#' Method1_k2 <- rep(seq_len(2), c(nc/2,nc/2))
#' Method1_k3 <- sample(seq_len(3), nc, replace = TRUE)
#' label_df <- data.frame(
#' 	Method1_k2 = Method1_k2,
#' 	Method1_k3 = Method1_k3)
#' 
#' ## select feature genes (see feature_gene_selection function)
#' selected_genes <- seq_len(30)
#' 
#' ## calculate size factor (see size_factor function)
#' size_factor_vec <- rep(1, nc)
#' 
#' calculate_CDI(
#'  X = X[selected_genes, ],
#' 	cand_lab_df = label_df,
#' 	cell_size_factor = size_factor_vec,
#' 	batch_label = Batches)
#' 
#' ## Input: SingleCellExperiment object
#' library(SingleCellExperiment)
#' sim_sce <- SingleCellExperiment(
#'   list(count = X),
#'   colData = data.frame(
#'     Cell_name = colnames(X),
#' 		 batch = Batches),
#' 	 rowData = data.frame(
#' 	   Gene_name = rownames(X)))
#' 
#' calculate_CDI(
#'  	X = sim_sce,
#'  	feature_gene_index = selected_genes, 
#'  	cand_lab_df = label_df,
#'  	cell_size_factor = size_factor_vec,
#'  	count_slot = "count",
#'  	batch_slot = "batch")
#' 
#' ## Input: Seurat object
#' library(Seurat)
#' library(SeuratObject)
#' sim_seurat <- CreateSeuratObject(counts = as.data.frame(X))
#' sim_seurat <- AddMetaData(sim_seurat, colnames(X), "Cell_name")
#' sim_seurat <- AddMetaData(sim_seurat, Batches, "batch")
#' 
#' calculate_CDI(
#' 	X = sim_seurat,
#' 	feature_gene_index = selected_genes,
#' 	cand_lab_df = label_df,
#' 	cell_size_factor = size_factor_vec,
#'  count_slot = "counts", 
#'  batch_slot = "batch")
#' 
#' ## parallel computing
#' library(BiocParallel)
#' ## single core
#' bp_object <- SerialParam()
#' ## multi-cores
#' ## bp_object  <- MulticoreParam(workers = 2)
#' calculate_CDI(
#' 	X = X[selected_genes, ],
#' 	cand_lab_df = label_df,
#' 	cell_size_factor = size_factor_vec,
#' 	batch_label = Batches,
#' 	lrt_pval_threshold = 0.01,
#' 	clustering_method = NULL,
#' 	BPPARAM = bp_object)
#'              
#' @references SMartin Morgan, Valerie Obenchain, Michel Lang, Ryan 
#' Thompson and Nitesh Turaga (2021). 
#' \doi{https://github.com/Bioconductor/BiocParallel}
#' @export
calculate_CDI <- function(
	X, 
	feature_gene_index = NULL,
	cand_lab_df, 
	cell_size_factor, 
	batch_label = NULL,
	count_slot = NULL, 
	batch_slot = NULL,
	lrt_pval_threshold = 0.01,
	clustering_method = NULL,
	BPPARAM = SerialParam()){
	# extract count and batch label
	gcmat <- extract_count(X, count_slot)
	batch_label <- extract_batch(X, batch_label, batch_slot)
	if(!is.null(batch_label) & length(batch_label) != ncol(gcmat)){
		stop("the length of batch_label does not match the number of cells (columns) in X.")
	}
	if(!is.null(feature_gene_index)){
			if(max(feature_gene_index) > nrow(gcmat)){
				stop("feature_gene_index exceed the number of genes in X.")
			} else{
				sub_gcmat <- as.matrix(gcmat[feature_gene_index, ])
			}
		# if feature genes are not given, use all genes
		} else{
			sub_gcmat <- as.matrix(gcmat)
		}
	rm(X, gcmat)
	## if cand_lab_df is a vector or a data frame with one column
	vec_1col_df <- ifelse(is.vector(cand_lab_df), TRUE, dim(cand_lab_df)[2] == 1)
  if(vec_1col_df){
    return(calculate_CDI_oneset(
    	sub_gcmat = sub_gcmat, 
    	candidate_label = unlist(cand_lab_df), 
    	batch_label = batch_label, 
    	cell_size_factor = cell_size_factor, 
    	BPPARAM = BPPARAM, 
    	lrt_pval_threshold = lrt_pval_threshold))
  ## if cand_lab_df is a a data frame with more than one column
  } else {
    lab_name <- colnames(cand_lab_df)
    cdi_return_df <- data.frame(Label_name = paste0("Label", seq_len(ncol(cand_lab_df))))
    if(!is.null(lab_name)){
      cdi_return_df["Label_name"] <- lab_name
      cdi_return_df["Cluster_method"] <- ifelse(
      	grepl(pattern = "^(\\w+)(_k)(\\d+)$", x = lab_name, ignore.case = TRUE), 
      	unlist(lapply(strsplit(lab_name, "_"), "[", 1)), 
      	NA)
     }
    if(!is.null(clustering_method)){
    	cdi_return_df["Cluster_method"] <- clustering_method
    }
    cdi_return <- apply(
    	X = cand_lab_df, 
    	MARGIN = 2, 
    	FUN = calculate_CDI_oneset,
    	sub_gcmat = sub_gcmat, 
    	batch_label = batch_label, 
    	cell_size_factor = cell_size_factor,
    	BPPARAM = BPPARAM, 
    	lrt_pval_threshold = lrt_pval_threshold) 
    tmp_df <- do.call(rbind.data.frame, cdi_return)
    cdi_return_df <- cbind(cdi_return_df, tmp_df)
    return(cdi_return_df)
  }
}


