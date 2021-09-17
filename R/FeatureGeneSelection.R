## ------------------------------------------------------------------------------------
##          Internal function -- OneBatchFeatureGeneSelection
## ------------------------------------------------------------------------------------

OneBatchFeatureGeneRank = function(gcmat,
                                   method = "wds", 
                                   nfeature = 500,
                                   zp_threshold = 0.95){
  
  ng = nrow(gcmat); nc = ncol(gcmat)
  if(method == "wds"){
    zp = rowMeans(gcmat == 0)
    zp_indx = c(seq_len(ng))[zp < zp_threshold]
    mu_g = rowMeans(gcmat); var_g = rowVars(gcmat)
    phi_est = (var_g-mu_g)/mu_g^2
    phi_df = data.frame(indx = seq_len(ng), 
                        phi_est = phi_est)
    phi_df[zp_indx, "phi_rank"] = rank(-phi_est[zp_indx])
    ## Assign genes with NA rank to a large value to remove NA
    phi_df$phi_rank[is.na(phi_df$phi_rank)] = ng + 1000
    return(phi_df$phi_rank)
  }
  if(method == "vst"){
    if(is.null(colnames(gcmat))){
      colnames(gcmat) = paste0("c", seq_len(nc))
    }
    if(is.null(rownames(gcmat))){
      rownames(gcmat) = paste0("g", seq_len(ng))
    }
    gene_indx_name = data.frame(indx = c(seq_len(ng)), grank = rep(ng, ng))
    rownames(gene_indx_name) = rownames(gcmat)
    seurat_obj = CreateSeuratObject(counts = gcmat, min.cells = 0, min.features = 0)
    seurat_obj = NormalizeData(seurat_obj, verbose = FALSE)
    seurat_obj = FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = nfeature, verbose = FALSE)
    # seu_sig_gname = seurat_obj@assays[["RNA"]]@var.features
    tmp_df = data.frame(gname = seurat_obj@assays[["RNA"]]@var.features, grank = seq_len(nfeature))
    gene_indx_name[as.character(tmp_df$gname), "grank"] = tmp_df$grank
    return(gene_indx_name$grank)
  }
  else{
    return("Error: Method not available!")
  }
}










## --------------------------------------------------------------------------
##        External function -- FeatureGeneSelection
## --------------------------------------------------------------------------

#' Select feature genes
#' 
#' This function selects a subset of feature genes that are expected to express 
#' differently across cell types before calculating CDI.
#'
#' @param gcmat A raw count matrix where each row represents a gene, and each column represents a cell.
#' @param method A character indicating the method used to select feature genes. "wds" 
#' (default) represent the working dispersion score proposed in our study; "vst" is the default method 
#' for feature gene selection in FindVariableFeatures function of Seurat package.
#' @param nfeature An integer indicating the number of features to select.
#' @param zp_threshold A number of zero proportion threshold. The range should be (0,1). Genes with 
#' zero proportion greater than this value will not be considered for feature selection.
#' @param batch_label A vector indicating the batch labels of the cells. The length of batch_label 
#' should be the same as the number of columns in gcmat.
#' @importFrom matrixStats rowVars rowMins
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures
#' @return A vector of indices of selected feature genes corresponding to the row indices of input count matrix.
#'
#' @examples
#' set.seed(100)
#' X = cbind(matrix(c(rnbinom(2500, size = 1, mu = 0.1), 
#'                    rnbinom(2500, size = 1, mu = 0.5)), 
#'                  nrow = 100, byrow = TRUE),
#'           matrix(c(rnbinom(2500, size = 1, mu = 1), 
#'                    rnbinom(2500, size = 1, mu = 0.5)), 
#'                  nrow = 100, byrow = TRUE))
#' batches = rep(c(1,2), ncol(X)/2)
#' FeatureGeneSelection(gcmat = X, 
#'                      batch_label = batches, 
#'                      nfeature = 20, 
#'                      zp_threshold = 0.95)
#'                      
#' @references Stuart and Butler et al. Comprehensive Integration of Single-Cell Data. Cell (2019) [Seurat V3]
#' @export
FeatureGeneSelection = function(gcmat, 
                                method = "wds", 
                                nfeature = 500, 
                                batch_label = NULL,
                                zp_threshold = 0.95){
  ng = nrow(gcmat)
  if(ng < nfeature){
    return("Error: The number of feature genes to select exceeds the number of genes.")
  }
  if(is.null(batch_label) | (length(unique(batch_label)) == 1)){
    return(c(seq_len(ng))[(OneBatchFeatureGeneRank(gcmat, method = method, nfeature = nfeature, zp_threshold = zp_threshold) <= nfeature)])
  } else{
    ## split gcmat columns according to the batch labels
    mat_list = lapply(split(seq_along(batch_label), batch_label), function(indx, mat) mat[,indx], mat = gcmat)
    ## rank genes in each matrix
    gene_rank_list = lapply(mat_list, OneBatchFeatureGeneRank, method = method, nfeature = nfeature, zp_threshold = zp_threshold)
    ## return top nfeature genes
    gene_rank_mat = do.call(cbind, gene_rank_list)
    gene_min_rank = rowMins(gene_rank_mat)
    names(gene_min_rank) = seq_len(nrow(gcmat))
    as.numeric(names(gene_min_rank)[order(gene_min_rank)[seq_len(nfeature)]])
    return(sort(as.numeric(names(gene_min_rank)[order(gene_min_rank)[seq_len(nfeature)]])))
  }
  
}

