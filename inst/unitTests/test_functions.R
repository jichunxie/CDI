
test_calculate_CDI <- function() {
	# correct input
	ng <- 100; nc <- 100
	set.seed(1)
	# count matrix
	X <- cbind(
			matrix(
				c(rnbinom(ng*nc/4, size = 1, mu = 0.1),
					rnbinom(ng*nc/4, size = 1, mu = 0.5)),
				nrow = ng,
				byrow = TRUE),
			matrix(
				c(rnbinom(ng*nc/4, size = 1, mu = 1),
					rnbinom(ng*nc/4, size = 1, mu = 0.5)),
				nrow = ng,
				byrow = TRUE))
	colnames(X) <- paste0('c', seq_len(nc))
	rownames(X) <- paste0('g', seq_len(ng))
	# batch label
	Batches <- rep(seq_len(2), nc/2)
	# cell clustering labels
	label_df <- data.frame(
			TrueLab = rep(c(1,2), c(nc/2, nc/2)),
			RandomLab = sample(c(1,2), size = nc, replace = TRUE))
			
	
	## select feature genes
	nf = 30
	selected_genes <- seq_len(nf)
	
	## calculate size factor
	size_factor_vec <- rep(1, nc)
	
	## Input: SingleCellExperiment object
	library(SingleCellExperiment)
	sim_sce <- SingleCellExperiment(
	  list(count = X),
	  colData = data.frame(
	    Cell_name = colnames(X),
			 batch = Batches),
		 rowData = data.frame(
		   Gene_name = rownames(X)))
	
	## Input: Seurat object
	library(Seurat)
	library(SeuratObject)
	sim_seurat <- CreateSeuratObject(counts = as.data.frame(X))
	sim_seurat <- AddMetaData(sim_seurat, colnames(X), "Cell_name")
	sim_seurat <- AddMetaData(sim_seurat, Batches, "batch")
	
	# correct input
	## matrix
	matrix_CDI <- calculate_CDI(X = X[selected_genes,], 
		cand_lab_df = label_df, 
		cell_size_factor = size_factor_vec)
	checkEquals(matrix_CDI[which.min(matrix_CDI$CDI_BIC), "Label_name"], "TrueLab")
	## sce
	sce_CDI <- calculate_CDI(X = sim_sce, 
		feature_gene_index = selected_genes, 
		cand_lab_df = label_df, 
		count_slot = "count",
		cell_size_factor = size_factor_vec)
	checkEquals(sce_CDI[which.min(sce_CDI$CDI_BIC), "Label_name"], "TrueLab")
	## seurat
	seurat_CDI <- calculate_CDI(X = sim_seurat, 
		feature_gene_index = selected_genes, 
		count_slot = "counts",
		cand_lab_df = label_df, 
		cell_size_factor = size_factor_vec)
	checkEquals(seurat_CDI[which.min(seurat_CDI$CDI_BIC), "Label_name"], "TrueLab")
	
	# incorrect input should generate error
	checkTrue(tryCatch(calculate_CDI(X = X, 
		feature_gene_index = seq_len(ng + 10),
		cand_lab_df = label_df, 
		cell_size_factor = size_factor_vec), error = function(e) return(TRUE)))
	checkTrue(tryCatch(calculate_CDI(X = sim_sce, 
		feature_gene_index = selected_genes, 
		cand_lab_df = label_df, 
		cell_size_factor = size_factor_vec, 
		batch_slot = "unknown"), error = function(e) return(TRUE)))
	checkTrue(tryCatch(calculate_CDI(X = sim_seurat, 
		feature_gene_index = selected_genes, 
		cand_lab_df = label_df,
		cell_size_factor = size_factor_vec, 
		count_slot = "unknown"), error = function(e) return(TRUE)))
	checkTrue(tryCatch(calculate_CDI(X = sim_seurat, 
		feature_gene_index = selected_genes, 
		cand_lab_df = label_df, 
		cell_size_factor = size_factor_vec, 
		batch_slot = "unknown"), error = function(e) return(TRUE)))
	
	
}



test_size_factor <- function() {
	# correct input
	ng <- 100; nc <- 100
	set.seed(1)
	# count matrix
	X <- cbind(
			matrix(
				c(rnbinom(ng*nc/4, size = 1, mu = 0.1),
					rnbinom(ng*nc/4, size = 1, mu = 0.5)),
				nrow = ng,
				byrow = TRUE),
			matrix(
				c(rnbinom(ng*nc/4, size = 1, mu = 1),
					rnbinom(ng*nc/4, size = 1, mu = 0.5)),
				nrow = ng,
				byrow = TRUE))
	colnames(X) <- paste0('c', seq_len(nc))
	rownames(X) <- paste0('g', seq_len(ng))
	
	library(SingleCellExperiment)
	sim_sce <- SingleCellExperiment(
	  list(count = X),
	  colData = data.frame(
	    Cell_name = colnames(X)),
		 rowData = data.frame(
		   Gene_name = rownames(X)))
	
	## Input: Seurat object
	library(Seurat)
	library(SeuratObject)
	sim_seurat <- CreateSeuratObject(counts = as.data.frame(X))
	sim_seurat <- AddMetaData(sim_seurat, colnames(X), "Cell_name")
	
	# Check 
	## matrix
	sf_return <- size_factor(X)
	checkEquals(sum(sf_return > 0), ncol(X))
	## sce
	sf_return <- size_factor(sim_sce, count_slot = "count")
	checkEquals(sum(sf_return > 0), ncol(X))
	checkTrue(tryCatch(size_factor(X = sim_sce, count_slot = "unknown"), error = function(e) return(TRUE)))
	## seurat
	sf_return <- size_factor(sim_seurat, count_slot = "counts")
	checkEquals(sum(sf_return > 0), ncol(X))
	checkTrue(tryCatch(size_factor(X = sim_seurat, count_slot = "unknown"), error = function(e) return(TRUE)))
}

test_feature_selection <- function() {
	ng <- 100; nc <- 100
	set.seed(1)
	# count matrix
	X <- cbind(
			matrix(
				c(rnbinom(ng*nc/4, size = 1, mu = 0.1),
					rnbinom(ng*nc/4, size = 1, mu = 0.5)),
				nrow = ng,
				byrow = TRUE),
			matrix(
				c(rnbinom(ng*nc/4, size = 1, mu = 1),
					rnbinom(ng*nc/4, size = 1, mu = 0.5)),
				nrow = ng,
				byrow = TRUE))
	colnames(X) <- paste0('c', seq_len(nc))
	rownames(X) <- paste0('g', seq_len(ng))
	# batch label
	Batches <- rep(seq_len(2), nc/2)
	
	## Input: SingleCellExperiment object
	library(SingleCellExperiment)
	sim_sce <- SingleCellExperiment(
	  list(count = X),
	  colData = data.frame(
	    Cell_name = colnames(X),
			 batch = Batches),
		 rowData = data.frame(
		   Gene_name = rownames(X)))
	
	## Input: Seurat object
	library(Seurat)
	library(SeuratObject)
	sim_seurat <- CreateSeuratObject(counts = as.data.frame(X))
	sim_seurat <- AddMetaData(sim_seurat, colnames(X), "Cell_name")
	sim_seurat <- AddMetaData(sim_seurat, Batches, "batch")
	
	nf = 20
	# correct input
	checkEquals(length(feature_gene_selection(X = X, 
		batch_label = Batches, 
		nfeature = nf)), nf)
	checkEquals(length(feature_gene_selection(X = sim_sce, 
		count_slot = "count", 
		batch_slot = "batch", 
		nfeature = nf)), nf)
	 checkEquals(length(feature_gene_selection(X = sim_seurat, 
		count_slot = "counts", 
		batch_slot = "batch", 
		nfeature = nf)), nf)
	# incorrect input
	checkTrue(tryCatch(feature_gene_selection(X = sim_sce, 
		count_slot = "unknown", 
		batch_slot = "batch",
		nfeature = nf), error = function(e) return(TRUE)))
	checkTrue(tryCatch(feature_gene_selection(X = sim_sce, 
		count_slot = "count",
		batch_slot = "unknown",
		nfeature = nf), error = function(e) return(TRUE)))
	checkTrue(tryCatch(feature_gene_selection(X = sim_seurat, 
		count_slot = "unknown", 
		batch_slot = "batch", 
		nfeature = nf), error = function(e) return(TRUE)))
	checkTrue(tryCatch(feature_gene_selection(X = sim_seurat, 
		count_slot = "counts", 
		batch_slot = "unknown", 
		nfeature = nf), error = function(e) return(TRUE)))
		
	
}
