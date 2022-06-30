library(Matrix)
library(matrixStats)
library(truncnorm)


rm(list = setdiff(ls(), lsf.str()))

# ------------------------------------------------------------------------------
#              Simulate one batch count matrix
# ------------------------------------------------------------------------------



# Simulation settings
set.seed(1)
set_5cl <- list(K = 5, 
                ngene = 3000, 
                nsig = 30,
                nfeature = 500,
                log2fc = 2.5,
                ncore = 20)

set_5cl$ncell <- rep(400, set_5cl$K)
set_5cl$clust.label.real <- rep(seq_len(set_5cl$K), set_5cl$ncell)
para.set <- set_5cl
set.indx <- "sim_5cl"
para.set$N0 <- sum(para.set$ncell)


# simulate mean and dispersion parameters
generic_mu <- 0.2
mu <- matrix(rep(round(rtruncnorm(para.set$ngene, a=0.001, b=1, 
	mean=generic_mu, sd=0.1), 4), para.set$K), 
	nrow=para.set$ngene)
r <- matrix(rep(round(rtruncnorm(para.set$ngene, a =0.2, b=5, 
	mean=0.5, sd=0.1), 4), para.set$K), 
	nrow=para.set$ngene)

ns <- para.set$nsig
nz <- round(para.set$nsig/3)
ls.sig <- split(seq_len((para.set$nsig * para.set$K)), 
               f = factor(rep(seq_len(para.set$K), rep(para.set$nsig,para.set$K))))

for(j in seq_len(para.set$K)){
  sig_g_indx <- ls.sig[[j]]
  non_sig_celltype <- setdiff(seq_len(para.set$K), j)
  zero_col <- sample(non_sig_celltype, round(length(non_sig_celltype)/2))
  zero_g_indx <- sample(ls.sig[[j]], nz)
  mu[zero_g_indx, zero_col] <- 0
}


for(j in seq_len(para.set$K)){
  mu[ls.sig[[j]],j] <- pmax(mu[ls.sig[[j]],j] *(2^para.set$log2fc), 
  	(generic_mu/2) * (2^para.set$log2fc)) 
}

for(j in seq_len(para.set$K)){
  r[ls.sig[[j]],j] <- r[ls.sig[[j]],j] + rnorm(ns, 0, sd = 0.05) 
}

# Simulate count matrix
X <- matrix(NA, nrow = para.set$ngene, ncol = para.set$N0)
for(g in seq_len(para.set$ngene)){
  tmp <- NULL
  for(cl in seq_len(para.set$K)){
    cur.mu <- mu[g,cl]
    if(cur.mu == 0){ 
    	tmp <- c(tmp, rep(0, para.set$ncell[cl]))
    	}
    else{ 
    	tmp <- c(tmp, rnbinom(para.set$ncell[cl], mu = cur.mu, size = r[g,cl]))
    	}
  }
  X[g,] <- tmp
}
colnames(X) <- paste0("rc",seq_len(para.set$N0))
rownames(X) <- paste0("g",seq_len(para.set$ngene))

one_batch_matrix <- X

## compress the count matrix
usethis::use_data(one_batch_matrix, compress = "bzip2", overwrite = TRUE)


one_batch_matrix_celltype <- paste0("type", set_5cl$clust.label.real)
usethis::use_data(one_batch_matrix_celltype, overwrite = TRUE)





# ------------------------------------------------------------------------------
#               Generate one batch dataset labels
# ------------------------------------------------------------------------------

## kmeans
set.seed(1)

cluster_number <- c(2,   3,    4,    5,    6,    7)
ncluster <- length(cluster_number)
method_name <- c("KMeans")
nmethod <- length(method_name)

gcmat_pc <- function(gcmat, npc){
  dat_pr <- prcomp(t(gcmat), scale = TRUE)[[5]][,seq_len(npc)]
  return(dat_pr)
}

X_pc <- gcmat_pc(X, npc = 200)

kmeans_df <- as.data.frame(
	matrix(NA, nrow = ncol(X), 
	ncol = ncluster*nmethod))
for(i in seq_len(ncluster)){
  for(m in seq_len(nmethod)){
    kmeans_df[,((i-1)*nmethod + m)] <- as.vector(
    	kmeans(X_pc, centers = cluster_number[i], nstart = 3)$cluster)
    colnames(kmeans_df)[((i-1)*nmethod + m)] <- paste0(method_name[m], "_k", cluster_number[i])
  }
}





## Seurat
set.seed(1)
library(Seurat)
Seurat.obj <- CreateSeuratObject(counts = X, project = "sc20a", 
	min.cells = 0, min.features = 0)
Seurat.obj <- NormalizeData(Seurat.obj, verbose = FALSE)
Seurat.obj <- FindVariableFeatures(Seurat.obj, selection.method = "vst", 
	nfeatures = 2000, verbose = FALSE)
Seurat.obj <- ScaleData(Seurat.obj, verbose = FALSE)
Seurat.obj <- RunPCA(Seurat.obj, npcs = 30, verbose = FALSE)
Seurat.obj <- FindNeighbors(Seurat.obj, reduction = "pca", dims = seq_len(20))



res_vec <- c(0.015, 0.018, 0.02, 0.5, 2.35, 2.5)
seurat_df <- matrix(NA, nrow = ncol(X), ncol = length(res_vec))
num_cl <- numeric(length(res_vec))
for(i in seq_len(length(res_vec))){
  Seurat.obj <- FindClusters(Seurat.obj, 
  	resolution = res_vec[i], 
  	verbose = FALSE)
  seurat_df[,i] <- as.vector(Seurat.obj@meta.data$seurat_clusters)
}

colnames(seurat_df) <- paste0("Seurat_k", cluster_number)


one_batch_matrix_label_df <- cbind(kmeans_df, seurat_df)
save(one_batch_matrix_label_df, file = "./data/one_batch_matrix_label_df.RData")


# ------------------------------------------------------------------------------
#               Session Information
# ------------------------------------------------------------------------------


# R version 4.1.1 (2021-08-10)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Big Sur 11.1
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
# 
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
# [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] CDI_0.99.2                  BiocParallel_1.26.1         roxygen2_7.1.1              devtools_2.4.2             
#  [5] usethis_2.0.1               SeuratObject_4.0.3          Seurat_4.0.5                SingleCellExperiment_1.14.1
#  [9] SummarizedExperiment_1.22.0 Biobase_2.52.0              GenomicRanges_1.44.0        GenomeInfoDb_1.28.4        
# [13] IRanges_2.26.0              S4Vectors_0.30.0            BiocGenerics_0.38.0         MatrixGenerics_1.4.3       
# [17] matrixStats_0.60.1    

