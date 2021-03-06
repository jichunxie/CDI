% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExtractInfo.R
\name{extract_seurat}
\alias{extract_seurat}
\title{Extract information from Seurat object}
\usage{
extract_seurat(seurat_obj, info_type, info_slot)
}
\arguments{
\item{seurat_obj}{Seurat object}

\item{info_type}{A character indicating the type of information to extract. 
"count" represents the count matrix; "label" represents cell labels; "batch"
represent batch labels.}

\item{info_slot}{Slot(s) in assays$RNA that stores the count matrix, or slot in 
meta.data that saves cluster / batch labels.}
}
\value{
A count matrix if info_type = "count"; A data frame if 
info_type = "label"; A vector of characters if info_type = "batch".
}
\description{
This function extracts count matrix, cell labels, and batch labels from 
Seurat object.
}
\examples{
## Simulate count matrix, batch, and cell type labels
ng = 100; nc = 100
set.seed(1)
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
Batches <- rep(seq_len(2), nc/2)
Method1_k2 <- rep(seq_len(2), c(nc/2,nc/2))
Method1_k3 <- sample(seq_len(3), nc, replace = TRUE)

library(Seurat)
library(SeuratObject)
sim_seurat <- CreateSeuratObject(counts = as.data.frame(X))
sim_seurat <- AddMetaData(sim_seurat, colnames(X), "Cell_name")
sim_seurat <- AddMetaData(sim_seurat, Batches, "Batch")
sim_seurat <- AddMetaData(sim_seurat, Method1_k2, "Method1_k2")
sim_seurat <- AddMetaData(sim_seurat, Method1_k3, "Method1_k3")
	
extract_seurat(sim_seurat, "count", "counts")
extract_seurat(sim_seurat, "label", c("Method1_k2", "Method1_k3"))
extract_seurat(sim_seurat, "batch", "Batch")
}
