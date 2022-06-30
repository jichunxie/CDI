test_calculate_CDI <- function() {
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
	
	label_df <- data.frame(
		TrueLab = rep(c(1,2), c(nc/2, nc/2)),
		RandomLab = sample(c(1,2), size = nc, replace = TRUE))
		CDI_return <- calculate_CDI(sub_gcmat = X,
			cand_lab_df = label_df,
			cell_size_factor = rep(1, 100))
	checkEquals(CDI_return[which.min(CDI_return$CDI_BIC), "Label_name"], "TrueLab")
}



test_size_factor <- function() {
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
	sf_return <- size_factor(X)
	checkEquals(sum(sf_return > 0), ncol(X))
}

test_feature_selection <- function() {
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
	nf = 50
	checkEquals(length(feature_gene_selection(X, nfeature = nf)), nf)
	
}