test_calculate_CDI <- function() {
	
	set.seed(100)
	X <- cbind(matrix(c(rnbinom(2500, size = 1, mu = 0.1), rnbinom(2500, size = 1, mu = 0.5)), 
	                  nrow = 100, byrow = TRUE),
	           matrix(c(rnbinom(2500, size = 1, mu = 1), rnbinom(2500, size = 1, mu = 0.5)), 
	                  nrow = 100, byrow = TRUE))
	labs <- data.frame(TrueLab = rep(c(1,2), c(50,50)), 
										 RandomLab = sample(c(1,2), size = 100, replace = TRUE))
	CDI_return <- calculate_CDI(sub_gcmat = X,
									            cand_lab_df = labs,
									            batch_label = rep(c(1,2), ncol(X)/2),
									            cell_size_factor = rep(1, 100))


    checkEquals(CDI_return[which.min(CDI_return$CDI_BIC), "Label_name"], "TrueLab")
}