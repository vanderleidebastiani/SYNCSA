#'@rdname cor.matrix
#'@export
pro.matrix<-function (m1, m2, x, y, permutations = 999, norm = FALSE, strata = NULL) 
{
	x <- cbind(x)
	y <- cbind(y)
	correlation<-procrustes.syncsa(x,y)
	value <- matrix(NA, nrow = permutations, ncol = 1)
	for (i in 1:permutations) {
		m2.permut <- permut.row.matrix(m2, strata = strata)
		x.permut <- m1 %*% m2.permut
		if (norm == "TRUE") {
			matrix.permut <- apply(x.permut^2, 2, sum)
			x.permut <- sweep(x.permut, 2, sqrt(matrix.permut), "/")
		}
		value[i,] <- procrustes.syncsa(x.permut,y) 
	}
	sig <- (sum(value >= correlation) + 1)/(permutations + 1)
return(list(Obs = correlation, p = sig))
}