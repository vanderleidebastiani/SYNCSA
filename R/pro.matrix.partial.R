#'@rdname cor.matrix
#'@export
pro.matrix.partial<-function (m1, m2, x, y, z, permutations = 999, norm = FALSE, strata = NULL)
{
	m1 <- as.matrix(m1)
	m2 <- as.matrix(m2)
	x <- as.matrix(x)
	y <- as.matrix(y)
	z <- as.matrix(z)
	statistic <- procrustes.partial(x,y,z)
	value <- matrix(NA, nrow = permutations, ncol = 1)
	for (i in 1: permutations) {
		m2.permut <- permut.row.matrix(m2, strata = strata)
		x.permut <- m1 %*% m2.permut
		if (norm == "TRUE") {
			matrix.permut <- apply(x.permut^2, 2, sum)
			x.permut <- sweep(x.permut, 2, sqrt(matrix.permut), "/")
		}
		value[i,] <- procrustes.partial(x.permut,y,z)
	}
	signif <- (sum(abs(value) >= abs(statistic)) + 1)/(permutations + 1)
return(list(Obs = statistic, p = signif))
}