procrustes.syncsa<-function (x, y)
{
	x <- cbind(x)
	y <- cbind(y)
	if (ncol(x) < ncol(y)){
		addcols <- ncol(y) - ncol(x)
		for (i in 1:addcols) x <- cbind(x, 0)
	}
	ctrace <- function(MAT) sum(diag(crossprod(MAT)))
	x <- scale(x, scale = FALSE)
	y <- scale(y, scale = FALSE)
	x <- x/sqrt(ctrace(x))
	y <- y/sqrt(ctrace(y))
	xy <- crossprod(x, y)
	sol <- svd(xy)
	c <- sum(sol$d)/ctrace(y)
	r2 <- ctrace(x) + c * c * ctrace(y) - 2 * c * sum(sol$d)
	ro<-sqrt(1 - r2)
return(ro)
}