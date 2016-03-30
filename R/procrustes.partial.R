#' @rdname procrustes.syncsa
#' @encoding UTF-8
#' @export
procrustes.partial<-function(x, y, z)
{
	x <- as.matrix(x)
	y <- as.matrix(y)
	z <- as.matrix(z)
	rxy <- procrustes.syncsa(x, y)
	rxz <- procrustes.syncsa(x, z)
	ryz <- procrustes.syncsa(y, z)
	statistic <- part.cor(rxy, rxz, ryz)
return(statistic)
}