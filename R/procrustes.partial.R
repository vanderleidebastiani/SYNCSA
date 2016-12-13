#' @rdname procrustes.syncsa
#' @encoding UTF-8
#' @export
#' @importFrom RcppArmadillo fastLm
procrustes.partial<-function(x, y, z)
{
	x <- as.matrix(x)
	y <- as.matrix(y)
	z <- as.matrix(z)

	pro.residuals <- function(Y,X){
	  X<-sweep(X,2,colMeans(X,na.rm=TRUE),check.margin = FALSE)
	  Y<-sweep(Y,2,colMeans(Y,na.rm=TRUE),check.margin = FALSE)
	  Yfit<-X%*%solve(t(X)%*%X)%*%(t(X))%*%Y
	  Yres<-Y-Yfit
	  return(Yres)
	}
	scoresofz<-stats::prcomp(z,scale = TRUE)$x
	nm<-round((dim(x)[1]-2)/2)
	if(nm<dim(scoresofz)[2]){
	  scoresofz<-scoresofz[,1:nm,drop=FALSE]
	}
	x.r<-pro.residuals(x,scoresofz)
	y.r<-pro.residuals(y,scoresofz)
	statistic<-procrustes.syncsa(x.r,y.r)
#	rxy <- procrustes.syncsa(x, y)
#	rxz <- procrustes.syncsa(x, z)
#	ryz <- procrustes.syncsa(y, z)
#	statistic <- part.cor(rxy, rxz, ryz)
return(statistic)
}
