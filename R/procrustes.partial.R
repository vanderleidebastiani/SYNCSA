#' @rdname procrustes.syncsa
#' @encoding UTF-8
#' @export
procrustes.partial<-function(x, y, z)
{
	x <- as.matrix(x)
	y <- as.matrix(y)
	z <- as.matrix(z)
#	pro.residuals<-function(x,scoresofz){
#		res<-matrix(NA,dim(x)[1],dim(x)[2])
#		for(i in 1:dim(x)[2]){ 
#			res[,i]<-cbind(stats::residuals(stats::lm(V1~.,data=as.data.frame(cbind(x[,i],scoresofz)))))
#		}
#	return(res)
#	}	
#	pro.residuals<-function(Y,X){
#		X<-base::scale(X, scale = FALSE)
#		Y<-base::scale(Y, scale = FALSE)
#		Yfit<-X%*%solve(t(X)%*%X)%*%(t(X))%*%Y
#		Yres<-Y-Yfit
#	return(Yres)
#	}
	pro.residuals <- function(xs, scores){
		res<-matrix(NA,dim(xs)[1],dim(xs)[2])
		for(i in 1:dim(xs)[2]){
			mat <- cbind(1, scores)
			fast_mod <- RcppArmadillo::fastLmPure(mat, x[,i])
			coeffs <- fast_mod$coefficients
			res[,i] <- x[,i] - mat %*% coeffs
		}
	return(res)
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