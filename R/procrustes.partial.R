#' @rdname procrustes.syncsa
#' @encoding UTF-8
#' @export
#' @importFrom RcppArmadillo fastLm
procrustes.partial<-function(x, y, z)
{
	x <- as.matrix(x)
	y <- as.matrix(y)
	z <- as.matrix(z)
	pro.residuals<-function(x,scoresofz){
		res<-matrix(NA,dim(x)[1],dim(x)[2])
		for(i in 1:dim(x)[2]){
		  dff <- data.frame(V1 = x[,i], scoresofz)
		  fmod <- fastLm(V1~.,data=dff)
			res[,i]<-cbind(fmod[["residuals"]])
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
