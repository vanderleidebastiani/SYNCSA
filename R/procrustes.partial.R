#' @rdname procrustes.syncsa
#' @encoding UTF-8
#' @export
procrustes.partial<-function(x, y, z)
{
	x <- as.matrix(x)
	y <- as.matrix(y)
	z <- as.matrix(z)
	pro.residuals<-function(x,scoresofz){
		res<-matrix(NA,dim(x)[1],dim(x)[2])
		for(i in 1:dim(x)[2]){ 
			res[,i]<-cbind(residuals(lm(V1~.,data=as.data.frame(cbind(x[,i],scoresofz)))))
		}
	return(res)
	}
	scoresofz<-prcomp(z,scale = TRUE)$x
	x.r<-pro.residuals(x,scoresofz)
	y.r<-pro.residuals(y,scoresofz)
#	x.r<-sweep(x.r, 2, sqrt(apply(x.r^2,2,sum)), "/")
#	y.r<-sweep(x.r, 2, sqrt(apply(x.r^2,2,sum)), "/")
	statistic<-procrustes.syncsa(x.r,y.r)
#	rxy <- procrustes.syncsa(x, y)
#	rxz <- procrustes.syncsa(x, z)
#	ryz <- procrustes.syncsa(y, z)
#	statistic <- part.cor(rxy, rxz, ryz)
return(statistic)
}