#' @rdname cor.matrix
#' @encoding UTF-8
#' @export
pro.matrix<-function (mx1, mx2, x, y, permutations = 999, norm = FALSE, strata = NULL, seqpermutation = NULL, parallel = NULL, newClusters=TRUE, CL =  NULL) 
{
	if(!is.null(seqpermutation)){
		if(dim(seqpermutation)[1]!=permutations){
			stop("\n The seqpermutation must be the dimension of permutations\n")
		}
	}
	x <- cbind(x)
	y <- cbind(y)
	correlation<-procrustes.syncsa(x,y)
	N <- dim(mx2)[1]
    if(is.null(seqpermutation)){
		seqpermutation <- permut.vector(N, strata = strata, nset = permutations)
	}
	if(!is.null(CL)){
		parallel<-length(CL)
	}
	ptest<-function(samp,mx1,mx2,norm,y){
		x.permut <- mx1 %*% mx2[samp,,drop=FALSE]
		if (norm) {
			matrix.permut <- apply(x.permut^2, 2, sum)
			x.permut <- sweep(x.permut, 2, sqrt(matrix.permut), "/")
		}
		res <- SYNCSA::procrustes.syncsa(x.permut,y)
		return(res)
	}
	if(is.null(parallel)){
    	value <- matrix(NA, nrow = permutations, ncol = 1)
	    for (i in 1: permutations) {
	        value[i,] <- ptest(samp = seqpermutation[i,],mx1=mx1,mx2=mx2, norm=norm,y=y)
    	}
	} else {
		if (newClusters) {
			CL <- parallel::makeCluster(parallel,type="PSOCK")
		}
		value <- cbind(parallel::parRapply(CL, seqpermutation, ptest,mx1=mx1,mx2=mx2, norm=norm,y=y))
		if (newClusters){
			parallel::stopCluster(CL)
		}
	}
	sig <- (sum(value >= correlation) + 1)/(permutations + 1)
return(list(Obs = correlation, p = sig))
}