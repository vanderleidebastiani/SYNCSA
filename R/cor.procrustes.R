#' @rdname cor.matrix
#' @encoding UTF-8
#' @export
cor.procrustes<-function (mx, my, permutations = 999, strata = NULL, na.rm = FALSE, seqpermutation = NULL, parallel = NULL, newClusters=TRUE, CL =  NULL) 
{
	if(!is.null(seqpermutation)){
		if(dim(seqpermutation)[1]!=permutations){
			stop("\n The seqpermutation must be the dimension of permutations\n")
		}
	}
    correlation <- procrustes.syncsa(mx,my)
	N<-dim(mx)[1]
    if(is.null(seqpermutation)){
		seqpermutation <- permut.vector(N, strata = strata, nset = permutations)
	}
	if(!is.null(CL)){
		parallel<-length(CL)
	}
	ptest <- function(samp,mx, my){ 
		res<-SYNCSA::procrustes.syncsa(mx, my[samp,,drop=FALSE])
		return(res)
	}
    if(is.null(parallel)){
    	value <- matrix(NA, nrow = permutations, ncol = 1)
	    for (i in 1: permutations) {
	        value[i,] <- ptest(samp = seqpermutation[i,],mx = mx, my=my)
    	}
	} else {
		if (newClusters) {
			CL <- parallel::makeCluster(parallel,type="PSOCK")
		}
		value <- cbind(parallel::parRapply(CL, seqpermutation, ptest, mx = mx,  my = my))
		if (newClusters){
			parallel::stopCluster(CL)
		}
	}
    signific <- (sum(abs(value) >= abs(correlation)) + 1)/(permutations + 1)
    return(list(Obs = correlation, p = signific))
}