#' @rdname cor.matrix
#' @encoding UTF-8
#' @export
pro.matrix.partial<-function (mx1, mx2, x, my1, my2, y, mz1, mz2, z, permute.my2 = FALSE, permute.mz2 = FALSE, permutations = 999, norm = FALSE, norm.y = FALSE, norm.z = FALSE, strata = NULL)
{
    mx1<-as.matrix(mx1)
    mx2<-as.matrix(mx2)
    x<-as.matrix(x)
    if(permute.my2){
	    my1<-as.matrix(my1)
    	my2<-as.matrix(my2)
    }
    y<-as.matrix(y)
    if(permute.mz2){
	    mz1<-as.matrix(mz1)
    	mz2<-as.matrix(mz2)
    }
    z<-as.matrix(z)
	statistic <- procrustes.partial(x,y,z)
	value <- matrix(NA, nrow = permutations, ncol = 1)
	for (i in 1: permutations) {
M.permut <- permut.row.matrix(mx2, strata = strata)
        mx2.permut <- M.permut$permut.matrix
        x.permut <- mx1 %*% mx2.permut
        if (norm == "TRUE") {
            matrix.permut <- apply(x.permut^2, 2, sum)
            x.permut <- sweep(x.permut, 2, sqrt(matrix.permut), "/")
        }
        if(permute.my2){
			y.permut <- my1 %*% my2[M.permut$samp,,drop=FALSE]
			if (norm.y == "TRUE") {
	            matrix.permut <- apply(y.permut^2, 2, sum)
    	        y.permut <- sweep(y.permut, 2, sqrt(matrix.permut), "/")
        	}
        }
        if(permute.mz2){
			z.permut <- mz1 %*% mz2[M.permut$samp,,drop=FALSE]
			if (norm.z == "TRUE") {
            	matrix.permut <- apply(z.permut^2, 2, sum)
	            z.permut <- sweep(z.permut, 2, sqrt(matrix.permut), "/")
    		}
		}
		if(!permute.my2 & !permute.mz2){
    	    value[i,] <- procrustes.partial(x.permut,y,z)
		}
		if(permute.my2 & !permute.mz2){
    	    value[i,] <- procrustes.partial(x.permut,y.permut,z)
		}
		if(!permute.my2 & permute.mz2){
    	    value[i,] <- procrustes.partial(x.permut,y,z.permut)
		}
		if(permute.my2 & permute.mz2){
    	    value[i,] <- procrustes.partial(x.permut,y.permut,z.permut)    	    
		}				
	}
	signif <- (sum(abs(value) >= abs(statistic)) + 1)/(permutations + 1)
return(list(Obs = statistic, p = signif))
}