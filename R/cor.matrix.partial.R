#' @rdname cor.matrix
#' @encoding UTF-8
#' @export
cor.matrix.partial<-function (mx1, mx2, x, my1, my2, y, mz1, mz2, z, method = "pearson", dist = "euclidean", permute.my2 = FALSE, permute.mz2 = FALSE, permutations = 999, norm = FALSE, norm.y = FALSE, norm.z = FALSE, strata = NULL, na.rm = FALSE) 
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
    dist.x <- vegan::vegdist(x, method = dist, na.rm = na.rm)
    dist.y <- vegan::vegdist(y, method = dist, na.rm = na.rm)
    dist.z <- vegan::vegdist(z, method = dist, na.rm = na.rm)
    rxy <- cor(dist.x, dist.y, method = method)
    rxz <- cor(dist.x, dist.z, method = method)
    ryz <- cor(dist.y, dist.z, method = method)
    statistic <- part.cor(rxy, rxz, ryz)
    value <- matrix(NA, nrow = permutations, ncol = 1)
    for (i in 1:permutations) {
    	M.permut <- permut.row.matrix(mx2, strata = strata)
        mx2.permut <- M.permut$permut.matrix
        x.permut <- mx1 %*% mx2.permut
        if (norm == "TRUE") {
            matrix.permut <- apply(x.permut^2, 2, sum)
            x.permut <- sweep(x.permut, 2, sqrt(matrix.permut), "/")
        }
        dist.x.permut <- vegan::vegdist(x.permut, method = dist, na.rm = na.rm)
        ####
        if(permute.my2){
			y.permut <- my1 %*% my2[M.permut$samp,,drop=FALSE]
			if (norm.y == "TRUE") {
	            matrix.permut <- apply(y.permut^2, 2, sum)
    	        y.permut <- sweep(y.permut, 2, sqrt(matrix.permut), "/")
        	}
	        dist.y.permut <- vegan::vegdist(y.permut, method = dist, na.rm = na.rm)	
        }
        if(permute.mz2){
			z.permut <- mz1 %*% mz2[M.permut$samp,,drop=FALSE]
			if (norm.z == "TRUE") {
            	matrix.permut <- apply(z.permut^2, 2, sum)
	            z.permut <- sweep(z.permut, 2, sqrt(matrix.permut), "/")
    		}
    		dist.z.permut <- vegan::vegdist(z.permut, method = dist, na.rm = na.rm)	
		}
		if(!permute.my2 & !permute.mz2){
	        rxy.temp <- cor(dist.x.permut, dist.y, method = method)
    	    rxz.temp <- cor(dist.x.permut, dist.z, method = method)	
    	    ryz.temp <- ryz
		}
		if(permute.my2 & !permute.mz2){
	        rxy.temp <- cor(dist.x.permut, dist.y.permut, method = method)
    	    rxz.temp <- cor(dist.x.permut, dist.z, method = method)	
    	    ryz.temp <- cor(dist.y.permut, dist.z, method = method)	
		}
		if(!permute.my2 & permute.mz2){
	        rxy.temp <- cor(dist.x.permut, dist.y, method = method)
    	    rxz.temp <- cor(dist.x.permut, dist.z.permut, method = method)	
    	    ryz.temp <- cor(dist.y, dist.z.permut, method = method)	
		}
		if(permute.my2 & permute.mz2){
	        rxy.temp <- cor(dist.x.permut, dist.y.permut, method = method)
    	    rxz.temp <- cor(dist.x.permut, dist.z.permut, method = method)	
    	    ryz.temp <- cor(dist.y.permut, dist.z.permut, method = method)	
		}
        value[i,] <- part.cor(rxy.temp, rxz.temp, ryz.temp)
    }
    signif <- (sum(abs(value) >= abs(statistic)) + 1)/(permutations + 1)
    return(list(Obs = statistic, p = signif))
}
