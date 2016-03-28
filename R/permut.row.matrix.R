permut.row.matrix <- function(matrix, strata = NULL){
	samp<-1:dim(matrix)[1]
	if(is.null(strata)){
		samp<-sample(samp)
	}
	if(!is.null(strata)){
		if(length(samp)!=length(strata)){
			stop("\n The strata must be the same length of number of rows\n")
		}
		inds<-levels(as.factor(strata))
		for(is in inds) {
			gr <- samp[strata == is]
			if (length(gr) > 1){
				samp[gr] <- sample(gr)
			}
		}
	}
	samp
	permut.matrix <- matrix[samp,]
	return(permut.matrix)
}