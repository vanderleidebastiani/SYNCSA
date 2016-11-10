#' Permutate rows in a matrix
#' 
#' Internal function to permutate rows in a matrix.
#' 
#' 
#' @encoding UTF-8
#' @param data A matrix
#' @param strata Argument to specify restricting permutations in rows 
#' (Default strata = NULL).
#' @param seqpermutation A pre set permutation vector (Default seqpermutation = NULL).
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @keywords SYNCSA
#' @export
permut.row.matrix <- function(data, strata = NULL, seqpermutation=NULL){
	N<-dim(data)[1]
	if(!is.null(strata) & N!=length(seqpermutation)){
		stop("\n The strata must be the length of number of row in the data\n")
	}
	if(!is.null(seqpermutation) & N!=length(seqpermutation)){
		stop("\n The seqpermutation must be the length of number of row in the data\n")
	}
	if(is.null(seqpermutation)){	
		samp<-permut.vector(N,strata=strata,nset=1)
	} else {
		samp<-seqpermutation
	}	
	permut.matrix <- data[samp,,drop=FALSE]
	res<-list(permut.matrix = permut.matrix,samp = samp)
	return(res)
}