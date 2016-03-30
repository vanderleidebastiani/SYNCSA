#' Permutate rows in a matrix
#' 
#' Internal function to permutate rows in a matrix.
#' 
#' 
#' @encoding UTF-8
#' @param matrix A matrix
#' @param strata Argument to specify restricting permutations in rows 
#' (Default strata = NULL).
#' @author Vanderlei Júlio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @keywords SYNCSA
#' @export
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
