#' Permutate rows in a matrix
#' 
#' Internal function to permutate rows in a matrix.
#' 
#' 
#' @author Vanderlei Júlio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @keywords SYNCSA
permut.row.matrix <-
function(matrix){
	row<-dim(matrix)[1]
	col<-dim(matrix)[2]
	samp<-sample(1:row,row)
	permut.matrix <- matrix[samp,1:col]
	return(permut.matrix)
	}

