#' @title Permutate rows in a matrix
#'
#' @description Internal function to permutate rows in a matrix.
#'
#' @encoding UTF-8
#' @param data A matrix.
#' @param strata Argument to specify restricting permutations in rows (Default strata = NULL).
#' @param seqpermutation A set of predefined permutation vector (Default seqpermutation = NULL).
#' @return \item{permut.matrix}{The matrix permuted.} \item{samp}{The sequence of permutation.}
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{syncsa}}, \code{\link{permut.vector}}
#' @keywords SYNCSA
#' @export
permut.row.matrix <- function(data, strata = NULL, seqpermutation = NULL)
{
  N <- dim(data)[1]
  if(!is.null(strata) & N!=length(seqpermutation)){
    stop("\n strata must be the length of number of row in the data\n")
  }
  if(!is.null(seqpermutation) & N!=length(seqpermutation)){
    stop("\n seqpermutation must be the length of number of row in the data\n")
  }
  if(is.null(seqpermutation)){
    samp <- permut.vector(N, strata = strata, nset = 1)
  } else {
    samp <- seqpermutation
  }
  permut.matrix <- data[samp, ,drop = FALSE]
  res <- list(permut.matrix = permut.matrix, samp = samp)
  return(res)
}
