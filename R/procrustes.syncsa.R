#' Procrustes and Partial Procrustes correlations.
#' 
#' 
#' Function to obtain the Procrustes correlation between two matrices and
#' partial Procrustes correlation between three matrices.
#' 
#' The function procrustes.syncsa is a small change in the function
#' \code{\link{procrustes}}. The function uses a correlation-like statistic
#' derived from the symmetric Procrustes sum of squares. Partial Procrustes
#' correlations is obtained by the Procrustes correlation of pairs of matrices
#' and applied the first-order partial correlation coefficient, in the same way
#' that the partial Mantel.
#' 
#' @aliases procrustes.syncsa procrustes.partial
#' @param x,y Matrix that will be correlated.
#' @param z Matrix whose effect will be removed from the correlation between x
#' and y.
#' @return Procrustes correlation between matrices.
#' @author Vanderlei Júlio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @references Peres-Neto, P.R. and Jackson, D.A. (2001). How well do
#' multivariate data sets match? The advantages of a Procrustean
#' superimposition approach over the Mantel test. Oecologia 129: 169-178.
#' @keywords SYNCSA
#' @export
procrustes.syncsa<-function (x, y)
{
	x <- cbind(x)
	y <- cbind(y)
	if (ncol(x) < ncol(y)){
		addcols <- ncol(y) - ncol(x)
		for (i in 1:addcols) x <- cbind(x, 0)
	}
	ctrace <- function(MAT) sum(diag(crossprod(MAT)))
	x <- scale(x, scale = FALSE)
	y <- scale(y, scale = FALSE)
	x <- x/sqrt(ctrace(x))
	y <- y/sqrt(ctrace(y))
	xy <- crossprod(x, y)
	sol <- svd(xy)
	c <- sum(sol$d)/ctrace(y)
	r2 <- ctrace(x) + c * c * ctrace(y) - 2 * c * sum(sol$d)
	ro<-sqrt(1 - r2)
return(ro)
}
