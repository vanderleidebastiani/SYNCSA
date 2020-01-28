#' @title Co-inertia and Partial Co-inertia correlations.
#'
#' @description Function to obtain the RV coefficient between two matrices and
#' partial RV coefficient between three matrices.
#'
#' @encoding UTF-8
#' @aliases coinertia.syncsa coinertia.partial.syncsa
#' @param x,y Matrix that will be correlated.
#' @param z Matrix whose effect will be removed from the correlation between x and y.
#' @param scale Standardized variables to unit variance.
#' @return RV coefficient (generalization of the Pearson correlation coefficient) between matrices.
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @references Legendre, P. and Legendre L. (2012). Numerical Ecology. 3nd English edition. Elsevier.
#' @seealso \code{\link{syncsa}}
#' @keywords SYNCSA
#' @export
coinertia.syncsa <- function(x, y, scale = FALSE){
  x <- as.matrix(scale(x, center = TRUE, scale = scale))
  y <- as.matrix(scale(y, center = TRUE, scale = scale))
  RV <- (sum((t(x)%*%y)^2))/sqrt(sum((t(x)%*%x)^2)*sum((t(y)%*%y)^2))
  return(RV)
}
