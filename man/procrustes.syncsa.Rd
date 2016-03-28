\encoding{UTF-8}
\name{procrustes.syncsa}
\alias{procrustes.syncsa}
\alias{procrustes.partial}

\title{Procrustes and Partial Procrustes correlations.}

\description{

Function to obtain the Procrustes correlation between two matrices and partial Procrustes correlation between three matrices.}

\usage{
procrustes.syncsa(x, y)
procrustes.partial(x, y, z)
}


\arguments{
\item{x, y}{Matrix that will be correlated.}
\item{z}{Matrix whose effect will be removed from the correlation between x and y.}	
}

\details{
The function procrustes.syncsa is a small change in the function \code{\link{procrustes}}. The function uses a correlation-like statistic derived from the symmetric Procrustes sum of squares. Partial Procrustes correlations is obtained by the Procrustes correlation of pairs of matrices and applied the first-order partial correlation coefficient, in the same way that the partial Mantel.
}

\value{
Procrustes correlation between matrices.
}

\references{
Peres-Neto, P.R. and Jackson, D.A. (2001). How well do multivariate data sets match? The advantages of a Procrustean superimposition approach over the Mantel test. Oecologia 129: 169-178.
}
		
\author{
Vanderlei Júlio Debastiani <vanderleidebastiani@yahoo.com.br>
}

\keyword{SYNCSA}