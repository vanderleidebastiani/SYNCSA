% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/permut.vector.R
\encoding{UTF-8}
\name{permut.vector}
\alias{permut.vector}
\title{Permutate a vector}
\usage{
permut.vector(n, strata = NULL, nset = 999)
}
\arguments{
\item{n}{The length of vector.}

\item{strata}{A vector to specify restricting permutations.}

\item{nset}{The number of permutations to generate for the set (Default strata = NULL).}
}
\value{
A matrix of permutations, where each row is a separate permutation.
}
\description{
Function to permutate a vector of size n using the function \code{\link{shuffleSet}}.
}
\seealso{
\code{\link{syncsa}}, \code{\link{permut.row.matrix}}
}
\author{
Vanderlei Julio Debastiani <vanderleidebastiani@yahoo.com.br>
}
\keyword{Auxiliary}
