\encoding{UTF-8}
\name{organize.syncsa}
\alias{organize.syncsa}

\title{Function for organize data for Package SYNCSA}

\description{
Package \strong{SYNCSA} requires that the species and community sequence in the data matrix must be the same for all matrices.

The Function organizes the data for the functions of the package \strong{SYNCSA}, placing the matrices of community, traits, phylogenetic distance and environmental varibles in the same order.
}

\usage{
organize.syncsa(comm, traits, dist.spp, envir)
}

\arguments{
\item{comm}{Community data, with species as columns and sampling units as rows.}
\item{traits}{Matrix data of species described by traits, with traits as columns and species as rows.}
\item{dist.spp}{Matrix containing phylogenetic distance between species. Must be a complete matrix (not a half diagonal matrix).This matrix can be larger than community data (more species) as long as it has at least all species that are in community data.}
\item{envir}{Environmental variables for each community, with variables as columns and sampling units as rows.}
}

\note{
The function organizes the matrices despite the absence of one of the matrices, provided that the community data had been entered. Unspecified matrices will appear as NULL.
}

\value{
The matrices of community, traits, phylogenetic distance and environmental variables. The objects returned belong to the matrix class.
}

\author{
Vanderlei Júlio Debastiani <vanderleidebastiani@yahoo.com.br>
}

\seealso{ \code{\link{matrix.t}}, \code{\link{matrix.x}}, \code{\link{matrix.p}}, \code{\link{syncsa}}
}

\keyword{SYNCSA}
