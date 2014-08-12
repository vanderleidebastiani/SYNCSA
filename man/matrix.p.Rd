\encoding{UTF-8}
\name{matrix.p}
\alias{matrix.p}

\title{Matrix P}

\description{
Function to obtain a matrix containing phylogeny-weighted species composition. For more details, see \code{\link{syncsa}}.
}

\usage{
matrix.p(comm, dist.spp, notification = TRUE)
}

\arguments{
\item{comm}{Community data, with species as columns and sampling units as rows. This matrix can contain either presence/absence or abundance data.}
\item{dist.spp}{Matrix containing phylogenetic distance between species. Must be a complete matrix (not a half diagonal matrix).}
\item{notification}{Logical argument (TRUE or FALSE) to specify if notifications of missing observations are shown (Default notification = TRUE).}
}

\value{
\item{matrix.w}{Standardized community matrix, where rows are communities and columns species. Row totals (communities) = 1.}
\item{matrix.q}{Standardized matrix containing the degree of belonging of species in relation to each other. Row totals (species) = 1.}
\item{matrix.P}{Phylogeny-weighted species composition matrix. Row totals (communities) = 1.}
}

\note{
\strong{IMPORTANT}: The sequence species show up in community data matrix MUST be the same as they show up in phylogenetic distance matrix. See \code{\link{organize.syncsa}}.
}

\references{
Pillar, V.D.; Duarte, L.d.S. (2010). A framework for metacommunity analysis of phylogenetic structure. Ecology Letters, 13, 587-596.
}

\author{
Vanderlei Júlio Debastiani <vanderleidebastiani@yahoo.com.br>
}

\seealso{ \code{\link{matrix.t}}, \code{\link{matrix.x}}, \code{\link{syncsa}}, \code{\link{organize.syncsa}}  }

\examples{
data(flona)
matrix.p(flona$community,flona$phylo)
}

\keyword{SYNCSA}
