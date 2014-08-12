\encoding{UTF-8}
\name{belonging}
\alias{belonging}

\title{Degree of belonging of species}

\description{
Function to obtain a matrix containing the degree of belonging of each species based on its ecological or phylogenetic resemblance to other species. For more details, see \code{\link{matrix.p}}, \code{\link{matrix.x}} and \code{\link{syncsa}}.
}

\usage{
belonging(dis, standardize = TRUE)
}

\arguments{
\item{dis}{Matrix containing distance between species.}
\item{standardize}{Logical argument (TRUE or FALSE) to specify if dis must be standardized in values ranging from 0 to 1 (Default standardize = TRUE).}
}

\value{
Standardized matrix containing the degree of belonging of species in relation to each other. Row totals (species) = 1.
}

\references{
Pillar, V.D.; Duarte, L.d.S. (2010). A framework for metacommunity analysis of phylogenetic structure. Ecology Letters, 13, 587-596.
}

\author{
Vanderlei Júlio Debastiani <vanderleidebastiani@yahoo.com.br>
}

\seealso{ \code{\link{matrix.p}}, \code{\link{matrix.x}}, \code{\link{syncsa}}  }

\examples{
data(flona)
belonging(flona$phylo)
}

\keyword{SYNCSA}