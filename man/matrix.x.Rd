\encoding{UTF-8}
\name{matrix.x}
\alias{matrix.x}

\title{Matrix X}

\description{
Function to obtain a matrix containing trait-weighted species composition. For more details, see \code{\link{syncsa}}.
}

\usage{
matrix.x(comm, traits, scale = TRUE, notification = TRUE)
}

\arguments{
\item{comm}{Community data, with species as columns and sampling units as rows. This matrix can contain either presence/absence or abundance data.}
\item{traits}{Matrix data of species described by traits, with traits as columns and species as rows.}
\item{scale}{Logical argument (TRUE or FALSE) to specify if the traits are measured on different scales (Default scale = TRUE). Scale = TRUE if traits are measured on different scales, the matrix of traits is subjected to standardization within each trait, and  Gower Index is used to calculate the degree of belonging to the species. Scale = FALSE if traits are measured on the same scale, the matrix of traits is not subjected to standardization, and Euclidean distance is calculated to determine the degree of belonging to the species.}
\item{notification}{Logical argument (TRUE or FALSE) to specify if notifications of missing observations are shown (Default notification = TRUE).}
}

\value{
\item{matriz.w}{Standardized community matrix, where rows are communities and columns species. Row totals (communities) = 1.}
\item{matriz.u}{Standardized matrix containing the degree of belonging of each species in relation to each other species. Row totals (species) = 1.}
\item{matriz.X}{Trait-weighted species composition matrix. Row totals (communities) = 1.}
}

\note{
\strong{IMPORTANT}: The sequence species show up in community data matrix MUST be the same as they show up in traits matrix. See \code{\link{organize.syncsa}}.
}

\references{
Pillar, V.D.; Duarte, L.d.S. (2010). A framework for metacommunity analysis of phylogenetic structure. Ecology Letters, 13, 587-596.

Pillar, V.D., Duarte, L.d.S., Sosinski, E.E. & Joner, F. (2009). Discriminating trait-convergence and trait-divergence assembly patterns in ecological community gradients. Journal of Vegetation Science, 20, 334–348.
}

\author{
Vanderlei Júlio Debastiani <vanderleidebastiani@yahoo.com.br>
}

\seealso{ \code{\link{matrix.t}}, \code{\link{matrix.p}}, \code{\link{syncsa}}, \code{\link{organize.syncsa}}  }

\examples{
data(flona)
matrix.x(flona$community,flona$traits,scale=TRUE)
}

\keyword{SYNCSA}