\encoding{UTF-8}
\name{rao.diversity}
\alias{rao.diversity}

\title{Rao's quadratic entropy}

\description{
Calculates Rao's quadratic entropy, functional and phylogenetic redundancy.
}

\usage{
rao.diversity(comm, traits = NULL, phylodist = NULL, checkdata = TRUE, ...)
}

\arguments{
\item{comm}{Community data, with species as columns and sampling units as rows. This matrix can contain either presence/absence or abundance data.}
\item{traits}{Matrix data of species described by traits, with traits as columns and species as rows (optional).}
\item{phylodist}{Matrix containing phylogenetic distance between species (optional).}
\item{checkdata}{Logical argument (TRUE or FALSE) to check if species sequence in the community data follows the same order as the one in the trait and in the phylodist matrices (Default checkdata = TRUE).}
\item{...}{Parameters for \code{\link{gowdis}} function.}
}

\details{
Rao's quadratic entropy is a measure of diversity of ecological communities defined by Rao (1982) and is based on the proportion of the abundance of species present in a community and some measure of dissimilarity among them. The dissimilarity range from 0 to 1 and is based on a set of specified functional traits or in the phylogenetic dissimilarity.

For the trait data , the function calculates the square root of the one-complement of Gower’s similarity index, in order to have a dissimilarity matrix with Euclidean metric properties. Gower’s index ranges from 0 to 1 and can handle traits measured indifferent scales. When the species are completely different in terms of their traits, Rao quadratic entropy is equivalent to the Gini–Simpson index.

Functional redundancy is defined purely as the difference between species diversity and Rao’s quadratic entropy based on their functional dissimilarity (de Bello et al. 2007). The same definition is used for phylogenetic redundancy.
}

\value{
\item{Simpson}{Gini-Simpson index within each community (equivalent to Rao quadratic entropy with null, crisp, similarities).}
\item{FunRao}{Rao quadratic entropy within each community, considering trait distance.}
\item{FunRedundancy}{Functional redundancy in each community.}
\item{PhyRao}{Rao quadratic entropy within each community, considering phylogenetic distance.}
\item{PhyRedundancy}{Phylogenetic redundancy in each community.}
}

\note{
\strong{IMPORTANT}: The sequence species show up in community data matrix MUST be the same as they show up in traits and phylodist matrices. See \code{\link{organize.syncsa}}.
}

\references{
de Bello, F.; Leps, J.; Lavorel, S. & Moretti, M. (2007). Importance of species abundance for assessment of trait composition: an example based on pollinator communities. Community Ecology, 8, 163–170.

Pillar, V.D.; Blanco, C.C.; Muler, S.C.; Sosinski, E.E.; Joner, F. & Duarte, L.d.S. (2013). Functional redundancy and stability in plant communities. Journal of Vegetation Science, 24, 963-974.

Rao, C.R. (1982). Diversity and dissimilarity coefficients: a unified approach. Theoretical Population Biology, 21, 24–43.
}

\author{
Vanderlei Júlio Debastiani <vanderleidebastiani@yahoo.com.br>
}

\seealso{\code{\link{organize.syncsa}}, \code{\link{gowdis}}, \code{\link{syncsa}}  }

\examples{
data(flona)
rao.diversity(flona$community)
rao.diversity(flona$community,traits=flona$traits)
}

\keyword{SYNCSA}