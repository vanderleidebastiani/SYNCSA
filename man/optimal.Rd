\encoding{UTF-8}
\name{optimal}
\alias{optimal}

\title{
Searching for optimal traits
}

\description{
Maximize trait-convergence assembly patterns (TCAP = roTE), trait-divergence assembly patterns (TDAP = roXE.T) or maximize both trait-divergence assembly patterns and trait-convergence assembly patterns (TCAP.TDAP = roXE). For more details, see \code{\link{syncsa}}.
}

\usage{
optimal(comm, envir, traits, subset.min = 2, subset.max = 3, pattern = "tcap" , 
		dist = "euclidean", method = "pearson", scale = TRUE, scale.envir = TRUE,
		na.rm = FALSE, notification = TRUE, progressbar = FALSE)
}

\arguments{
\item{comm}{Community data, with species as columns and sampling units as rows. This matrix can contain either presence/absence or abundance data.}
\item{traits}{Matrix data of species described by traits, with traits as columns and species as rows.}
\item{envir}{Environmental variables for each community, with variables as columns and sampling units as rows.}
\item{subset.min}{Minimum of traits in each subset (Default subset.min=2).}
\item{subset.max}{Maximum of traits in each subset (Default subset.max=3).}
\item{pattern}{Patterns for maximize correlation, "tcap","tdap" or "tcap.tdap" (Default pattern="tcap").}
\item{method}{Correlation method, as accepted by cor: "pearson", "spearman" or "kendall".}
\item{dist}{Dissimilarity index, as accepted by vegdist: "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial" or "chao".}
\item{scale}{Logical argument (TRUE or FALSE) to specify if the traits are measured on different scales (Default Scale = TRUE). Scale = TRUE if traits are measured on different scales, the matrix T is subjected to standardization within each trait. Scale = FALSE if traits are measured on the same scale, the matrix T is not subjected to standardization. Furthermore, if Scale = TRUE the matrix of traits is subjected to standardization within each trait, and  Gower Index is used to calculate the degree of belonging to the species, and if Scale = FALSE the matrix of traits is not subjected to standardization, and Euclidean distance is calculated to determine the degree of belonging to the species.}
\item{scale.envir}{Logical argument (TRUE or FALSE) to specify if the environmental variables are measured on different scales (Default Scale = TRUE). If the enviromental variables are measured on different scales, the matrix is subjected to centralization and standardization within each variable.}
\item{na.rm}{Logical argument (TRUE or FALSE) to specify if pairwise deletion of missing observations when computing dissimilarities (Default na.rm = FALSE).}
\item{notification}{Logical argument (TRUE or FALSE) to specify if notifications of missing observations are shown (Default notification = TRUE).}
\item{progressbar}{Logical argument (TRUE or FALSE) to specify if display a progress bar on the R console (Default progressbar = FALSE).}
}

\value{
\item{Subset}{Subset of traits that maximizes the correlation.}
\item{ro}{Correlation for the subset of traits.}
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

\seealso{ \code{\link{syncsa}}, \code{\link{organize.syncsa}}  }

\examples{
data(flona)
optimal(flona$community,flona$environment,flona$traits,subset.min=3,subset.max=5,pattern="tcap")
optimal(flona$community,flona$environment,flona$traits,subset.min=3,subset.max=5,pattern="tdap")
optimal(flona$community,flona$environment,flona$traits,
	subset.min=3,subset.max=5,pattern="tcap.tdap")
}

\keyword{SYNCSA}