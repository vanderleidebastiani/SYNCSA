\encoding{UTF-8}
\name{cor.matrix.partial}
\alias{cor.matrix.partial}

\title{
Function to obtain the partial matrix correlation between three matrices.
}

\description{
Function to obtain the partial matrix correlation between three matrices, similar to the function mantel.partial. The significance of the statistics is evaluated differently from mantel.partial. For more details, see \code{\link{syncsa}}.
}

\usage{
cor.matrix.partial(m1, m2, x, y, z, method = "pearson", dist = "euclidean",
		permutations = 999, norm = FALSE, na.rm = FALSE)
}

\arguments{
\item{m1}{Matrix that multiplied by m2 results in the matrix x.}
\item{m2}{Matrix that when multiplied by m1 results in the matrix x. See `details` below.}
\item{x}{Matrix obtained by multiplication of m1 and m2.}
\item{y}{Matrix that will be correlated with the matrix x.}
\item{z}{Matrix whose effect will be removed from the correlation between x and y.}
\item{method}{Correlation method, as accepted by cor: "pearson", "spearman" or "kendall".}
\item{dist}{Dissimilarity index, as accepted by vegdist: "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial" or "chao". However, some of these will not make sense in this case.}
\item{permutations}{Number of permutations in assessing significance.}
\item{norm}{Logical argument (TRUE or FALSE) to specify if x is standardized within variables (Default norm = FALSE).}
\item{na.rm}{Logical argument (TRUE or FALSE) to specify if pairwise deletion of missing observations when computing dissimilarities (Default na.rm = FALSE).}
}

\details{
The null model is based on permutations in the matrix m2, typically the matrices B, U and Q. 

Null model described by Pillar et al. (2009) and Pillar & Duarte (2010). For more details on the matrices and the null model, see \code{\link{syncsa}}.
}

\value{
\item{Obs}{Correlation between matrices.}
\item{p}{Significance level based on permutations.}
}

\references{
Pillar, V.D.; Duarte, L.d.S. (2010). A framework for metacommunity analysis of phylogenetic structure. Ecology Letters, 13, 587-596.

Pillar, V.D., Duarte, L.d.S., Sosinski, E.E. & Joner, F. (2009). Discriminating trait-convergence and trait-divergence assembly patterns in ecological community gradients. Journal of Vegetation Science, 20, 334–348.
}

\author{
Vanderlei Júlio Debastiani <vanderleidebastiani@yahoo.com.br>
}

\seealso{ \code{\link{cor.matrix}}, \code{\link{syncsa}}  }

\keyword{SYNCSA}