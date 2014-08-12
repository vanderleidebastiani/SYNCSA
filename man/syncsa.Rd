\encoding{UTF-8}
\name{syncsa}
\alias{syncsa}

\title{Syncsa
}

\description{
This function integrates several steps for the analysis of phylogenetic assembly patterns and their links to traits and ecological processes in a metacommunity (Pillar et al. 2009, Pillar & Duarte 2010). It requires data organized into the following matrices: (1) the presences or abundances of species in a set of communities (\strong{W}); (2) the phylogenetic pairwise dissimilarities of these species (\strong{DF}, in the range 0 to 1); (3) a set of functional traits describing the species (\strong{B}), which may be a mixture of binary and quantitative traits (ordinal, interval, ratio scales), but not nominal ones (these should be expanded into binary traits); and (4) the ecological gradient of interest, which may be one or more factors to which the communities respond or ecosystem effects of the communities (\strong{E}). The function computes several matrix correlations that express trait-convergence assembly patterns (TCAP), trait-divergence assembly patterns (TDAP), and phylogenetic signal in functional traits at the species pool level and at the metacomunity level. This function also generates P-values by permutation testing based on null models (Pillar et al. 2009, Pillar & Duarte 2010).

The function implement methods that have been available in the SYNCSA application written in C++ (by Valério Pillar, available at http://ecoqua.ecologia.ufrgs.br/ecoqua/SYNCSA.html).
}

\usage{
syncsa(comm, traits, dist.spp, envir, method = "pearson", dist = "euclidean",
		scale = TRUE, scale.envir = TRUE, permutations = 999, na.rm = FALSE, 
		notification = TRUE)
}

\arguments{
\item{comm}{Community data, with species as columns and sampling units as rows. This matrix can contain either presence/absence or abundance data.}
\item{traits}{Matrix data of species described by traits, with traits as columns and species as rows.}
\item{dist.spp}{Matrix containing phylogenetic distance between species. Must be a complete matrix (not a half diagonal matrix).}
\item{envir}{Environmental variables for each community, with variables as columns and sampling units as rows.}
\item{method}{Correlation method, as accepted by cor: "pearson", "spearman" or "kendall".}
\item{dist}{Dissimilarity index, as accepted by vegdist: "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial" or "chao". However, some of these will not make sense in this case.}
\item{scale}{Logical argument (TRUE or FALSE) to specify if the traits are measured on different scales (Default Scale = TRUE). Scale = TRUE if traits are measured on different scales, the matrix T is subjected to standardization within each trait. Scale = FALSE if traits are measured on the same scale, the matrix T is not subjected to standardization. Furthermore, if Scale = TRUE the matrix of traits is subjected to standardization within each trait, and  Gower Index is used to calculate the degree of belonging to the species, and if Scale = FALSE the matrix of traits is not subjected to standardization, and Euclidean distance is calculated to determine the degree of belonging to the species.}
\item{scale.envir}{Logical argument (TRUE or FALSE) to specify if the environmental variables are measured on different scales (Default Scale = TRUE). If the enviromental variables are measured on different scales, the matrix is subjected to centralization and standardization within each variable.}
\item{permutations}{Number of permutations in assessing significance.}
\item{na.rm}{Logical argument (TRUE or FALSE) to specify if pairwise distances should be deleted in cases of missing observations (Default na.rm = FALSE).}
\item{notification}{Logical argument (TRUE or FALSE) to specify if notification of missing observations should to be shown (Default notification = TRUE).}
}

\details{
\strong{ro(TE)}

This matrix correlation refers to trait-convergence assembly patterns related to the ecological gradient (TCAP, Pillar et al. 2009). For evaluating TCAP, by matrix multiplication we define \strong{T = WB}, which with previous standardization of \strong{W} to unit column totals will contain the trait averages in each community. The elements in \strong{T} are community weighted mean values or community functional parameters (Violle et al. 2007). Standardization of the traits (rows) in \strong{T} is needed if the trait set contains traits measured with different scales. By using matrix correlation, we evaluate how the trait patterns in \strong{T} are associated to ecological gradients in \strong{E}. For relating \strong{T} to \strong{E}, we define a distance matrix of the communities (\strong{DT}) using \strong{T}, and another distance matrix of the community sites (\strong{DE}) using \strong{E}. The matrix correlation ro(\strong{TE}) = ro(\strong{DT};\strong{DE}) measures the level of congruence between TCAP and \strong{E}. A strong correlation ro(\strong{TE}) indicates the factors directly or indirectly represented in \strong{E} are involved in ecological filtering of species that, at least for the traits considered in the analysis, consistently produce trait-convergence assembly patterns along the gradient comprising the metacommunity.

\strong{ro(XE) and ro(XE.T)}

These matrix correlations refer to trait-divergence assembly patterns related to the ecological gradient (TDAP, Pillar et al. 2009). For the identification of TDAP, in a first step the species pairwise similarities (in the range 0 to 1) in matrix \strong{SB} based on traits in \strong{B} are used to define matrix \strong{U} with degrees of belonging of species to fuzzy sets. By matrix multiplication \strong{X = WU} will contain the species composition of the communities after fuzzy-weighting by their trait similarities (each row in \strong{X} will refer to a species). Matrix \strong{X} expresses both TCAP and TDAP (Pillar et al. 2009). By using matrix correlation, we evaluate how the trait patterns in \strong{X} (TCAP and TDAP) are associated to ecological gradients in \strong{E}. For relating \strong{X} to \strong{E}, we define a distance matrix of the communities (\strong{DX}) using \strong{X}, and another distance matrix of the community sites (\strong{DE}) using \strong{E}. The matrix correlation ro(\strong{XE}) = ro(\strong{DX};\strong{DE}) between \strong{X} and \strong{E} is defined. We then remove the trait-convergence component ro(\strong{TE}) from ro(\strong{XE}) by computing the partial matrix (Mantel) correlation ro(\strong{XE.T}), which measures the level of congruence between TDAP and \strong{E}. Trait-divergence assembly patterns (TDAP, Pillar et al. 2009) may result from community assembly processes related to biotic interactions (Stubbs & Wilson 2004; Wilson 2007). 

\strong{ro(PE)}

This matrix correlation refers to the phylogenetic structure  related to the ecological gradient comprising the metacommunity. The phylogenetic pairwise dissimilarities in \strong{DF} are transformed into similarities and used to define degrees of belonging qij to fuzzy sets. This is analogous to the definition of functional fuzzy sets (Pillar & Orlóci 1991; Pillar et al. 2009). Based on the phylogenetic similarities, every species i among s species in the metacommunity specifies a fuzzy set to which every species j (j = 1 to s species, including species i) belongs with a certain degree of belonging in the interval [0, 1]. In our definition, each row in matrix \strong{Q} with the degrees of belonging must add to unit, i.e., the degrees of belonging of a given species across the fuzzy sets are standardized to unit total. By matrix multiplication \strong{P = WQ} will contain the composition of the communities after fuzzy-weighting of species presences or abundances by the species' phylogenetic similarities. Each column in matrix \strong{P} holds the phylogenetic structure of a community. The standardization of \strong{Q} is essential for the community totals in each column in \strong{W} remaining the same in \strong{P}. Further, matrix \strong{W} is adjusted to unit column totals prior to multiplication, so that the total richness or abundance within each community in \strong{W} will be standardized. Matrix correlation ro(\strong{PE}) = ro(\strong{DP};\strong{DE}) measures the strength of the association between community distances based on their phylogenetic structure in \strong{DP} and  distances based on their ecological conditions (\strong{DE}). Further, \strong{P} can be explored for phylogenetic patterns at the metacommunity level by using, e.g., ordination techniques. 

\strong{ro(PT) and ro(PX.T)} 

These matrix correlations measure phylogenetic signal at the metacommunity level related to TCAP and to TDAP. We define phylogenetic signal at the metacommunity level related to TCAP (PSMT) as the correlation between the phylogenetic structure described in matrix \strong{P} and the trait-convergence structure described in matrix \strong{T}. For this, a proper distance matrix (e.g. Euclidean distances) of communities (\strong{DP}) is computed using \strong{P} and another distance matrix of the same communities (\strong{DT}) is computed using \strong{T}. Then matrix correlation ro(\strong{PT}) = ro(\strong{DP};\strong{DT}) will measure the level of congruence between variation in \strong{P} and \strong{T}, which is a measure of PSMT. A strong phylogenetic signal at the metacommunity level is expected when communities that are more similar in terms of phylogenetic structure are also similar regarding their average trait values.
We also define phylogenetic signal at the metacommunity level related to TDAP (PSMX.T) as the partial matrix correlation ro(\strong{PX.T}) = ro(\strong{DP};\strong{DX.DT}) between community distances DP computed on phylogenetic structure and community distances \strong{DX} computed on species composition after fuzzy-weighting by the species‚ or trait similarities, removing the effect of TCAP (\strong{DT}). This is analogous to TDAP.

\strong{ro(BF)}

This matrix correlation measures phylogenetic signal at the species pool level (PSS, Pillar & Duarte 2010). We define PSS as the matrix correlation ro(\strong{FB}) = ro(\strong{DF};\strong{DB}) between species phylogenetic dissimilarities (already defined as matrix \strong{DF}) and species trait dissimilarities (derived from already defined matrix \strong{SB}) computed on any number of traits from matrix \strong{B}. The species pool refers to the species present in the metacommunity.

\strong{Additional matrix correlations} 

The matrix correlations ro(\strong{TE.P}) and ro(\strong{XE.P}) are also computed, which may be useful for evaluating causal models in path analysis.

\strong{Testing against null models}

All the matrix correlations are tested against null models. The null model is defined accoding to the correlation being tested. For ro(\strong{TE}), each permutation generates a random matrix \strong{T} calculated after the permutation among the species vectors in matrix \strong{B}. For ro(\strong{XE}) and ro(\strong{XE.T}), each permutation generates a random matrix \strong{X} after the permutation among species fuzzy sets (rows) in matrix \strong{U}. For ro(\strong{PE}), ro(\strong{PT}), and ro(\strong{PX.T}), each permutation generates a random matrix \strong{P} after the permutation among species fuzzy sets (rows) in matrix \strong{Q}. For ro(\strong{BF}), a conventional Mantel test is performed with dissimilarity matrices \strong{DF} and \strong{DB}. Analogous null models are used for testing the additional matrix correlations; that is, the same null model for ro(\strong{TE}) is used for ro(\strong{TE.P}), the same model for ro(\strong{XE}) is used for ro(\strong{XE.P}).
}

\value{
Correlations roTE, roXE, roPE, roPT, roPX.T, roXE.T, roTE.P, roXE.P and roBF, and their significance levels based on permutations. 
}

\note{
The function calculates the correlations despite the lack of one of the matrices, provided that community data had been entered. Correlations including unspecified matrices will appear with ro = 0.

\strong{IMPORTANT}: The sequence of species in the community data matrix MUST be the same as that in the phylogenetic distance matrix and in traits matrix. Similarly, the sequence of communities in the community data matrix MUST be the same as that in the environmental data matrix. See \code{\link{organize.syncsa}}.

The functions ignore missing data when specified. In the case of direct multiplication of matrices (matrices \strong{W} and matrix \strong{B}) the missing data are replaced by 0, ignoring the cell with missing value. Result matrices are shown without missing values. Where the matrices are calculated using a dissimilarity index (matrix \strong{U} and correlations between matrices) the missing data are ignored as in \code{\link{vegdist}} function. In some cases the dissimilarity matrices obtained by the function \code{\link{vegdist}} still contain some missing values. In these cases the rest of the procedure will be affected. In these cases you can find solutions in the package \code{\link{mice}}.

}

\references{

Pillar, V.D.; Duarte, L.d.S. (2010). A framework for metacommunity analysis of phylogenetic structure. Ecology Letters, 13, 587-596.

Pillar, V.D., Duarte, L.d.S., Sosinski, E.E. & Joner, F. (2009). Discriminating trait-convergence and trait-divergence assembly patterns in ecological community gradients. Journal of Vegetation Science, 20, 334–348.

Pillar, V.D. & Orlóci, L. (1991). Fuzzy components in community level comparisons. In: Computer Assisted Vegetation Analysis (eds Feoli, E. & Orlóci, L.). Kluwer, Dordrecht, pp. 87–93.

Stubbs, W.J. & Wilson, J.B. (2004). Evidence for limiting similarity in a sand dune community. Journal of Ecology, 92, 557–567.

Violle, C., Navas, M.L., Vile, D., Kazakou, E., Fortunel, C., Hummel, I. & Garnier, E. (2007). Let the concept of trait be functional! Oikos, 116, 882–892.

Wilson, J.B. (2007). Trait-divergence assembly rules have been demonstrated: limiting similarity lives! A reply to Grime. Journal of Vegetation Science, 18, 451–452.

}

\author{
Vanderlei Júlio Debastiani <vanderleidebastiani@yahoo.com.br>
}

\seealso{ \code{\link{matrix.t}}, \code{\link{matrix.x}}, \code{\link{matrix.p}}, \code{\link{optimal}}, \code{\link{belonging}}, \code{\link{organize.syncsa}}, \code{\link{rao.diversity}} 
}

\examples{
data(flona)
syncsa(comm=flona$community,traits=flona$traits,dist.spp=flona$phylo,envir=flona$environment)
syncsa(flona$community,traits=flona$traits,envir=flona$environment)
}

\keyword{SYNCSA}
