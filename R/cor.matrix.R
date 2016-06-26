#' Function to obtain the correlation between two matrices and partial matrix
#' correlation between three matrices.
#' 
#' Function obtains the correlation between two matrices or to obtain the
#' partial matrix correlation between three matrices. The functions cor.matrix
#' and cor.matrix.partial are similar the function mantel and partial mantel,
#' although the significance of the statistics is evaluated differently from
#' mantel. The functions pro.matrix and pro.matrix.partial uses symmetric
#' Procrustes as a measure of concordance between data sets. For more details,
#' see \code{\link{procrustes}} and \code{\link{syncsa}}.
#' 
#' The null model is based on permutations in the matrix m2, typically the
#' matrices B, U and Q.
#' 
#' Null model described by Pillar et al. (2009) and Pillar & Duarte (2010). For
#' more details on the matrices and the null model, see \code{\link{syncsa}}.
#' 
#' @encoding UTF-8
#' @importFrom vegan vegdist
#' @aliases cor.matrix cor.matrix.partial pro.matrix pro.matrix.partial
#' @param m1 Matrix that multiplied by m2 results in the matrix x.
#' @param m2 Matrix that when multiplied by m1 results in the matrix x. See
#' `details` below.
#' @param x Matrix obtained by multiplication of m1 and m2.
#' @param y Matrix that will be correlated with the matrix x.
#' @param z Matrix whose effect will be removed from the correlation between x
#' and y.
#' @param method Correlation method, as accepted by cor: "pearson", "spearman"
#' or "kendall".
#' @param dist Dissimilarity index, as accepted by vegdist: "manhattan",
#' "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower",
#' "altGower", "morisita", "horn", "mountford", "raup" , "binomial" or "chao".
#' However, some of these will not make sense in this case.
#' @param permutations Number of permutations in assessing significance.
#' @param norm Logical argument (TRUE or FALSE) to specify if x is standardized
#' within variables (Default norm = FALSE).
#' @param strata Argument to specify restricting permutations within species
#' groups (Default strata = NULL).
#' @param na.rm Logical argument (TRUE or FALSE) to specify if pairwise
#' deletion of missing observations when computing dissimilarities (Default
#' na.rm = FALSE).
#' @return \item{Obs}{Correlation between matrices.} \item{p}{Significance
#' level based on permutations.}
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{procrustes}}, \code{\link{mantel}},
#' \code{\link{syncsa}}
#' @references Pillar, V.D.; Duarte, L.d.S. (2010). A framework for
#' metacommunity analysis of phylogenetic structure. Ecology Letters, 13,
#' 587-596.
#' 
#' Pillar, V.D., Duarte, L.d.S., Sosinski, E.E. & Joner, F. (2009).
#' Discriminating trait-convergence and trait-divergence assembly patterns in
#' ecological community gradients. Journal of Vegetation Science, 20, 334:348.
#' @keywords SYNCSA
#' @export
cor.matrix<-function (m1, m2, x, y, method = "pearson", dist = "euclidean", permutations = 999, norm = FALSE, strata = NULL, na.rm = FALSE) 
{
    m1<-as.matrix(m1)
    m2<-as.matrix(m2)
    x<-as.matrix(x)
    y<-as.matrix(y)
    dist.y <- vegan::vegdist(y, method = dist, na.rm = na.rm)
    dist.x <- vegan::vegdist(x, method = dist, na.rm = na.rm)
    correlation <- cor(dist.x, dist.y, method = method)
    value <- matrix(NA, nrow = permutations, ncol = 1)
    for (i in 1: permutations) {
        m2.permut <- permut.row.matrix(m2, strata = strata)$permut.matrix
        x.permut <- m1 %*% m2.permut
        if (norm == "TRUE") {
            matrix.permut <- apply(x.permut^2, 2, sum)
            x.permut <- sweep(x.permut, 2, sqrt(matrix.permut), "/")
        }
        dist.x.permut <- vegan::vegdist(x.permut, method = dist , na.rm = na.rm)
        cor.x.permut <- cor(dist.x.permut, dist.y, method = method)
        value[i,] <- cor.x.permut
    }
    signific <- (sum(abs(value) >= abs(correlation)) + 1)/(permutations + 1)
    return(list(Obs = correlation, p = signific))
}
