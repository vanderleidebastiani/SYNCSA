#' Function to obtain the partial matrix correlation between three matrices.
#' 
#' Function to obtain the partial matrix correlation between three matrices,
#' similar to the function mantel.partial. The significance of the statistics
#' is evaluated differently from mantel.partial. For more details, see
#' \code{\link{syncsa}}.
#' 
#' The null model is based on permutations in the matrix m2, typically the
#' matrices B, U and Q.
#' 
#' Null model described by Pillar et al. (2009) and Pillar & Duarte (2010). For
#' more details on the matrices and the null model, see \code{\link{syncsa}}.
#' 
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
#' @param na.rm Logical argument (TRUE or FALSE) to specify if pairwise
#' deletion of missing observations when computing dissimilarities (Default
#' na.rm = FALSE).
#' @return \item{Obs}{Correlation between matrices.} \item{p}{Significance
#' level based on permutations.}
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{cor.matrix}}, \code{\link{syncsa}}
#' @references Pillar, V.D.; Duarte, L.d.S. (2010). A framework for
#' metacommunity analysis of phylogenetic structure. Ecology Letters, 13,
#' 587-596.
#' 
#' Pillar, V.D., Duarte, L.d.S., Sosinski, E.E. & Joner, F. (2009).
#' Discriminating trait-convergence and trait-divergence assembly patterns in
#' ecological community gradients. Journal of Vegetation Science, 20, 334?348.
#' @keywords SYNCSA
#' @export
cor.matrix.partial<-function (m1, m2, x, y, z, method = "pearson", dist = "euclidean", permutations = 999, norm = FALSE, na.rm = FALSE) 
{
    permut <- permutations
    m1<-as.matrix(m1)
    m2<-as.matrix(m2)
    x<-as.matrix(x)
    y<-as.matrix(y)
    z<-as.matrix(z)
    part.cor <- function(rxy, rxz, ryz) {
        (rxy - rxz * ryz)/sqrt(1 - rxz * rxz)/sqrt(1 - ryz * ryz)
    }
    dist.x <- vegdist(x, method = dist, na.rm = na.rm)
    dist.y <- vegdist(y, method = dist, na.rm = na.rm)
    dist.z <- vegdist(z, method = dist, na.rm = na.rm)
    rxy <- cor(dist.x, dist.y, method = method)
    rxz <- cor(dist.x, dist.z, method = method)
    ryz <- cor(dist.y, dist.z, method = method)
    statistic <- part.cor(rxy, rxz, ryz)
    if ((rxz == 1 | ryz == 1) == TRUE) {
        statistic <- 0
    }
    value <- matrix(NA, nrow = 1, ncol = permut)
    for (i in 1:permut) {
        m2.permut <- permut.row.matrix(m2)
        x.permut <- m1 %*% m2.permut
        if (norm == "TRUE") {
            matrix.permut <- apply(x.permut^2, 2, sum)
            x.permut <- sweep(x.permut, 2, sqrt(matrix.permut), "/")
        }
        dist.x.permut <- vegdist(x.permut, method = dist, na.rm = na.rm)
        rxy <- cor(dist.x.permut, dist.y, method = method)
        rxz <- cor(dist.x.permut, dist.z, method = method)
        value[i] <- part.cor(rxy, rxz, ryz)
        if ((rxz == 1 | ryz == 1) == TRUE) {
            value[i] <- 0
        }
    }
    signif <- (sum(value >= statistic) + 1)/(permut + 1)
    return(list(Obs = statistic, p = signif))
}
