#' @rdname cor.matrix
#' @encoding UTF-8
#' @export
cor.matrix.partial<-function (m1, m2, x, y, z, method = "pearson", dist = "euclidean", permutations = 999, norm = FALSE, strata = NULL, na.rm = FALSE) 
{
    m1<-as.matrix(m1)
    m2<-as.matrix(m2)
    x<-as.matrix(x)
    y<-as.matrix(y)
    z<-as.matrix(z)
    dist.x <- vegdist(x, method = dist, na.rm = na.rm)
    dist.y <- vegdist(y, method = dist, na.rm = na.rm)
    dist.z <- vegdist(z, method = dist, na.rm = na.rm)
    rxy <- cor(dist.x, dist.y, method = method)
    rxz <- cor(dist.x, dist.z, method = method)
    ryz <- cor(dist.y, dist.z, method = method)
    statistic <- part.cor(rxy, rxz, ryz)
    value <- matrix(NA, nrow = permutations, ncol = 1)
    for (i in 1:permutations) {
        m2.permut <- permut.row.matrix(m2, strata = strata)
        x.permut <- m1 %*% m2.permut
        if (norm == "TRUE") {
            matrix.permut <- apply(x.permut^2, 2, sum)
            x.permut <- sweep(x.permut, 2, sqrt(matrix.permut), "/")
        }
        dist.x.permut <- vegdist(x.permut, method = dist, na.rm = na.rm)
        rxy <- cor(dist.x.permut, dist.y, method = method)
        rxz <- cor(dist.x.permut, dist.z, method = method)
        value[i,] <- part.cor(rxy, rxz, ryz)
    }
    signif <- (sum(abs(value) >= abs(statistic)) + 1)/(permutations + 1)
    return(list(Obs = statistic, p = signif))
}
