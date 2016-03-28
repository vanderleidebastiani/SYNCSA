cor.matrix<-function (m1, m2, x, y, method = "pearson", dist = "euclidean", permutations = 999, norm = FALSE, strata = NULL, na.rm = FALSE) 
{
    m1<-as.matrix(m1)
    m2<-as.matrix(m2)
    x<-as.matrix(x)
    y<-as.matrix(y)
    dist.y <- vegdist(y, method = dist, na.rm = na.rm)
    dist.x <- vegdist(x, method = dist, na.rm = na.rm)
    correlation <- cor(dist.x, dist.y, method = method)
    value <- matrix(NA, nrow = permutations, ncol = 1)
    for (i in 1: permutations) {
        m2.permut <- permut.row.matrix(m2, strata = strata)
        x.permut <- m1 %*% m2.permut
        if (norm == "TRUE") {
            matrix.permut <- apply(x.permut^2, 2, sum)
            x.permut <- sweep(x.permut, 2, sqrt(matrix.permut), "/")
        }
        dist.x.permut <- vegdist(x.permut, method = dist , na.rm = na.rm)
        cor.x.permut <- cor(dist.x.permut, dist.y, method = method)
        value[i,] <- cor.x.permut
    }
    signific <- (sum(abs(value) >= abs(correlation)) + 1)/(permutations + 1)
    return(list(Obs = correlation, p = signific))
}
