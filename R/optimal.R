#' Searching for optimal traits
#' 
#' Maximize trait-convergence assembly patterns (TCAP = roTE), trait-divergence
#' assembly patterns (TDAP = roXE.T) or maximize both trait-divergence assembly
#' patterns and trait-convergence assembly patterns (TCAP.TDAP = roXE). For
#' more details, see \code{\link{syncsa}}.
#' 
#' 
#' @param comm Community data, with species as columns and sampling units as
#' rows. This matrix can contain either presence/absence or abundance data.
#' @param traits Matrix data of species described by traits, with traits as
#' columns and species as rows.
#' @param envir Environmental variables for each community, with variables as
#' columns and sampling units as rows.
#' @param subset.min Minimum of traits in each subset (Default subset.min=2).
#' @param subset.max Maximum of traits in each subset (Default subset.max=3).
#' @param pattern Patterns for maximize correlation, "tcap","tdap" or
#' "tcap.tdap" (Default pattern="tcap").
#' @param method Correlation method, as accepted by cor: "pearson", "spearman"
#' or "kendall".
#' @param dist Dissimilarity index, as accepted by vegdist: "manhattan",
#' "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower",
#' "altGower", "morisita", "horn", "mountford", "raup" , "binomial" or "chao".
#' @param scale Logical argument (TRUE or FALSE) to specify if the traits are
#' measured on different scales (Default Scale = TRUE). Scale = TRUE if traits
#' are measured on different scales, the matrix T is subjected to
#' standardization within each trait. Scale = FALSE if traits are measured on
#' the same scale, the matrix T is not subjected to standardization.
#' Furthermore, if Scale = TRUE the matrix of traits is subjected to
#' standardization within each trait, and Gower Index is used to calculate the
#' degree of belonging to the species, and if Scale = FALSE the matrix of
#' traits is not subjected to standardization, and Euclidean distance is
#' calculated to determine the degree of belonging to the species.
#' @param scale.envir Logical argument (TRUE or FALSE) to specify if the
#' environmental variables are measured on different scales (Default Scale =
#' TRUE). If the enviromental variables are measured on different scales, the
#' matrix is subjected to centralization and standardization within each
#' variable.
#' @param na.rm Logical argument (TRUE or FALSE) to specify if pairwise
#' deletion of missing observations when computing dissimilarities (Default
#' na.rm = FALSE).
#' @param notification Logical argument (TRUE or FALSE) to specify if
#' notifications of missing observations are shown (Default notification =
#' TRUE).
#' @param progressbar Logical argument (TRUE or FALSE) to specify if display a
#' progress bar on the R console (Default progressbar = FALSE).
#' @return \item{Subset}{Subset of traits that maximizes the correlation.}
#' \item{ro}{Correlation for the subset of traits.}
#' @note \strong{IMPORTANT}: The sequence species show up in community data
#' matrix MUST be the same as they show up in traits matrix. See
#' \code{\link{organize.syncsa}}.
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{syncsa}}, \code{\link{organize.syncsa}}
#' @references Pillar, V.D.; Duarte, L.d.S. (2010). A framework for
#' metacommunity analysis of phylogenetic structure. Ecology Letters, 13,
#' 587-596.
#' 
#' Pillar, V.D., Duarte, L.d.S., Sosinski, E.E. & Joner, F. (2009).
#' Discriminating trait-convergence and trait-divergence assembly patterns in
#' ecological community gradients. Journal of Vegetation Science, 20, 334?348.
#' @keywords SYNCSA
#' @examples
#' 
#' data(flona)
#' optimal(flona$community,flona$environment,flona$traits,subset.min=3,subset.max=5,pattern="tcap")
#' optimal(flona$community,flona$environment,flona$traits,subset.min=3,subset.max=5,pattern="tdap")
#' optimal(flona$community,flona$environment,flona$traits,
#' 	subset.min=3,subset.max=5,pattern="tcap.tdap")
#' 
#' @importFrom vegan vegdist
#' @export
optimal<-function (comm, envir, traits, subset.min = 2, subset.max = 3, pattern = "tcap", dist = "euclidean", method = "pearson", scale = TRUE, scale.envir = TRUE , na.rm = FALSE, notification = TRUE, progressbar=FALSE) 
{
    part.cor <- function(rxy, rxz, ryz) {
        (rxy - rxz * ryz)/sqrt(1 - rxz * rxz)/sqrt(1 - ryz * ryz)
    }
    comm<-as.matrix(comm)
    envir<-as.matrix(envir)
    traits<-as.matrix(traits)
    if(notification==TRUE){
    	c.NA <- apply(comm, 2, is.na)
    	if(length(which(unique(as.vector(c.NA))==TRUE))>0){
			warning("Warning: NA in community data",call.=FALSE)	
    	}
		t.NA <- apply(traits, 2, is.na)
    	if(length(which(unique(as.vector(t.NA))==TRUE))>0){
			warning("Warning: NA in traits matrix",call.=FALSE)	
		}
		e.NA <- apply(envir, 2, is.na)
		if(length(which(unique(as.vector(e.NA))==TRUE))>0){
			warning("Warning: NA in environmental data",call.=FALSE)	
    	}
    }
    colnames(traits) <- colnames(traits, do.NULL = FALSE, prefix = "T")
    if (scale.envir == "TRUE") {
        envir <- cent.norm(envir,na.rm = na.rm)
    }
    dist.y <- vegdist(envir, method = dist, na.rm = na.rm)
    m <- dim(traits)[2]
    if (subset.max > m) {
        stop("\n Subset must be lower than the number of traits\n")
    }
    PATTERNS <- c("tcap", "tdap", "tcap.tdap")
    pattern <- pmatch(pattern, PATTERNS)
    if (length(pattern) > 1) {
        stop("\n Only one argument is accepted in pattern \n")
    }
    if (is.na(pattern)) {
        stop("\n Invalid pattern \n")
    }
    p <- 1:subset.max
    bin <- factorial(m)/(factorial(p) * factorial(m - p))
    nT<-sum(bin[subset.min:subset.max])
    comb <- matrix(NA, nrow = sum(bin[subset.min:subset.max]), ncol = 1)
    n=0
    for (i in subset.min:subset.max) {
        combinations <- combn(colnames(traits), i, simplify = TRUE)
        for (j in 1:bin[i]) {
        	n=n+1
            comb[n, 1] <- paste(combinations[,j], collapse = " ")
        }
    }
    n=0
    correlation <- matrix(NA, nrow = sum(bin[subset.min:subset.max]), ncol = 1)
    for (i in subset.min:subset.max) {
        combinations1 <- combn(colnames(traits), i, simplify = TRUE)
        for (j in 1:bin[i]) {
            if (pattern == 1) {
            	n=n+1
                T <- matrix.t(comm, as.matrix(traits[, combinations1[, j]]), scale = scale, notification = FALSE)
                correlation[n, 1] <- cor(vegdist(as.matrix(T$matrix.T), method = dist, na.rm = na.rm), dist.y, method = method)
                if(progressbar){
					ProgressBAR(n,nT,style=3)
				}
            }
            if (pattern == 2) {
            	n=n+1
                T <- matrix.t(comm, as.matrix(traits[, combinations1[, j]]), scale = scale, notification = FALSE)
                X <- matrix.x(comm, as.matrix(traits[, combinations1[, j]]), scale = scale, notification = FALSE)
                dist.x <- vegdist(X$matrix.X, method = dist, na.rm = na.rm)
                dist.z <- vegdist(T$matrix.T, method = dist, na.rm = na.rm)
                rxy <- cor(dist.x, dist.y, method = method)
                rxz <- cor(dist.x, dist.z, method = method)
                ryz <- cor(dist.y, dist.z, method = method)
                correlation[n, 1] <- part.cor(rxy, rxz, ryz)
                if(!is.na(rxz)){
                	if(!is.na(ryz)){
                		if(!is.na(rxy)){
                			if ((rxz == 1 | ryz == 1) == TRUE) {
	                  			correlation[n, 1] <- 0
	                  		}
        	        	}
    	            }
                }
				if(progressbar){
					ProgressBAR(n,nT,style=3)
				}
            }
            if (pattern == 3) {
            	n=n+1
                X <- matrix.x(comm, as.matrix(traits[, combinations1[, j]]), scale = scale, notification = FALSE)
                correlation[n, 1] <- cor(vegdist(as.matrix(X$matrix.X), method = dist, na.rm = na.rm), dist.y, method = method)
				if(progressbar){
					ProgressBAR(n,nT,style=3)
				}
            }
        }
    }
    result <- data.frame(Subset = comb, ro = correlation, stringsAsFactors = FALSE)
    return(result[order(result[, 2], decreasing = TRUE), ])
}
