#' Matrix P
#' 
#' Function to obtain a matrix containing phylogeny-weighted species
#' composition. For more details, see \code{\link{syncsa}}.
#' 
#' 
#' @param comm Community data, with species as columns and sampling units as
#' rows. This matrix can contain either presence/absence or abundance data.
#' @param dist.spp Matrix containing phylogenetic distance between species.
#' Must be a complete matrix (not a half diagonal matrix).
#' @param notification Logical argument (TRUE or FALSE) to specify if
#' notifications of missing observations are shown (Default notification =
#' TRUE).
#' @return \item{matrix.w}{Standardized community matrix, where rows are
#' communities and columns species. Row totals (communities) = 1.}
#' \item{matrix.q}{Standardized matrix containing the degree of belonging of
#' species in relation to each other. Row totals (species) = 1.}
#' \item{matrix.P}{Phylogeny-weighted species composition matrix. Row totals
#' (communities) = 1.}
#' @note \strong{IMPORTANT}: The sequence species show up in community data
#' matrix MUST be the same as they show up in phylogenetic distance matrix. See
#' \code{\link{organize.syncsa}}.
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{matrix.t}}, \code{\link{matrix.x}},
#' \code{\link{syncsa}}, \code{\link{organize.syncsa}}
#' @references Pillar, V.D.; Duarte, L.d.S. (2010). A framework for
#' metacommunity analysis of phylogenetic structure. Ecology Letters, 13,
#' 587-596.
#' @keywords SYNCSA
#' @examples
#' 
#' data(flona)
#' matrix.p(flona$community,flona$phylo)
#' 
#' @export
matrix.p<-function (comm, dist.spp, notification = TRUE)
{
	comm<-as.matrix(comm)
	dist.spp<-as.matrix(dist.spp)
    matrix.w <- sweep(comm, 1, rowSums(comm, na.rm=TRUE), "/")
	w.NA <- apply(matrix.w, 2, is.na)
    matrix.w[w.NA] <-0
    if(notification==TRUE){
    	if(length(which(unique(as.vector(w.NA))==TRUE))>0)
    	{
			warning("Warning: NA in community data",call.=FALSE)		
    	}  	 
    }
    similar.phy <- 1 - (dist.spp/max(dist.spp))
    matrix.phy <- 1/colSums(similar.phy)
    matrix.q <- sweep(similar.phy, 1, matrix.phy, "*")
    matrix.P <- matrix.w %*% matrix.q
    return(list(matrix.w = matrix.w, matrix.q = matrix.q, matrix.P = matrix.P))
}
