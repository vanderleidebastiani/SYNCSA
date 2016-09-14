#' Matrix T
#' 
#' Function to obtain a matrix containing trait averages at community level.
#' For more details, see \code{\link{syncsa}}.
#' 
#' 
#' @encoding UTF-8
#' @param comm Community data, with species as columns and sampling units as
#' rows. This matrix can contain either presence/absence or abundance data.
#' @param traits Matrix data of species described by traits, with traits as
#' columns and species as rows.
#' @param scale Logical argument (TRUE or FALSE) to specify if the traits are
#' measured on different scales (Default Scale = TRUE). Scale = TRUE if traits
#' are measured on different scales, the matrix T is subjected to
#' standardization within each trait. Scale = FALSE is traits are measured on
#' the same scale, the matrix T is not subjected to standardization.
#' @param notification Logical argument (TRUE or FALSE) to specify if
#' notifications of missing observations are shown (Default notification =
#' TRUE).
#' @param ord Method to be used for ordinal variables, see \code{\link{gowdis}}
#' (Default ord = "metric").
#' @return \item{matriz.w}{Standardized community matrix, where rows are
#' communities and columns species. Row totals (communities) = 1.}
#' \item{matriz.b}{Matrix of traits, exactly the same data input.}
#' \item{matriz.T}{Matrix containing trait averages at community level. If
#' Scale = TRUE the matrix T is standardized within the traits.}
#' @note \strong{IMPORTANT}: The sequence species show up in community data
#' matrix MUST be the same as they show up in traits matrix. See
#' \code{\link{organize.syncsa}}.
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{matrix.p}}, \code{\link{matrix.x}},
#' \code{\link{syncsa}}, \code{\link{organize.syncsa}}
#' @references Pillar, V.D.; Duarte, L.d.S. (2010). A framework for
#' metacommunity analysis of phylogenetic structure. Ecology Letters, 13,
#' 587-596.
#' 
#' Pillar, V.D., Duarte, L.d.S., Sosinski, E.E. & Joner, F. (2009).
#' Discriminating trait-convergence and trait-divergence assembly patterns in
#' ecological community gradients. Journal of Vegetation Science, 20, 334:348.
#' @keywords SYNCSA
#' @examples
#' 
#' data(flona)
#' matrix.t(flona$community,flona$traits,scale=TRUE)
#' 
#' @export
matrix.t<-function (comm, traits, scale = TRUE, notification = TRUE, ord = "metric") 
{
	comm<-as.matrix(comm)
	vartype<-vartype(traits)
	if(any(vartype=="n")){
		stop("\n trait must contain only numeric, binary or ordinal variables \n")
	}
	for(i in 1:length(vartype)){
		if(ord != "classic" & vartype[i]=="o"){
			traits[, i] <- rank(traits[, i], na.last = "keep")
		}
		traits[, i] <- as.numeric(traits[, i])
	}
	traits<-as.matrix(traits)
	matrix.w <- sweep(comm, 1, rowSums(comm, na.rm=TRUE), "/")
	w.NA <- apply(matrix.w, 2, is.na)
    matrix.w[w.NA] <-0
    if(notification){
    	if(any(w.NA)){
			warning("Warning: NA in community data",call.=FALSE)		
    	}  	 
    }
    b.NA <- apply(traits, 2, is.na)
    traits[b.NA] <-0
    if(notification){
    	if(any(b.NA)){
			warning("Warning: NA in traits matrix",call.=FALSE)	
    	}
    }
	matrix.b <- traits
    matrix.T <- matrix.w %*% traits
    if (scale == "TRUE") {
        matrix.traits <- apply(matrix.T^2, 2, sum)
        matrix.T <- sweep(matrix.T, 2, sqrt(matrix.traits), "/")
    }
    return(list(matrix.w = matrix.w, matrix.b = matrix.b, matrix.T = matrix.T))
}