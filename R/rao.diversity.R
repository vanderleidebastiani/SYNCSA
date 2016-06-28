#' Rao's quadratic entropy
#'
#' Calculates Rao's quadratic entropy, functional and phylogenetic redundancy.
#'
#' Rao's quadratic entropy is a measure of diversity of ecological communities
#' defined by Rao (1982) and is based on the proportion of the abundance of
#' species present in a community and some measure of dissimilarity among them.
#' The dissimilarity range from 0 to 1 and is based on a set of specified
#' functional traits or in the phylogenetic dissimilarity.
#'
#' For the trait data , the function calculates the square root of the
#' one-complement of Gower`s similarity index, in order to have a dissimilarity
#' matrix with Euclidean metric properties. Gower`s index ranges from 0 to 1
#' and can handle traits measured indifferent scales. When the species are
#' completely different in terms of their traits, Rao quadratic entropy is
#' equivalent to the Gini-Simpson index.
#'
#' Functional redundancy is defined purely as the difference between species
#' diversity and Rao`s quadratic entropy based on their functional
#' dissimilarity (de Bello et al. 2007). The same definition is used for
#' phylogenetic redundancy.
#'
#' @encoding UTF-8
#' @importFrom FD gowdis
#' @param comm Community data, with species as columns and sampling units as
#' rows. This matrix can contain either presence/absence or abundance data.
#' @param traits Matrix data of species described by traits, with traits as
#' columns and species as rows (optional).
#' @param phylodist Matrix containing phylogenetic distance between species
#' (optional).
#' @param checkdata Logical argument (TRUE or FALSE) to check if species
#' sequence in the community data follows the same order as the one in the
#' trait and in the phylodist matrices (Default checkdata = TRUE).
#' @param ord Method to be used for ordinal variables, see \code{\link{gowdis}}
#' (Default ord = "metric").
#' @param put.together List to specify group traits that are added or removed
#' together (Default put.together = NULL). This argument must be a list, see
#' examples.
#' @param ... Parameters for \code{\link{gowdis}} function.
#' @return \item{Simpson}{Gini-Simpson index within each community (equivalent
#' to Rao quadratic entropy with null, crisp, similarities).} \item{FunRao}{Rao
#' quadratic entropy within each community, considering trait distance.}
#' \item{FunRedundancy}{Functional redundancy in each community.}
#' \item{PhyRao}{Rao quadratic entropy within each community, considering
#' phylogenetic distance.} \item{PhyRedundancy}{Phylogenetic redundancy in each
#' community.}
#' @note \strong{IMPORTANT}: The sequence species show up in community data
#' matrix MUST be the same as they show up in traits and phylodist matrices.
#' See \code{\link{organize.syncsa}}.
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{organize.syncsa}}, \code{\link{gowdis}},
#' \code{\link{syncsa}}
#' @references de Bello, F.; Leps, J.; Lavorel, S. & Moretti, M. (2007).
#' Importance of species abundance for assessment of trait composition: an
#' example based on pollinator communities. Community Ecology, 8, 163:170.
#'
#' Pillar, V.D.; Blanco, C.C.; Muler, S.C.; Sosinski, E.E.; Joner, F. & Duarte,
#' L.d.S. (2013). Functional redundancy and stability in plant communities.
#' Journal of Vegetation Science, 24, 963:974.
#'
#' Rao, C.R. (1982). Diversity and dissimilarity coefficients: a unified
#' approach. Theoretical Population Biology, 21, 24:43.
#' @keywords SYNCSA
#' @examples
#'
#' data(flona)
#' rao.diversity(flona$community)
#' rao.diversity(flona$community,traits=flona$traits)
#'
#' @export
rao.diversity<-function(comm, traits = NULL, phylodist = NULL, checkdata = TRUE, ord = "metric", put.together = NULL, ...){
	comm <- as.matrix(comm)
    N <- dim(comm)[1]
    S <- dim(comm)[2]
    tij2 <- 1 - diag(x = rep(1, S))
    if (!is.null(traits)) {
		if (checkdata) {
        if (is.null(rownames(traits))) {
            stop("\n Erro in row names of traits\n")
        }
        if (is.null(colnames(comm))) {
            stop("\n Erro in row names of comm\n")
        }
		match.names <- match(colnames(comm), rownames(traits))
        if (sum(is.na(match.names)) > 0) {
            stop("\n There are species from community data that are not on traits matrix\n")
        }
    		traits<-as.data.frame(traits[match.names,])
    	}
       	m <- dim(traits)[2]
		weights<-rep(1,m)
		make.names<-is.null(colnames(traits))
		colnames(traits) <- colnames(traits, do.NULL = FALSE, prefix = "T")
		names(weights)<-colnames(traits)
		if(!is.null(put.together)){
			if(class(put.together)!="list"){
				stop("\n The put.together must be a object of class list\n")
			}
			if(make.names){
				for(k in 1:length(put.together)){
					put.together[[k]]<-paste("T", put.together[[k]],sep="")
				}
			}
			if(max(table(unlist(put.together)))>1){
				stop("\n The same trait appears more than once in put.together\n")
			}
			if(length(setdiff(unlist(put.together),colnames(traits)))>0){
				stop("\n Check traits names in put.together\n")
			}
			for(k in 1:length(put.together)){
				weights[put.together[[k]]]<-1/length(put.together[[k]])
			}
		}    	
    	D1<-as.matrix(FD::gowdis(x=traits, asym.bin = NULL, ord = ord, w = weights))
    	#S1<-1-D1
    	#tij<-sqrt(1-S1)
    	tij<-sqrt(D1)
    }
    if (!is.null(phylodist)) {
		if (checkdata) {
        if (is.null(rownames(phylodist))) {
            stop("\n Erro in row names of phylodist\n")
        }
		if (is.null(colnames(phylodist))) {
            stop("\n Erro in column names of phylodist\n")
        }
        if (is.null(colnames(comm))) {
            stop("\n Erro in row names of comm\n")
        }
        match.names <- match(colnames(comm), colnames(phylodist))
        if (sum(is.na(match.names)) > 0) {
            stop("\n There are species from community data that are not on phylogenetic distance matrix\n")
        }
        phylodist <- as.matrix(phylodist[match.names, match.names])
    	}
    	D1<-as.matrix(phylodist)
    	tij3<-D1/max(D1)
    }
	comm <- sweep(comm, 1, rowSums(comm, na.rm = TRUE), "/")
	inter<-comm%*%tij2
	SD<-rowSums(sweep(comm,1,inter,"*",check.margin=F))
	if (!is.null(traits)){
		inter<-comm%*%tij
		RD<-rowSums(sweep(comm,1,inter,"*",check.margin=F))
	}
	if (!is.null(phylodist)){
		inter<-comm%*%tij3
		FRD<-rowSums(sweep(comm,1,inter,"*",check.margin=F))
	}
    Res<-list(Simpson=SD)
    if (!is.null(traits)){
		Res<-list(Simpson=SD,FunRao=RD,FunRedundancy=SD-RD)
    }
    if (!is.null(phylodist)){
		Res<-list(Simpson=SD,PhyRao=FRD,PhyRedundancy=SD-FRD)
    }
    if (!is.null(phylodist) & !is.null(traits)){
		Res<-list(Simpson=SD,FunRao=RD,FunRedundancy=SD-RD,PhyRao=FRD,PhyRedundancy=SD-FRD)
    }
return(Res)
}