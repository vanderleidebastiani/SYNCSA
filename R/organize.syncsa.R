#' Function for organize data for Package SYNCSA
#' 
#' Package \strong{SYNCSA} requires that the species and community sequence in
#' the dataframe or matrix must be the same for all dataframe/matrices.
#' 
#' The function organizes the data for the functions of the package
#' \strong{SYNCSA}, placing the matrices of community, traits, phylogenetic
#' distance and environmental varibles in the same order.
#' 
#' 
#' @encoding UTF-8
#' @param comm Community data, with species as columns and sampling units as
#' rows.
#' @param traits Matrix data of species described by traits, with traits as
#' columns and species as rows.
#' @param dist.spp Matrix containing phylogenetic distance between species.
#' Must be a complete matrix (not a half diagonal matrix).This matrix can be
#' larger than community data (more species) as long as it has at least all
#' species that are in community data.
#' @param envir Environmental variables for each community, with variables as
#' columns and sampling units as rows.
#' @param strata Strata nomed vector to specify restricting permutations within
#' species groups.
#' @param check.comm Logical argument (TRUE or FALSE) to remove sampling units and 
#' species with total sums equal or less than zero (Default check.comm = "TRUE").
#' @return The dataframes or matrices of community, traits, phylogenetic distance and
#' environmental variables. The strata vector for permutations.
#' @note The function organizes the data despite the absence of one of the
#' dataframes or matrices, provided that the community data had been entered. Unspecified
#' data will appear as NULL.
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{matrix.t}}, \code{\link{matrix.x}},
#' \code{\link{matrix.p}}, \code{\link{syncsa}}
#' @keywords SYNCSA
#' @export
organize.syncsa <- function (comm, traits, dist.spp, envir, strata, check.comm = TRUE){
	if (missing(comm)=="TRUE"){
		stop("\n Community not fount\n")
	}
	if (is.null(colnames(comm))){
		stop("\n Column names of comm are null\n")
	}
	if (is.null(rownames(comm))){
		stop("\n Row names of comm are null\n")
	}
	if(check.comm){
		col.rm<-colnames(comm)[!colSums(comm,na.rm=TRUE)>0]
		row.rm<-rownames(comm)[!rowSums(comm,na.rm=TRUE)>0]
		if(length(col.rm)>0){
			print("Species removed from community data:",quote=FALSE)
			print(col.rm)
		}
		if(length(row.rm)>0){
			print("Communities removed from community data:",quote=FALSE)
			print(row.rm)	
		}		
		comm<-comm[,colSums(comm,na.rm=TRUE)>0,drop=FALSE]
		comm<-comm[rowSums(comm,na.rm=TRUE)>0,,drop=FALSE]
	}
	if(any(is.na(comm))){
		warning("Warning: NA in community data",call.=FALSE)
	}
	commvartype<-vartype(comm)
	if(any(commvartype=="N")){
		stop("\n comm must contain only numeric, binary or ordinal variables \n")
	}
    if (!missing(traits) == "TRUE") {
    	if (is.null(colnames(traits))){
		stop("\n Column names of traits are null\n")
		}
		if (is.null(rownames(traits))){
		stop("\n Row names of traits are null\n")
		}
    	if(any(is.na(traits))){
			warning("Warning: NA in traits matrix",call.=FALSE)
    	}
		match.names<-match(colnames(comm), rownames(traits))
		if(sum(is.na(match.names))>0){
			print("There are species from community data that are not on traits matrix:",quote=FALSE)
			print(setdiff(colnames(comm), rownames(traits)))
			stop("\n Species not found on traits matrix\n")
		}
		traits<-as.data.frame(traits[match.names,,drop=FALSE])
		traitsvartype<-vartype(traits)
		if(any(traitsvartype=="N")){
			stop("\n trait must contain only numeric, binary or ordinal variables \n")
		}
	}
	if (!missing(dist.spp) == "TRUE") {
		if (is.null(colnames(dist.spp))){
			stop("\n Column names of dist.spp are null\n")
		}
		if (is.null(rownames(dist.spp))){
			stop("\n Row names of dist.spp are null\n")
		}
		if(any(is.na(dist.spp))){
			warning("Warning: NA in phylogenetic distance matrix",call.=FALSE)
		}
		match.names<-match(colnames(comm), colnames(dist.spp))
		if(sum(is.na(match.names))>0){
			print("There are species from community data that are not on phylogenetic distance matrix:",quote=FALSE)
			print(setdiff(colnames(comm), colnames(dist.spp)))
			stop("\n Species not found on phylogenetic distance matrix\n")
		}
		dist.spp<-dist.spp[match.names,match.names,drop=FALSE]
		dist.sppvartype<-vartype(dist.spp)
		if(any(dist.sppvartype=="N")){
			stop("\n dist.spp must contain only numeric, binary or ordinal variables \n")
		}
	}
	if (!missing(envir) == "TRUE") {
		if (is.null(colnames(envir))){
			stop("\n Column names of envir are null\n")
		}		
		if (is.null(rownames(envir))){
			stop("\n Row names of envir are null\n")
		}
		if(any(is.na(envir))){
			warning("Warning: NA in environmental data",call.=FALSE)
		}
		match.names<-match(rownames(comm), rownames(envir))
		if(sum(is.na(match.names))>0){
			print("The are community that are not on environmental matrix:",quote=FALSE)
			print(setdiff(rownames(comm), rownames(envir)))
			stop("\n Sampling units not found on environmental matrix\n")
		}
		envir<-envir[match.names,,drop=FALSE]
		envirvartype<-vartype(envir)
		if(any(envirvartype=="N")){
			stop("\n envir must contain only numeric, binary or ordinal variables \n")
		}
    }
    if (!missing(strata)) {
		if (is.null(names(strata))){
			stop("\n Names of strata factor are null\n")
		}		
		match.names<-match(colnames(comm), names(strata))
		if(sum(is.na(match.names))>0){
			print("There are species from community data that are not on strata vector :",quote=FALSE)
			print(setdiff(colnames(comm), names(strata)))
			stop("\n Species not found on strata\n")
		}
		strata<-strata[match.names]
    }
    if (missing(traits)){
		traits<-NULL
		traitsvartype=NULL
	}
	if (missing(dist.spp)){
		dist.spp<-NULL
		dist.sppvartype =NULL
	}
	if (missing(envir)){
		envir<-NULL
		envirvartype=NULL
	}
	if (missing(strata)){
		strata<-NULL
	}
	res<-list(community=comm,traits=traits,dist.spp=dist.spp,environmental=envir,community.var.type= commvartype,traits.var.type= traitsvartype, dist.spp.var.type = dist.sppvartype, environmental.var.type = envirvartype, strata = strata)		
    return(res)
}