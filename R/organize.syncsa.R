#' Function for organize data for Package SYNCSA
#' 
#' Package \strong{SYNCSA} requires that the species and community sequence in
#' the data matrix must be the same for all matrices.
#' 
#' The Function organizes the data for the functions of the package
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
#' @return The matrices of community, traits, phylogenetic distance and
#' environmental variables. The objects returned belong to the matrix class.
#' @note The function organizes the matrices despite the absence of one of the
#' matrices, provided that the community data had been entered. Unspecified
#' matrices will appear as NULL.
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{matrix.t}}, \code{\link{matrix.x}},
#' \code{\link{matrix.p}}, \code{\link{syncsa}}
#' @keywords SYNCSA
#' @export
organize.syncsa <-
function (comm, traits, dist.spp, envir) 
	{
	if (missing(comm)=="TRUE"){
		stop("\n Community not fount\n")
		}
	if (is.null(colnames(comm))){
		stop("\n Erro in column names of comm\n")
		}
	if (is.null(rownames(comm))){
		stop("\n Erro in row names of comm\n")
		}	
    if (!missing(traits) == "TRUE") {
    	if (is.null(colnames(traits))){
		stop("\n Erro in column names of traits\n")
		}
		if (is.null(rownames(traits))){
		stop("\n Erro in row names of traits\n")
		}
		match.names<-match(colnames(comm), rownames(traits))
		if(sum(is.na(match.names))>0)
		{
		print("There are species from community data that are not on traits matrix:",quote=FALSE)
		print(setdiff(colnames(comm), rownames(traits)))
		stop("\n Species not found on traits matrix\n")
		}
		traits<-as.matrix(traits[match.names,])
		}
	if (!missing(dist.spp) == "TRUE") {
		if (is.null(colnames(dist.spp))){
		stop("\n Erro in column names of dist.spp\n")
		}
		if (is.null(rownames(dist.spp))){
		stop("\n Erro in row names of dist.spp\n")
		}
		match.names<-match(colnames(comm), colnames(dist.spp))
		if(sum(is.na(match.names))>0)
		{
		print("There are species from community data that are not on phylogenetic distance matrix:",quote=FALSE)
		print(setdiff(colnames(comm), colnames(dist.spp)))
		stop("\n Species not found on phylogenetic distance matrix\n")
		}
	dist.spp<-as.matrix(dist.spp[match.names,match.names])
	}
	if (!missing(envir) == "TRUE") {
		if (is.null(colnames(envir))){
		stop("\n Erro in column names of envir\n")
		}		
		if (is.null(rownames(envir))){
		stop("\n Erro in row names of envir\n")
		}
		match.names<-match(rownames(comm), rownames(envir))
		if(sum(is.na(match.names))>0)
		{
		print("The are community that are not on environmental matrix:",quote=FALSE)
		print(setdiff(rownames(comm), rownames(envir)))
		stop("\n Sampling units not found on environmental matrix\n")
		}
		envir<-as.matrix(envir[match.names,])   
    }
    if (missing(traits)=="TRUE"){
	traits<-NULL
		}
	if (missing(dist.spp)=="TRUE"){
	dist.spp<-NULL
		}
	if (missing(envir)=="TRUE"){
	envir<-NULL
		}		
    return(list(community=as.matrix(comm),traits=traits,dist.spp=dist.spp,environmental=envir))
}

