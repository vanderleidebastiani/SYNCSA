#' @title Function for organize data for Package SYNCSA
#'
#' @description Package \strong{SYNCSA} requires that the species and community sequence in
#' the dataframe or matrix must be the same for all dataframe/matrices.
#'
#' @details The function organizes the data for the functions of the package
#' SYNCSA, placing the matrices of community, traits, phylogenetic
#' distance, environmental varibles and strata vector in the same order.
#' The function use as reference the community data for organize all data frame or matrices
#' in the same order that the sampling units names and species names found in community
#' data set. For this all data sets entered must be correctly named, with rows and columns
#' named. The matrices phylodist, traits, envir can be larger than community
#' data (more species and/or more sampling units) as long as it has at least
#' all species and/or sampling units that are in community data. The function
#' organizes the data despite the absence of one of the dataframes or matrices,
#' provided that the community data had been entered. Unspecified data will
#' appear as NULL.
#'
#' When trait is a data.frame with differents types of variables correctly identified and
#' convert.traits is TRUE factor traits are expanded in dummy traits and ordinal variables
#' are converted in numeric according to ranks argument.
#'
#' The strata must be a named vector. The strata vector is a vector for restrict
#' permutation within species groups, insofar as the SYNCSA package the null
#' models are based in permutation of species rather than permutation of sample units.
#'
#' @encoding UTF-8
#' @param comm Community data, with species as columns and sampling units as
#' rows.
#' @param traits Matrix or dataframe of species described by traits, with traits as
#' columns and species as rows (Default traits = NULL).
#' @param phylodist Matrix containing phylogenetic distance between species.
#' Must be a complete matrix (not a half diagonal matrix).This matrix can be
#' larger than community data (more species) as long as it has at least all
#' species that are in community data (Default phylodist = NULL).
#' @param envir Environmental variables for each community, with variables as
#' columns and sampling units as rows (Default envir = NULL).
#' @param strata Strata named vector to specify restricting permutations within
#' species groups (Default strata = NULL).
#' @param check.comm Logical argument (TRUE or FALSE) to remove sampling units and
#' species with total sums equal or less than zero (Default check.comm = TRUE).
#' @param convert.traits Logical argument (TRUE or FALSE) to convert factor traits in
#' dummy traits and/or convert ordinal variables in numeric (see ranks argument)
#' (Default convert.traits = FALSE).
#' @param ranks Logical argument (TRUE or FALSE) to specify if ordinal variables are
#' convert to ranks. If ranks = TRUE all ordered variable are replaced by their ranks
#' and if ranks = FALSE all ordinal variables are simply treated as continuous variables
#' (Default ranks = TRUE).
#' @return The dataframes or matrices of community, traits, phylogenetic distance and
#' environmental variables, as well as the type of each varible in each dataframe or matrix
#' (where 'c' to continuous/numeric, 'o' to ordinal, 'b' to binary and 'f' to factor, nominal
#' are not allowed). The strata vector for permutations. If convert.traits is TRUE and some
#' traits is of class factor the function returm the traits with expand traits in dummy variables
#' and return a list with sugestion to group of traits that are analysed together. Warning messenger
#' are returned as a list of warning.
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{syncsa}}, \code{\link{var.dummy}}, \code{\link{var.type}}
#' @keywords SYNCSA
#' @examples
#' data(ADRS)
#' organize.syncsa(ADRS$community, ADRS$traits, ADRS$phylo, ADRS$envir)
#' @export
organize.syncsa <- function (comm, traits = NULL, phylodist = NULL, envir = NULL, strata = NULL,
                             check.comm = TRUE, convert.traits = FALSE, ranks = TRUE)
{
  res <- list(call = match.call())
  if (missing(comm)){
    stop("\n comm not fount\n")
  }
  if(class(comm) != "data.frame" & class(comm) != "matrix"){
    stop("comm must be a matrix or a data.frame")
  }
  if (is.null(colnames(comm))){
    stop("\n Column names of comm are null\n")
  }
  if (is.null(rownames(comm))){
    stop("\n Row names of comm are null\n")
  }
  commvartype <- var.type(comm)
  if(any(commvartype == "n") | any(commvartype == "f") | any(commvartype == "o")){
    stop("\n comm must contain only numeric or binary variables \n")
  }
  list.warning <- list()
  put.together <- NULL
  if(check.comm){
    col.rm <- colnames(comm)[!colSums(comm, na.rm = TRUE)>0]
    row.rm <- rownames(comm)[!rowSums(comm, na.rm = TRUE)>0]
    if(length(col.rm)>0){
      warning("Species removed from community data - Check list of warning in list.warning", call. = FALSE)
      list.warning$comm$spp <- data.frame(species.removed = col.rm)
    }
    if(length(row.rm)>0){
      warning("Communities removed from community data - Check list of warning in list.warning", call. = FALSE)
      list.warning$comm$comm <- data.frame(communities.removed=row.rm)
    }
    comm <- comm[, colSums(comm, na.rm = TRUE)>0, drop = FALSE]
    comm <- comm[rowSums(comm,na.rm=TRUE)>0, , drop = FALSE]
  }
  if(any(is.na(comm))){
    warning("Warning: NA in community data", call. = FALSE)
  }
  if (!is.null(traits)) {
    if(class(traits) != "data.frame" & class(traits) != "matrix"){
      stop("traits must be a matrix or a data.frame")
    }
    if (is.null(colnames(traits))){
      stop("\n Column names of traits are null\n")
    }
    if (is.null(rownames(traits))){
      stop("\n Row names of traits are null\n")
    }
    traitsvartype <- var.type(traits)
    if(any(traitsvartype == "n")){
      stop("\n trait must contain only numeric, binary, factor or ordinal variables \n")
    }
    match.names <- match(colnames(comm), rownames(traits))
    if(sum(is.na(match.names))>0){
      list.warning$traits$spp <- data.frame(species.not.on.traits = setdiff(colnames(comm), rownames(traits)))
      warning("ERROR - Check list of warning in list.warning$traits", call. = FALSE)
      return(list(list.warning = list.warning))
    }
    traits <- as.data.frame(traits[match.names, , drop = FALSE])
    if(convert.traits){
      if(any(traitsvartype == "f")){
        warning("Factor variables expanded in dummy variables", call. = FALSE)
      }
      traits.dummy.temp <- var.dummy(traits)
      traits <- traits.dummy.temp$data
      put.together <- traits.dummy.temp$together
      traitsvartype <- var.type(traits)
      if(any(traitsvartype == "o")){
        warning("Ordinal variables transformed in continual variables", call. = FALSE)
      }
      traits <- data.matrix(traits)
      for (i in 1:length(traitsvartype)){
        if (traitsvartype[i] == "o"){
          if (ranks){
            traits[,i] <- rank(traits[,i], na.last = "keep")
          } else {
            traits[,i] <- as.numeric(traits[, i])
          }
        }
      }
      traitsvartype <- var.type(traits)
    }
    if(any(is.na(traits))){
      warning("Warning: NA in traits matrix", call. = FALSE)
    }
  }
  if (!is.null(phylodist)) {
    if(class(phylodist) != "data.frame" & class(phylodist) != "matrix"){
      stop("phylodist must be a matrix or a data.frame")
    }
    if (is.null(colnames(phylodist))){
      stop("\n Column names of phylodist are null\n")
    }
    if (is.null(rownames(phylodist))){
      stop("\n Row names of phylodist are null\n")
    }
    phylodistvartype <- var.type(phylodist)
    if(any(phylodistvartype == "n") | any(phylodistvartype == "f") | any(phylodistvartype == "o")){
      stop("\n phylodist must contain only numeric or binary variables \n")
    }
    match.names <- match(colnames(comm), colnames(phylodist))
    if(sum(is.na(match.names))>0){
      list.warning$phylodist$spp <- data.frame(species.not.on.phylodist = setdiff(colnames(comm), colnames(phylodist)))
      warning("ERROR - Check list of warning in list.warning$phylodist", call. = FALSE)
      return(list.warning)
    }
    phylodist <- phylodist[match.names, match.names, drop = FALSE]
    if(any(is.na(phylodist))){
      warning("Warning: NA in phylogenetic distance matrix",call.=FALSE)
    }
  }
  if (!is.null(envir)) {
    if(class(envir) != "data.frame" & class(envir) != "matrix"){
      stop("envir must be a matrix or a data.frame")
    }
    if (is.null(colnames(envir))){
      stop("\n Column names of envir are null\n")
    }
    if (is.null(rownames(envir))){
      stop("\n Row names of envir are null\n")
    }
    envirvartype <- var.type(envir)
    if(any(envirvartype == "n") | any(envirvartype == "f") | any(envirvartype == "o")){
      stop("\n envir must contain only numeric or binary variables \n")
    }
    match.names <- match(rownames(comm), rownames(envir))
    if(sum(is.na(match.names))>0){
      list.warning$envir$comm <- data.frame(comm.not.on.envir = setdiff(rownames(comm), rownames(envir)))
      warning("ERROR - Check list of warning in list.warning$envir",call.=FALSE)
      return(list.warning)
    }
    envir <- envir[match.names,,drop=FALSE]
    if(any(is.na(envir))){
      warning("Warning: NA in environmental data", call. = FALSE)
    }
  }
  if (!is.null(strata)) {
    if (is.null(names(strata))){
      stop("\n Names of strata factor are null\n")
    }
    match.names <- match(colnames(comm), names(strata))
    if(sum(is.na(match.names))>0){
      list.warning$strata <- data.frame(species.not.on.strata = setdiff(colnames(comm), names(strata)))
      warning("ERROR - Check list of warning in list.warning$strata",call.=FALSE)
      return(list.warning)
    }
    strata <- strata[match.names]
  }
  if (is.null(traits)){
    traits <- NULL
    traitsvartype <- NULL
  }
  if (is.null(phylodist)){
    phylodist <- NULL
    phylodistvartype <- NULL
  }
  if (is.null(envir)){
    envir <- NULL
    envirvartype <- NULL
  }
  if (is.null(strata)){
    strata <- NULL
  }
  if(length(list.warning)>0){
    res$list.warning <- list.warning
  }
  res$community <- comm
  res$traits <- traits
  res$phylodist <- phylodist
  res$environmental <- envir
  res$community.var.type <- commvartype
  res$traits.var.type <- traitsvartype
  res$phylodist.var.type <- phylodistvartype
  res$environmental.var.type <- envirvartype
  res$strata <- strata
  res$put.together <- put.together
  return(res)
}
