#' @title Matrix X
#'
#' @description Function to obtain a matrix containing trait-weighted species composition.
#' For more details, see \code{\link{syncsa}}.
#'
#' @encoding UTF-8
#' @importFrom vegan vegdist
#' @importFrom FD gowdis
#' @param comm Community data, with species as columns and sampling units as
#' rows. This matrix can contain either presence/absence or abundance data.
#' @param traits Matrix or data frame of species described by traits, with traits as
#' columns and species as rows.
#' @param scale Logical argument (TRUE or FALSE) to specify if the traits are
#' measured on different scales (Default scale = TRUE). When scale = TRUE traits
#' are measured on different scales the matrix of traits is subjected to
#' standardization within each trait, and Gower Index is used to calculate the
#' degree of belonging to the species. When scale = FALSE traits are measured on
#' the same scale the matrix of traits is not subjected to standardization,
#' and Euclidean distance is calculated to determine the degree of belonging to
#' the species.
#' @param ranks Logical argument (TRUE or FALSE) to specify if ordinal variables are
#' convert to ranks (Default ranks = TRUE).
#' @param ord Method to be used for ordinal variables, see \code{\link{gowdis}}, if any
#' method is provided the rank parameter is not apply.
#' @param transformation Method to community data transformation, "none", "standardized"
#' or "weights" (Default transformation = "standardized").
#' @param spp.weights Vector with 0 or 1 to specify individual species weights (Default
#' spp.weights = NULL).
#' @param notification Logical argument (TRUE or FALSE) to specify if
#' notifications of missing observations are shown (Default notification =
#' TRUE).
#' @param ... Parameters for \code{\link{gowdis}} function.
#' @return \item{matriz.w}{Standardized community matrix, where rows are
#' communities and columns species. If default transformation, row totals (communities) = 1.}
#' \item{matriz.u}{Standardized matrix containing the degree of belonging of
#' each species in relation to each other species. Row totals (species) = 1.}
#' \item{matriz.X}{Trait-weighted species composition matrix. If default transformation,
#' row totals (communities) = 1.}
#' @note \strong{IMPORTANT}: The sequence species show up in community data
#' matrix MUST be the same as they show up in traits matrix or in the spp.weights vector. See
#' \code{\link{organize.syncsa}}.
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{syncsa}}, \code{\link{organize.syncsa}}, \code{\link{belonging}},
#' \code{\link{matrix.t}}, \code{\link{matrix.p}}, \code{\link{gowdis}}
#' @references Pillar, V.D.; Duarte, L.d.S. (2010). A framework for
#' metacommunity analysis of phylogenetic structure. Ecology Letters, 13,
#' 587-596.
#'
#' Pillar, V.D., Duarte, L.d.S., Sosinski, E.E. & Joner, F. (2009).
#' Discriminating trait-convergence and trait-divergence assembly patterns in
#' ecological community gradients. Journal of Vegetation Science, 20, 334:348.
#' @keywords SYNCSA
#' @examples
#' data(ADRS)
#' matrix.x(ADRS$community, ADRS$traits)
#' @export
matrix.x <- function (comm, traits, scale = TRUE, ranks = TRUE, ord,
                      transformation = "standardized", spp.weights = NULL, notification = TRUE, ...)
{
  vartype <- var.type(traits)
  if(any(vartype == "n")){
    stop("\n trait must contain only numeric, binary or ordinal variables \n")
  }
  if(missing(ord)){
    for(i in 1:length(vartype)){
      if(ranks & vartype[i] == "o"){
        traits[, i] <- rank(traits[, i], na.last = "keep")
      }
      traits[, i] <- as.numeric(traits[, i])
    }
    traits <- as.matrix(traits)
  }
  TRANS <- c("none", "standardized", "weights", "beals")
  trans <- pmatch(transformation, TRANS)
  if (length(trans) > 1) {
    stop("\n Only one argument is accepted in transformation \n")
  }
  if (is.na(trans) | trans == 4) {
    stop("\n Invalid transformation \n")
  }
  matrix.w <- matrix.w.transformation(comm, transformation = transformation,
                                      spp.weights = spp.weights, notification = notification)
  x.NA <- apply(traits, 2, is.na)
  if(notification){
    if(any(x.NA)){
      warning("Warning: NA in traits matrix", call. = FALSE)
    }
  }
  if (scale) {
    dist.traits <- as.matrix(FD::gowdis(traits, ...))
    ## In case a species has NO trait info, replace the diag with NA
    diag(dist.traits)[which(rowSums(is.na(dist.traits))==(ncol(dist.traits)-1))] <- NA
    similar.traits <- 1 - dist.traits
    matrix.traits <- 1/colSums(similar.traits, na.rm = TRUE)
    matrix.u <- sweep(similar.traits, 1, matrix.traits, "*")
  }
  else{
    dist.traits <- as.matrix(vegan::vegdist(traits, method = "euclidean", diag = TRUE, upper = TRUE, na.rm = TRUE))
    similar.traits <- 1 - (dist.traits/max(dist.traits, na.rm = TRUE))
    matrix.traits <- 1/colSums(similar.traits, na.rm = TRUE)
    matrix.u <- sweep(similar.traits, 1, matrix.traits, "*")
  }
  u.NA <- apply(matrix.u, 2, is.na)
  if(notification){
    if(any(u.NA)){
      warning("Warning: NA in matrix U", call. = FALSE)
    }
  }
  matrix.u[u.NA] <- 0
  matrix.X <- matrix.w %*% matrix.u
  # matrix.X <- matmult.syncsa(matrix.w, matrix.u)
  if(any(u.NA)){
    matrix.X <- sweep(matrix.X, 1, rowSums(matrix.X), "/")
  }
  return(list(matrix.w = matrix.w, matrix.u = matrix.u, matrix.X = matrix.X))
}
