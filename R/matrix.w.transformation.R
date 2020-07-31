#' @title Community data transformation
#'
#' @description Function to transformation community data and replace missing data. See details.
#'
#' @details The function applies standardization/transformation for community data. The options
#' are: "none" no transformation is applied; "standardized" the community data is standardized
#' to row totals will be 1; "weights" community data is first "standardized" and when, individual
#' species weights are multiplied to each species entries; and "beals" option applies Beals smoothing
#' using the function \code{\link{beals}}. The arguments "reference", "type" and "include" are the same
#' used in \code{\link{beals}} function. Missing data are replaced by 0 after the standardization/transformation.
#'
#' @encoding UTF-8
#' @importFrom vegan beals
#' @param comm Community data, with species as columns and sampling units as
#' rows. This matrix can contain either presence/absence or abundance data.
#' @param transformation Method to community data transformation, "none", "standardized", "weights"
#' or "beals" (Default transformation = "standardized").
#' @param spp.weights Vector with 0 or 1 to specify individual species weights (Default spp.weights = NULL).
#' @param reference Data to compute joint occurrences in \code{\link{beals}} function. If NULL, comm is
#' used as reference to compute the joint occurrences (Default reference = NULL).
#' @param type Method to specifies how abundance values are used in \code{\link{beals}} function, type = 0
#' presence/absence, type = 1 abundances are used to compute conditioned probabilities, type = 2 abundances
#' are used to compute weighted averages of conditioned probabilities or type = 3 abundances are used to
#' compute both conditioned probabilities and weighted averages (Default type = 0).
#' @param include Logical argument (TRUE or FALSE) to specify if target species are included when computing
#' the mean of the conditioned probabilities in \code{\link{beals}} (Default include = TRUE).
#' @param notification Logical argument (TRUE or FALSE) to specify if notifications for missing observations
#' are to be shown (Default notification = TRUE).
#' @return Transformed community matrix, where rows are communities and columns species.
#' @note \strong{IMPORTANT}: Species sequence in the community data
#' matrix MUST be the same as the one in the reference matrix or in spp.weights vector. See
#' \code{\link{organize.syncsa}}.
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{syncsa}}, \code{\link{organize.syncsa}}, \code{\link{matrix.p}},
#' \code{\link{matrix.t}}, \code{\link{matrix.x}}
#' @references Pillar, V.D.; Duarte, L.d.S. (2010). A framework for
#' metacommunity analysis of phylogenetic structure. Ecology Letters, 13,
#' 587-596.
#' De Cáceres, M.; Legendre, P. 2008. Beals smoothing revisited. Oecologia 156: 657–669.
#' @keywords Auxiliary
#' @examples
#' data(ADRS)
#' matrix.w.transformation(ADRS$community)
#' @export
matrix.w.transformation <- function(comm, transformation = "standardized", spp.weights = NULL,
                                     reference = NULL, type = 0, include = TRUE, notification = TRUE)
{
  matrix.w <- as.matrix(comm)
  TRANS <- c("none", "standardized", "weights", "beals")
  trans <- pmatch(transformation, TRANS)
  if (length(trans) > 1) {
    stop("\n Only one argument is accepted in transformation \n")
  }
  if (is.na(trans)) {
    stop("\n Invalid transformation \n")
  }
  w.NA <- apply(matrix.w, 2, is.na)
  if(notification){
    if(any(w.NA)){
      warning("Warning: NA in community data", call.= FALSE)
    }
  }
  if(trans == 2){
    matrix.w <- sweep(matrix.w, 1, rowSums(matrix.w, na.rm = TRUE), "/")
  }
  if(trans == 3){
    is.bin.weights <- all(spp.weights %in% c(0, 1))
    if(!is.bin.weights | is.null(spp.weights)){
      stop("\n spp.weights must be 0 or 1\n")
    }
    if((ncol(matrix.w) != length(spp.weights))){
      stop("\n spp.weights must be the same length of number of species \n")
    }
    matrix.w <- sweep(matrix.w, 1, rowSums(matrix.w, na.rm = TRUE), "/")
    spp.weights <- sapply(spp.weights, function(x, n) rep(x, n), n = nrow(matrix.w))
    matrix.w <- matrix.w*spp.weights
  }
  if(trans == 4){
    if(is.null(reference)){
      reference <- matrix.w
    }
    if((ncol(matrix.w) != ncol(reference)) | (nrow(matrix.w) != nrow(reference))){
      stop("\n comm and reference data must be the same dimensions \n")
    }
    matrix.w <- vegan::beals(matrix.w, reference = reference, type = type, include = include)
  }
  w.NA <- apply(matrix.w, 2, is.na)
  matrix.w[w.NA] <- 0
  return(matrix.w)
}
