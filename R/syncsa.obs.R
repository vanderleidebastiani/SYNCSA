#' @rdname syncsa
#' @encoding UTF-8
#' @export
syncsa.obs <- function (comm, traits = NULL, phylodist = NULL, envir = NULL,
                        scale = TRUE, scale.envir = TRUE, ranks = TRUE,
                        asym.bin = NULL, ord = "metric", put.together = NULL,
                        na.rm = FALSE, transformation = "standardized",
                        spp.weights = NULL, notification = TRUE)
  # syncsa.obs <- function (comm, traits = NULL, phylodist = NULL, envir = NULL,
  #                         scale = TRUE, scale.envir = TRUE, ranks = TRUE, ord,
  #                         put.together = NULL, na.rm = FALSE, transformation = "standardized",
  #                         spp.weights = NULL, notification = TRUE)
{
  res <- list(W = NULL, Y = NULL,
              B = NULL, T = NULL,
              U = NULL, X = NULL,
              Q = NULL, P = NULL, traits.weights = NULL)
  if (is.null(traits) & is.null(phylodist)){
    res$W <- matrix.w.transformation(comm, transformation = transformation,
                                     spp.weights = spp.weights, notification = FALSE)
  }
  if (!is.null(traits)) {
    m <- ncol(traits)
    traits.weights <- rep(1, m)
    make.names <- is.null(colnames(traits))
    colnames(traits) <- colnames(traits, do.NULL = FALSE, prefix = "T")
    names(traits.weights) <- colnames(traits)
    if(!is.null(put.together)){
      if(!inherits(put.together, "list")){
        stop("\n put.together must be a object of class list\n")
      }
      if(make.names){
        for(k in 1:length(put.together)){
          put.together[[k]] <- paste("T", put.together[[k]], sep = "")
        }
      }
      if(max(table(unlist(put.together)))>1){
        stop("\n The same trait appears more than once in put.together\n")
      }
      if(length(setdiff(unlist(put.together),colnames(traits)))>0){
        stop("\n Check traits names in put.together\n")
      }
      for(k in 1:length(put.together)){
        traits.weights[put.together[[k]]] <- 1/length(put.together[[k]])
      }
    }
    matrixT <- matrix.t(comm, traits, scale = scale, ranks = ranks,
                        transformation = transformation, spp.weights = spp.weights, notification = FALSE)
    check.U <- function(traits, scale, ranks, ...){
      # check.U <- function(traits, scale, ranks, ord, ...){
      vartype <- var.type(traits)
      # if(missing(ord)){
      #   for(i in 1:length(vartype)){
      #     if(ranks & vartype[i] == "o"){
      #       traits[, i] <- rank(traits[, i], na.last = "keep")
      #     }
      #     traits[, i] <- as.numeric(traits[, i])
      #   }
      #   traits <- as.matrix(traits)
      # }
      for(i in 1:length(vartype)){
        if(ranks & vartype[i] == "o"){
          traits[, i] <- rank(traits[, i], na.last = "keep")
        }
      }
      if (scale) {
        # dist.traits <- FD::gowdis(traits, asym.bin = NULL, ...)
        dist.traits <- FD::gowdis(traits, ...)
      }
      else{
        dist.traits <- as.matrix(vegan::vegdist(traits, method = "euclidean", diag = TRUE, upper = TRUE, na.rm = TRUE))
      }
      res <- any(is.na(dist.traits))
      return(res)
    }
    if(notification){
      # if(check.U(traits, scale = scale, ranks = ranks, ord, w = weights)){
      if(check.U(traits, scale = scale, ranks = ranks, asym.bin = asym.bin, ord = ord, w = traits.weights)){
        warning("Warning: NA in distance matrix between species based in traits", call. = FALSE)
      }
    }
    matrixX <- matrix.x(comm, traits, scale = scale, ranks = ranks, notification = FALSE,
                        asym.bin = asym.bin, ord = ord, w = traits.weights,
                        transformation = transformation, spp.weights = spp.weights)
    res$W <- matrixT$matrix.w
    res$B <- matrixT$matrix.b
    res$T <- matrixT$matrix.T
    res$U <- matrixX$matrix.u
    res$X <- matrixX$matrix.X
    res$traits.weights <- traits.weights
  }
  if (!is.null(phylodist)) {
    matrixP <- matrix.p(comm, phylodist, transformation = transformation, spp.weights = spp.weights, notification = FALSE)
    res$W <- matrixP$matrix.w
    res$Q <- matrixP$matrix.q
    res$P <- matrixP$matrix.P
  }
  if (!is.null(envir)) {
    E <- envir
    if (scale.envir) {
      E <- cent.norm(envir, na.rm = na.rm)
    }
    res$E <- E
  }
  res$Y <- matrix.w.transformation(res$W, transformation = "beals", reference = NULL,
                                   type = 3, include = TRUE, notification = FALSE)
  return(res)
}
