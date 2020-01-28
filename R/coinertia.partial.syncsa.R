#' @rdname coinertia.syncsa
#' @encoding UTF-8
#' @export
coinertia.partial.syncsa <- function(x, y, z, scale = FALSE){
  rxy.temp <- SYNCSA::coinertia.syncsa(x, y, scale = scale)
  rxz.temp <- SYNCSA::coinertia.syncsa(x, z, scale = scale)
  ryz.temp <- SYNCSA::coinertia.syncsa(y, z, scale = scale)
  RV <- SYNCSA::part.cor(rxy.temp, rxz.temp, ryz.temp)
  return(RV)
}
