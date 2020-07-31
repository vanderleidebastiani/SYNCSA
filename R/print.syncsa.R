#' @rdname syncsa
#' @include syncsa.R
#' @encoding UTF-8
#' @export
print.syncsa <- function(x, ...)
{
  cat("Call:\n")
  cat(deparse(x$call), "\n\n")
  if(length(x$list.warning)>0){
    cat("List of warning:\n")
    namestemp <- CollectNames(x$list.warning, prefix = "$list.warning$")
    cat(unlist(namestemp, use.names = FALSE), sep = "\n")
  }
  cat("\nList of matrices:\n")
  namestemp <- CollectNames(x$matrices, prefix = "$matrices$")
  cat(unlist(namestemp[!sapply(x$matrices, is.null)], use.names = FALSE), sep = "\n")
  if(length(x$notes)>0){
    cat("\n$notes:\n")
    print(x$notes)
    cat("\n$statistics:\n")
    print(as.matrix(x$statistics), ...)
  }
  if(length(x$matrices.null)>0){
    cat("\nList of null matrices:\n")
    namestemp <- paste("$matrices.null$", names(x$matrices.null), sep = "")
    cat(unlist(namestemp, use.names = FALSE), sep = "\n")
  }
  invisible(x)
}
