#' @rdname optimal
#' @encoding UTF-8
#' @export
print.optimal <- function(x, ...)
{
  n <- 10
  if (dim(x$optimization)[1] < n) {
    n <- dim(x$optimization)[1]
  }
  cat("Call:\n")
  cat(deparse(x$call), "\n\n")
  cat("Number of subset:\n")
  cat(deparse(x$N_subset), "\n\n")
  cat("Trait subset optimization:\n")
  print(x$optimization[1:n,])
  invisible(x)
}