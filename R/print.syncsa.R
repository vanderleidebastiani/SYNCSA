#' @rdname syncsa
#' @encoding UTF-8
#' @export
print.syncsa <- function(x, ...)
{
  cat("Call:\n")
  cat(deparse(x$call), "\n\n")
  cat("Notes:\n")
  print(x$notes)
  cat("\n")
  cat("Statistics:\n")
  print(x$statistics)
  invisible(x)
}