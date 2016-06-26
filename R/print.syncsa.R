#' @rdname syncsa
#' @encoding UTF-8
#' @export
print.syncsa<-function(x , ...){
	cat("Call:\n")
    cat(deparse(x$call), "\n\n")
    cat("notes:\n")
	print(x$notes)
	cat("\n")
	cat("statistics:\n")
	print(x$statistics)
	invisible(x)
}