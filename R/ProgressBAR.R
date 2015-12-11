#' Text Progress Bar
#' 
#' Text progress bar in the R console. See \code{\link{txtProgressBar}}.
#' 
#' @param n is a number
#' @param N is the total number
#' @param ... additional arguments to txtProgressBar
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @keywords SYNCSA
#' @export
ProgressBAR<-function(n,N,...){
	n=n/N
	for(i in c(0, n)){
		A<-txtProgressBar(...)	
		setTxtProgressBar(A,i)
	}
}
