#' Text Progress Bar
#' 
#' Text progress bar in the R console. See \code{\link{txtProgressBar}}.
#' 
#' 
#' @encoding UTF-8
#' @param n Number of current progress
#' @param N Total number of cases
#' @param ... Other parameters for the txtProgressBar function
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
