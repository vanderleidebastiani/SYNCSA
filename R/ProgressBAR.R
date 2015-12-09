#' Text Progress Bar
#' 
#' Text progress bar in the R console. See \code{\link{txtProgressBar}}.
#' 
#' 
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @keywords SYNCSA
ProgressBAR<-function(n,N,...){
	n=n/N
	for(i in c(0, n)){
		A<-txtProgressBar(...)	
		setTxtProgressBar(A,i)
	}
}
