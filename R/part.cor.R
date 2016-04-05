#' First-order partial correlation coefficient
#' 
#' Internal function to obtain the first-order partial correlation coefficient.
#' 
#' 
#' @encoding UTF-8
#' @param rxy Correlation between xy
#' @param rxz Correlation between xz
#' @param ryz Correlation between yz
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @keywords SYNCSA
#' @export
part.cor <- function(rxy, rxz, ryz)
{
	rxy<-round(rxy,10)
	rxz<-round(rxz,10)
	ryz<-round(ryz,10)
	if (rxz == 1 | ryz == 1) {
		res <- 0
	}   	
	else {
		res<-(rxy - rxz * ryz)/(sqrt(1 - rxz * rxz)*sqrt(1 - ryz * ryz))
	}
return(res)
}
