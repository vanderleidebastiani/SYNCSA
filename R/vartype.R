#' Check the type of variables 
#' 
#' Function extracted (with small changes) of the function \code{\link{gowdis}} 
#'to check the type of variables in a dataframe or matrix.
#'
#' The variable types, where 'C' is continuous/numeric, 'O' is ordinal, 'B' is
#' symmetric binary and 'N' is nominal. 
#'
#'
#' @encoding UTF-8
#' @param data A dataframe or matrix.
#' @return \item{type}{The variable type.}
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso  \code{\link{syncsa}}, \code{\link{organize.syncsa}}
#' @keywords SYNCSA
#' @export
vartype<-function(data){
	if(class(data) != "data.frame" & class(data) != "matrix"){
		stop("data must be a matrix or a data.frame")
	}
	colnames(data)<-colnames(data,do.NULL = FALSE, "var")
	is.bin <- function(k) all(k[!is.na(k)] %in% c(0, 1))
	nc <- ncol(data)
	if (is.data.frame(data)) {
		type <- sapply(data, data.class)
		type2 <- type
		bin.var <- rep(NA, nc)
		for (i in 1:nc) {
			bin.var[i] <- is.bin(data[, i])
		}	
		type[type %in% c("numeric", "integer")] <- "C"
		type[type == "ordered"] <- "O"
		type[type %in% c("factor", "character")] <- "N"
		type[bin.var] <- "B" # rever
		type[type2 %in% c("factor", "character")] <- "N"
		#
		
		names(type)<-NULL
	}
    else {
    	if(any(sapply(data, data.class)=="character")){
    		stop("\n If data are a matrix class it must be entirely numeric \n")
    	}
        type <- rep("C", nc)
    }
return(type)
}