#' Generate dunny variable
#' 
#' Function to expand factor variables in dummy variables in a dataframe.
#'
#' Variables in the dataframe of class factor is expanded in dummy variables. Each
#' level of factor produce an new column in the dataframe, with presence (1) or
#' absense (0) of level. The name of columns is a combination of orginal variable name 
#' plus the level separate by underscore ( _ ). Ordered factor and character class are
#' not expanded.
#'
#'
#' @encoding UTF-8
#' @param data A dataframe.
#' @return \item{data}{The data with all variables of class factor expanded} \item{together}{A list
#' with sugestion to group of traits that are analysed together.}
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso  \code{\link{syncsa}}, \code{\link{organize.syncsa}}
#' @keywords SYNCSA
#' @export
var.dummy<-function(data){
	if(class(data) != "data.frame"){
		stop("data must be a data.frame")
	}
	colnames(data)<-colnames(data,do.NULL = FALSE, "var")
	type<-var.type(data)
	n<-dim(data)[1]
	RES <- as.data.frame(rep(NA,n))
	rownames(RES)<-rownames(data)
	m<-1
	l<-1
	together<-list()
	for(i in 1:dim(data)[2]){
		if(type[i]=="f"){
			xx_temp<-table(1:n,as.factor(data[,i]))
			colnames(xx_temp)<-paste(colnames(data)[i],colnames(xx_temp),sep="_")
			together[[l]]<-colnames(xx_temp)
			l<-l+1
			for(j in 1:dim(xx_temp)[2]){
				m<-m+1
				RES[,m]<-as.numeric(xx_temp[,j,drop=FALSE])
				colnames(RES)[m]<-colnames(xx_temp)[j]
			}
		}else{
			m<-m+1
			RES[,m]<-data[,i,drop=FALSE]
		}	
	}
	RES<-RES[,-1,drop=FALSE]
return(list(data=RES,together=together))
}