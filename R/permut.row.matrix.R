#' Permutate rows in a matrix
#' 
#' Internal function to permutate rows in a matrix.
#' 
#' 
#' @encoding UTF-8
#' @importFrom permute how Within Plots shuffleSet
#' @param data A matrix
#' @param strata Argument to specify restricting permutations in rows 
#' (Default strata = NULL).
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @keywords SYNCSA
#' @export
permut.row.matrix <- function(data, strata = NULL){
	N<-dim(data)[1]
	if(is.null(strata)){
		CTRL<-permute::how(within = permute::Within(type = "free"),plots = permute::Plots(type = "free"))
		samp<-as.vector(permute::shuffleSet(n = N, nset = 1, control = CTRL, check = FALSE))
	}
	if(!is.null(strata)){
		if(N!=length(strata)){
			stop("\n The strata must be the same length of number of rows\n")
		}
		CTRL <- permute::how(within = permute::Within(type = "free"), plots = permute::Plots(strata = strata, type = "none"))
		samp<-as.vector(permute::shuffleSet(n = N, nset = 1, control = CTRL, check = FALSE))
	}
	permut.matrix <- data[samp,,drop=FALSE]
	res<-list(permut.matrix = permut.matrix,samp = samp)
	return(res)
}


# Mudar blocos de lugar (somente se os blocos tiverem o mesmo tamanho)
#CTRL <- how(within = Within(type = "none"),
#            plots = Plots(strata = strata, type = "free"))
#shuffleSet(N, 1, control = CTRL, check = FALSE)
#numPerms(N, control = CTRL)
#
# Mudar blocos e parcelas dentro de blocos (somente se os blocos tiverem o mesmo tamanho)
#CTRL <- how(within = Within(type = "free"),
#            plots = Plots(strata = strata, type = "free"))
#shuffleSet(N, 1, control = CTRL, check = FALSE)
#numPerms(N, control = CTRL)