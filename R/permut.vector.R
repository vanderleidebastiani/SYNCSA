#' Permutate a vector
#' 
#' Internal function to permutate a vector of size n.
#' 
#' 
#' @encoding UTF-8
#' @importFrom permute how Within Plots shuffleSet
#' @param n The length of vector
#' @param strata Argument to specify restricting permutations 
#' @param nset The number of permutations to generate for the set
#' (Default strata = NULL).
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @keywords SYNCSA
#' @export
permut.vector <- function(n, strata = NULL, nset = 999){
	if(is.null(strata)){
		CTRL<-permute::how(within = permute::Within(type = "free"),plots = permute::Plots(type = "free"))
		samp<-permute::shuffleSet(n = n, nset = nset, control = CTRL, check = FALSE)
	}
	if(!is.null(strata)){
		if(n!=length(strata)){
			stop("\n The strata must be the length of n\n")
		}
		CTRL <- permute::how(within = permute::Within(type = "free"), plots = permute::Plots(strata = strata, type = "none"))
		samp<-permute::shuffleSet(n = n, nset = nset, control = CTRL, check = FALSE)
	}
	return(as.matrix(samp))
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