#' Syncsa
#'
#' This function integrates several steps for the analysis of phylogenetic
#' assembly patterns and their links to traits and ecological processes in a
#' metacommunity (Pillar et al. 2009, Pillar & Duarte 2010). It requires data
#' organized into the following matrices: (1) the presences or abundances of
#' species in a set of communities (\strong{W}); (2) the phylogenetic pairwise
#' dissimilarities of these species (\strong{DF}, in the range 0 to 1); (3) a
#' set of functional traits describing the species (\strong{B}), which may be a
#' mixture of binary and quantitative traits (ordinal, interval, ratio scales),
#' but not nominal ones (these should be expanded into binary traits); and (4)
#' the ecological gradient of interest, which may be one or more factors to
#' which the communities respond or ecosystem effects of the communities
#' (\strong{E}). The function computes several correlations (Mantel or
#' Procrustes) that express trait-convergence assembly patterns (TCAP),
#' trait-divergence assembly patterns (TDAP), and phylogenetic signal in
#' functional traits at the species pool level and at the metacomunity level.
#' This function also generates P-values by permutation testing based on null
#' models (Pillar et al. 2009, Pillar & Duarte 2010).
#'
#' The function implement methods that have been available in the SYNCSA
#' application written in C++ (by Valerio Pillar, available at
#' http://ecoqua.ecologia.ufrgs.br/ecoqua/SYNCSA.html).
#'
#' \strong{ro(TE)}
#'
#' This correlation refers to trait-convergence assembly patterns related to
#' the ecological gradient (TCAP, Pillar et al. 2009). For evaluating TCAP, by
#' matrix multiplication we define \strong{T = WB}, which with previous
#' standardization of \strong{W} to unit column totals will contain the trait
#' averages in each community. The elements in \strong{T} are community
#' weighted mean values or community functional parameters (Violle et al.
#' 2007). Standardization of the traits (rows) in \strong{T} is needed if the
#' trait set contains traits measured with different scales. By using matrix
#' correlation, we evaluate how the trait patterns in \strong{T} are associated
#' to ecological gradients in \strong{E}. For relating \strong{T} to
#' \strong{E}, using Mantel correlation we define a distance matrix of the
#' communities (\strong{DT}) using \strong{T}, and another distance matrix of
#' the community sites (\strong{DE}) using \strong{E}. The matrix correlation
#' ro(\strong{TE}) = ro(\strong{DT};\strong{DE}) measures the level of
#' congruence between TCAP and \strong{E}. A strong correlation ro(\strong{TE})
#' indicates the factors directly or indirectly represented in \strong{E} are
#' involved in ecological filtering of species that, at least for the traits
#' considered in the analysis, consistently produce trait-convergence assembly
#' patterns along the gradient comprising the metacommunity.
#'
#' \strong{ro(XE) and ro(XE.T)}
#'
#' These matrix correlations refer to trait-divergence assembly patterns
#' related to the ecological gradient (TDAP, Pillar et al. 2009). For the
#' identification of TDAP, in a first step the species pairwise similarities
#' (in the range 0 to 1) in matrix \strong{SB} based on traits in \strong{B}
#' are used to define matrix \strong{U} with degrees of belonging of species to
#' fuzzy sets. By matrix multiplication \strong{X = WU} will contain the
#' species composition of the communities after fuzzy-weighting by their trait
#' similarities (each row in \strong{X} will refer to a species). Matrix
#' \strong{X} expresses both TCAP and TDAP (Pillar et al. 2009). By using
#' matrix correlation, we evaluate how the trait patterns in \strong{X} (TCAP
#' and TDAP) are associated to ecological gradients in \strong{E}. For relating
#' \strong{X} to \strong{E}, we define a distance matrix of the communities
#' (\strong{DX}) using \strong{X}, and another distance matrix of the community
#' sites (\strong{DE}) using \strong{E}. The matrix correlation ro(\strong{XE})
#' = ro(\strong{DX};\strong{DE}) between \strong{X} and \strong{E} is defined.
#' We then remove the trait-convergence component ro(\strong{TE}) from
#' ro(\strong{XE}) by computing the partial matrix correlation
#' ro(\strong{XE.T}), which measures the level of congruence between TDAP and
#' \strong{E}. Trait-divergence assembly patterns (TDAP, Pillar et al. 2009)
#' may result from community assembly processes related to biotic interactions
#' (Stubbs & Wilson 2004; Wilson 2007).
#'
#' \strong{ro(PE)}
#'
#' This matrix correlation refers to the phylogenetic structure related to the
#' ecological gradient comprising the metacommunity. The phylogenetic pairwise
#' dissimilarities in \strong{DF} are transformed into similarities and used to
#' define degrees of belonging qij to fuzzy sets. This is analogous to the
#' definition of functional fuzzy sets (Pillar & Orloci 1991; Pillar et al.
#' 2009). Based on the phylogenetic similarities, every species i among s
#' species in the metacommunity specifies a fuzzy set to which every species j
#' (j = 1 to s species, including species i) belongs with a certain degree of
#' belonging in the interval [0, 1]. In our definition, each row in matrix
#' \strong{Q} with the degrees of belonging must add to unit, i.e., the degrees
#' of belonging of a given species across the fuzzy sets are standardized to
#' unit total. By matrix multiplication \strong{P = WQ} will contain the
#' composition of the communities after fuzzy-weighting of species presences or
#' abundances by the species` phylogenetic similarities. Each column in matrix
#' \strong{P} holds the phylogenetic structure of a community. The
#' standardization of \strong{Q} is essential for the community totals in each
#' column in \strong{W} remaining the same in \strong{P}. Further, matrix
#' \strong{W} is adjusted to unit column totals prior to multiplication, so
#' that the total richness or abundance within each community in \strong{W}
#' will be standardized. Matrix correlation ro(\strong{PE}) =
#' ro(\strong{DP};\strong{DE}) measures the strength of the association between
#' community distances based on their phylogenetic structure in \strong{DP} and
#' distances based on their ecological conditions (\strong{DE}). Further,
#' \strong{P} can be explored for phylogenetic patterns at the metacommunity
#' level by using, e.g., ordination techniques.
#'
#' \strong{ro(PT) and ro(PX.T)}
#'
#' These matrix correlations measure phylogenetic signal at the metacommunity
#' level related to TCAP and to TDAP. We define phylogenetic signal at the
#' metacommunity level related to TCAP (PSMT) as the correlation between the
#' phylogenetic structure described in matrix \strong{P} and the
#' trait-convergence structure described in matrix \strong{T}. For this, a
#' proper distance matrix (e.g. Euclidean distances) of communities
#' (\strong{DP}) is computed using \strong{P} and another distance matrix of
#' the same communities (\strong{DT}) is computed using \strong{T}. Then matrix
#' correlation ro(\strong{PT}) = ro(\strong{DP};\strong{DT}) will measure the
#' level of congruence between variation in \strong{P} and \strong{T}, which is
#' a measure of PSMT. A strong phylogenetic signal at the metacommunity level
#' is expected when communities that are more similar in terms of phylogenetic
#' structure are also similar regarding their average trait values. We also
#' define phylogenetic signal at the metacommunity level related to TDAP
#' (PSMX.T) as the partial matrix correlation ro(\strong{PX.T}) =
#' ro(\strong{DP};\strong{DX.DT}) between community distances DP computed on
#' phylogenetic structure and community distances \strong{DX} computed on
#' species composition after fuzzy-weighting by the species, or trait
#' similarities, removing the effect of TCAP (\strong{DT}). This is analogous
#' to TDAP.
#'
#' \strong{ro(BF)}
#'
#' This matrix correlation measures phylogenetic signal at the species pool
#' level (PSS, Pillar & Duarte 2010). We define PSS as the matrix correlation
#' ro(\strong{FB}) = ro(\strong{DF};\strong{DB}) between species phylogenetic
#' dissimilarities (already defined as matrix \strong{DF}) and species trait
#' dissimilarities (derived from already defined matrix \strong{SB}) computed
#' on any number of traits from matrix \strong{B}. The species pool refers to
#' the species present in the metacommunity.
#'
#' \strong{Additional matrix correlations}
#'
#' The matrix correlations ro(\strong{TE.P}) and ro(\strong{XE.P}) are also
#' computed, which may be useful for evaluating causal models in path analysis.
#'
#' \strong{Mantel correlations}
#'
#' \strong{Procrustes correlations}
#'
#' \strong{Partial correlations}
#'
#' \strong{Testing against null models}
#'
#' All the matrix correlations are tested against null models. The null model
#' is defined accoding to the correlation being tested. For ro(\strong{TE}),
#' each permutation generates a random matrix \strong{T} calculated after the
#' permutation among the species vectors in matrix \strong{B}. For
#' ro(\strong{XE}) and ro(\strong{XE.T}), each permutation generates a random
#' matrix \strong{X} after the permutation among species fuzzy sets (rows) in
#' matrix \strong{U}. For ro(\strong{PE}), ro(\strong{PT}), and
#' ro(\strong{PX.T}), each permutation generates a random matrix \strong{P}
#' after the permutation among species fuzzy sets (rows) in matrix \strong{Q}.
#' For ro(\strong{BF}), a conventional Mantel test is performed with
#' dissimilarity matrices \strong{DF} and \strong{DB}. Analogous null models
#' are used for testing the additional matrix correlations; that is, the same
#' null model for ro(\strong{TE}) is used for ro(\strong{TE.P}), the same model
#' for ro(\strong{XE}) is used for ro(\strong{XE.P}).
#'
#' @encoding UTF-8
#' @importFrom stats cor as.dist
#' @importFrom vegan wcmdscale protest vegdist
#' @importFrom permute how
#' @importFrom parallel makeCluster stopCluster
#' @param comm Community data, with species as columns and sampling units as
#' rows. This matrix can contain either presence/absence or abundance data.
#' @param traits Matrix data of species described by traits, with traits as
#' columns and species as rows.
#' @param dist.spp Matrix containing phylogenetic distance between species.
#' Must be a complete matrix (not a half diagonal matrix).
#' @param envir Environmental variables for each community, with variables as
#' columns and sampling units as rows.
#' @param ro.method Method to obtain the correlation, "mantel" or "procrustes"
#' (Default ro.method="mantel").
#' @param method Mantel correlation method, as accepted by cor: "pearson",
#' "spearman" or "kendall".
#' @param dist Dissimilarity index used for Mantel correlation, as accepted by
#' vegdist: "manhattan", "euclidean", "canberra", "bray", "kulczynski",
#' "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" ,
#' "binomial" or "chao". However, some of these will not make sense in this
#' case.
#' @param scale Logical argument (TRUE or FALSE) to specify if the traits are
#' measured on different scales (Default Scale = TRUE). Scale = TRUE if traits
#' are measured on different scales, the matrix T is subjected to
#' standardization within each trait. Scale = FALSE if traits are measured on
#' the same scale, the matrix T is not subjected to standardization.
#' Furthermore, if Scale = TRUE the matrix of traits is subjected to
#' standardization within each trait, and Gower Index is used to calculate the
#' degree of belonging to the species, and if Scale = FALSE the matrix of
#' traits is not subjected to standardization, and Euclidean distance is
#' calculated to determine the degree of belonging to the species.
#' @param scale.envir Logical argument (TRUE or FALSE) to specify if the
#' environmental variables are measured on different scales (Default Scale =
#' TRUE). If the enviromental variables are measured on different scales, the
#' matrix is subjected to centralization and standardization within each
#' variable.
#' @param permutations Number of permutations in assessing significance.
#' @param strata Argument to specify restricting permutations within species
#' groups (Default strata = NULL).
#' @param put.together List to specify group traits that are analysed
#' together (Default put.together = NULL). This argument must be a list, see
#' examples.
#' @param ord Method to be used for ordinal variables, see \code{\link{gowdis}}
#' (Default ord = "metric").
#' @param na.rm Logical argument (TRUE or FALSE) to specify if pairwise
#' distances should be deleted in cases of missing observations (Default na.rm
#' = FALSE).
#' @param notification Logical argument (TRUE or FALSE) to specify if
#' notification of missing observations should to be shown (Default
#' notification = TRUE).
#' @param parallel Number of parallel processes.  Tip: use detectCores() (Default parallel = NULL).
#' @param x An object of class syncsa.
#' @param ... Other parameters for the respective functions.
#' @return Correlations roTE, roXE, roPE, roPT, roPX.T, roXE.T, roTE.P, roXE.P
#' and roBF, and their significance levels based on permutations.
#' @note The function calculates the correlations despite the lack of one of
#' the matrices, provided that community data had been entered. Correlations
#' including unspecified matrices will appear with ro = 0.
#'
#' \strong{IMPORTANT}: The sequence of species in the community data matrix
#' MUST be the same as that in the phylogenetic distance matrix and in traits
#' matrix. Similarly, the sequence of communities in the community data matrix
#' MUST be the same as that in the environmental data matrix. See
#' \code{\link{organize.syncsa}}.
#'
#' The functions ignore missing data when specified. In the case of direct
#' multiplication of matrices (matrices \strong{W} and matrix \strong{B}) the
#' missing data are replaced by 0, ignoring the cell with missing value. Result
#' matrices are shown without missing values. Where the matrices are calculated
#' using a dissimilarity index (matrix \strong{U} and correlations between
#' matrices) the missing data are ignored as in \code{\link{vegdist}} function.
#' In some cases the dissimilarity matrices obtained by the function
#' \code{\link{vegdist}} still contain some missing values. In these cases the
#' rest of the procedure will be affected. In these cases you can find
#' solutions in the package mice.
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{matrix.t}}, \code{\link{matrix.x}},
#' \code{\link{matrix.p}}, \code{\link{optimal}}, \code{\link{belonging}},
#' \code{\link{organize.syncsa}}, \code{\link{rao.diversity}}
#' @references
#'
#' Pillar, V.D.; Duarte, L.d.S. (2010). A framework for metacommunity analysis
#' of phylogenetic structure. Ecology Letters, 13, 587:596.
#'
#' Pillar, V.D., Duarte, L.d.S., Sosinski, E.E. & Joner, F. (2009).
#' Discriminating trait-convergence and trait-divergence assembly patterns in
#' ecological community gradients. Journal of Vegetation Science, 20, 334:348.
#'
#' Pillar, V.D. & Orloci, L. (1991). Fuzzy components in community level
#' comparisons. In: Computer Assisted Vegetation Analysis (eds Feoli, E. &
#' Orloci, L.). Kluwer, Dordrecht, pp. 87:93.
#'
#' Stubbs, W.J. & Wilson, J.B. (2004). Evidence for limiting similarity in a
#' sand dune community. Journal of Ecology, 92, 557:567.
#'
#' Violle, C., Navas, M.L., Vile, D., Kazakou, E., Fortunel, C., Hummel, I. &
#' Garnier, E. (2007). Let the concept of trait be functional! Oikos, 116,
#' 882:892.
#'
#' Wilson, J.B. (2007). Trait-divergence assembly rules have been demonstrated:
#' limiting similarity lives! A reply to Grime. Journal of Vegetation Science,
#' 18, 451:452.
#' @keywords SYNCSA
#' @examples
#'
#' data(flona)
#' syncsa(comm=flona$community,traits=flona$traits,dist.spp=flona$phylo,envir=flona$environment)
#' syncsa(flona$community,traits=flona$traits,envir=flona$environment)
#'
#' @export
syncsa<-function (comm, traits, dist.spp, envir, ro.method = "mantel", method = "pearson", dist = "euclidean", scale = TRUE, scale.envir = TRUE, permutations = 999, strata = NULL, put.together = NULL, ord = "metric", na.rm = FALSE, notification = TRUE, parallel = NULL)
{
    N <- permutations
    roTE <- 0
    roXE <- 0
    roPE <- 0
    roPT <- 0
    roPX.T <- 0
    roXE.T <- 0
    roTE.P <- 0
    roXE.P <- 0
    roBF <- 0
    note.roTE <- paste("Trait-convergence assembly patterns (TCAP): roTE")
    note.roXE <- paste("Both trait-convergence assembly patterns and trait-divergence assembly patterns: roXE")
    note.roXE.T <- paste("Trait-divergence assembly patterns (TDAP): roXE.T")
    note.roBF <- paste("Phylogenetic signal at species level: roBF")
    note.roPE <- paste("Correlation of phylogenetically structured assembly patterns to ecological variables: roPE")
    note.roPT <- paste("Correlation of phylogenetically structured assembly patterns to trait-convergence assembly patterns: roPT")
    note.roPX.T <- paste("Correlation of phylogenetically structured assembly patterns to trait-divergence assembly patterns: roPX.T")
    note.roTE.P <- paste("Removing phylogeny from trait-convergence assembly patterns: roTE.P")
    note.roXE.P <- paste("Removing phylogeny from both trait-convergence assembly patterns and trait-divergence assembly patterns: roXE.P")
    note <- rbind(note.roTE, note.roXE, note.roXE.T, note.roBF, note.roPE, note.roPT, note.roPX.T, note.roTE.P, note.roXE.P)
    colnames(note) = "Correlation meanings"
    RES<-list()
    RES.matrices<-list()
	RES$notes<-note
    roMETHOD <- c("mantel", "procrustes")
    romethod <- pmatch(ro.method, roMETHOD)
    if (length(romethod) > 1) {
        stop("\n Only one argument is accepted in ro.method \n")
    }
    if (is.na(romethod)) {
        stop("\n Invalid ro.method \n")
    }
    if (!missing(comm)=="TRUE"){
		commvartype<-var.type(comm)
		if(any(commvartype=="n")){
			stop("\n comm must contain only numeric, binary or ordinal variables \n")
		}
	}
    if(notification==TRUE){
    	if (!missing(comm) == "TRUE") {
    		c.NA <- apply(comm, 2, is.na)
    		if(length(which(unique(as.vector(c.NA))==TRUE))>0)
    		{
				warning("Warning: NA in community data",call.=FALSE)
    		}
    	}
    	if (!missing(traits) == "TRUE") {
			t.NA <- apply(traits, 2, is.na)
    		if(length(which(unique(as.vector(t.NA))==TRUE))>0)
    		{
				warning("Warning: NA in traits matrix",call.=FALSE)
    		}
		}
		if (!missing(envir) == "TRUE") {
			e.NA <- apply(envir, 2, is.na)
			if(length(which(unique(as.vector(e.NA))==TRUE))>0)
    		{
				warning("Warning: NA in environmental data",call.=FALSE)
    		}
		}
    }
    if(!is.null(strata)){
    	if(length(strata)!=dim(comm)[2]){
    		stop("\n  The strata must be the same length of number of species \n")
    	}
    }
    seqpermutation<-permut.vector(dim(comm)[2],strata = strata, nset = permutations)
    if (!missing(traits) == "TRUE") {
		traitsvartype<-var.type(traits)
		if(any(traitsvartype=="n")){
			stop("\n trait must contain only numeric, binary or ordinal variables \n")
		}
	}
	if (!missing(dist.spp) == "TRUE") {
		dist.sppvartype<-var.type(dist.spp)	
		if(any(dist.sppvartype=="n")){
			stop("\n dist.spp must contain only numeric, binary or ordinal variables \n")
		}	
	}
	if (!missing(envir) == "TRUE") {
		envirvartype<-var.type(envir)
		if(any(envirvartype=="n")){
			stop("\n envir must contain only numeric, binary or ordinal variables \n")
		}
    }
    if(!is.null(parallel)){
    	CL <- parallel::makeCluster(parallel,type="PSOCK")	
    } else {
    	CL <- NULL
    }
    if (!missing(traits) == "TRUE") {	
    	m <- dim(traits)[2]
		weights<-rep(1,m)
		make.names<-is.null(colnames(traits))
		colnames(traits) <- colnames(traits, do.NULL = FALSE, prefix = "T")
		names(weights)<-colnames(traits)
		if(!is.null(put.together)){
			if(class(put.together)!="list"){
				stop("\n The put.together must be a object of class list\n")
			}
			if(make.names){
				for(k in 1:length(put.together)){
					put.together[[k]]<-paste("T", put.together[[k]],sep="")
				}
			}
			if(max(table(unlist(put.together)))>1){
				stop("\n The same trait appears more than once in put.together\n")
			}
			if(length(setdiff(unlist(put.together),colnames(traits)))>0){
				stop("\n Check traits names in put.together\n")
			}
			for(k in 1:length(put.together)){
				weights[put.together[[k]]]<-1/length(put.together[[k]])
			}
		}
        matrixT <- matrix.t(comm, traits, scale = scale, notification = FALSE)
        matrixX <- matrix.x(comm, traits, scale = scale, notification = FALSE, ord = ord, w = weights)
        W <- matrixT$matrix.w
        B <- matrixT$matrix.b
        T <- matrixT$matrix.T
        U <- matrixX$matrix.u
        X <- matrixX$matrix.X
        RES.matrices$W<-W
        RES.matrices$B<-B
        RES.matrices$T<-T
        RES.matrices$U<-U
        RES.matrices$X<-X
        RES$weights<-weights
        if (!missing(envir) == "TRUE") {
            E <- envir
            if (scale.envir == "TRUE") {
                E <- cent.norm(envir,na.rm = na.rm)
            }
            RES.matrices$E<-E
            if(romethod == 1){
            	roTE <- cor.matrix(mx1 = W, mx2 = B, x = T, y = E, method = method, dist = dist, permutations = N, norm = scale, strata = strata, na.rm = na.rm, seqpermutation = seqpermutation, parallel = parallel, newClusters = FALSE, CL = CL)
	            roXE <- cor.matrix(mx1 = W, mx2 = U, x = X, y = E, method = method, dist = dist, permutations = N, strata = strata, na.rm = na.rm, seqpermutation = seqpermutation, parallel = parallel, newClusters = FALSE, CL = CL)
    	        roXE.T <- cor.matrix.partial(mx1 = W, mx2 = U, x = X, y = E, mz1 = W, mz2 = B, z = T, permute.my2 = FALSE, permute.mz2 = TRUE, method = method, dist = dist, permutations = N, strata = strata, na.rm = na.rm, norm.z = scale, seqpermutation = seqpermutation, parallel = parallel, newClusters = FALSE, CL = CL)
        	}
        	if(romethod == 2){
            	roTE <- pro.matrix(mx1 = W, mx2 = B, x = T, y = E, permutations = N, norm = scale, strata = strata, seqpermutation = seqpermutation, parallel = parallel, newClusters = FALSE, CL = CL)
	            roXE <- pro.matrix(mx1 = W, mx2 = U, x = X, y = E, permutations = N, strata = strata, seqpermutation = seqpermutation, parallel = parallel, newClusters = FALSE, CL = CL)
    	        roXE.T <- pro.matrix.partial(mx1 = W, mx2 = U, x = X, y = E, mz1 = W, mz2 = B, z = T, permute.my2 = FALSE, permute.mz2 = TRUE, permutations = N, strata = strata, norm.z = scale, seqpermutation = seqpermutation, parallel = parallel, newClusters = FALSE, CL = CL)
        	}
        }
    }
    if (!missing(dist.spp) == "TRUE") {
        matrixP <- matrix.p(comm, dist.spp, notification = FALSE)
        W <- matrixP$matrix.w
        Q <- matrixP$matrix.q
        P <- matrixP$matrix.P
        RES.matrices$W<-W
        RES.matrices$Q<-Q
        RES.matrices$P<-P
        if (!missing(envir) == "TRUE") {
            E <- envir
            if (scale.envir == "TRUE") {
                E <- cent.norm(envir,na.rm = na.rm)
            }
            RES.matrices$E<-E
            if(romethod == 1){
            	roPE <- cor.matrix(mx1 = W, mx2 = Q, x = P, y = E, method = method, dist = dist, permutations = N, strata = strata, na.rm = na.rm, seqpermutation = seqpermutation, parallel = parallel, newClusters = FALSE, CL = CL)
            }
            if(romethod == 2){
            	roPE <- pro.matrix(mx1 = W, mx2 = Q, x = P, y = E, permutations = N, strata = strata, seqpermutation = seqpermutation, parallel = parallel, newClusters = FALSE, CL = CL)
            }
            if (!missing(traits) == "TRUE") {
            	if(romethod == 1){
                	roTE.P <- cor.matrix.partial(mx1 = W, mx2 = B, x = T, y = E, mz1 = W, mz2 = Q, z = P, permute.my2 = FALSE, permute.mz2 = FALSE, method = method, dist = dist, permutations = N, norm = scale, strata = strata, na.rm = na.rm, seqpermutation = seqpermutation, parallel = parallel, newClusters = FALSE, CL = CL)
	                roXE.P <- cor.matrix.partial(mx1 = W, mx2 = U, x = X, y = E, mz1 = W, mz2 = Q, z = P, permute.my2 = FALSE, permute.mz2 = FALSE, method = method, dist = dist, permutations = N, strata = strata, na.rm = na.rm, seqpermutation = seqpermutation, parallel = parallel, newClusters = FALSE, CL = CL)
                }
                if(romethod == 2){
                	roTE.P <- pro.matrix.partial(mx1 = W, mx2 = B, x = T, y = E, mz1 = W, mz2 = Q, z = P, permute.my2 = FALSE, permute.mz2 = FALSE, permutations = N, norm = scale, strata = strata, seqpermutation = seqpermutation, parallel = parallel, newClusters = FALSE, CL = CL)
	                roXE.P <- pro.matrix.partial(mx1 = W, mx2 = U, x = X, y = E, mz1 = W, mz2 = Q, z = P, permute.my2 = FALSE, permute.mz2 = FALSE, permutations = N, strata = strata, seqpermutation = seqpermutation, parallel = parallel, newClusters = FALSE, CL = CL)
                }
            }
        }
        if (!missing(traits) == "TRUE") {
        	if(romethod == 1){
	            roPT <- cor.matrix(mx1 = W, mx2 = Q, x = P, y = T, method = method, dist = dist, permutations = N, strata = strata, na.rm = na.rm, seqpermutation = seqpermutation, parallel = parallel, newClusters = FALSE, CL = CL)
    	        roPX.T <- cor.matrix.partial(mx1 = W, mx2 = Q, x = P, my1 = W, my2 = U, y = X, mz1 = W, mz2 = B, z = T, permute.my2 = FALSE, permute.mz2 = FALSE, method = method, dist = dist, permutations = N, strata = strata, na.rm = na.rm, norm.z = scale, seqpermutation = seqpermutation, parallel = parallel, newClusters = FALSE, CL = CL)
    	    }
    	    if(romethod == 2){
	            roPT <- pro.matrix(mx1 = W, mx2 = Q, x = P, y = T, permutations = N, strata = strata, seqpermutation = seqpermutation, parallel = parallel, newClusters = FALSE, CL = CL)
    	        roPX.T <- pro.matrix.partial(mx1 = W, mx2 = Q, x = P, my1 = W, my2 = U, y = X, mz1 = W, mz2 = B, z = T, permute.my2 = FALSE, permute.mz2 = FALSE, permutations = N, strata = strata, norm.z = scale, seqpermutation = seqpermutation, parallel = parallel, newClusters = FALSE, CL = CL)
    	    }
    	    if(romethod == 1){
            	dist.traits <- vegan::vegdist(traits, method = "euclidean", diag = TRUE, upper = TRUE, na.rm = na.rm) # VER TUDO DA BF, CLUSTER, PROCRUSTES
	            if (scale == "TRUE") {
    	            dist.traits <- vegan::vegdist(traits, method = "gower", diag = TRUE, upper = TRUE, na.rm = na.rm)
        	    }
            	roBF <- cor.mantel(dist.traits, stats::as.dist(dist.spp), method = method, permutations = N, strata = strata, na.rm = na.rm, seqpermutation = seqpermutation, parallel = parallel, newClusters = FALSE, CL = CL) #VER
	            #roBF <- c(BF$statistic, BF$signif)
            }
            if(romethod == 2){
			    #dist.spp.t <- sweep(dist.spp, 2, sqrt(apply(dist.spp^2,2,sum)), "/")
			    dist.spp.t<-dist.spp
			    vectors <- vegan::wcmdscale(dist.spp.t,eig = TRUE)$points
			    #vectors <- prcomp(dist.spp.t,scale = TRUE)$x
			    traits.t <- sweep(traits, 2, sqrt(apply(traits^2,2,sum)), "/")
				#BF <- suppressWarnings(vegan::protest(vectors,traits.t,permutations = permute::how(nperm = N, blocks = strata))) #VER
				roBF <- cor.procrustes(vectors,traits.t,permutations = N, strata = strata, na.rm = na.rm,seqpermutation = seqpermutation, parallel = parallel, newClusters = FALSE, CL = CL)
	            #roBF <- c(sqrt(1 - BF$ss), BF$signif)
            }
        }
    }
    SYNCSA <- rbind(roTE, roXE, roPE, roPT, roPX.T, roXE.T, roTE.P, roXE.P, roBF)
    if(!is.null(parallel)){
    	parallel::stopCluster(CL)
    }
    RES$statistics<-SYNCSA
    RES$call<-match.call()
    RES$matrices<-RES.matrices
    class(RES) <- "syncsa"
    return(RES)
}