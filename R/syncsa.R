syncsa<-function (comm, traits, dist.spp, envir, ro.method = "mantel", method = "pearson", dist = "euclidean", scale = TRUE, scale.envir = TRUE, permutations = 999, strata = NULL, na.rm = FALSE, notification = TRUE) 
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
    roMETHOD <- c("mantel", "procrustes")
    romethod <- pmatch(ro.method, roMETHOD)
    if (length(romethod) > 1) {
        stop("\n Only one argument is accepted in ro.method \n")
    }
    if (is.na(romethod)) {
        stop("\n Invalid ro.method \n")
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
    if (!missing(traits) == "TRUE") {
        matrixT <- matrix.t(comm, traits, scale = scale, notification = FALSE)
        matrixX <- matrix.x(comm, traits, scale = scale, notification = FALSE)
        W <- matrixT$matrix.w
        B <- matrixT$matrix.b
        T <- matrixT$matrix.T
        U <- matrixX$matrix.u
        X <- matrixX$matrix.X
        if (!missing(envir) == "TRUE") {
            E <- envir
            if (scale.envir == "TRUE") {
                E <- cent.norm(envir,na.rm = na.rm)
            }
            if(romethod == 1){
            	roTE <- cor.matrix(W, B, T, E, method = method, dist = dist, permutations = N, norm = scale, strata = strata, na.rm = na.rm)
	            roXE <- cor.matrix(W, U, X, E, method = method, dist = dist, permutations = N, strata = strata, na.rm = na.rm)
    	        roXE.T <- cor.matrix.partial(W, U, X, E, T, method = method, dist = dist, permutations = N, strata = strata, na.rm = na.rm)
        	}
        	if(romethod == 2){
            	roTE <- pro.matrix(W, B, T, E, permutations = N, norm = scale, strata = strata)
	            roXE <- pro.matrix(W, U, X, E, permutations = N, strata = strata)
    	        roXE.T <- pro.matrix.partial(W, U, X, E, T, permutations = N, strata = strata)
        	}
        }
    }
    if (!missing(dist.spp) == "TRUE") {
        matrixP <- matrix.p(comm, dist.spp, notification = FALSE)
        W <- matrixP$matrix.w
        Q <- matrixP$matrix.q
        P <- matrixP$matrix.P
        if (!missing(envir) == "TRUE") {
            E <- envir
            if (scale.envir == "TRUE") {
                E <- cent.norm(envir,na.rm = na.rm)
            }
            if(romethod == 1){
            	roPE <- cor.matrix(W, Q, P, E, method = method, dist = dist, permutations = N, strata = strata, na.rm = na.rm)
            }
            if(romethod == 2){
            	roPE <- pro.matrix(W, Q, P, E, permutations = N, strata = strata)
            }
            if (!missing(traits) == "TRUE") {
            	if(romethod == 1){
                	roTE.P <- cor.matrix.partial(W, B, T, E, P, method = method, dist = dist, permutations = N, norm = scale, strata = strata, na.rm = na.rm)
	                roXE.P <- cor.matrix.partial(W, U, X, E, P, method = method, dist = dist, permutations = N, strata = strata, na.rm = na.rm)
                }
                if(romethod == 2){
                	roTE.P <- pro.matrix.partial(W, B, T, E, P, permutations = N, norm = scale, strata = strata)
	                roXE.P <- pro.matrix.partial(W, U, X, E, P, permutations = N, strata = strata)
                }
            }
        }
        if (!missing(traits) == "TRUE") {
        	if(romethod == 1){
	            roPT <- cor.matrix(W, Q, P, T, method = method, dist = dist, permutations = N, strata = strata, na.rm = na.rm)
    	        roPX.T <- cor.matrix.partial(W, Q, P, X, T, method = method, dist = dist, permutations = N, strata = strata, na.rm = na.rm)
    	    }
    	    if(romethod == 2){
	            roPT <- pro.matrix(W, Q, P, T, permutations = N, strata = strata)
    	        roPX.T <- pro.matrix.partial(W, Q, P, X, T, permutations = N, strata = strata)
    	    }
    	    if(romethod == 1){    
            	dist.traits <- vegdist(traits, method = "euclidean", diag = TRUE, upper = TRUE, na.rm = na.rm)
	            if (scale == "TRUE") {
    	            dist.traits <- vegdist(traits, method = "gower", diag = TRUE, upper = TRUE, na.rm = na.rm)
        	    }
            	BF <- mantel(dist.traits, dist.spp, method = method, permutations = N, strata = strata, na.rm = na.rm)
	            roBF <- c(BF$statistic, BF$signif)
            }
            if(romethod == 2){
			    dist.spp.t <- sweep(dist.spp, 2, sqrt(apply(dist.spp^2,2,sum)), "/")
			    vectors <- wcmdscale(dist.spp.t,eig = TRUE)$points
				BF <- suppressWarnings(protest(vectors,traits,permutations = how(nperm = N, blocks = strata)))
	            roBF <- c(sqrt(1 - BF$ss), BF$signif)
            }
        }
    }
    SYNCSA <- rbind(roTE, roXE, roPE, roPT, roPX.T, roXE.T, roTE.P, roXE.P, roBF)
    return(list(Notes = note, Statistics = SYNCSA))
}