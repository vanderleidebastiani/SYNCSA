optimal<-function (comm, envir, traits, subset.min = 2, subset.max = 3, 
    pattern = "tcap", ro.method = "mantel", dist = "euclidean", 
    method = "pearson", scale = TRUE, scale.envir = TRUE, na.rm = FALSE, 
    notification = TRUE, put.together = NULL, progressbar = FALSE) 
{
	comm <- as.matrix(comm)
	envir <- as.matrix(envir)
	traits <- as.matrix(traits)
	roMETHOD <- c("mantel", "procrustes")
	romethod <- pmatch(ro.method, roMETHOD)
	if (length(romethod) > 1) {
		stop("\n Only one argument is accepted in ro.method \n")
	}
	if (is.na(romethod)) {
		stop("\n Invalid ro.method \n")
	}
	PATTERNS <- c("tcap", "tdap", "tcap.tdap")
	pattern <- pmatch(pattern, PATTERNS)
	if (length(pattern) > 1) {
		stop("\n Only one argument is accepted in pattern \n")
	}
	if (is.na(pattern)) {
		stop("\n Invalid pattern \n")
	}
	if (notification == TRUE) {
		c.NA <- apply(comm, 2, is.na)
		if (length(which(unique(as.vector(c.NA)) == TRUE)) > 0) {
			warning("Warning: NA in community data", call. = FALSE)
		}
		t.NA <- apply(traits, 2, is.na)
		if (length(which(unique(as.vector(t.NA)) == TRUE)) > 0) {
			warning("Warning: NA in traits matrix", call. = FALSE)
		}
		e.NA <- apply(envir, 2, is.na)
		if (length(which(unique(as.vector(e.NA)) == TRUE)) > 0) {
			warning("Warning: NA in environmental data", call. = FALSE)
		}
	}
	make.names<-is.null(colnames(traits))
	colnames(traits) <- colnames(traits, do.NULL = FALSE, prefix = "T")
	if (scale.envir == "TRUE") {
		envir <- cent.norm(envir, na.rm = na.rm)
	}
	if (romethod == 1) {
		dist.y <- vegdist(envir, method = dist, na.rm = na.rm)
	}
	m <- dim(traits)[2]
	p <- 1:subset.max
	names.traits<-colnames(traits)
	if(!missing(put.together)){
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
			names.traits[which(names.traits==put.together[[k]][1])]<-paste(put.together[[k]], collapse = " ")
			names.traits<-setdiff(names.traits,put.together[[k]][-1])
		}
		m<-length(names.traits)
		put.together2<-list()
		for(k in 1:length(put.together)){
			put.together2[[k]]<-paste(put.together[[k]], collapse = " ")
		}
	}
	if (subset.max > m) {
		stop("\n Subset must be lower than the number of traits\n")
	}
	bin <- factorial(m)/(factorial(p) * factorial(m - p))
	nT <- sum(bin[subset.min:subset.max])
	comb <- matrix(NA, nrow = sum(bin[subset.min:subset.max]), ncol = 1)
	n = 0
	for (i in subset.min:subset.max) {
		combinations <- combn(names.traits, i, simplify = TRUE) #mudei
		for (j in 1:bin[i]) {
			n = n + 1
			comb[n, 1] <- paste(combinations[, j], collapse = " ")
		}
	}
	n = 0
	correlation <- matrix(NA, nrow = sum(bin[subset.min:subset.max]), ncol = 1)
	for (i in subset.min:subset.max) {
		combinations1 <- combn(names.traits, i, simplify = TRUE)
		for (j in 1:bin[i]) {
			if (pattern == 1) {
				n = n + 1
				choose.traits<-combinations1[, j]
				if(!missing(put.together)){
					if(sum(match(choose.traits,unlist(put.together2)),na.rm=TRUE)>0){
						choose.traits2<-intersect(choose.traits,unlist(put.together2))
						choose.traits3<-c()
						for(k in 1:length(choose.traits2)){
							for(l in 1:length(put.together)){			        	
								if(choose.traits2[k]==put.together2[[l]]){
									choose.traits3<-c(choose.traits3,put.together[[l]])
								}
							}
						}
						choose.traits<-c(choose.traits3,setdiff(choose.traits,unlist(put.together2)))
					}
				}
				T <- matrix.t(comm, as.matrix(traits[, choose.traits]), scale = scale, notification = FALSE)
				if (romethod == 1) {
					correlation[n, 1] <- cor(vegdist(as.matrix(T$matrix.T), method = dist, na.rm = na.rm), dist.y, method = method)
					if (progressbar) {
						ProgressBAR(n, nT, style = 3)
					}
				}
				if (romethod == 2) {
					correlation[n, 1] <- procrustes.syncsa(T$matrix.T, envir)
					if (progressbar) {
						ProgressBAR(n, nT, style = 3)
					}
				}
			}
			if (pattern == 2) {
				n = n + 1
				choose.traits<-combinations1[, j]
				if(!missing(put.together)){
					if(sum(match(choose.traits,unlist(put.together2)),na.rm=TRUE)>0){
						choose.traits2<-intersect(choose.traits,unlist(put.together2))
						choose.traits3<-c()
						for(k in 1:length(choose.traits2)){
							for(l in 1:length(put.together)){			        	
								if(choose.traits2[k]==put.together2[[l]]){
									choose.traits3<-c(choose.traits3,put.together[[l]])
								}
							}
						}
						choose.traits<-c(choose.traits3,setdiff(choose.traits,unlist(put.together2)))
					}
				}
				T <- matrix.t(comm, as.matrix(traits[, choose.traits]), scale = scale, notification = FALSE)
				X <- matrix.x(comm, as.matrix(traits[, choose.traits]), scale = scale, notification = FALSE)
				if (romethod == 1) {
					dist.x <- vegdist(X$matrix.X, method = dist, na.rm = na.rm)
					dist.z <- vegdist(T$matrix.T, method = dist, na.rm = na.rm)
					rxy <- cor(dist.x, dist.y, method = method)
					rxz <- cor(dist.x, dist.z, method = method)
					ryz <- cor(dist.y, dist.z, method = method)
					correlation[n, 1] <- part.cor(rxy, rxz, ryz)
					if (progressbar) {
						ProgressBAR(n, nT, style = 3)
					}
				}
				if (romethod == 2) {
					correlation[n, 1] <- procrustes.partial(X$matrix.X, envir, T$matrix.T)
					if (progressbar) {
						ProgressBAR(n, nT, style = 3)
					}
				}
			}
			if (pattern == 3) {
				n = n + 1
				choose.traits<-combinations1[, j]
				if(!missing(put.together)){
					if(sum(match(choose.traits,unlist(put.together2)),na.rm=TRUE)>0){
						choose.traits2<-intersect(choose.traits,unlist(put.together2))
						choose.traits3<-c()
						for(k in 1:length(choose.traits2)){
							for(l in 1:length(put.together)){			        	
								if(choose.traits2[k]==put.together2[[l]]){
									choose.traits3<-c(choose.traits3,put.together[[l]])
								}
							}
						}
						choose.traits<-c(choose.traits3,setdiff(choose.traits,unlist(put.together2)))
					}
				}
				X <- matrix.x(comm, as.matrix(traits[, choose.traits]), scale = scale, notification = FALSE)
				if (romethod == 1) {
					correlation[n, 1] <- cor(vegdist(as.matrix(X$matrix.X), method = dist, na.rm = na.rm), dist.y, method = method)
					if (progressbar) {
						ProgressBAR(n, nT, style = 3)
					}
				}
				if (romethod == 2) {
					correlation[n, 1] <- procrustes.syncsa(X$matrix.X, envir)
					if (progressbar) {
						ProgressBAR(n, nT, style = 3)
					}
				}
			}
		}
	}
	result <- data.frame(Subset = comb, ro = correlation, stringsAsFactors = FALSE)
	result <- result[order(abs(result[, 2]), decreasing = TRUE), ]
	Res <- list(call = match.call(), N_subset = nT, optimization = result)
	class(Res) <- "optimal"
	return(Res)
}