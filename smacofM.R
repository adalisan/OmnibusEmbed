#
# Embedding in-sample with weighted raw-stress criterion
#

smacofM <- function(D,
                    ndim    = 2,
                    W       = NULL,
                    init    = NULL,
                    verbose = FALSE,
                    itmax   = 1000,
                    eps     = 1e-6)
{

    n <- nrow(D)
    rnames <- rownames(D)
    rownames(D) <- NULL
    colnames(D) <- NULL

    if (is.null(W)) {
        W <- 1 - diag(n)
    }
    distW <- as.dist(W)
    distD <- as.dist(D)
    require(MASS)
    if (is.null(init)) {
        X <- cmdscale(D, ndim)
    } else {
        X <- init
    }
    rownames(X) <- NULL
    V <- -W
    diag(V) <- rowSums(W)
	if (verbose){
		
		print(V)
		print(W)
	}
		
    Vinv <- ginv(V)
    distE <- dist(X)
    stressOld <- sum(distW * (distD - distE)^2)

    ##--------------- begin majorization --------------------
    for (itel in 1:itmax) {
		
        B <- - as.matrix(distW * distD / ifelse(distE > 1e-08, distE, Inf))
        diag(B) <- - rowSums(B)
        Z <- Vinv %*% B %*% X
        distE <- dist(Z)
        stress <- sum(distW * (distD - distE)^2)

        if (verbose)
            cat("Iteration:", formatC(itel,width=3, format="d"),
                " Stress:", formatC(c(stressOld,stress),digits=5,width=7,format="f"),"\n")
        if (stressOld - stress < eps){
			stop.crit<-"eps"
            break()
		}
        
        X <- Z
        stressOld <- stress
		if (itel==itmax)
			stop.crit<-"itmax"
    }
	
    ##------------------ end majorization ---------------
    rownames(X) <- rnames
    X
}


smacofOos <- function(D,
					X,
					
					oos.flag,
					W       = NULL,
					
                    init    = NULL,
                    verbose = FALSE,
                    itmax   = 1000,
                    eps     = 1e-6)
{
	n <- sum(oos.flag)
	m<- sum(!oos.flag)
	n.dim<-ncol(X)
	rnames <- rownames(D)
	rownames(D) <- NULL
	colnames(D) <- NULL
	
	if (is.null(init))
		Z<-cmdscale(D[!oos.flag,!oos.flag],n.dim)
	else Z<- matrix(init,nrow=m)
	if (is.null(W)) {
		W <- 1 - diag(m+n)
	}
	distW <- as.dist(W)
	distD <- as.dist(D)
	rownames(X) <- NULL
	V <- -W
	diag(V) <- rowSums(W)
	V.m<-V[!oos.flag,!oos.flag]
	V.mn<-V[!oos.flag,oos.flag]
    Vinv<-ginv(V.m)
	distE <- dist(rbind(X,Z))
	stressOld <- sum(distW * (distD - distE)^2)
	
	
	for (itel in 1:itmax) {
		B <- - as.matrix(distW * distD / ifelse(distE > 1e-08, distE, Inf))
		diag(B) <- - rowSums(B)
		B.mn<-B[!oos.flag,oos.flag]
		B.m <-B[!oos.flag,!oos.flag]
		Y <- Vinv %*% (B.mn-V.mn) %*% X+Vinv%*%B.m%*%Z
		distE <- dist(rbind(X,Y))
		stress <- sum(distW * (distD - distE)^2)
		
		if (verbose)
			cat("Iteration:", formatC(itel,width=3, format="d"),
					" Stress:", formatC(c(stressOld,stress),digits=5,width=7,format="f"),"\n")
		if (stressOld - stress < eps)
			break()
		
		Z <- Y
		
		stressOld <- stress
	}
	##------------------ end majorization ---------------
	rownames(X) <- rnames[1:n]
	rownames(Y) <- rnames[(n+1):(n+m)]
	Y
}
