
#'  Embedding in-sample dissimilarities with weighted raw-stress criterion using SMACOF
#'
#' @param D
#' @param ndim
#' @param W
#' @param init
#' @param verbose
#' @param itmax
#' @param eps
#' @param debug.mode
#'
#' @return
#' @export
#'
#' @examples
smacofM <- function(D,
                    ndim    = 2,
                    W       = NULL,
                    init    = NULL,
                    verbose = FALSE,
                    itmax   = 1000,
                    eps     = 1e-6,
					debug.mode =FALSE)
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
    if ((sum(is.null(distW))>0) |(sum(is.na(distW))>0)) {
   stop("invalid distW: stress can't be computed")}
    if ((sum(is.null(distD))>0) |(sum(is.na(distD))>0) ) {
   stop("invalid distD: strees can't be computed")}
    require(MASS)
    require(MCMCpack)
    if (is.null(init)) {
        X <- cmdscale(D, ndim)
    } else {
        X <- init
    }


   if ((sum(is.null(X))>0) |(sum(is.na(X))>0)| (sum(is.infinite(X))>0)| (ncol(X)==0)) {
    X <- mvrnorm(n, mu = rep(0,ndim), Sigma = max(D,na.rm=TRUE)*diag(ndim))
       X <- matrix(X, nrow=n, ncol=ndim)
    if (debug.mode)  {sink("debug.X.n.txt")
    print(D)
   print(ndim)
    print(X)
    print(str(X))

        sink()
  }
  }
  if (debug.mode)  {
 sink("debug.X.txt")
   print(ndim)
    print(X)
    print(str(X))

        sink()
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
 if ((sum(is.null(distE))>0) |(sum(is.na(distE))>0) ) {
   sink("debug.distE.txt")
   print(X)
   print(distE)
        sink()
   stop("invalid distE: strees can't be computed")}

    stressOld <- sum(distW * (distD - distE)^2)
   if ((is.null(stressOld))|(is.na(stressOld))) {
   stop("Old Stress can't be computed")}


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


#'
#'  Embedding out-of-sample dissimilarities with weighted raw-stress criterion using SMACOF
#'
#' @param D
#' @param X
#' @param oos.flag
#' @param W
#' @param init
#' @param verbose
#' @param itmax
#' @param eps
#'
#' @return
#' @export
#'
#' @examples
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
