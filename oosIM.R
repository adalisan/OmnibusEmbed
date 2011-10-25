oosIM <- function(D, X,
                  init     = "random",
                  verbose  = FALSE,
                  itmax    = 1000,
                  eps      = 1e-8,
                  W        = NULL,
                  isWithin = NULL,
                  bwOos    = TRUE) {
  ## Input:
  ##       D   : the full (n+m)x(n+m) distance matrix.
  ##       X   : the within-sample embeddings.
  ##    init   : "random", "gower" or a given numerical matrix.
  ##   isWithin: NULL if the first n rows and colums in D correspond to
  ##             within-sample observations. Otherwise, specify a vector with
  ##             {0,1}-elements of length n+m. 1 indicates the corresponding
  ##             observation is within-sample and 0 means it is an out-of-sample
  ##             observation.
  ##          W: NULL if using equal weights. Otherwise, specify a (n+m)x(n+m)
  ##             weight matrix, whose up-left corner is the nxn 0 matrix. 
  ##   bwOos   : whether to use distances between oos observations.
  ##
  ## Return: a k x ncol(X) matrix

  N <- nrow(D)
  n <- nrow(X)
  m <- N - n
  p <- ncol(X)
  rnames <- rownames(D)
  unname(D)
  unname(X)

  if (!is.null(isWithin)) {
    i1 <- which(isWithin == 1)
    i0 <- which(isWithin == 0)
    i10 <- c(i1, i0)
    D <- D[i10, i10]
  }
  require(MASS)
  if (is.character(init) && init == "random") {
    ## Y <- mvrnorm(m, mu = colMeans(X), Sigma = diag(apply(X, 2, var)))
    Y <- matrix(mvrnorm(m, colMeans(X), diag(apply(X,2,var),ncol(X),ncol(X))))
    Y <- matrix(Y, m, p)
  } else if (is.character(init) && init == "gower") {
    ## initializing by Gower's formula
    ## more likely leads to a local minima
    diSq <- rowSums(X^2)
    dMat <- t(matrix(diSq, m, n) - D[-(1:n), 1:n]^2)
    Y <- t(ginv(t(X) %*% X) %*% t(X) %*% dMat)/2
  } else if (is.numeric(init)) {
    Y <- matrix(init, m, p)
  } else {
    stop("init shuld be 'random', 'gower', or a matrix!")
  }

  ## calculate V
  if (is.null(W)) {
    W = matrix(1, n+m, n+m)
    W[1:n, 1:n] <- 0
  }
  if (bwOos != TRUE) {
    W[-(1:n), -(1:n)] <- 0
  }

  V <- -W
  diag(V) <- colSums(W) - diag(W)
	#if(sum(is.na(V)>0)) print("is.na V ")

  Vmn <- V[-(1:n), 1:n]
  Vm <- V[-(1:n), -(1:n)]
  VmInv <- ginv(Vm)
  #if(sum(is.na(V)>0)) print("is.na VmInv ")
  ## distance *between* in-samples and out-of-samples
  dissBw <- D[-(1:n), 1:n]
  dEucBw <- t((outer(rowSums(X^2),rowSums(Y^2),"+") - 2*X%*%t(Y)))

  dEucBw <- sqrt(dEucBw+abs(dEucBw)/2)
#if(sum(is.na(dEucBw)>0)) print("is.na dEucBw ")
  if (bwOos == TRUE) {
    dissWn <- as.dist(D[-(1:n), -(1:n)])
    dEucWn <- dist(Y)
  } else {
    dissWn <- 0
    dEucWn <- 0
  }
#if(sum(is.na(dEucWn)>0)) {
#	print("is.na dEucWn ")
#	print(Y)
#	}
  stressOld <- sum((dissBw - dEucBw)^2) + sum((dissWn - dEucWn)^2)
   #if (is.na(stressOld)|is.nan(stressOld)) {
	#print("init matrices")
	#if (sum(is.na(dissWn))>0) print("dissWn")
	#print(dissWn)
	#if (sum(is.na(dEucBw))>0) print("dEucBw")
	#print(dEucBw)
	#
	#if (sum(is.na(dEucWn))>0) print("dEucWn")
	#print(dEucWn)
	#print("init matrices end ")
	#}
  for (itel in 1:itmax) {
    ## Bmn is an off-diagnoal blcok of B
	
	
	dEucBw[dEucBw<1e-8]<-Inf

	
    Bmn <- - W[-(1:n), 1:n]*dissBw / dEucBw 
    ## calculate Bm
    ## notice that the off-diagnoal entries of Bn should be all 0's
	#if (sum(is.na((dissBw)))) print("NA in dissBw")
	#if (sum(is.na((dEucBw)))) print("NA in dEucBw")
	#if (sum(is.na((dissWn)))) print("NA in dissWn")
	#if (sum(is.na((Bmn))))    print("NA in Bmn")
    if(m == 1 || bwOos == FALSE) {
      Bm <- - diag(rowSums(Bmn), m, m)
    } else {
	  dEucWn[dEucWn<1e-8]<-Inf
      Bm <- -W[-(1:n),-(1:n)]*as.matrix(dissWn / dEucWn)
      diag(Bm) <- - (rowSums(Bm) + rowSums(Bmn))
    }
	
	#if (sum(is.na((Bmn)))) print("NA in Bmn")
	#if (sum(is.na((Bm)))) print("NA in Bm")
	#if (sum(is.na((VmInv)))) print("NA in VmInv")
	#if (sum(is.na((Vmn)))) print("NA in Vmn")
    Z <- VmInv %*% ((Bmn - Vmn) %*% X + Bm %*% Y)

    dEucBw <- t((outer(rowSums(X^2),rowSums(Z^2),"+") - 2*X%*%t(Z)))
    dEucBw <- sqrt((dEucBw + abs(dEucBw) )/2)
    if (bwOos == TRUE) {
      dEucWn <- dist(Z) 
    } else {
      dEucWn <- 0
    }
    stress <- sum((dissBw - dEucBw)^2) + sum((dissWn - dEucWn)^2)
	#if (is.na(stress)|is.nan(stress)){
		#print(paste("Iter",itel,"begin"))
		#print("dissBw")
	#print(dissBw[1,1])
	#print("dEucBw")
	#print(dEucBw[1])
##
	#print("dEucWn")
	#print(dEucWn[1])
	#print(paste("Iter",itel,"end"))
	#}
    if (verbose == TRUE) {
      cat("Iteration:",
          formatC(itel, width=3, format="d"),
          " Stress:",
          formatC(c(stressOld,stress), digits=5, width=7, format="f"),"\n")
    }
    
    if (stressOld - stress < eps) {
      break()
    }

    Y <- Z
    stressOld <- stress
  }

  rownames(Y) <- rnames[-(1:n)]
  Y
}
