
oosMDS <- function(D, X, w=c(rep(1, nrow(X)), rep(0, nrow(as.matrix(D))-nrow(X))),
                   init = "gower", itmax = 100) {
  ## Args:
  ##   D   : the full (n+k)x(n+k) distance matrix
  ##   X   : the within-sample (CMDS) embeddings
  ##   init: "random", "gower" or a numerical matrix
  ##   w   : a vector with {0,1}-elements of length n+k. 1 indicates the
  ##         corresponding observation is within-sample and 0 means it is
  ##         an out-of-sample observation.
  ## 
  n <- nrow(X)
  d <- ncol(X)
  k <- length(w) - sum(w)
  D <- as.matrix(D)
  if (nrow(D) != n+k || ncol(D) != n+k || length(w) != n+k) {
	print("n,d,k,len.w")
print(c(n,d,k,length(w)))
    stop("non-match length!")
  }
  rnames <- rownames(D)[which(w == 0)]  # rownames for out-of-samples
  unname(D)
  unname(X)
  B <- tau(D, w)

  if (init == "random") {
    require(MASS)
   # set.seed(12345)
    y0 <- matrix(mvrnorm(k, colMeans(X), diag(apply(X,2,var),ncol(X),ncol(X))))
  } else if (init == "gower") {
    diSq <- rowSums(X^2)
    dMat <- t(matrix(diSq, k, n) - D[-(1:n), 1:n]^2)
    y0 <- t(ginv(t(X) %*% X) %*% t(X) %*% dMat)/2
  } else {
    y0 <- matrix(init, k, p)
  }
  y0 <- as.vector(y0)
    
  nlm.out <- try(nlm(f, p=y0, B=B, X=X, iterlim=itmax))
  if (!is.list(nlm.out)){
  if (attr(nlm.out,"class")=="try-error"){
	print("oos Embedding failed,trying CMDS")
	tmp<-cmdscale(D[w==0,w==0],k=d)
	rownames(tmp)<-rnames
	return(tmp)
	}
	}
  Y <- matrix(nlm.out$estimate,ncol=d)
  rownames(Y) <- rnames
  Y
}


##* some needed functions

f <- function(y, B, X) {
  ## nlm optimization
  n <- nrow(X)
  d <- ncol(X)
  Y <- matrix(y,ncol=d)
  i <- 1:n
  res <- 2*sum((B[i,-i]- X %*% t(Y))^2) + sum((B[-i,-i]- Y %*% t(Y))^2)
  G <- -4*t(B[i,-i]- X %*% t(Y)) %*% X - 4*(B[-i,-i]- Y %*% t(Y)) %*% Y
  attr(res,"gradient") <- as.vector(G)
  res
}

mds.tau.w <- function(H, w) {
  ## This function computes tau_w for specified w
  n <- length(w)
  w <- matrix(w,ncol=1)
  e <- matrix(1,nrow=n,ncol=1)
  s <- sum(w)
  P <- diag(n) - e %*% t(w)/s
  Q <- diag(n) - w %*% t(e)/s
  -0.5 * P %*% H %*% Q
}

tau.e <- function(M) {
  ## double-centering
  n <- nrow(M)
  M <- M - matrix(colMeans(M), n, n, byrow=TRUE)
  -0.5*(M - matrix(rowMeans(M), n, n))
}

tau.w <- function(tau.e.Delta2, Delta2, a2) {
  n <- nrow(Delta2)
  P.Delta2 <- Delta2 - matrix(colMeans(Delta2),n,n,byrow=T)
  ## P %*% Delta2 %*% q      +       P %*% a2
  b <- -0.5 * ( - rowMeans(P.Delta2) + a2 - mean(a2) )
  ## t(q) %*% Delta2 %*% q
  beta <- -0.5 * (mean(Delta2) - 2 * mean(a2))
  rbind(cbind(tau.e.Delta2, b), cbind(t(b), beta))
}

tau <- function(fullD, w) {
  ## alternative to mds.tau.w
  ## difference: mds.tau.w requires fullD^2
  N <- length(w)
  n <- sum(w)
  k <- N - n
  if (k == 0)
    return(tau.e(A2))
  i1 <- which(w == 1)
  i0 <- which(w == 0)
  i10 <- c(i1, i0)
  fullD <- fullD[i10, i10]
  w <- w[i10]
  A2 <- as.matrix(fullD^2)
  D2 <- A2[1:n, 1:n]
  a2 <- A2[-(1:n), 1:n, drop=FALSE]
  m2 <- A2[-(1:n), -(1:n)]
  Bn <- tau.e(D2)
  Bxy <- -0.5 * scale(-matrix(rowSums(D2), n, k)/n + t(a2), center=TRUE, scale=FALSE)
  if (nrow(a2) == 1) {
    a2Jt <- - sum(a2)/n
  } else {
    a2Jt <- -matrix(rowSums(a2), k, k)/n
  }
  Byy <- -0.5 * (matrix(sum(D2), k, k)/n^2 + a2Jt + t(a2Jt) + m2)
  B <- rbind(cbind(Bn, Bxy), cbind(t(Bxy), Byy))
  B
}

## compare tau.w vs mds.tau.w
## Trosset and Priebe think that the oos observations are the last k ones
##
## N <- 1382
## oos <- 1300:1382
## k <- length(oos)
## n <- N - k
## w <- c(rep(1, n), rep(0, k))
## A2 <- DE^2
## t <- proc.time()
## B.old <- mds.tau.w(A2, w)
## proc.time() - t                         # 2.920   0.010   2.935
## t <- proc.time()
## B.new <- tau(A2,w)
## proc.time() - t                         # 0.230   0.000   0.229
## all.equal(B.new, B.old)                 # TRUE

## Delta2 <- (DE[-oos, -oos])^2
## a2 <- (DE[oos, -oos])^2
## t <- proc.time()
## tau.e.Delta2 <- tau.e(Delta2)
## B.new <- tau.w(tau.e.Delta2, Delta2, a2)
## proc.time() - t
## ## ##  user  system elapsed
## ## ## 0.536   0.321   1.008
## sum(round(B.new,6) == round(B.old,6))   # 1909924 == 1382^2 ==> all equal

##
## compute the embeddings for out-of-sample observations
##
## Input:
##   fullD: dissimilarity matrix (full): (n+k) x (n+k)
##       X: train matrix: n x d
##       w: the indices for out-of-sample observations
##
## Output:
##       Y: oos matrix of test: (N-n) x d
## 

