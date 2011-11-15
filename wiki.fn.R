## Time-stamp: <wiki.fn.R zma 2010-03-03 18:40>

get.a2 <- function(M, newObs, method="average") {
  ## get the vector a2 used for out-of-sample embedding
  N <- nrow(M)/2
  k <- length(newObs)
  n <- N - k/2

  ## for each pair of e_new and f_new (no matter they match or not),
  ## use (u1 + v2)/2 as the estimate of v1 and u2.
  get.a2.average <- function() {
    pair <- matrix(newObs[c(1,3,2,4,1,4,2,3)], 4, 2, byrow=T)
    a2 <- matrix(0, 2*k, 2*n)
    rownames(a2) <- c("11.e","11.f","22.e","22.f", "12.e","12.f","21.e","21.f")
    for (i in 1:k) {
      u1 <- M[pair[i,1], -newObs][1:n]             # dist(e, Se)
      v2 <- M[pair[i,2], -newObs][(n+1):(2*n)]     # dist(f, Sf)
      u2 <- v1 <- (u1 + v2)/2                      # dist(f, Se), dist(e, Sf)
      a2[(2*i-1):(2*i),] <- rbind(c(u1,v1), c(u2,v2))^2
    }
    a2
  }

  get.a2.sqrtSumSq <- function() {
    pair <- matrix(newObs[c(1,3,2,4,1,4,2,3)], 4, 2, byrow=T)
    a2 <- matrix(0, 2*k, 2*n)
    rownames(a2) <- c("11.e","11.f","22.e","22.f", "12.e","12.f","21.e","21.f")
    for (i in 1:k) {
      u1 <- M[pair[i,1], -newObs][1:n]             # dist(e, Se)
      v2 <- M[pair[i,2], -newObs][(n+1):(2*n)]     # dist(f, Sf)
      u2 <- v1 <- sqrt(u1^2 + v2^2)                # dist(f, Se), dist(e, Sf)
      a2[(2*i-1):(2*i),] <- rbind(c(u1,v1), c(u2,v2))^2
    }
    a2
  }

  get.a2.sqrtSumSq2 <- function() {
    pair <- matrix(newObs[c(1,3,2,4,1,4,2,3)], 4, 2, byrow=T)
    a2 <- matrix(0, 2*k, 2*n)
    rownames(a2) <- c("11.e","11.f","22.e","22.f", "12.e","12.f","21.e","21.f")
    for (i in 1:k) {
      u1 <- M[pair[i,1], -newObs][1:n]             # dist(e, Se)
      v2 <- M[pair[i,2], -newObs][(n+1):(2*n)]     # dist(f, Sf)
      u2 <- v1 <- sqrt(u1^2 + v2^2)/2              # dist(f, Se), dist(e, Sf)
      a2[(2*i-1):(2*i),] <- rbind(c(u1,v1), c(u2,v2))^2
    }
    a2
  }

  ## W == 1 - I ==> u2 = v1 = (1,...,1)
  get.a2.oneMinusI <- function() {
    u1 <- M[newObs[1:2], -newObs][,1:n]           # dist(e, Se)
    v2 <- M[newObs[3:4], -newObs][,(n+1):(2*n)]   # dist(f, Sf)
    a2 <- rbind(cbind(u1, matrix(1,2,n)), cbind(matrix(1,2,n), v2))^2
    rownames(a2) <- c("1.e","2.e","1.f", "2.f")
    a2
  }

  ## use u1 (or v2) itself as v1 (or u2)
  get.a2.self <- function() {
    u1 <- M[newObs[1:2], -newObs][,1:n]           # dist(e, Se)
    v2 <- M[newObs[3:4], -newObs][,(n+1):(2*n)]   # dist(f, Sf)
    a2 <- rbind(cbind(u1, u1), cbind(v2, v2))^2
    rownames(a2) <- c("1.e","2.e","1.f", "2.f")
    a2
  }

  ## find the k nearest neighbors of e_new (or f_new) in the orignal
  ## dissmilarity space, then use the colMeans of corresponding row of
  ## W as the estimate of v1 (or u2).
  get.a2.tunnel.nnDissmW <- function() {
    u1 <- M[newObs[1:2], -newObs][,1:n]           # dist(e, Se)
    v2 <- M[newObs[3:4], -newObs][,(n+1):(2*n)]   # dist(f, Sf)
    u2 <- matrix(0, 2, n); v1 <- matrix(0, 2, n)
    for (i in 1:2) {
      r.e <- order(u1[i,] + 10e-6*rnorm(n))[1:k.nn] # break tie randomly
      v1[i,] <- colMeans(M[-newObs,-newObs][r.e, (n+1):(2*n)])
      r.f <- order(v2[i,] + 10e-6*rnorm(n))[1:k.nn] # break tie randomly
      u2[i,] <- colMeans(M[-newObs,-newObs][r.f+n, 1:n])
    }
    a2 <- rbind(cbind(u1, v1), cbind(u2, v2))^2
    rownames(a2) <- c("1.e","2.e","1.f", "2.f")
    a2
  }

  ## find the k nearest neighbors of e_new (or f_new) in the orignal
  ## dissmilarity space, then use the colMeans of corresponding row of
  ## DF (or DE) as the estimate of u1 (or v2)'s pair, then use their
  ## average as the estimate of v1 (or u2).
  get.a2.tunnel.nnDissmD <- function() {
    u1 <- M[newObs[1:2], -newObs][,1:n]           # dist(e, Se)
    v2 <- M[newObs[3:4], -newObs][,(n+1):(2*n)]   # dist(f, Sf)
    u2 <- matrix(0, 2, n); v1 <- matrix(0, 2, n)
    for (i in 1:2) {
      r.e <- order(u1[i,] + 10e-6*rnorm(n))[1:k.nn] # break tie randomly
      v1[i,] <- (colMeans(M[-newObs,-newObs][r.e+n, (n+1):(2*n)]) + u1[i,]) / 2
      r.f <- order(v2[i,] + 10e-6*rnorm(n))[1:k.nn] # break tie randomly
      u2[i,] <- (colMeans(M[-newObs,-newObs][r.f, 1:n]) + v2[i,]) / 2
    }
    a2 <- rbind(cbind(u1, v1), cbind(u2, v2))^2
    rownames(a2) <- c("1.e","2.e","1.f", "2.f")
    a2
  }

  ## find k nearest neighbors of u1 (or v2), use the colMeans of
  ## corresponding rows of W as the estimate of v1 (or u2).
  get.a2.tunnel.nnW <- function() {
    u1 <- M[newObs[1:2], -newObs][,1:n]           # dist(e, Se)
    v2 <- M[newObs[3:4], -newObs][,(n+1):(2*n)]   # dist(f, Sf)
    DE.n <- M[-newObs,-newObs][1:n,1:n]
    DF.n <- M[-newObs,-newObs][(n+1):(2*n),(n+1):(2*n)]
    u2 <- matrix(0, 2, n); v1 <- matrix(0, 2, n)
    for (i in 1:2) {
      r.e <- order(rowSums((DE.n - matrix(u1[i,], n, n, byrow=T))^2))[1:k.nn]
      v1[i,] <- colMeans(M[-newObs,-newObs][r.e, (n+1):(2*n)])
      r.f <- order(rowSums((DF.n - matrix(v2[i,], n, n, byrow=T))^2))[1:k.nn]
      u2[i,] <- colMeans(M[-newObs,-newObs][r.f+n, 1:n])
    }
    a2 <- rbind(cbind(u1, v1), cbind(u2, v2))^2
    rownames(a2) <- c("1.e","2.e","1.f", "2.f")
    a2
  }

  ## find k nearest neighbors of u1 (or v2), use the colMeans of
  ## corresponding rows of DF (or DE) as the estimate of u1 (or v2)'s
  ## pair, then use their average as the estimate of v1 (or u2)
  get.a2.tunnel.nnD <- function() {
    u1 <- M[newObs[1:2], -newObs][,1:n]           # dist(e, Se)
    v2 <- M[newObs[3:4], -newObs][,(n+1):(2*n)]   # dist(f, Sf)
    DE.n <- M[-newObs,-newObs][1:n,1:n]
    DF.n <- M[-newObs,-newObs][(n+1):(2*n),(n+1):(2*n)]
    u2 <- matrix(0, 2, n); v1 <- matrix(0, 2, n)
    for (i in 1:2) {
      r.e <- order(rowSums((DE.n - matrix(u1[i,], n, n, byrow=T))^2))[1:k.nn]
      v1[i,] <- (colMeans(M[-newObs,-newObs][r.e+n, (n+1):(2*n)]) + u1[i,]) / 2
      r.f <- order(rowSums((DF.n - matrix(v2[i,], n, n, byrow=T))^2))[1:k.nn]
      u2[i,] <- (colMeans(M[-newObs,-newObs][r.f, 1:n]) + v2[i,]) / 2
    }
    a2 <- rbind(cbind(u1, v1), cbind(u2, v2))^2
    rownames(a2) <- c("1.e","2.e","1.f", "2.f")
    a2
  }
  ## end-of-define functions

  ## extract k.nn
  strs <- unlist(strsplit(method, "\\."))
  if (1 %in% match(strs, "nnW")      ||
      1 %in% match(strs, "nnD")      ||
      1 %in% match(strs, "nnDissmW") ||
      1 %in% match(strs, "nnDissmD")) {
    method <- paste(strs[1], strs[3], sep=".")
    k.nn <- as.integer(strs[2])
  }

  switch(method,
         average          = get.a2.average(),
         oneMinusI        = get.a2.oneMinusI(),
         self             = get.a2.self(),
         tunnel.nnDissmW  = get.a2.tunnel.nnDissmW(),
         tunnel.nnDissmD  = get.a2.tunnel.nnDissmD(),
         tunnel.nnW       = get.a2.tunnel.nnW(),
         tunnel.nnD       = get.a2.tunnel.nnD(),
         sqrtSumSq        = get.a2.sqrtSumSq(),
         sqrtSumSq2       = get.a2.sqrtSumSq2())
}

wOutOfSample <- function(M, d = NDIM, newObs, method = "average")
{
    X <- cmdscale(M[-newObs, -newObs], d)
    a2 <- get.a2(M, newObs, method = method)
    Delta2 <- M[-newObs,-newObs]^2
    tau.e.Delta2 <- tau.e(Delta2)
    d <- ncol(X)
    k <- nrow(a2)
    Y <- matrix(0, k, d)
    rownames(Y) <- rownames(a2)
    set.seed(12345)
    for (i in 1:k) {
        B <- tau.w(tau.e.Delta2, Delta2, a2[i,])
        y0 <- rnorm(d)*10
        nlm.out <- nlm(f,p=y0,B=B,X=X)
        Y[i,] <- matrix(nlm.out$estimate,ncol=d)
    }
    Y
}

## 11.e 11.f can be newObs embedded simultaneously, because the distance between
## them are 0 (i.e., a = 0). I tried and the results are almost the same as
## embedding them separately.
## ## separate
## 11.e  0.3960020 0.1645827 0.6892661 1.177684 0.3097395 -0.2239567
## 11.f -0.7106901 0.4715190 0.8072221 1.343949 0.5543344 -0.3660075
## ## simultaneous
## 11.e  0.3954620 0.1652384 0.6906816 1.180742 0.3111866 -0.2252196
## 11.f -0.7103875 0.4717460 0.8084286 1.346626 0.5551401 -0.3667754

wInSample <- function(M, d = NDIM, newObs, method = "average")
{
    D <- M[-newObs, -newObs]
    n <- nrow(D)
    a2 <- get.a2(M, newObs, method = method)
    k <- nrow(a2)
    Y <- matrix(0, k, d)
    rownames(Y) <- rownames(a2)
    d2List <- vector("list", k/2)
    names(d2List) <- c("d11", "d22", "d12", "d21")
    for (i in 1:(k/2)) {
        uv <- sqrt(a2[(2*i-1):(2*i), , drop=FALSE]) # uv must be matrix
        Dnew <- rbind(cbind(D,  t(uv)), cbind(uv, matrix(0,2,2)))
        X <- cmdscale(Dnew, d)
        Y[(2*i-1):(2*i), ] <- X[-(1:n), ]
        d2match <- rowSums((X[1:(n/2), ] - X[(n/2+1):(n), ])^2)
        d2unmatch <- rowSums((X[2:(n/2), ] - X[(n/2+1):(n-1), ])^2)
        d2new <- sum((Y[2*i-1, ] - Y[2*i,])^2)
        d2List[[i]] <- list(d2match=d2match, d2unmatch=d2unmatch, d2new=d2new)
    }
    list(Y = Y, d2List = d2List)
}

pOutOfSample <- function(DE, DF, d = NDIM, newObs)
{
    suppressMessages(require("vegan"))
    XE <- cmdscale(DE[-newObs, -newObs], d)
    XF <- cmdscale(DF[-newObs, -newObs], d)
    pro <- procrustes(XF, XE, scale=FALSE)
    XErot <- pro$Yrot
    Q <- pro$rotation    # no need for translation, because XE
                         # and XF have column means equal to 0

    ## out-of-sample embedding
    YE <- YF <- matrix(0, 2, NDIM)
    rownames(YE) <- c("e1", "e2")
    rownames(YF) <- c("f1", "f2")

    a2E <- (DE[newObs, -newObs])^2
    a2F <- (DF[newObs, -newObs])^2
    Delta2E <- (DE[-newObs, -newObs])^2
    Delta2F <- (DF[-newObs, -newObs])^2
    tau.Delta2E <- tau.e(Delta2E)
    tau.Delta2F <- tau.e(Delta2F)
    set.seed(12345)
    for (i in 1:2) {
        B  <- tau.w(tau.Delta2E, Delta2E, a2E[i,])
        y0 <- rnorm(d)*10
        nlm.out <- nlm(f, p=y0, B=B, X=XE)
        YE[i,] <- matrix(nlm.out$estimate, ncol = d)

        B  <- tau.w(tau.Delta2F, Delta2F, a2F[i,])
        y0 <- rnorm(d)*10
        nlm.out <- nlm(f, p=y0, B=B, X=XF)
        YF[i,] <- matrix(nlm.out$estimate, ncol = d)
    }
    YErot <- YE %*% Q
    rbind(YErot, YF)
}


pInSample <- function(DE, DF, d = NDIM, newObs)
{
    suppressMessages(require("vegan"))
    N <- nrow(DE)
    YE <- YF <- matrix(0, 2, NDIM)
    rownames(YE) <- c("e1", "e2")
    rownames(YF) <- c("f1", "f2")
    for (i in 1:2) {
        idx <- c( (1:N)[-newObs], newObs[i] )
        XE <- cmdscale(DE[idx, idx], d)
        XF <- cmdscale(DF[idx, idx], d)
        XErot <- procrustes(XF, XE, scale=FALSE)$Yrot
        YE[i,] <- XErot[N-1,]   # length(newObs) = 2
        YF[i,] <- XF[N-1,]
    }
    rbind(YE, YF)
}


getDist2 <- function(Y)
{
    n <- nrow(Y)
    D2 <- as.matrix(dist(Y)^2)
    if (n == 4) {
        rownames(Y) <- c("e1", "e2", "f1", "f2")
        d2 <- c(D2["e1","f1"], D2["e2","f2"], D2["e1","f2"], D2["e2","f1"])
    } else if (n == 8) {
          d2 <- c(D2["11.e", "11.f"], D2["22.e", "22.f"],
                  D2["12.e", "12.f"], D2["21.e", "21.f"])
    }
    d2
}


getPower <- function(d2)
{
    n <- nrow(d2)
    ecdf.A <- ecdf(d2[,3:4])
    alpha <- seq(0, 1, by=0.05)
    power <- rep(0, length(alpha))
    for (i in 1:20) {
        power[i] <- 1 - ecdf.A(sort(d2[,1:2])[round(2*n*(1-alpha[i]))])
    }
    power[21] <- 1
    power
}


getTrueRates <- function(d0, dA, alternative = "greater")
{
    n0 <- length(d0)
    n1 <- length(dA)
    calpha <- c(-Inf, sort(unique(c(d0, dA))), Inf)
    if (alternative == "greater") {
        notRejectH0 <- t(sapply(calpha, function(xx) d0 < xx))
        acceptHA <- t(sapply(calpha, function(xx) dA >= xx))
    } else if (alternative == "less") {
        notRejectH0 <- t(sapply(calpha, function(xx) d0 > xx))
        acceptHA <- t(sapply(calpha, function(xx) dA <= xx))
    } else {
        stop("alternative must be either 'greater' or 'less'!")
    }
    specificity <- rowSums(notRejectH0)/n0
    sensitivity <- rowSums(acceptHA)/n1
    cbind(specificity, sensitivity)
}


## getAlphaBeta <- function(d0, dA)
## {
##     calpha <- c(-Inf, seq(xrange[1], xrange[2], length=length), Inf)
##     cdfHA <- ecdf(d2[,3:4])
##     cdfH0 <- ecdf(d2[,1:2])
##     alpha <- 1 - cdfH0(calpha)
##     power <- 1 - cdfHA(calpha)
##     cbind(alpha = alpha, beta = power)
## }


getAlphaBeta <- function(d2, xrange = range(d2), length = 1000)
{
    calpha <- c(-Inf, seq(xrange[1], xrange[2], length=length), Inf)
    cdfHA <- ecdf(d2[,3:4])
    cdfH0 <- ecdf(d2[,1:2])
    alpha <- 1 - cdfH0(calpha)
    power <- 1 - cdfHA(calpha)
    cbind(alpha = alpha, beta = power)
}


##
## getTrueRates calculates the true positive and true negative rates for
## different critical values. (tpr == sensitivity; tnr == specificity )
##
## Args:
##   cdf0, cdfA  : the cdfs corresponding to H0 and HA
##   calpha      : the set of critical values
##   alternative : indicates the alternative hypothesis and must be '"less"' or
##                 '"greater"'
##
## Returns:
##   a matrix with length(calpha) rows and 2 columns [specificity, sensitivity]
##
## getTrueRates <- function(cdf0, cdfA, alternative="greater", calpha)
## {
##     if (alternative == "greater") {
##         specificity <- cdf0(calpha)
##         sensitivity <- 1 - cdfA(calpha)
##     } else if (alternative == "less") {
##         specificity <- 1 - cdf0(calpha)
##         sensitivity <- cdfA(calpha)
##     }
##     cbind(specificity = specificity, sensitivity = sensitivity)
## }


get.d2 <- function(d2List)
{
    d2 <- sapply(d2List, function(xx) sapply(xx, function(yy) yy$d2new))
    t(d2)
}


##
## determine the percentage of matched pairs among the k nearest neighbors of a
## pair of new documents (e_new, f_new)
##
## Args:
##   DE, DF: the full distance matrices (including the new observations)
##   newObs: the indices of the *two* new observations
##           newObs = 1:2 ==> oos = c(1:2, 1:2 + 1382)
##        k: the # of nearest neighbors in consideration
##
## Returns:
##   a vector of length 4, indicating the percentage of matched pairs among the
##   kNN of the two new documents. e.g., if the 3NN of e_new are e_i, e_j and
##   e_k, while the 3NN of f_new are f_j, f_p and f_q, then the percentage is
##   1/k = 1/3 = 0.33
##
matchedKnn <- function(DE, DF, newObs, k)
{
    n <- nrow(DE) - length(newObs)
    eknn <- t(apply(DE[newObs, -newObs], 1, order, rnorm(n)))[, 1:k]
    fknn <- t(apply(DF[newObs, -newObs], 1, order, rnorm(n)))[, 1:k]
    pairs <- matrix(c(1,1, 2,2, 1,2, 2,1), ncol=2, byrow=T)
    nMatched <- apply(pairs, 1, function(r)
                      length(intersect(eknn[r[1], ], fknn[r[2], ])))
    nMatched / k
}

##
## the last one is the new observation
##
matchedKnn1 <- function(DE, DF, k)
{
    n <- nrow(DE)
    eknn <- order(DE[n, -n], rnorm(n-1))[1:k]
    fknn <- order(DF[n, -n], rnorm(n-1))[1:k]

    length(intersect(eknn, fknn))/k
}


##*
##* matchedKnn + w-approach
##*
##       D1     D2     A1 B2
##     +------+------+
##     |      |      |
##     |      |      |
##     +------+------+
##     |      |      |
##     |      |      |
##     +------+------+
## A1  --r----   ?      0  a
## B2    ?     --s----  a  0
##
## - D1 and D2 are a known matched pair
## - A1 and B2 are the two new docs
## - d(A1, D1) = r
## - d(B2, D2) = s
## - need to determine the value of
##   1) a
##   2) d(A1, D2)
##   3) d(B2, D1)
##
## IDEA:
## - a should be large when the percentage of matchedKnn (p) is small
## - if A1 and B2 are indeed a matched pair, then
##    - d(A1, D2) == d(B2, D1) == (r + s)/2 =: t
## - if A1 and B2 are not a matched pair
##   . let H1, I1 and J1 be the 3 nearest neighbors of A1
##   . let L2, M2 and N2 be the 3 nearest neighbors of B2
##   . let the corresponding matched pairs of H1, I1, J1, L2, M2 and N2
##     are H2, I2, J2, L1, M1 and N1
##   then
##   - d(A1, D2) can be approximated by the average of d(H1, D2), d(I1, D2)
##     and d(J1, D2)
##   - d(B2, D1) can be approximated by the average of d(L2, D1), d(M2, D1)
##     and d(N2, D1)
##
## APPROACH:
## - let u1 and v2 defined as before
## - a = p * 0 + (1 - p) * mean( (u1 + v2)/2 )
## - d(A1, D2) = p * t + (1 - p) * [d(H1, D2) + d(I1, D2) + d(J1, D2)]/3
## - d(B2, D1) = p * t + (1 - p) * [d(L2, D1) + d(M2, D1) + d(N2, D1)]/3
##
## OUTPUT:
##   X  : the within-sample embeddings obtained by w-approach
##   Y  : the out-of-sample ebeedings obtained by matchedKnn+w-approach
##   Td : the squared distances between the four pairs of new documents
wOutOfSampleMatchedKnn <- function(M, d=NDIM, newObs, returnX = FALSE)
{
    N <- nrow(M)/2
    n <- N - length(newObs)/2
    M0 <- M[-newObs, -newObs]
    X <- cmdscale(M0, d)
    p <- matchedKnn(DE = M[1:N, 1:N],
                    DF = M[(N+1):(2*N), (N+1):(2*N)],
                    newObs[1:2], 10)
    pair <- matrix(newObs[c(1,3,2,4,1,4,2,3)], 4, 2, byrow=T)

    Td <- rep(0, 4)
    Y <- matrix(0, 8, d)
    rownames(Y) <- c("11.e","11.f","22.e","22.f", "12.e","12.f","21.e","21.f")
    for(i in 1:4) {
        u1 <- M[pair[i,1], -newObs][1:n]
        v2 <- M[pair[i,2], -newObs][(n+1):(2*n)]
        u2Avg <- v1Avg <- (u1 + v2)/2
        e3nn <- order(u1, rnorm(n))[1:3]
        f3nn <- order(v2, rnorm(n))[1:3]
        u2Knn <- colMeans(M0[e3nn, (n+1):(2*n)])
        v1Knn <- colMeans(M0[f3nn, 1:n])
        u2 <- p[i] * u2Avg + (1 - p[i]) * u2Knn
        v1 <- p[i] * v1Avg + (1 - p[i]) * v1Knn
        a <- p[i] * 0 + (1 - p[i]) * mean((u1 + v2)/2)
        uv <- rbind(c(u1, u2), c(v1, v2))
        A <- matrix(c(0, a, a, 0), nrow=2, byrow=T)
        fullD <- rbind(cbind(M0, t(uv)), cbind(uv, A))
        w <- c(rep(1,2*n), 0, 0)
        Y[(2*i-1):(2*i), ] <- oosMDS(fullD, X, w)
        Td[i] <- as.vector(dist(Y[(2*i-1):(2*i), ]))
    }

    if (returnX == TRUE) {
        list(X=X, Y=Y, Td=Td^2)
    } else {
        list(Y=Y, Td=Td^2)
    }

}

matchManifolds <- function(D0, D1, label,
                           J.red   = 0:2,
                           NDIM    = 6,
                           method  = "P",
                           block   = FALSE,
                           crit    = "strain",
                           use.knn = TRUE,
                           k       = 3) {
  ## Xi_0, represented by D0, has both red and blue objects in the training data
  ## Xi_1, represented by D1, has only red objects in the trainin data. The goal
  ## is to build a classifier based on the blue objects in Xi_0 to classify the
  ## blue objects in Xi_1
  ##
  ## Args:
  ##   D0, D1  : the two Distance matrices in manifold mathcing figure
  ##             correspding to Xi_0 and Xi_1
  ##   label   : labels for training data
  ##   J.red   : classes available in the training data for both Xi_0 and Xi_1
  ##   NDIM    : embedding dimension
  ##   method  : "P" or "W" approach
  ##   block   : TRUE means not assuming the 1 - 1 correspondence between the
  ##             red objects in Xi_0 and Xi_1
  ##   crit    : "strain" implies the use of cmdscale and TP's out-of-sample
  ##             embedding; "stress" implies the use of smacof procedures to do
  ##             within-sample and out-of-sample embedding
  ##   use.knn : whether to use knn together with TP's oos in W-approach
  ##   k       : k-nn
  ##
  ## Returns:
  ##   X0.red  : the embeddings for the red objects in Xi_0
  ##   X1.red  : the embeddings for the red objects in Xi_1
  ##   X0.blue : the embeddings for the blue objects in Xi_0
  ##   X1.blue : the embeddings for the blue objects in Xi_1
  ##   rho     : the posteriors for the blue objects in Xi_1
  ##   class   : the predicted class labels for the blue objects in Xi_1

  suppressMessages(require("MASS"))
  if (crit == "strain") {
    inMds  <- "cmdscale"
    oosMds <- "oosMDS"
  } else {
    inMds  <- "smacofM"
    oosMds <- "smacofOos"
  }
  N          <- length(label)
  r.red      <- (1:N)[label %in% J.red]
  r.blue     <- (1:N)[-r.red]
  n.red      <- length(r.red)
  n.blue     <- length(r.blue)
  D0.red     <- D0[r.red, r.red]
  D1.red     <- D1[r.red, r.red]
  D0.blue    <- D0[r.blue, r.blue]
  D1.blue    <- D1[r.blue, r.blue]
  label.red  <- label[r.red]
  label.blue <- label[r.blue]

  if (method == "P") {
    suppressMessages(require("vegan"))
    X0.red <- match.fun(inMds)(D0.red, NDIM)
    X1.red <- match.fun(inMds)(D1.red, NDIM)
    w <- as.integer(label %in% J.red)

    if (block == TRUE) {
      m0.red <- apply(X0.red, 2, function(xx) tapply(xx, label.red, mean))
      m1.red <- apply(X1.red, 2, function(xx) tapply(xx, label.red, mean))
      ## Q <- procrustes(m0.red, m1.red, scale=FALSE)$rotation
      proc <- procrustes(m0.red, m1.red, scale=TRUE)
      Q <- proc$rotation * proc$scale
      ## m0.red and m1.red may not be column-centered. Therefore it possibly
      ## involves translation to Procrustean transform m1.red to m0.red.
      ## Then should the translatin be incorporated in transforming X1.blue
      ## to X0.blue? The answer to this question should be the same as the
      ## answer to the question: should the translatin be incorporated in
      ## transforming X1.red to X0.red? I think the answer is NO, because
      ## both X0.red and X1.red are column-centered.
    } else {
      ## Q <- procrustes(X0.red, X1.red, scale=FALSE)$rotation
      proc <- procrustes(X0.red, X1.red, scale=TRUE)
      Q <- proc$rotation * proc$scale
    }
    X0.blue <- match.fun(oosMds)(D0, X0.red, w=w)
    X1.blue <- match.fun(oosMds)(D1, X1.red, w=w) %*% Q
  } else if (method == "W") {
    W.red <- (D0.red + D1.red)/2
    if (block == TRUE) {
      for (i in J.red) {
        block.rows <- (label.red == i)
        for (j in J.red) {
          block.cols <- (label.red == j)
          if (i == j) {
            W.red[block.rows, block.cols] <- 0
          } else {
            W.red[block.rows, block.cols] <- mean(W.red[block.rows, block.cols])
          }
        }
      }
    }
    M.red <- rbind(cbind(D0.red, W.red), cbind(W.red, D1.red))
    X.red <- match.fun(inMds)(M.red, NDIM)
    X0.red <- X.red[1:n.red, ]
    X1.red <- X.red[-(1:n.red), ]

    if (use.knn == FALSE) { # then TP's oosMDS solutions aren't optimal
      w <- as.integer(label %in% J.red)
      X0.red <- scale(X.red[1:n.red, ], center=TRUE, scale=FALSE)
      X1.red <- scale(X.red[-(1:n.red), ], center=TRUE, scale=FALSE)
      ## mX0.red <- apply(X.red[1:n.red, ], 2, mean)
      ## mX1.red <- apply(X.red[-(1:n.red), ], 2, mean)
      X0.blue <- match.fun(oosMds)(D0, X0.red, w=w) # + matrix(mX0.red, k, NDIM, byrow=T)
      X1.blue <- match.fun(oosMds)(D1, X1.red, w=w) # + matrix(mX1.red, k, NDIM, byrow=T)
      ## X0.blue <- scale(X0.blue, scale=FALSE)        # column center
      ## X1.blue <- scale(X1.blue, scale=FALSE)        # column center
    } else {
      ## use knn to construct an omnibus matrix for out-of-sample embedding
      knn0 <- t(apply(D0[r.blue, r.red], 1, order, rnorm(n.red)))[, 1:k]
      knn1 <- t(apply(D1[r.blue, r.red], 1, order, rnorm(n.red)))[, 1:k]

      W0   <- t(apply(knn0, 1, function(r) colMeans(W.red[r, ])))
      W1   <- t(apply(knn1, 1, function(r) colMeans(W.red[r, ])))
      mat0 <- cbind(D0[r.blue, r.red], W0)
      mat1 <- cbind(D1[r.blue, r.red], W1)
      M0   <- rbind(cbind(M.red, t(mat0)), cbind(mat0, D0.blue))
      M1   <- rbind(cbind(M.red, t(mat1)), cbind(mat1, D0.blue))
      ## start oos embedding
      w <- c(rep(1, 2*n.red), rep(0, n.blue))
      X0.blue <- match.fun(oosMds)(M0, X.red, w=w)
      X1.blue <- match.fun(oosMds)(M1, X.red, w=w)
    }
  } else {
    stop("method must be either 'P' or 'W'!")
  }
  z <- lda(X0.blue, label.blue)
  pred <- predict(z, X1.blue)
  rho <- pred$posterior
  class <- pred$class

  list(X0.red  = X0.red,
       X1.red  = X1.red,
       X0.blue = X0.blue,
       X1.blue = X1.blue,
       rho     = rho,
       class   = class)
}


getEpsilon <- function(D0, D1, label,
                       J.red  = 0:2,
                       NDIM   = 6,
                       method = "P",
                       block  = FALSE) {
  suppressMessages(require("MASS"))
  N          <- length(label)
  r.red      <- (1:N)[label %in% J.red]
  r.blue     <- (1:N)[-r.red]
  n.red      <- length(r.red)
  n.blue     <- length(r.blue)
  D0.red     <- D0[r.red, r.red]
  D1.red     <- D1[r.red, r.red]
  label.red  <- label[r.red]
  label.blue <- label[r.blue]

  if (method == "P") {
    suppressMessages(require("vegan"))
    X0.red <- cmdscale(D0.red, NDIM)
    X1.red <- cmdscale(D1.red, NDIM)
    if (block == TRUE) {
      m0.red <- apply(X0.red, 2, function(xx) tapply(xx, label.red, mean))
      m1.red <- apply(X1.red, 2, function(xx) tapply(xx, label.red, mean))
      Q <- procrustes(m0.red, m1.red, scale=FALSE)$rotation
    } else {
      Q <- procrustes(X0.red, X1.red, scale=FALSE)$rotation
    }
    X1.red <- X1.red %*% Q
  } else if (method == "W") {
    W.red <- (D0.red + D1.red)/2
    if (block == TRUE) {
      for (i in J.red) {
        block.rows <- (label.red == i)
        for (j in J.red) {
          block.cols <- (label.red == j)
          if (i == j) {
            W.red[block.rows, block.cols] <- 0
          } else {
            W.red[block.rows, block.cols] <- mean(W.red[block.rows, block.cols])
          }
        }
      }
    }
    M.red <- rbind(cbind(D0.red, W.red), cbind(W.red, D1.red))
    X.red <- cmdscale(M.red, NDIM)
    X0.red <- X.red[1:n.red, ]
    X1.red <- X.red[-(1:n.red), ]
  }

  DX0.red <- dist(X0.red)
  DX1.red <- dist(X1.red)

  ## ## sum_{i<j} [delta_ij - d_ij(X)]^2 / sum_{i<j} [delta_ij^2 + d_ij(X)^2]
  ## e0 <- sum((as.dist(D0.red) - DX0.red)^2) /
  ##   (sum(as.dist(D0.red)^2) + sum(DX0.red^2))
  ## e1 <- sum((as.dist(D1.red) - DX1.red)^2) /
  ##   (sum(as.dist(D1.red)^2) + sum(DX1.red^2))
  ## ec <- sum((DX0.red - DX1.red)^2) /
  ##   (sum(DX0.red^2) + sum(DX1.red^2))

  ## strain
  ## source("~/work/r_functions/oosMDS.R")
  ## XX0.red <- tcrossprod(X0.red)
  ## XX1.red <- tcrossprod(X1.red)
  ## B0.red <- tau.e(D0.red^2)
  ## B1.red <- tau.e(D1.red^2)

  ## e0 <- sum((XX0.red - B0.red)^2) / sum(XX0.red^2 + B0.red^2)
  ## e1 <- sum((XX1.red - B1.red)^2) / sum(XX1.red^2 + B1.red^2)
  ## ## ec <- sum((XX0.red - XX1.red)^2) / sum(XX0.red^2 + XX1.red^2)
  ## ec <- sum((DX0.red - DX1.red)^2) /
  ##   (sum(DX0.red^2) + sum(DX1.red^2))

  ## 1 - congruence
  source("~/work/r_functions/zmR/R/CoxC.R")
  e0 <- 1 - CoxC(D0.red, as.matrix(DX0.red))
  e1 <- 1 - CoxC(D1.red, as.matrix(DX1.red))
  ec <- 1- CoxC(DX0.red, DX1.red)

  c(e0 = e0, e1 = e1, ec = ec, e = e0 + e1 + ec)
}


## old

manifoldMatching <- function(DL, DR, label, J0=0:2, NDIM=6, method="P", block=FALSE)
{
    n <- length(label)
    r0 <- (1:n)[label %in% J0]
    r1 <- (1:n)[-r0]
    n0 <- length(r0)
    n1 <- length(r1)

    if (method == "P") {
        suppressMessages(require("vegan"))
        XLJ0 <- cmdscale(DL[r0, r0], NDIM)
        XRJ0 <- cmdscale(DR[r0, r0], NDIM)
        if (block == TRUE) {
            mXLJ0 <- apply(XLJ0, 2, function(xx) tapply(xx, label[r0], mean))
            mXRJ0 <- apply(XRJ0, 2, function(xx) tapply(xx, label[r0], mean))
            translation <- apply(mXRJ0, 2, mean)
            Q <- procrustes(mXLJ0, mXRJ0, scale=FALSE)$rotation
        } else {
            translation <- apply(XRJ0, 2, mean)
            Q <- procrustes(XLJ0, XRJ0, scale=FALSE)$rotation
        }
        XLJ1 <- oosMDS(DL, XLJ0, w=as.integer(label %in% J0))
        XRJ1 <- oosMDS(DR, XRJ0, w=as.integer(label %in% J0))
        XRJ1 <- (XRJ1 - matrix(translation, n1, NDIM, byrow=TRUE)) %*% Q
    } else if (method == "W") {
        W0 <- (DL[r0, r0] + DR[r0, r0])/2
        if (block == TRUE) {
            for (i in J0) {
                for (j in J0) {
                    W0[label[r0]==i, label[r0]==j] <- mean(W0[label[r0]==i, label[r0]==j])
                }
            }
        }
        M0 <- rbind(cbind(DL[r0, r0], W0), cbind(W0, DR[r0, r0]))
        XJ0 <- cmdscale(M0, NDIM)
        ## XLJ0 <- scale(XJ0[1:n0, ], center=TRUE, scale=FALSE)
        ## XRJ0 <- scale(XJ0[-(1:n0), ], center=TRUE, scale=FALSE)
        XLJ0 <- XJ0[1:n0, ]
        XRJ0 <- XJ0[-(1:n0), ]
        XLJ1 <- oosMDS(DL, XLJ0, w=as.integer(label %in% J0))
        XRJ1 <- oosMDS(DR, XRJ0, w=as.integer(label %in% J0))
    } else {
        stop("method must be either 'P' or 'W'!")
    }
    suppressMessages(require("MASS"))
    z <- lda(XLJ1, label[r1])
    pred <- predict(z, XRJ1)
    rho <- pred$posterior
    class <- pred$class

    list(XLJ0=XLJ0, XRJ0=XRJ0, XLJ1=XLJ1, XRJ1=XRJ1, rho=rho, class=class)
}


## ##
## ## given full matrix D and the leaved-out samples, return a2
## ##
## get.a2.4x4 <- function(D, newObs)
## {
##     u1 <- D[newObs[1:4], -newObs][,1:(2*n)]
##     v2 <- D[newObs[5:8], -newObs][,(2*n+1):(4*n)]
##     ## c(1,3,2,4,1,3,2,4)
##     v1 <- u2 <- array(0, dim=c(8, ncol(v2)))
##     rownames(u1) <-  c("1.A","2.A","1.T","2.T")
##     rownames(v2) <- rownames(u1)
##     rownames(v1) <- c("11.A","11.T","22.A","22.T","12.A","12.T","21.A","21.T")
##     rownames(u2) <- rownames(v1)
##     ## 12.A: (a) the 1 refers to the 1st English doc
##     ##       (b) the 2 refers to the 2nd French doc
##     ##       (c) the A refers to rows corresponding to A (link)

##     ## pair <- c("11","22","12","21")
##     pair <- matrix(c(1,2,1,2,1,2,2,1), 4, 2)
##     col.AE <- col.AF <- 1:n
##     col.TE <- col.TF <- (n+1):(2*n)
##     for (k in 1:4) {
##         i.e <- pair[k,1]
##         i.f <- pair[k,2]
##         ru1 <- rv2 <- list(AE = paste(i.e,"A",sep="."),
##                            TE = paste(i.e,"T",sep="."),
##                            AF = paste(i.f,"A",sep="."),
##                            TF = paste(i.f,"T",sep="."))
##         rv1 <- ru2 <- list(A = paste(i.e,i.f,".A",sep=""),
##                            T = paste(i.e,i.f,".T",sep=""))

##         v1.AEAF <- (u1[ru1$AE, col.AE] + v2[rv2$AF, col.AF])/2
##         v1.AETF <- (u1[ru1$AE, col.AE] + v2[rv2$TF, col.TF])/2
##         v1.TEAF <- (u1[ru1$TE, col.TE] + v2[rv2$AF, col.AF])/2
##         v1.TETF <- (u1[ru1$TE, col.TE] + v2[rv2$TF, col.TF])/2

##         v1[rv1$A,] <- c(v1.AEAF, v1.AETF)
##         v1[rv1$T,] <- c(v1.TEAF, v1.TETF)
##         u2[ru2$A,] <- c(v1.AEAF, v1.TEAF)
##         u2[ru2$T,] <- c(v1.AETF, v1.TETF)
##     }
##     a2 <- rbind(cbind(u1[c("1.A","1.T","2.A","2.T","1.A","1.T","2.A","2.T"),],
##                       v1),
##                 cbind(u2,
##                       v2[c("1.A","1.T","2.A","2.T","2.A","2.T","1.A","1.T"),]))^2
##     rownames(a2) <- c(paste(rownames(v1),"E", sep=""),
##                        paste(rownames(u2), "F", sep=""))
##     ## 11.AE
##     ## 11.TE
##     ## 22.AE
##     ## 22.TE
##     ## 12.AE
##     ## 12.TE
##     ## 21.AE
##     ## 21.TE
##     ## 11.AF
##     ## 11.TF
##     ## 22.AF
##     ## 22.TF
##     ## 12.AF
##     ## 12.TF
##     ## 21.AF
##     ## 21.TF
##     a2
## }

## ##
## ## OOS.4x4
## ##
## OOS.4x4 <- function(X, D, newObs)
## {
##     a2 <- get.a2.4x4(D, newObs)
##     Delta2 <- D[-newObs,-newObs]^2
##     tau.e.Delta2 <- tau.e( Delta2 )
##     d <- ncol(X)
##     k <- nrow(a2)
##     Y <- matrix(0, k, d)
##     rownames(Y) <- rownames(a2)
##     set.seed(12345)
##     for (i in 1:k) {
##         B <- tau.w(tau.e.Delta2, Delta2, a2[i,])
##         y0 <- rnorm(d)*10
##         nlm.out <- nlm(f,p=y0,B=B,X=X)
##         Y[i,] <- matrix(nlm.out$estimate,ncol=d)
##     }
##     Y
## }

## ##
## ## given Y (output of OOS.4x4), return various (squared) distances
## ##
## dist.4x4 <- function(Y)
## {
##     Td <- c(sum((Y["11.AE",] - Y["11.AF",])^2),
##             sum((Y["22.AE",] - Y["22.AF",])^2),
##             sum((Y["12.AE",] - Y["12.AF",])^2),
##             sum((Y["21.AE",] - Y["21.AF",])^2),
##             sum((Y["11.TE",] - Y["11.TF",])^2),
##             sum((Y["22.TE",] - Y["22.TF",])^2),
##             sum((Y["12.TE",] - Y["12.TF",])^2),
##             sum((Y["21.TE",] - Y["21.TF",])^2))
##     Td
## }

## ##
## ##
## ## critical.values: given a distance vector, return critical values
## ## for different alpha \in ALPHA=seq(0, 1, by=0.05)
## ##
## critical.values <- function(d0, ALPHA=seq(0, 1, by=0.05))
## {
##     len <- length(ALPHA)
##     calpha <- rep(0, len)
##     for (k in 1:(len-1)) {
##         calpha[k] <- sort(d0)[round(length(d0)*(1-ALPHA[k]))]
##     }
##     calpha
## }


## ##
## ## LOO.critical.values: for each alpha, and each i \in {1, 2, ..., 200},
## ## return the critical values resulted from Td[-i,1:2]
## ##
## LOO.critical.values <- function(Td, ALPHA=seq(0, 1, by=0.05))
## {
##     nr <- nrow(Td)
##     nc <- length(ALPHA)
##     calpha <- matrix(0, nr, nc)
##     for (k in 1:(nc-1)) {
##         alpha <- ALPHA[k]
##         for (i in 1:nr) {
##             calpha[i,k] <- sort(as.vector(Td[-i,1:2]))[round((nr-1)*2*(1-alpha))]
##         }
##     }
##     calpha
## }


## ##
## ## decision.error: for each alpha, and each i \in {1, 2, ..., 200},
## ## compare Td[i,] to the critical value resulting from Td[-i,1:2], and
## ## record the decision error.
## ##

## decision.error <- function(d2)
## {
##     calpha <- LOO.critical.values(d2)
##     error <- array(0, dim(calpha))
##     for (i in 1:nrow(d2)) {
##         error[i,] <- (d2[i,1]  > calpha[i,]) + (d2[i,2]  > calpha[i,]) +
##                      (d2[i,3] <= calpha[i,]) + (d2[i,4] <= calpha[i,])
##     }
##     error
## }

## ##
## ## Find the rotation matrix R such that ||AR - B||_F is minimized.
## ## NOTE:
## ## 1. A and B is centered first.
## ## 2. return R, AR and a function f.
## ## 3. f takes as input a vector (or matrix) x which has same dimensions
## ##    as A and B, and outputs its procruste transformation xR (if from
## ##    =="A"), or x - colMeans(B) (if from == "B").
## ##
## ## procrustes <- function(A, B)
## ## {
## ##     A <- as.matrix(A)
## ##     B <- as.matrix(B)
## ##     A.colmean <- colMeans(A)
## ##     B.colmean <- colMeans(B)
## ##     A <- scale(A, scale=FALSE)
## ##     B <- scale(B, scale=FALSE)
## ##     PS <- t(A) %*% B
## ##     svd1 <- svd(PS)
## ##     R <- svd1$u %*% t(svd1$v)
## ##     f <- function(x, from="A") {
## ##         if (from == "A") {
## ##             x <- scale(as.matrix(x), center=A.colmean, scale=FALSE)
## ##             return (x %*% R)
## ##         } else if (from == "B") {
## ##             return (scale(as.matrix(x), center=B.colmean, scale=FALSE))
## ##         } else
## ##         stop ("from must be either A or B")
## ##     }
## ##     list(R = R, rot = A %*% R, f = f)
## ## }


##* functions needed for oos

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
  M <- M - matrix(colMeans(M), n, n, byrow=T)
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


ThreewayMDS<- function(D1.all.obj,D2.all.obj,model,ndim,oos)
{
	D1.base <-
	D2.base <-
	Threeway.Embed<- smacofIndDiff(delta=list(D1.all.obj,D2.all.obj), ndim = ndim, weightmat = NULL, init = NULL, metric = TRUE,
              ties = "primary", constraint = model, verbose = FALSE, modulus = 1,
              itmax = 1000, eps = 1e-6)
	Embed.info <-list(common.config = Threeway.Embed$gspace,conditional.config=Threeway.Embed$conf, cond.space.maps=Threeway.Embed$cweights)
}