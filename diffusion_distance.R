
# Return the transition matrix given a similarity measure
transition.matrix <- function(W, dissimilarity = TRUE){
if(dissimilarity == TRUE)
  W[W != 0] <- 1/W[W != 0]
  
S <- apply(W,1,sum)
P <- W/S

P[is.nan(P)]<-0
P[is.na(P)]<-0
return(P)
}

diffusion.distance <- function(P, t, directed = FALSE){
  
  n <- nrow(P)
  w <- stationary.dist(P)
  Q <- outer(seq(1,1,length.out = n), w)

  if(directed == FALSE){
      P.2t <- mtx.exp(P,2*t)
      D.squared <- kappa(P.2t * (1/Q))
  }

  if(directed == TRUE){
      P.t <- mtx.exp(P,t)
      D.squared <- kappa(P.t %*% diag(1/w) %*% t(P.t))
  }
  
  D <- sqrt(D.squared)

  return(D)
}

stationary.dist <- function(P){
    
    n <- nrow(P)
    decomposition <- eigen(t(P),symmetric=FALSE)
    w <- abs(decomposition$vectors[,1])
    w <- w/sum(w)
#   Q <- outer(seq(1,1,length.out=n), w)
    return(w)
}

diffusion.map <- function(P, t, dim=2){

    n <- nrow(P)
    w <- stationary.dist(P)
    
    Q <- outer(seq(1,1,length.out = n), w) 
    Q.sqrt <- sqrt(Q)
    
    A <- t(Q.sqrt) * P * (1/Q.sqrt)
    
    decomp <- eigen(A)
    eigen.vals <- decomp$values[2:(dim+1)]
    eigen.vects <- decomp$vectors[,2:(dim+1)]
    eigen.vals.powered <- eigen.vals^t

    Psi <- outer(1/sqrt(w), seq(1,1,length.out = dim)) * eigen.vects * outer(seq(1,1,length.out = n), eigen.vals.powered)
    return(Psi)
}

kappa <- function(X){
    n <- nrow(X)
    D <- outer(diag(X),seq(1,1,length.out = n)) - X - t(X) + outer(seq(1,1,length.out = n), diag(X))
    return(D)
}

# Compute the exponent of matrix
mtx.exp <- function(X,n){
    if(n != round(n)){
        n <- round(n)
        warning("rounding exponent `n` to", n)
    }
    phi <- diag(nrow = nrow(X))
    pot <- X
    while(n > 0){
        if(n %% 2)
            phi <- phi %*% pot
        n <- n %/% 2
        pot <- pot %*% pot
    }
    return(phi)
}

