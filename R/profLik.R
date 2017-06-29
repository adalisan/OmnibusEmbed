#################################################################
#
# Main Functions for Profile Likelihood Dimension Selection
# by Mu Zhu
# First version = Jan 2005
# This version  = Jun 2007
#
#################################################################

plot1<-function(x, n, top)
{
  barplot(x, col=c(rep("black", top), rep("grey", n-top)), 
          border=c(rep("black", top), rep("grey", n-top)),
          xlab="Dimension",
          ylab="",
          main="Scree Plot",
          cex.main=1.75,
          cex.lab=1.75,
          cex.axis=1.5); 
}

plot2<-function(lik, top, legendx, legendy, preface)
{
  if (missing(legendx) || missing(legendy))
  {
    legendx<-top;
    legendy<-(max(lik)+min(lik))/2
  }
  if (missing(preface)) preface<-"q ="
  plot(1:length(lik),lik, type="l", 
       xlab="Dimension",
       ylab="",
       main="Profile Log-Likelihood",
       cex.main=1.75,
       cex.lab=1.75,
       cex.axis=1.5);
  abline(v=top, lty=2);
  if (preface != "none")
    text(legendx, legendy, paste(preface, top), cex=1.5)
}

scree.thresh <- function(x, trace=F, filename=NULL, ...)
{
  
  myvar <- function(x)
  {
    if (length(x) == 1) 
      return (0)
    else 
      return(var(x))
  }
  
  x <- sort(x, decr=T);
  n <- length(x);
  lik <- rep(0, n);
  for (i in 1:(n-1))
  {
    u <- mean(x[1:i]); s1 <- myvar(x[1:i]);
    v <- mean(x[(i+1):n]); s2 <- myvar(x[(i+1):n]);
    s <- sqrt(((i-1)*s1+(n-i-1)*s2)/(n-2))
    lik[i] <- sum(dnorm(x[1:i], u, s, log=T))+sum(dnorm(x[(i+1):n], v, s, log=T)) 
  }
  
  u <- mean(x); s <- sqrt(var(x));
  lik[n] <- sum(dnorm(x, u, s, log=T)) 
  
  top<-which(lik==max(lik));
  
  if (is.null(filename)) 
  {
    par(mfrow=c(2,2))
    plot1(x, n, top)
    plot2(lik, top, ...)
  }
  else 
  {
    postscript(file=paste(filename, "-bar", ".ps", sep=""), horiz=F, height=8, width=8)
    plot1(x, n, top)
    dev.off()
    postscript(file=paste(filename, "-pro", ".ps", sep=""), horiz=F, height=8, width=8)
    plot2(lik, top, ...)
    dev.off()
  }
  
  if (!trace) return(top)
  else return(list(top=top, lik=lik))
}