
require(gdata)
require(Hmisc)
require(rgl)
require(animation)
if (!exists("smacofM")) source("./lib/smacofM.R")
if (!exists("oosIM")) source("./lib/oosIM.R")
verbose <- FALSE
run.in.linux<-FALSE
results.dir<-"./graphs"
plot.in.3d<-TRUE
create.ani <-FALSE

fps=10

ani.options(outdir=file.path(Sys.getenv("PROJECT_DIR"),"graphs"))

meshgrid <- function(a,b) {
  list(
    x=outer(b*0,a,FUN="+"),
    y=outer(b,a*0,FUN="+")
  )
}

raw.stress.at <-function(config){
  sum(as.dist(W)*((dist(config)-d.X)^2))
}

hessian.mat <- function (X.embed.2.norm,n) {
  hess.mat.size<-4*n^2
  hess.mat<-matrix(0,2*n,2*n)
  for (i.r in 1:n){
    for (i.d in 1:2){
      hess.idx.r<-(2*(i.r-1)+i.d)
      X.embed.2.norm.minus <- X.embed.2.norm.plus<-X.embed.2.norm
      X.embed.2.norm.plus[i.r,i.d] <- X.embed.2.norm.plus[i.r,i.d] + epsilon
      X.embed.2.norm.minus[i.r,i.d] <- X.embed.2.norm.minus[i.r,i.d] - epsilon        
      dir.deriv[w.i,i,j,]<-(raw.stress.at(X.embed.2.norm.plus)- raw.stress.at(X.embed.2.norm.minus))/(2*epsilon)
      for (j.r in i.r:n){
        for (j.d in 1:2){
          hess.idx.c<-(2*(j.r-1)+j.d)
          X.embed.2.norm.back.forw<-X.embed.2.norm.forw.forw <- X.embed.2.norm.forw.back <-  X.embed.2.norm.back.back <-X.embed.2.norm
          X.embed.2.norm.forw.forw[i.r,i.d] <- X.embed.2.norm.forw.forw[i.r,i.d] + epsilon
          X.embed.2.norm.forw.forw[j.r,j.d] <- X.embed.2.norm.forw.forw[j.r,j.d] + epsilon
          
          X.embed.2.norm.forw.back[i.r,i.d] <- X.embed.2.norm.forw.back[i.r,i.d] + epsilon
          X.embed.2.norm.forw.back[j.r,j.d] <- X.embed.2.norm.forw.back[j.r,j.d] - epsilon
          
          X.embed.2.norm.back.back[i.r,i.d] <- X.embed.2.norm.back.back[i.r,i.d] - epsilon
          X.embed.2.norm.back.back[j.r,j.d] <- X.embed.2.norm.back.back[j.r,j.d] - epsilon
          
          X.embed.2.norm.back.forw[i.r,i.d] <- X.embed.2.norm.back.forw[i.r,i.d] - epsilon
          X.embed.2.norm.back.forw[j.r,j.d] <- X.embed.2.norm.back.forw[j.r,j.d] + epsilon
          
          # Approximate entry of hessian matrix by finite difference
          hess.mat[hess.idx.r,hess.idx.c]<- (raw.stress.at(X.embed.2.norm.forw.forw) - raw.stress.at(X.embed.2.norm.forw.back)
                                             - raw.stress.at(X.embed.2.norm.back.forw) + raw.stress.at(X.embed.2.norm.back.back))/(4*epsilon*epsilon)
        }
      }
      
    }
  }
  return(hess.mat)
}

stress.plot3d <- function (time,sign.hessian.at.pt,x.coords,y.coords,
                           stress.at.loc.w,grid.seq.x,grid.seq.y,w.vals,
                           rotate.z.angle=0) {
  w.i <- min(floor(time/1.5)+1,length(w.vals))
  print('time and w.i')
  print(time)
  print(w.i)
  
  
  col.matrix = sign.hessian.at.pt[w.i,,]
  col.matrix[sign.hessian.at.pt[w.i,,]==1]="orange"
  col.matrix[sign.hessian.at.pt[w.i,,]==2]="red"
  col.matrix[sign.hessian.at.pt[w.i,,]==0]="green"
  col.matrix[sign.hessian.at.pt[w.i,,]==-1]="blue"
  col.matrix[sign.hessian.at.pt[w.i,,]==-2]="pink"
  clear3d(type="shapes")
  
  plot3d(x=x.coords, y=y.coords,
         #plot3d(unmatrix(mesh.grid.coords$x,byrow=FALSE),
         #   unmatrix(mesh.grid.coords$y,byrow=FALSE),[]
         z=stress.at.loc.w[w.i,,],col=unmatrix(col.matrix),byrow=FALSE,
         xlab="x",ylab="y",zlab="Stress" ,box=FALSE,axes=FALSE,
         xlim=c(min(x.coords),max(x.coords)),ylim=c(min(y.coords),max(y.coords)),
         zlim=c(0,max(stress.at.loc.w)))
  decorate3d(xlim=c(min(x.coords),max(x.coords)),ylim=c(min(y.coords),max(y.coords)),
             zlim=c(0,max(stress.at.loc.w)))
  surface3d(grid.seq.x,grid.seq.y,stress.at.loc.w[w.i,,])
  title3d()
  title3d(paste("Stress w=",eval(expression(w.vals[w.i]))),col="red")
  if (time==0 || w.i==1 && time<(2.0/fps)){
    new_mat=transform3d(par3d("userMatrix"),angle=rotate.z.angle,x=0,y=0,z=1)
  } else{
    new_mat=par3d("userMatrix")
  }
  return(list(userMatrix=new_mat ))
}


animate.config.w <- function () {
  
  oopt = ani.options(interval = 1, nmax = length(w.vals))
  ## use a loop to create images one by one
  for (w.i in 1:length(w.vals)) {
    par(pch=1)
    if (is.null(far.to.init.X5.for.w[[w.i]]) || nrow(far.to.init.X5.for.w[[w.i]])>0)
      plot(x=far.to.init.X5.for.w[[w.i]][,1],y=far.to.init.X5.for.w[[w.i]][,2],col="red",
           xlim=c(min(grid.seq.x),max(grid.seq.x))
           ,ylim=c(min(grid.seq.y),max(grid.seq.y)),xlab="x",ylab="y")
    else{
      
      plot.new()
      plot.window(xlim=c(min(grid.seq.x),max(grid.seq.x))
                  ,ylim=c(min(grid.seq.y),max(grid.seq.y)))
      
      axis(1)
      axis(2)
      
      title(xlab="x",ylab="y")
      
      box()
    }
    
    
    par(pch=3)
    if (is.null(close.to.init.X5.for.w[[w.i]]) || nrow(close.to.init.X5.for.w[[w.i]])>0)
      points(x=close.to.init.X5.for.w[[w.i]][,1],y=close.to.init.X5.for.w[[w.i]][,2],col="red",
             xlim=c(min(grid.seq.x),max(grid.seq.x))
             ,ylim=c(min(grid.seq.y),max(grid.seq.y)),xlab="x",ylab="y")
    par(pch=1)
    if (is.null(far.to.init.X6.for.w[[w.i]]) || (nrow(far.to.init.X6.for.w[[w.i]])>0))
      points(x=far.to.init.X6.for.w[[w.i]][,1],y=far.to.init.X6.for.w[[w.i]][,2],col="blue",
             xlim=c(min(grid.seq.x),max(grid.seq.x))
             ,ylim=c(min(grid.seq.y),max(grid.seq.y)),xlab="x",ylab="y")
    par(pch=3)
    if (is.null(close.to.init.X6.for.w[[w.i]]) || (nrow(close.to.init.X6.for.w[[w.i]])>0))
      points(x=close.to.init.X6.for.w[[w.i]][,1],y=close.to.init.X6.for.w[[w.i]][,2],col="blue",
             xlim=c(min(grid.seq.x),max(grid.seq.x))
             ,ylim=c(min(grid.seq.y),max(grid.seq.y)),xlab="x",ylab="y")
    title(paste("Final point config. of X_5 and X_6 for w=",eval(expression(w.vals[w.i]))))
    
    ani.pause() ## pause for a while ('interval')
  }
  ## restore the options
  ani.options(oopt)
  
}


animate.final.loc.stress.w <- function () {
  
  oopt = ani.options(interval = 1, nmax = length(w.vals))
  ## use a loop to create images one by one
  for (w.i in 1:length(w.vals)) {
    par(pch=1)
    stress.vals<-as.vector(stress.at.loc.w[w.i,,])
    
    
    if (((max(stress.vals))-min(stress.vals))>5E-3){
      value.intervals <- hist(stress.vals,breaks=3,
                              plot=FALSE)
      print(levels(value.intervals))
      pt.at.level.1 <- (stress.vals<=value.intervals$breaks[2])
      pt.at.level.2 <- (stress.vals>value.intervals$breaks[2]) & (stress.vals<=value.intervals$breaks[3])
      pt.at.level.3 <- (stress.vals>value.intervals$breaks[3])
      
    }else{
      #value.intervals==levels(value.intervals)[1]
      pt.at.level.1 = 1:length(stress.vals)
      pt.at.level.2 = c()
      pt.at.level.3 = c()
      
    }
    
    plot  (x=final.coords.x.5.w[w.i,,][ pt.at.level.1],
           y=final.coords.y.5.w[w.i,,][ pt.at.level.1],col="red",
           xlim=c(min(grid.seq.x),max(grid.seq.x)),
           ylim=c(min(grid.seq.y),max(grid.seq.y)),
           xlab="x",ylab="y"
    )
    points(x=final.coords.x.5.w[w.i,,][ pt.at.level.2],
           y=final.coords.y.5.w[w.i,,][ pt.at.level.2],col="purple")
    
    points  (x=final.coords.x.5.w[w.i,,][ pt.at.level.3],
             y=final.coords.y.5.w[w.i,,][ pt.at.level.3],col="blue")        
    
    par(pch=3)
    points  (x=final.coords.x.6.w[w.i,,][ pt.at.level.1],
             y=final.coords.y.6.w[w.i,,][ pt.at.level.1],col="red",
             xlim=c(min(grid.seq.x),max(grid.seq.x)),
             ylim=c(min(grid.seq.y),max(grid.seq.y))
    )
    points(x=final.coords.x.6.w[w.i,,][ pt.at.level.2],
           y=final.coords.y.6.w[w.i,,][ pt.at.level.2],col="purple")
    
    points  (x=final.coords.x.6.w[w.i,,][ pt.at.level.3],
             y=final.coords.y.6.w[w.i,,][ pt.at.level.3],col="blue")
    X[5,]<-c(1,0)
    X[6,]<-c(0,1)
    points(x=c(1,0),y=c(0,1),col="black",pch=c(1,3))
    
    legends.txt<-c("Lowest Stress","Medium","Highest","True Coords","X_5","X_6")
    legend(legend=legends.txt,x="topright",col=c("red","purple","blue","black","red","red"),
           pch=c(1,1,1,1,1,3))
    
    ## restore the options
    title(paste("Final point config. of X_5 and X_6 for w=",eval(expression(w.vals[w.i]))))
    
    title(eval(expression(w.vals[w.i])))
    ani.pause() ## pause for a while ('interval')
  }
  
  
  ani.options(oopt)
  
}


