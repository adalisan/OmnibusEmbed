
draw.plots<-function(sim.res,model,params,plot.w.vals,compare.pom.cca,cca.reg,compute.bound){
  
  #plot.title<-"Power Curves for Match Testing: JOFC vs Others"
  #
  # Plotting of ROC curves
  #
  sim.res.pow<-sim.res$power
  #w.val.len <-length(params$w.vals)
  w.val.len <-length(plot.w.vals)
  num.w.val<- length(params$w.vals)
  
  
  
  require(RColorBrewer)
  
  colors.vec.alt <- c("red","green","aquamarine","purple",
                      "darkblue","salmon","rosybrown","magenta","orange","gold4","darkblue","darkorange4")
  
  colors.vec.len <- w.val.len
  if (w.val.len<=8){
    colors.vec<- brewer.pal(w.val.len+3,"Spectral")
  } else if (w.val.len<=11){
    colors.vec<- brewer.pal(w.val.len,"Spectral")
    colors.vec[(w.val.len+1):(w.val.len+3)]<-  brewer.pal(3,"Set3")
  } else{
    colors.vec<- brewer.pal(11,"Spectral")
    colors.vec[12:(w.val.len+3)]<-  brewer.pal(length(12:(w.val.len+3),"Set3"))
  }
  
  
  # colors.vec <- colors.vec[-(1:3)]
  #colors.vec[colors.vec.len+1]<-"cornflowerblue"
  #colors.vec[colors.vec.len+2]<-"red"
  #colors.vec[colors.vec.len+3]<-"cyan"
  
  colors.vec.len <- length(colors.vec)
  
  #palette("YlOrRd")
  
  
  
  linewd<-3
  lty.i.vec<-c()
  par(cex=1.4)
#   par(cex.lab=1.8,cex.axis=1.8,cex.main=1.8,cex.sub=1.8)
  for (color.i in 1:w.val.len){
     i<- plot.w.vals[color.i]
      
      lty.i <- 1+((color.i-1)%%10)
      
      lty.i.vec <- c(lty.i.vec,lty.i)
      par(lty=lty.i)
      plot.ROC.with.CI(sim.res.pow[i,,],plot.title="",plot.col = colors.vec[color.i],
                       conf.int=FALSE,add=(color.i>1),ylim=1,linewd=linewd)
    
  }
  
  if (compare.pom.cca) {
    i<- 1+((w.val.len)%%10)
    par(lty=i)
    plot.ROC.with.CI(sim.res$power.cmp$pom,plot.title="",plot.col = colors.vec[w.val.len+1],
                     conf.int=FALSE,add=TRUE,ylim=1,linewd=linewd)
    lty.i.vec <- c(lty.i.vec,i)  		
    i<- 1+((w.val.len)%%10)
    par(lty=i)
    plot.ROC.with.CI(sim.res$power.cmp$cca,plot.title="",plot.col = colors.vec[w.val.len+2],
                     conf.int=FALSE,add=TRUE,ylim=1,linewd=linewd)
    lty.i.vec <- c(lty.i.vec,i)	
    if (cca.reg){
      i<- 1+((w.val.len+2)%%10)
      par(lty=i) 
      plot.ROC.with.CI(sim.res$power.cmp$reg.cca,plot.title="",plot.col = colors.vec[w.val.len+3],
                       conf.int=FALSE,add=TRUE,ylim=1,linewd=linewd)
      lty.i.vec <- c(lty.i.vec,i)	
      #lwd.i.vec <- c(lwd.i.vec,line.width)
    }
    
    
  }
  if (compute.bound) {
    i<- 1+((w.val.len+3)%%10)
    par(lty=3)
    lines(sim.res$optim.power,col="black")
    lty.i.vec <- c(lty.i.vec,i)
  }
  
  legend.txt <- (params$w.vals)[plot.w.vals]
  legend.txt.prefix <- c("w =",rep(c("   "),length(plot.w.vals)-1))
  legend.txt<- paste(legend.txt.prefix,legend.txt)
  if (compare.pom.cca)
    legend.txt <-c(legend.txt ,"pom","cca")
  if (compute.bound) legend.txt <-c(legend.txt ,"bound")
  if (cca.reg)
    legend.txt <-c(legend.txt ,"cca.reg")
  if (compute.bound) {legend.txt <-c(legend.txt ,"bound")}
  
  legend("bottomright",legend=legend.txt,
         col=colors.vec,#c(colors.vec[plot.w.vals],colors.vec[(w.val.len+1:3)]),
         lty=lty.i.vec,lwd=linewd) #,text.width=0.15
  #title(plot.title)
  par(lty=1)
  fname<- file.path('graphs',paste(c(model,"-FC-Tradeoff-",ifelse(oos,"OOS","noOOS"),"c",params$c.val),collapse=""))
  
  
  if(!run.in.linux&(!run.for.Sweave))  savePlot(paste("./",fname,".pdf",sep="",collapse=""),type="pdf")
  if(!run.in.linux&(!run.for.Sweave))  savePlot(filename=paste("./",fname,".ps",sep="",collapse=""),type="ps")
  if(!run.in.linux&(!run.for.Sweave)) savePlot(filename=paste("./",fname,".png",collapse="",sep=""),type="png")
  
  #dev.off()
  
  
  
  
 # par(ps=20)
  
  #
  #  Power vs w plotting
  
  #
  if( run.in.linux) X11() else {windows()}
  
  
  #temp
  par(lty=1)
  lwd.old <- par()$lwd
  par(lwd=3)
  par (cex=1.4)
  w.val.len <-length(params$w.vals)
  #
  # Power vs w plots for fixed alpha
  avg.power.w <-rep(0,length(w.val.len))
  for (i in 1:num.w.val){
    beta.w<-sim.res$power[i,,2]
    avg.power.w[i]<-mean(beta.w)		
  }
  
  
  #print(which.max(avg.power.w))
  #print(params$w.vals)	
  
  #w
  max.w <- which.max(avg.power.w)
  sim.res$wstar.estim <- params$w.vals [max.w]
  if (verbose) print("Estimate of wstar for average power curve")
  if (verbose) print(sim.res$wstar.estim)
  #sink("debug-wstar.txt")
  #print(sim.res$power[,,6])
  
  #Which w value was the best w for the highest number of mc replicates
  #Note that multiple w.values might have the best power
  sim.res$wstar.idx.estim.mc<- apply(sim.res$power[,,6],2,function(x) which(x==max(x)))
  
  print(sim.res$wstar.idx.estim.mc)
  #sink()
  avg.power.w.2 <-rep(0,length(w.val.len))
  for (i in 1:num.w.val){
    beta.w<-sim.res$power[i,,6]
    avg.power.w.2[i]<-mean(beta.w)
  }
  
  avg.power.w.3 <-rep(0,length(w.val.len))
  for (i in 1:num.w.val){
    beta.w<-sim.res$power[i,,11]
    avg.power.w.3[i]<-mean(beta.w)			
  }
  
  x.vals<- 1:num.w.val
  par(lwd=3)
  par(cex.lab=1.2)
  plot(x=x.vals,y=avg.power.w,col="blue",type="l",ylim=c(0,1),log='x',xlab=expression(w),ylab=expression(beta),xaxt="n",ps=10)
  lines(x=x.vals,y=avg.power.w.2,col="red",xlog=TRUE)
  lines(x=x.vals,y=avg.power.w.3,col="green",xlog=TRUE)
  par(cex.axis=0.9)
  axis(side=1,at=1:length(params$w.vals), labels = params$w.vals)
  legend("topright",expression(alpha==0.01,alpha==0.05,
                               alpha==0.1),col=c("blue","red","green"),lty=rep(1,3),lwd=3)
  #,adj=c(0.5,0)
  
  #title(paste(model," model: power vs w plot"))
  fname <- file.path(paste(results.dir,ifelse(oos,"OOS","noOOS"),model,"-power-w-c",params$c.val,sep="", collapse=""))
  if(!run.in.linux&(!run.for.Sweave))	savePlot(paste(fname,".pdf",sep="",collapse=""),"pdf")
  if(!run.in.linux&(!run.for.Sweave)) savePlot(paste(fname,".png",sep="",collapse=""),"png")
  
  if(!run.in.linux&(!run.for.Sweave)) savePlot(paste(fname,".ps",sep="",collapse=""),"ps")
  
  par(lwd=lwd.old)
  
  #dev.print(png,file=fname,height=600,width=600)
  #dev.off()
  
  
  
}

draw.auc.error.bar.argmax.bar <- function (w.vals) {
  require(ggplot2)
  require(reshape2)
  
  auc.meas<- read.csv("./cache/auc_meas.txt",header=TRUE,row.names=1)
  
  auc.avg <-apply(auc.meas,2,mean)
  
  
  w.val.len <- length(w.vals)
  #w.vals <-  c(0.1,0.4,0.5,0.8,0.85,0.9,0.925,0.95,0.99,0.999)
  auc.mat<- matrix(auc.avg,1,length(w.vals))
  
  colnames(auc.mat) <- params$w.vals # c(0.1,0.4,0.5,0.8,0.85,0.9,0.925,0.95,0.99,0.999)
  
  
  argmax.w<-apply(auc.meas,1,which.max)
  
  qplot(x=argmax.w, geom="histogram")
  
  
  argmax.w.fac<-factor(argmax.w,levels=1:length(w.vals))
  argmax.data<- data.frame(argmax.w=argmax.w)
  ggplot(argmax.data,aes(argmax.w.fac))+geom_bar(width=0.2)+scale_x_discrete(breaks=1:length(w.vals),labels=w.vals)+
    theme_minimal()+theme(text=element_text(size=32))+labs(x="w*")
  hist.plot<-ggplot(argmax.data,aes(argmax.w.fac))+geom_bar(width=0.2)+scale_x_discrete(breaks=1:length(w.vals),labels=w.vals)+
    theme_minimal()+theme(text=element_text(size=32))+labs(title="AUC",x="w*")
  
  ggsave(hist.plot,filename="./graphs/auc_argmax_hist.pdf")
  
  
  auc.meas.melt<-melt(auc.meas)
  corr.summ<-summarySE(data=auc.meas.melt,measurevar="value",groupvars=c("variable"))
  
  
  auc.errorbar.plot <- ggplot(corr.summ, aes(x=variable, y=value,group=1)) + 
    geom_line(aes(x=variable, y=value),size=1.2) +
    geom_errorbar(aes(ymin=value-ci, ymax=value+ci),size=1.2, width=.1) +
    geom_point(size=2)+
    scale_x_discrete(labels=w.vals)+
    theme_minimal()+theme(text=element_text(size=16)) +
    labs(title="AUC",x=expression(w),y=((expression(AUC(w)))))
  
  ggsave(auc.errorbar.plot,filename="./graphs/auc_errorbar_plot.pdf")
  
  save.image("./cache/argmax.RData")
}

