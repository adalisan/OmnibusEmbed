simulate.generate.test.model.plot<-function(model,params,run.parallel.sf){
	wstar.List<-c()
	wstar.estim<-0
	sim.result<-with(params,{
		#		if ((model=="MVN") && run.parallel)              call.func<-gaussian_simulation_jofc_tradeoff_par
				if ((model=="MVN") && !run.parallel.sf)          call.func<-gaussian_simulation_jofc_tradeoff
			#	if ((model=="Dirichlet") && run.parallel)        call.func<-dirichlet_simulation_jofc_tradeoff_par
				if ((model=="Dirichlet") && !run.parallel.sf)    call.func<-dirichlet_simulation_jofc_tradeoff
        if ((model=="MVN") && run.parallel.sf)           call.func<-gaussian_simulation_jofc_tradeoff_sf
        if ((model=="Dirichlet") && run.parallel.sf)     call.func<-dirichlet_simulation_jofc_tradeoff_sf
				sim.res<-list()
				print("c")
				print(c.val)
				print(p)
				print("w values")
				print(w.vals)
        if (model=="MVN") real.dim<- p+q
            if (model=="Dirichlet") real.dim <- p+q+2
				begin.time.g <-Sys.time()
				args.for.func.call<-list(p=p, r=r, q=q, c.val=c.val,d=d,
          
						Wchoice     = "avg", 
						pre.scaling = TRUE,
						oos         = oos,
						alpha       = NULL,
						n = n, m = s, nmc = nmc,
					
						old.gauss.model.param=old.gauss.model,
						separability.entries.w=separability.entries.w,
						compare.pom.cca=compare.pom.cca,
						oos.use.imputed=oos.use.imputed,
						
						
						rival.w=rival.w,
						proc.dilation=FALSE,
						assume.matched.for.oos =assume.matched.for.oos,
						w.vals=w.vals,
						wt.equalize=wt.equalize,
						verbose=verbose,  power.comparison.test=power.comparison.test)
        
        
        if (!run.parallel.sf)
        args.for.func.call<- c(args.for.func.call,list(pprime1= real.dim, pprime2= real.dim))
        
        ############################################
        #  Running simulation function 
        ############################################
				sim.res <- do.call(call.func,args=args.for.func.call)
        
        
        
				if (verbose) print("Struct of sim.res")
				if (verbose)  print(str(sim.res))
				print("Simulation completed.Starting Plotting functions")
				#sim.res<- call.func(p=p, r=r, q=q, c.val=c.val,d=d,
				#Wchoice     = "avg", 
				#pre.scaling = TRUE,
				#oos         = oos,
				#alpha       = NULL,
				#n = n, m = s, nmc = nmc,
				#sim.grass=grassmannian.dist,
				#old.gauss.model.param=old.gauss.model,
				#separability.entries.w=separability.entries.w,
				#compare.pom.cca=compare.pom.cca,
				#oos.use.imputed=oos.use.imputed,
				#
				#
				#rival.w=rival.w,
				#proc.dilation=FALSE,
				#assume.matched.for.oos =assume.matched.for.oos,
				#w.vals=w.vals,verbose=verbose)
				try({	
                w.val.len<- length(w.vals)
							run.time.g <-Sys.time()-begin.time.g
							print("Total sim time:")
							print(run.time.g)
							print("Sim Time per MC replicate")
							print(run.time.g/nmc)
							
							
							sim.res.pow <- sim.res$power
							
							print(dim(sim.res.pow))
							
							#
							# Sign.rank.sum and sign.tests for  test of power comparison
							#
							
							#Binom.test logging start
							sink(paste(results.dir,model.letter,"c",params$c.val,"binom-tests.txt",collapse=""))
							sign.test <- try(sign.test.cont.table(sim.res$conting.table.list))
							
							if (inherits(sign.test,"try-error")) {
								print(paste("error in ",model,collapse=""))
								print(sim.res$conting.table.list)
							}
							else{
							print("Cont Table List: sign test p-val")
							print(sign.test$p.value)
							}
		
		
							sign.rank.sum.test<-sign.rank.sum.test.cont.table(sim.res$conting.table.list)
							print("Cont Table List: signed rank sum test p-val")
							print(sign.rank.sum.test$p.value)
							sink()
							#Binom.test logging end
							
							
							#McNemar's test for comparing equal weighted JOFC with weighted JOFC
							#
							# McNemar's test for comparing classifiers(corresp. to different JOFC weights)
							#
							if (run.mcnemars.test) {
								
								
								#mcnemars.test logging start
							sink(paste(model,c.val,"McNemars.test.R",collapse=""))
							avg.cont.table<- (sim.res$conting.table+0.001)/nmc
							print("Aggregate Cont Table: ")
							print(avg.cont.table)
							p.val <-mcnemar.test(avg.cont.table)$p.value
							print("Aggregate Cont Table: McNemar's p-val")
							print(p.val)
							
							sim.res$conting.table.list <- lapply(sim.res$conting.table.list,function(x) (x)+0.001)
							p.vals.list.two.sided <- lapply(sim.res$conting.table.list,mcnemar.test)
							p.vals.list <- lapply(sim.res$conting.table.list,exact2x2 ,
												paired=TRUE, alternative="greater",tsmethod="central")
							p.vals<-c()
							for (t in p.vals.list)
								p.vals<- c(p.vals,t$p.value)
							p.vals <- sort(p.vals)
							print("MC Replicate Cont Tables: McNemar's p-val")
							print(p.vals)
							sink()
							#mcnemars.test logging end
							
							if( run.in.linux) X11() else {windows()}
							hist(p.vals,15,main="Histogram of p-values")
							
							p.vals.two.sided<-c()
							for (t in p.vals.list)
								p.vals.two.sided<- c(p.vals.two.sided,t$p.value)
							p.vals.two.sided <- sort(p.vals.two.sided)
							
							#two-sided mcnemars.test logging end
							sink(paste(model,c.val,"McNemars.two.sided.test.R",collapse=""))
							print("MC Replicate Cont Tables: McNemar's p-val")
							print(p.vals.two.sided)
							sink()
							#two-sided mcnemars.test logging end

							if( run.in.linux) X11() else {windows()}
							hist(p.vals.two.sided,15,main="Histogram of p-values")
							
							
							
							#
							if(!run.in.linux&(!run.for.Sweave))	savePlot(paste(model,c.val,"p-values-McNemars-hist",".pdf",collapse="",sep=""),"pdf")
							if(!run.in.linux&(!run.for.Sweave)) savePlot(paste(model,c.val,"p-values-McNemars-hist",".png",collapse="",sep=""),"png")
							
							if(!run.in.linux&(!run.for.Sweave)) savePlot(paste(model,c.val,"p-values-McNemars-hist",".ps",sep="",collapse=""),"ps")
						
              
              
              
              if( run.in.linux) X11() else {windows()}
							
							 p.vals.density<-density(p.vals,from=0,to=1)
							 plot.density(p.vals.density,main="Kernel Density Estimate of p-values")
							
							#
							if(!run.in.linux&(!run.for.Sweave))	savePlot(paste(model,c.val,"p-values-McNemars-kde",".pdf",collapse="",sep=""),"pdf")
							if(!run.in.linux&(!run.for.Sweave)) savePlot(paste(model,c.val,"p-values-McNemars-kde",".png",collapse="",sep=""),"png")
							
							if(!run.in.linux&(!run.for.Sweave)) savePlot(paste(model,c.val,"p-values-McNemars-kde",".ps",collapse="",sep=""),"ps")
							
              
              
              if( run.in.linux) X11() else {windows()}
							
#	
# Error checking :Convergence of JOFC config to PoM config as w->0
#
						#	plot.graph.with.CI(sim.res$config.dist[,,1],plot.title="Config. Mismatch of wMDS with PoM",plot.col="black",conf.int=FALSE,fp.points=1:length(w.vals))
							
							
							if( run.in.linux) X11() else {windows()}
							
						
						#	plot.graph.with.CI(sim.res$min.stress[,1:w.val.len],plot.title="Minimum Stress",plot.col="black",conf.int=FALSE,fp.points=1:length(w.vals))
						#	lines(x=1:w.val.len, y=rep(mean(sim.res$min.stress[,w.val.len+1]),w.val.len),col="red")
							
							if( run.in.linux) X11() else {windows()}
							
#
# Plotting of ROC curves
#
							
							
							lty.i.vec<-c()
							for (i in 1:w.val.len){
								lty.i <- 1+((i-1)%%10)
								
								lty.i.vec <- c(lty.i.vec,lty.i)
								par(lty=lty.i)
								plot.ROC.with.CI(sim.res.pow[i,,],plot.title="",plot.col = colors.vec[i],
										conf.int=FALSE,add=(i>1),ylim=1)
								
							}
							
							if (compare.pom.cca) {
								i<- 1+((w.val.len)%%10)
								par(lty=i)
								plot.ROC.with.CI(sim.res$power.cmp$pom,plot.title="",plot.col = colors.vec[w.val.len+1],
										conf.int=FALSE,add=TRUE,ylim=1)
								lty.i.vec <- c(lty.i.vec,i)			
								i<- 1+((w.val.len+1)%%10)
								par(lty=i)
								plot.ROC.with.CI(sim.res$power.cmp$cca,plot.title="",plot.col = colors.vec[w.val.len+2],
										conf.int=FALSE,add=TRUE,ylim=1)
								lty.i.vec <- c(lty.i.vec,i)	
							}
							i<- 1+((w.val.len+3)%%10)
							par(lty=i)
						    lines(sim.res$optim.power,col=colors.vec[w.val.len+3])
							lty.i.vec <- c(lty.i.vec,i)	
						
							legend.txt <- w.vals
							if (compare.pom.cca)
								legend.txt <-c(legend.txt ,"pom","cca")
							legend.txt <-c(legend.txt ,"pom","cca","bound")
							legend("bottomright",legend=legend.txt,
									col=colors.vec,lty=lty.i.vec)
							title(plot.title)
							par(lty=1)
							fname<- file.path('graphs',paste(c(model,"-FC-Tradeoff-",ifelse(oos,"OOS","noOOS"),"c",c.val),collapse=""))
							
              
              if(!run.in.linux&(!run.for.Sweave))  savePlot(paste(fname,".pdf",sep="",collapse=""),type="pdf")
							if(!run.in.linux&(!run.for.Sweave))  savePlot(filename=paste(fname,".ps",sep="",collapse=""),type="ps")
							if(!run.in.linux&(!run.for.Sweave)) savePlot(filename=paste(fname,".png",collapse="",sep=""),type="png")
							
							dev.off()
	
							}

					
							if (verbose) print("F's and C")
							Fid1<-sim.res$FidComm.Terms$F1
							Fid2<-sim.res$FidComm.Terms$F2
							Comm<-sim.res$FidComm.Terms$C
							
              if (verbose) print(c("Fid1","Fid2","Comm"))
              if (verbose) print(cbind(Fid1,Fid2,Comm))
						if( run.in.linux) X11() else {windows()}
							
							plot.graph.with.CI(Fid1,"Fid and Comm Terms","red",conf.int=TRUE,   add=FALSE,fp.points=w.vals,
									customx.labels=NULL,customy.labels=NULL,ispowercurve=FALSE)
								traceback()
							plot.graph.with.CI(Fid2,"Fid and Comm Terms","orange",conf.int=TRUE,add=TRUE,fp.points=w.vals,
									customx.labels=NULL,customy.labels=NULL,ispowercurve=FALSE)
								traceback()
							plot.graph.with.CI(Comm,"Fid and Comm Terms","blue",conf.int=TRUE,  add=TRUE,fp.points=w.vals,
									customx.labels=NULL,customy.labels=NULL,ispowercurve=FALSE)
							traceback()
							fname<- file.path('graphs',paste(c(model,"-FidCommTerms-n",n,"-",ifelse(oos,"OOS","noOOS"),"c",c.val),collapse="",sep=""))
							if(!run.in.linux&(!run.for.Sweave))  savePlot(paste(fname,".pdf",sep="",collapse=""),type="pdf")
							if(!run.in.linux&(!run.for.Sweave))  savePlot(filename=paste(fname,".ps",sep="",collapse=""),type="ps")
							if(!run.in.linux&(!run.for.Sweave)) savePlot(filename=paste(fname,".png",collapse="",sep=""),type="png")
														
                
                
							if( run.in.linux) X11() else {windows()}
							plot.default(x=params$w.vals,y=rep(0,length(params$w.vals)),ylim=c(0,1.5))
							for (n.plot in 1:nmc){
							
								points(x=params$w.vals,Fid1[n.plot,],main="Fid and Comm Terms",col="red",   )
								traceback()
							points(x=params$w.vals,Fid2[n.plot,],main="Fid and Comm Terms",col="orange"	)
								traceback()
							points(x=params$w.vals,Comm[n.plot,],main="Fid and Comm Terms",col="blue")
							traceback()
							}
							
							
							fname<- file.path('graphs',paste(c(model,"-points-FidCommTerms-n",n,"-",ifelse(oos,"OOS","noOOS"),"c",c.val),collapse="",sep=""))
							if(!run.in.linux&(!run.for.Sweave))  savePlot(paste(fname,".pdf",sep="",collapse=""),type="pdf")
							if(!run.in.linux&(!run.for.Sweave))  savePlot(filename=paste(fname,".ps",sep="",collapse=""),type="ps")
							if(!run.in.linux&(!run.for.Sweave)) savePlot(filename=paste(fname,".png",collapse="",sep=""),type="png")
							
							
							avgFid1.n<-colMeans(Fid1)
							avgFid2.n<-colMeans(Fid2)
							avgComm.n<-colMeans(Comm)
							
							print("FC.ratio.n for MC replicates")
							if (verbose) print(sim.res$FC.ratios$f.c)
							print("wtFC.ratio.n  for MC replicates")
							if (verbose) print(sim.res$FC.ratios$wtf.c)
							
							FC.ratio.n  <-colMeans(sim.res$FC.ratios$f.c)
							wtFC.ratio.n<-colMeans(sim.res$FC.ratios$wtf.c)
							
							print("FC.ratio: average of MC replicates")
							if (verbose) print(FC.ratio.n)
							print("wtFC.ratio.n: average of MC replicates")
							if (verbose) print(wtFC.ratio.n)
							
							
							avg.Fid1<-rbind(avg.Fid1,avgFid1.n)
							avg.Fid2<-rbind(avg.Fid2,avgFid2.n)
							avg.Comm<-rbind(avg.Comm,avgComm.n)
							
							
							avg.FCratios.List<-rbind(avg.FCratios.List,FC.ratio.n)
							avg.wtFCratios.List<-rbind(avg.wtFCratios.List,wtFC.ratio.n)
							#FCRatios.List<-c(FCRatios.List,list(sim.res$FC.ratios$f.c))
							#wtFCRatios.List<-c(wtFCRatios.List,list(sim.res$FC.ratios$wtf.c))
							
							
							
#
#  Power vs w plotting

#
							  if( run.in.linux) X11() else {windows()}

							
				#temp
							par(lty=1)
							lwd.old <- par()$lwd
							par(lwd=3)
							#
							# Power vs w plots for fixed alpha
							avg.power.w <-rep(0,length(w.val.len))
							for (i in 1:w.val.len){
								beta.w<-sim.res$power[i,,2]
								avg.power.w[i]<-mean(beta.w)		
							}
						
						
						print(which.max(avg.power.w))
						print(params$w.vals)	
						
						max.w <- which.max(avg.power.w)
						sim.res$wstar.estim <- params$w.vals [max.w]
						print(wstar.estim)
						#wstar.List<-c(wstar.List,wstar.estim)
						
							avg.power.w.2 <-rep(0,length(w.val.len))
							for (i in 1:w.val.len){
								beta.w<-sim.res$power[i,,6]
								avg.power.w.2[i]<-mean(beta.w)
							}
							
							avg.power.w.3 <-rep(0,length(w.val.len))
							for (i in 1:w.val.len){
								beta.w<-sim.res$power[i,,11]
								avg.power.w.3[i]<-mean(beta.w)			
							}
							
							x.vals<- 1:w.val.len
							
							plot(x=x.vals,y=avg.power.w,col="blue",type="l",ylim=c(0,1),log='x',xlab=expression(w),ylab=expression(beta),xaxt="n",ps=10)
							lines(x=x.vals,y=avg.power.w.2,col="red",xlog=TRUE)
							lines(x=x.vals,y=avg.power.w.3,col="green",xlog=TRUE)
							par(cex.axis=0.9)
							axis(side=1,at=1:length(params$w.vals), labels = params$w.vals)
							legend("bottomright",expression(alpha==0.01,alpha==0.05,
											alpha==0.1),col=c("blue","red","green"),lty=rep(1,3))
						
										
							#title(paste(model," model: power vs w plot"))
							fname <- file.path(paste(results.dir,ifelse(oos,"OOS","noOOS"),model,"-power-w-c",c.val,sep="", collapse=""))
							if(!run.in.linux&(!run.for.Sweave))	savePlot(paste(fname,".pdf",sep="",collapse=""),"pdf")
							if(!run.in.linux&(!run.for.Sweave)) savePlot(paste(fname,".png",sep="",collapse=""),"png")
							
							if(!run.in.linux&(!run.for.Sweave)) savePlot(paste(fname,".ps",sep="",collapse=""),"ps")
							
							par(lwd=lwd.old)
							
							#dev.print(png,file=fname,height=600,width=600)
							dev.off()
						
						}
				) #end of try
				sim.res	
			} #end of argument of with()
	)	
	print("wstar.estim")
	print(sim.result$wstar.estim)
	#sim.result$wstar<- wstar.estim
	return(sim.result)
}
