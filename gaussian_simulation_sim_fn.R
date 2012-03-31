

gaussian_simulation_jofc_tradeoff <- function(p, r, q, c.val,
		d           = p-1,
		pprime1     = p+q,   # cca arguments
		pprime2     = p+q,   # cca arguments
		Wchoice     = "avg", 
		pre.scaling = TRUE,
		oos         = TRUE,
		alpha       = NULL,
		n = 100, m = 100, nmc = 100,
		
		old.gauss.model.param=FALSE,
		separability.entries.w,
		compare.pom.cca=TRUE,
		oos.use.imputed,
		level.mcnemar=0.01,
		def.w=0.5,
		rival.w,
		proc.dilation=FALSE,
		assume.matched.for.oos,
		w.vals,
		wt.equalize,
		verbose,
		power.comparison.test=TRUE) {
	## p: draw observations (signal) on Delta^p \in R^{p+1}
	## r: determine the matchedness between matched pairs --
	##    r = 0, independent(?); r = infinity, exact match
	## q: noise from Delta^q \in R^{q+1}
	## c.val: weight on noise
	##    X_{ik} ~ [(1-c.val)*signal, c.val*noise], i = 1, ..., n; k = 1, 2
	## d: dimension of projection space
	## alpha: Normal on Delta^p (NULL) or given
	##        if given, alpha is a list of length nmc, in which each
	##        component is a (n+2*m)x(p+1) matrix
	## Wchoice:
	##    "avg" => (D1 + D2)/2
	##    "sqrt" => sqrt((D1^2 + D2^2)/2)
	##    
	## n  : number of matched pairs
	## m  : number of testing pairs -- m matched pairs and m unmatched pairs
	## nmc: number of MC replicates
	##
	## Return: powers
	
	w.max.index <- length(w.vals)
	
	
	Fid.Terms.1<-c()
	Fid.Terms.2<-c()
	Comm.Terms <-c()
	
	F.to.C.ratio <-  c()
	wtF.to.C.ratio <- c()
	F.bar.to.C.bar.ratio <-c()
	
	
	size <- seq(0, 1, 0.01)
	len <- length(size)
	power <- array(0,dim=c(w.max.index,nmc,len))
	power.cmp<-list(pom= array(0,dim=c(nmc,len)), cca= array(0,dim=c(nmc,len)) )
	config.dist<- array(0,dim=c(nmc,w.max.index,3))
	
	agg.cont.table <- matrix(0,2,2)
	empty.cont.tab<- list(matrix(0,2,2))
	cont.tables<-rep(empty.cont.tab,nmc)
	min.stress<- array(0,dim=c(nmc,w.max.index+1))
	optim.power.vec <- array(0,dim=c(nmc,len))
	best.w.vals <-rep (0,nmc)
	
	for(mc in 1:nmc) {
		
		set.seed(mc)
		print(mc)
		
		mc.run<-try(run.mc.replicate("gaussian",p, r, q, c.val,
						d           = d,
						pprime1     = pprime1,   # cca arguments
						pprime2     = pprime2,   # cca arguments
						Wchoice     = Wchoice, 
						pre.scaling = pre.scaling,
						oos         = oos,
						alpha       = alpha,
						n = n, m = m,
						
						old.gauss.model.param=old.gauss.model.param,
						separability.entries.w=separability.entries.w,
						compare.pom.cca=compare.pom.cca,
						oos.use.imputed=oos.use.imputed,
						level.mcnemar=level.mcnemar,
						def.w=def.w,
						rival.w=rival.w,
						proc.dilation=proc.dilation,
						assume.matched.for.oos=assume.matched.for.oos,
						w.vals=w.vals,
						wt.equalize=wt.equalize,
						verbose=verbose,
						power.comparison.test=power.comparison.test))
		if (inherits(mc.run,"try-error")){
			print(paste("error in iter",mc,collapse="")	)
			power.mc= array(NA,dim=c(w.max.index,len))
			power.cca.mc = array(NA,dim=c(len))
			power.pom.mc = array(NA,dim=c(len))
			config.mismatch <-  list(frob.norm=array(NA,dim=c(w.max.index)))
			min.stress.mc = array(NA,dim=c(w.max.index+1))
			means <- array(NA , dim=c(w.max.index,2*d))			
			
			cont.table<-matrix(c(3,0,0,3),nrow=2,ncol=2)
			cont.tables[[mc]]<-(cont.table)
			
			if (compare.pom.cca) {
				power.cmp$cca[mc,]<-power.cca.mc
				power.cmp$pom[mc,]<-power.pom.mc
				config.dist[mc,,1]<-rep(NA,w.max.index)
				
				min.stress[mc,]   <-min.stress.mc
				
			}
			
			for (l in 1:w.max.index){
				power[l,mc, ] <- power.mc[l,]
			}
			
			
			next
			
		}
		cont.tables[[mc]]<-(mc.run$cont.tables)
		agg.cont.table <- agg.cont.table + cont.tables[[mc]]
		for (l in 1:w.max.index){
			power[l,mc, ] <- mc.run$power.mc[l,]
		}
		
		if (compare.pom.cca) {
			power.cmp$cca[mc,]<-mc.run$power.cmp$cca
			power.cmp$pom[mc,]<-mc.run$power.cmp$pom
			config.dist[mc,,1]<-mc.run$config.dist$frob.norm
			
			min.stress[mc,]   <-mc.run$min.stress
			
		}
		
		
		Fid.Terms.1 <- rbind(Fid.Terms.1, mc.run$FidComm.Terms$F1)
		Fid.Terms.2 <- rbind(Fid.Terms.2, mc.run$FidComm.Terms$F2)
		Comm.Terms <- rbind(Comm.Terms, mc.run$FidComm.Terms$C)
		F.to.C.ratio <-  rbind(F.to.C.ratio,mc.run$F.to.C.ratio)
		wtF.to.C.ratio <- rbind(wtF.to.C.ratio,mc.run$wtF.to.C.ratio)
		F.bar.to.C.bar.ratio <- rbind(F.bar.to.C.bar.ratio,mc.run$F.bar.to.C.bar.ratio)
		
		
		optim.power.vec[mc,]<- mc.run$optim.power
		best.w.vals[mc] <- mc.run$best.w
		
		
		
		
	}
	FC.terms<-list(F1=Fid.Terms.1, F2=Fid.Terms.2, C=Comm.Terms)
	FC.ratios<-list(f.c=F.to.C.ratio,wtf.c=wtF.to.C.ratio,f.c.bar=F.bar.to.C.bar.ratio)
	
	
	
	return (list(power=power,power.cmp=power.cmp, conting.table=agg.cont.table,conting.table.list=cont.tables,
					config.dist=config.dist,min.stress=min.stress,FidComm.Terms=FC.terms,
					FC.ratios=FC.ratios,optim.power =  optim.power.vec,wstar.estim= The.mode( best.w.vals) ))
	
	
}



gaussian_simulation_jofc_tradeoff_par <- function(p, r, q, c.val,
		d           = p-1,
		pprime1     = p+q,   # cca arguments
		pprime2     = p+q,   # cca arguments
		Wchoice     = "avg", 
		pre.scaling = TRUE,        
		oos         = TRUE,
		alpha       = NULL,
		n = 100, m = 100, nmc = 100,
		
		old.gauss.model.param=FALSE,
		separability.entries.w,
		compare.pom.cca=TRUE,
		oos.use.imputed,
		level.mcnemar=0.01,
		def.w=0.5,
		rival.w,		
		proc.dilation=FALSE,
		assume.matched.for.oos,
		w.vals,
		wt.equalize,
		verbose=FALSE,
		power.comparison.test=TRUE)  {
	## p: draw observations (signal) on Delta^p \in R^{p+1}
	## r: determine the matchedness between matched pairs --
	##    r = 0, independent(?); r = infinity, exact match
	## q: noise from Delta^q \in R^{q+1}
	## c.val: weight on noise
	##    X_{ik} ~ [(1-c.val)*signal, c.val*noise], i = 1, ..., n; k = 1, 2
	## d: dimension of projection space
	## alpha: Normal on Delta^p (NULL) or given
	##        if given, alpha is a list of length nmc, in which each
	##        component is a (n+2*m)x(p+1) matrix
	## Wchoice:
	##    "avg" => (D1 + D2)/2
	##    "sqrt" => sqrt((D1^2 + D2^2)/2)
	##    
	## n  : number of matched pairs
	## m  : number of testing pairs -- m matched pairs and m unmatched pairs
	## nmc: number of MC replicates
	##
	## Return: powers
	
	w.max.index <- length(w.vals)
	
	size <- seq(0, 1, 0.01)
	len <- length(size)
	power <- array(0,dim=c(w.max.index,nmc,len))
	power.cmp<-list(pom= array(0,dim=c(nmc,len)), cca= array(0,dim=c(nmc,len)) )
	optim.power <- array(0,dim=c(nmc,len))
	agg.cont.table <- matrix(0,2,2)
	empty.cont.tab<- list(matrix(0,2,2))
	cont.tables<-rep(empty.cont.tab,nmc)
	config.dist<- array(0,dim=c(nmc,w.max.index,3))
	min.stress<- array(0,dim=c(nmc,w.max.index+1))
	
	Fid.Terms.1<-c()
	Fid.Terms.2<-c()
	Comm.Terms <-c()
	
	F.to.C.ratio <-  c()
	wtF.to.C.ratio <- c()
	F.bar.to.C.bar.ratio <-c()
	
	seeds<-rep(list(),nmc)
	par.mc.result<- foreach(mc =1:nmc,.packages=c("MASS","MCMCpack","smacof") ) %dopar% { #,.export=c("power.comparison.test")
		
		source("./src/simulation_math_util_fn.R")
		
		source("./src/oosMDS.R")
		source("./src/smacofM.R")
		source("./src/oosIM.R")
		sink(file=file.path('logs',paste("debug-G-",mc,".txt",collapse="")))
		set.seed(mc)
		#if(mc<=4) {for(i in 1:mc) print(mvrnorm(4,mu=rep(0,4),Sigma=diag(4)))}
		print(runif(2))
		print(rnorm(2))
		
		print(mc)
		#seeds<-c(seeds,list(.Random.seed))
		
		seeds[[mc]]<-list(.Random.seed)
		tmp<- run.mc.replicate("gaussian",p, r, q, c.val,  ##try(
				d           = d,
				pprime1     = pprime1,   # cca arguments
				pprime2     = pprime2,   # cca arguments
				Wchoice     = Wchoice, 
				pre.scaling = pre.scaling,
				oos         = oos,
				alpha       = alpha,
				n = n, m = m,
				old.gauss.model.param=old.gauss.model.param,
				separability.entries.w=separability.entries.w,
				compare.pom.cca=compare.pom.cca,
				oos.use.imputed=oos.use.imputed,
				level.mcnemar=level.mcnemar,
				def.w=def.w,
				rival.w=rival.w,
				proc.dilation=proc.dilation,
				assume.matched.for.oos=assume.matched.for.oos,
				w.vals=w.vals,
				wt.equalize=wt.equalize,
				verbose=verbose,
				power.comparison.test=power.comparison.test) 
		#)
		if (verbose) sink(file=file.path('logs',paste("traceback-debug-G-",mc,".txt",collapse="")))
		if (verbose) traceback()
		if (verbose) sink()
		#If  simulation fails, generate a results list of NAs
		if (inherits(tmp,"try-error")) {
			print(paste("error in iter ",mc,collapse=""))
			
			power.mc= array(NA,dim=c(w.max.index,len))
			power.cca.mc = array(NA,dim=c(len))
			power.pom.mc = array(NA,dim=c(len))
			config.mismatch <-  list(frob.norm=array(NA,dim=c(w.max.index)))
			min.stress.mc = array(NA,dim=c(w.max.index+1))
			means <- array(NA , dim=c(w.max.index,2*d))			
			
			
			cont.table<-matrix(c(3,0,0,3),nrow=2,ncol=2)
			FC.terms <- list(F1=Fid.Terms.1, F2=Fid.Terms.2, C=Comm.Terms)
			
			
			tmp<- list(power.mc=power.mc,power.cmp=list(cca = power.cca.mc,pom = power.pom.mc), cont.tables=cont.table,
					config.dist= config.mismatch, min.stress=min.stress.mc,means=means,FidComm.Terms=FC.terms)
			
		}
		print(paste("Iteration complete:",mc,collapse=""))				
		sink()
		tmp
		
		
	}
	
	
	
	# Number of elements in the par.mc.result list for each mc replicate	
	num.value.per.mc.rep <-12
	if (verbose) print(str(par.mc.result))
	for (i in 1:nmc){
		mc.res.i<-par.mc.result[[i]]
		power[,i,]<-mc.res.i[[1]]
		#power.cmp is has values of 0  if compare.pom.cca is FALSE
		if (compare.pom.cca) {
			power.cmp$cca[i,] <-mc.res.i[[2]]$cca
			power.cmp$pom[i,] <-mc.res.i[[2]]$pom
		}
		cont.tables[[i]]<-mc.res.i[[3]]
		
		agg.cont.table<-agg.cont.table+cont.tables[[i]]
		print (paste(i,"th contingency table"))
		print(cont.tables[[i]])
		config.dist[i,,1]<- mc.res.i[[4]]$frob.norm
		
		min.stress[i,] <- mc.res.i[[5]]
#		means[i,]      <- par.mc.result[[(num.value.per.mc.rep*(i-1))+6]]
		Fid.Term.1.i <- mc.res.i[[7]]$F1
		Fid.Term.2.i <- mc.res.i[[7]]$F2
		Comm.Term.i <- mc.res.i[[7]]$C
		Fid.Terms.1 <- rbind(Fid.Terms.1,Fid.Term.1.i)
		Fid.Terms.2 <- rbind(Fid.Terms.2,Fid.Term.2.i)
		Comm.Terms <- rbind(Comm.Terms,Comm.Term.i)
		F.to.C.ratio.i <-mc.res.i[[9]]
		wtF.to.C.ratio.i <- mc.res.i[[10]]
		F.bar.to.C.bar.ratio.i <- mc.res.i[[11]]
		F.to.C.ratio <-  rbind(F.to.C.ratio,F.to.C.ratio.i)
		wtF.to.C.ratio <- rbind(wtF.to.C.ratio,wtF.to.C.ratio.i)
		F.bar.to.C.bar.ratio <- rbind(F.bar.to.C.bar.ratio,F.bar.to.C.bar.ratio.i)
		optim.power[i,]<-mc.res.i[[12]]
		
		
	}
	print(agg.cont.table)
	
	FC.terms <- list(F1=Fid.Terms.1, F2=Fid.Terms.2, C=Comm.Terms)
	#FC.sum.terms<-
	FC.ratios<-list(f.c=F.to.C.ratio,wtf.c=wtF.to.C.ratio,f.c.bar=F.bar.to.C.bar.ratio)
	if (verbose) print(str(FC.ratios))
	return (list(power=power,power.cmp = power.cmp, conting.table=agg.cont.table,conting.table.list=cont.tables,
					config.dist=config.dist,min.stress=min.stress,seeds=seeds, FidComm.Terms=FC.terms,
					FC.ratios=FC.ratios,optim.power=optim.power))
}




gaussian_simulation_jofc_tradeoff_sf <- function(p, r, q, c.val,
		d           = p-1,
		pprime1     = p+q,   # cca arguments
		pprime2     = p+q,   # cca arguments
		Wchoice     = "avg", 
		pre.scaling = TRUE,        
		oos         = TRUE,
		alpha       = NULL,
		n = 100, m = 100, nmc = 100,
		hardest.alt = TRUE,
		old.gauss.model.param=FALSE,
		separability.entries.w,
		compare.pom.cca=TRUE,
		oos.use.imputed,
		level.mcnemar=0.01,
		def.w=0.5,
		rival.w,		
		proc.dilation=FALSE,
		assume.matched.for.oos,
		w.vals,
		wt.equalize,
		power.comparison.test,
		verbose=FALSE) {
	## p: draw observations (signal) on Delta^p \in R^{p+1}
	## r: determine the matchedness between matched pairs --
	##    r = 0, independent(?); r = infinity, exact match
	## q: noise from Delta^q \in R^{q+1}
	## c.val: weight on noise
	##    X_{ik} ~ [(1-c.val)*signal, c.val*noise], i = 1, ..., n; k = 1, 2
	## d: dimension of projection space
	## alpha: Normal on Delta^p (NULL) or given
	##        if given, alpha is a list of length nmc, in which each
	##        component is a (n+2*m)x(p+1) matrix
	## Wchoice:
	##    "avg" => (D1 + D2)/2
	##    "sqrt" => sqrt((D1^2 + D2^2)/2)
	##    
	## n  : number of matched pairs
	## m  : number of testing pairs -- m matched pairs and m unmatched pairs
	## nmc: number of MC replicates
	##
	## Return: powers
	
	w.max.index <- length(w.vals)
	
	size <- seq(0, 1, 0.01)
	len <- length(size)
	power <- array(0,dim=c(w.max.index,nmc,len))
	power.cmp<-list(pom= array(0,dim=c(nmc,len)), cca= array(0,dim=c(nmc,len)) )
	optim.power <- array(0,dim=c(nmc,len))
	agg.cont.table <- matrix(0,2,2)
	empty.cont.tab<- list(matrix(0,2,2))
	cont.tables<-rep(empty.cont.tab,nmc)
	config.dist<- array(0,dim=c(nmc,w.max.index,3))
	min.stress<- array(0,dim=c(nmc,w.max.index+1))
	
	Fid.Terms.1<-c()
	Fid.Terms.2<-c()
	Comm.Terms <-c()
	
	Fid.Sum.Terms.1<-c()
	Fid.Sum.Terms.2<-c()
	Comm.Sum.Terms <-c()
	
	F.to.C.ratio <-  c()
	wtF.to.C.ratio <- c()
	F.bar.to.C.bar.ratio <-c()
	
	seeds<-rep(list(),nmc)
	sfInit( parallel=TRUE, cpus=num.cpus )
	
	
	
	#p.prime.cond = p+q
	sfExport( "p", "r", "q", "c.val",  ##try(
			"d","pprime1","pprime2",#"p.prime.cond",
			"Wchoice" ,
			"pre.scaling",
			"oos",
			"alpha",
			"n",
			"m",
			
			"old.gauss.model.param",
			"separability.entries.w",
			"compare.pom.cca",
			"oos.use.imputed",
			"level.mcnemar",
			"def.w",
			"rival.w",
			"proc.dilation",
			"assume.matched.for.oos",
			"w.vals",
			"wt.equalize",
			"verbose",
			"power.comparison.test" )
	
	
	print("Starting parallelization in gaussian_simulation_jofc_tradeoff_sf") 
	par.mc.result <- sfLapply( 1:nmc, run.mc.rep.with.seed)
	sfStop()
	
	
	
	
	#sink(file=file.path('logs',paste("traceback-debug-G-",mc,".txt",collapse="")))
	#traceback()
	#sink()
	#If  simulation fails, generate a results list of NAs
	
	
	
	
	
	
	
	# Number of elements in the par.mc.result list for each mc replicate	
	num.value.per.mc.rep <-13
	#if (verbose) print("str of par.mc.result")
	#if (verbose) print(str(par.mc.result))
	for (i in 1:nmc){
		mc.res.i <- par.mc.result[[i]]
		sink("debug.mc.txt")
		print(str(mc.res.i))
		print(mc.res.i)
		sink()
		power[,i,] <- mc.res.i[[1]]
		#power.cmp is has values of 0  if compare.pom.cca is FALSE
		if (compare.pom.cca) {
			power.cmp$cca[i,] <-mc.res.i[[2]]$cca
			power.cmp$pom[i,] <-mc.res.i[[2]]$pom
		}
		cont.tables[[i]]<-mc.res.i[[3]]
		
		agg.cont.table<-agg.cont.table+cont.tables[[i]]
		
		if (verbose) print(str(par.mc.result[[i]]))
		if (verbose) print (paste(i,"th contingency table"))
		if (verbose) print(cont.tables[[i]])
		config.dist[i,,1]<- mc.res.i[[4]]$frob.norm
		
		min.stress[i,] <- mc.res.i[[5]]
#		means[i,]      <- par.mc.result[[(num.value.per.mc.rep*(i-1))+6]]
		Fid.Term.1.i <- mc.res.i[[7]]$F1
		Fid.Term.2.i <- mc.res.i[[7]]$F2
		Comm.Term.i <- mc.res.i[[7]]$C
		Fid.Terms.1 <- rbind(Fid.Terms.1,Fid.Term.1.i)
		Fid.Terms.2 <- rbind(Fid.Terms.2,Fid.Term.2.i)
		Comm.Terms <- rbind(Comm.Terms,Comm.Term.i)
		
		Fid.Sum.Term.1.i <- mc.res.i[[8]]$F1
		Fid.Sum.Term.2.i <- mc.res.i[[8]]$F2
		Comm.Sum.Term.i <- mc.res.i[[8]]$C
		Fid.Sum.Terms.1 <- rbind(Fid.Sum.Terms.1,Fid.Sum.Term.1.i)
		Fid.Sum.Terms.2 <- rbind(Fid.Sum.Terms.2,Fid.Sum.Term.2.i)
		Comm.Sum.Terms <- rbind(Comm.Sum.Terms,Comm.Sum.Term.i)
		
		
		
		
		F.to.C.ratio.i <-mc.res.i[[9]]
		wtF.to.C.ratio.i <- mc.res.i[[10]]
		F.bar.to.C.bar.ratio.i <- mc.res.i[[11]]
		F.to.C.ratio <-  rbind(F.to.C.ratio,F.to.C.ratio.i)
		wtF.to.C.ratio <- rbind(wtF.to.C.ratio,wtF.to.C.ratio.i)
		F.bar.to.C.bar.ratio <- rbind(F.bar.to.C.bar.ratio,F.bar.to.C.bar.ratio.i)
		optim.power[i,]<-mc.res.i[[12]]
	}	
	
	
	print(agg.cont.table)
	
	FC.terms <- list(F1=Fid.Terms.1, F2=Fid.Terms.2, C=Comm.Terms)
	#FC.sum.terms<-
	FC.ratios<-list(f.c=F.to.C.ratio,wtf.c=wtF.to.C.ratio,f.c.bar=F.bar.to.C.bar.ratio)
	if (verbose) print(str(FC.ratios))
	return (list(power=power,power.cmp = power.cmp, conting.table=agg.cont.table,conting.table.list=cont.tables,
					config.dist=config.dist,min.stress=min.stress,seeds=seeds, FidComm.Terms=FC.terms,
					FC.ratios=FC.ratios,optim.power=optim.power))
}

run.mc.rep.with.seed <-function(seed){
	
	source(file.path("lib","simulation_math_util_fn.R"))
	source(file.path("lib","oosMDS.R"))
	source(file.path("lib","smacofM.R"))
	source(file.path("lib","oosIM.R"))
	library(MASS)
	print("Lib functions loaded")
	
	
	
	sink(file=file.path('logs',paste("debug-G-",seed,".txt",collapse="")))
	set.seed(seed)
	#if(mc<=4) {for(i in 1:mc) print(mvrnorm(4,mu=rep(0,4),Sigma=diag(4)))}
	print(runif(2))
	print(rnorm(2))
	
	print(seed)
	#seeds<-c(seeds,list(.Random.seed))
	
	#seeds[[mc]]<-list(.Random.seed)
	sink()
	
	sink(file=file.path('logs',paste("debug-G-mc-rep-",seed,".txt",collapse="")))
	print("Running run.mc.replicate function")
	tmp<- run.mc.replicate("gaussian",p, r, q, c.val,  ##try(
			d           = d,
			pprime1=pprime1,   pprime2=pprime2,
			Wchoice     = Wchoice, 
			pre.scaling = pre.scaling,
			oos         = oos,
			alpha       = alpha,
			n = n, m = m,
			
			old.gauss.model.param=old.gauss.model.param,
			separability.entries.w=separability.entries.w,
			compare.pom.cca=compare.pom.cca,
			oos.use.imputed=oos.use.imputed,
			level.mcnemar=level.mcnemar,
			def.w=def.w,
			rival.w=rival.w,
			proc.dilation=proc.dilation,
			assume.matched.for.oos=assume.matched.for.oos,
			w.vals=w.vals,
			wt.equalize=wt.equalize,
			verbose=verbose,
			power.comparison.test=power.comparison.test) 
	#)
	
	sink()
	return(tmp)	
	
}
