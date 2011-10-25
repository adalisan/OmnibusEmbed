


gaussian_simulation <- function(p, r, q, c.val,
		d           = p-1,
		pprime1     = p+q,   # cca arguments
		pprime2     = p+q,   # cca arguments
		Wchoice     = "avg", 
		pre.scaling = TRUE,
		oos         = TRUE,
		alpha       = NULL,
		n = 100, m = 100, nmc = 100,
		old.gauss.model.param=FALSE,
		sim.grass=FALSE,
		eigen.spectrum=FALSE) {
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
	
	size <- seq(0, 1, 0.01)
	len <- length(size)
	power <- list(cca  = matrix(0, nmc, len),
			pom  = matrix(0, nmc, len),
			jofc = matrix(0, nmc, len),
			threeway.iden = matrix(0, nmc, len),
			threeway.diag = matrix(0, nmc, len),
			threeway.idio = matrix(0, nmc, len)
	)
	epsilon.pom <- matrix(0,nmc,3)
	initial.vectors.x1 <- matrix(0,p+q,p+q)
	initial.vectors.x2 <- matrix(0,p+q,p+q)
	eigenvalue.sample.x1 <- matrix (NA,p+q,nmc)
	eigenvalue.sample.x2 <- matrix (NA,p+q,nmc)
	grass.metric <- matrix(0,nmc,5)
	geo.dist <- matrix(0,nmc,5)
	Haursdorf <- matrix(0,nmc,5)
	
	for(mc in 1:nmc) {
		set.seed(mc)
		print(mc)
		sigma <- matrix(0,p,p)
		
		if (is.null(alpha)) {
			sigma<- diag(p)
			if (old.gauss.model) sigma <-Posdef(p,r)
			alpha.mc <- mvrnorm(n+(2*m), rep(0,p),sigma)
		} else {
			alpha.mc <- alpha[[mc]]
		}
		
		## n pairs of matched points
		xlist <- matched_rnorm(n, p, q, c.val, r, alpha=alpha.mc[1:n, ],sigma.alpha=sigma,old.gauss.model.param=old.gauss.model.param)
		X1 <- xlist$X1
		X2 <- xlist$X2
		D1 <- dist(X1)
		D2 <- dist(X2)
		if (pre.scaling) {
			s <- lm(as.vector(D1) ~ as.vector(D2) + 0)$coefficients
		} else {
			s <- 1
		}
		
		if (eigen.spectrum){
			
			e.decomp.x<-eigen(cov(X1))
			e.decomp.y<-eigen(cov(X2))
			
			
			if (mc==1){
				initial.vectors.x1 <-	e.decomp.x$vectors
				initial.vectors.x2 <-	e.decomp.y$vectors
				eigenvalue.sample.x1[,1] <- e.decomp.x$values
				eigenvalue.sample.x2[,1] <- e.decomp.y$values
			}
			else{
				for (cur.evector.ind in 1:(p+q)){
					cur.evector.match.1<-F
					cur.evector.match.2<-F
					for (init.evector.ind in 1:(p+q)){
						dot.pr.1 <- sum(e.decomp.x$vectors[,cur.evector.ind]*initial.vectors.x1[,init.evector.ind])
						dot.pr.2 <- sum(e.decomp.y$vectors[,cur.evector.ind]*initial.vectors.x2[,init.evector.ind])
						if ((!cur.evector.match.1)&dot.pr.1>0.95){
							eigenvalue.sample.x1[init.evector.ind,mc] <- 
									e.decomp.x$values[cur.evector.ind]
							cur.evector.match.1<-T
						}
						if ((!cur.evector.match.2)&dot.pr.2>0.95){
							eigenvalue.sample.x2[init.evector.ind,mc] <- 
									e.decomp.y$values[cur.evector.ind]
							cur.evector.match.2<-T
							
						}
						if (cur.evector.match.1&cur.evector.match.2)  break
					}
				}
			}
		}
		## test observations -- m pairs of matched and m pairs of unmatched
		ylist <- matched_rnorm(m, p, q, c.val, r, alpha=alpha.mc[(n+1):(n+m), ],sigma.alpha=sigma,old.gauss.model.param=old.gauss.model.param)
		Y1 <- ylist$X1
		Y20 <- ylist$X2
		Y2A <- matched_rnorm(m, p, q, c.val, r, alpha=alpha.mc[(n+m+1):(n+m+m), ],sigma.alpha=sigma,old.gauss.model.param=old.gauss.model.param)$X2
		D10A <- as.matrix(dist(rbind(X1, Y1)))
		D20 <- as.matrix(dist(rbind(X2, Y20))) * s
		D2A <- as.matrix(dist(rbind(X2, Y2A))) * s
		
		## ==== cca ====
		if (oos == TRUE) {
			if (c.val==0){
				X1t <- cmdscale(D1, p)
				X2t <- cmdscale(D2, p)
				
			} else{
				X1t <- cmdscale(D1, pprime1)
				X2t <- cmdscale(D2, pprime2)
			}
			
			
			xcca <- cancor(X1t, X2t)
			
#			if (use.Euc.points.y){
#				if (!use.Euc.points.x & c.val==0){
#				Y1t  <- (cbind(Y1,matrix(0,dim(Y1)[1],q))   %*% xcca$xcoef)[, 1:d]
#				Y20t <- (cbind(Y20,matrix(0,dim(Y20)[1],q)) %*% xcca$ycoef)[, 1:d]
#				Y2At <- (cbind(Y2A,matrix(0,dim(Y2A)[1],q)) %*% xcca$ycoef)[, 1:d]
#				} else{
#				Y1t  <- (Y1  %*% xcca$xcoef)[, 1:d]
#				Y20t <- (Y20 %*% xcca$ycoef)[, 1:d]
#				Y2At <- (Y2A %*% xcca$ycoef)[, 1:d]				
#				}
#			
#			} else{
			Y1t <- (oosMDS(D10A, X1t) %*% xcca$xcoef)[, 1:d]
			Y20t <- (oosMDS(D20, X2t) %*% xcca$ycoef)[, 1:d]
			Y2At <- (oosMDS(D2A, X2t) %*% xcca$ycoef)[, 1:d]
#		    }
		} else {
			X1t <- cmdscale(D10A, pprime1)
			X2t <- cmdscale(dist(rbind(X2, Y20, Y2A)), pprime2)
			center1 <- colMeans(X1t[1:n, ])   # column means of training obs
			center2 <- colMeans(X2t[1:n, ])
			X1t <- X1t - matrix(center1, n+m, pprime1, byrow=TRUE) # column-center training only
			X2t <- X2t - matrix(center2, n+2*m, pprime2, byrow=TRUE)
			cca <- cancor(X1t[1:n, ], X2t[1:n, ])
			Y1t <- (X1t[(n+1):(n+m), ] %*% cca$xcoef )[, 1:d]
			Y20t <- (X2t[(n+1):(n+m), ] %*% cca$ycoef)[, 1:d]
			Y2At <- (X2t[(n+m+1):(n+2*m), ] %*% cca$ycoef)[, 1:d]
		}
		T0 <- rowSums((Y1t - Y20t)^2)
		TA <- rowSums((Y1t - Y2At)^2)
		power$cca[mc, ] <- get_power(T0, TA, size)
		
		## ==== pom = procrustes o mds ====
		if (oos == TRUE) {
			X1t <- cmdscale(D1, d)
			X2t <- cmdscale(D2, d)
			proc <- procrustes(X2t, X1t, dilation=TRUE)
			Y1t <- oosMDS(D10A, X1t)
			Y20t <- oosMDS(D20, X2t) %*% proc$R * proc$s
			Y2At <- oosMDS(D2A, X2t) %*% proc$R * proc$s
		} else {
			X1t <- cmdscale(D10A, d)
			X2t <- cmdscale(dist(rbind(X2, Y20, Y2A)), d)
			center1 <- colMeans(X1t[1:n, ])
			center2 <- colMeans(X2t[1:n, ])
			X1t <- X1t - matrix(center1, n+m, d, byrow=TRUE) # column-center training only
			X2t <- X2t - matrix(center2, n+2*m, d, byrow=TRUE)
			proc <- procrustes(X2t[1:n, ], X1t[1:n, ], dilation=TRUE)
			Y1t <- X1t[(n+1):(n+m), ]
			Y20t <- X2t[(n+1):(n+m), ] %*% proc$R * proc$s
			Y2At <- X2t[(n+m+1):(n+2*m), ] %*% proc$R * proc$s
		}
		T0 <- rowSums((Y1t - Y20t)^2)
		TA <- rowSums((Y1t - Y2At)^2)
		power$pom[mc, ] <- get_power(T0, TA, size)
		
		
		#
		# Grassmannian metric
		if (sim.grass){
			if (oos == TRUE) {
				pc1 <- prcomp(X1t)
				pc2 <- prcomp(X2t)
				X1t.sub <- pc1$x[, 1:d]
				X2t.sub <- pc2$x[, 1:d]
				proc <- procrustes(X2t.sub, X1t.sub, dilation=TRUE)}
			else{
				pc1 <- prcomp(X1t)
				pc2 <- prcomp(X2t)
				proc <- procrustes(X2t[1:n, ], X1t[1:n, ], dilation=TRUE)
				X1t.sub <- pc1$x[1:n, 1:d]
				X2t.sub <- pc2$x[1:n, 1:d] %*% proc$R * proc$s
				
			}
			
			epsilon.pom[mc, ] <- get_epsilon(D1, dist(X1t[1:n,]), D2, dist(X2t[1:n,]),
					X1t[1:n, ], proc$X.new[1:n, ])
			
			
			R1 <- pc1$rotation[, 1:d]
			R2 <- pc2$rotation[, 1:d]
			
			#Q1<- R1%*%t(R1)
			#Q2<- R2%*%t(R2)
			Q1 <- R1
			Q2 <- R2
			grass.metric[mc, 1] <- grassmannian(Q1, Q2) # plain
			geo.dist[mc, 1] <- geo_dist(Q1, Q2)
			Haursdorf [mc, 1] <- Haursdorf_dist (Q1, Q2)
			
			sa <- lm(as.vector(X1t) ~ as.vector(X2t) + 0)$coefficients # scale after embedding
			
			#Q2<- R2%*%t(R2)
			Q2 <- R2 * sa
			grass.metric[mc, 2] <- grassmannian(Q1, Q2 )           # scale after embedding
			geo.dist[mc, 2] <- geo_dist(Q1, Q2)
			Haursdorf [mc, 2] <- Haursdorf_dist(Q1, Q2 )
			
			
			sb <- lm(as.vector(D1) ~ as.vector(D2) + 0)$coefficients   # scale before embedding
			Q2 <- R2 * sb
			#Q2<- R2%*%t(R2)
			
			grass.metric[mc, 3] <- grassmannian(Q1, Q2)           # scale before embedding
			geo.dist[mc, 3] <- geo_dist(Q1, Q2 )
			Haursdorf [mc, 3] <- Haursdorf_dist(Q1, Q2)
			
			P <- polarity(X2t, X1t)
			Q2 <- R2 %*% P * proc$s
			#Q2<- R2%*%t(R2)
			
			
			grass.metric[mc, 4] <- grassmannian(Q1, Q2)     # scale after embedding & polarity
			geo.dist[mc, 4] <- geo_dist(Q1, Q2 %*% P * proc$s)
			Haursdorf [mc, 4] <- Haursdorf_dist(Q1, Q2)
			
			
			Q2 <- R2%*% proc$R * proc$s
			#Q2<- R2%*%t(R2)
			#Q2 <- R2
			grass.metric[mc, 5] <- grassmannian(Q1, Q2 ) # procrustes
			geo.dist[mc, 5] <- geo_dist(Q1, Q2 )
			Haursdorf [mc, 5] <- Haursdorf_dist(Q1, Q2  )
		}
		
		
		
		
		## ==== jofc ====
		if (Wchoice == "avg") {
			W <- (D1 + D2)/2
		} else if (Wchoice == "sqrt") {
			W <- sqrt((D1^2 + D2^2)/2)
		}
		if (oos == TRUE) {
			M <- omnibusM(D1, D2, W)
			X <- cmdscale(M, d)
			X1t <- X[1:n, ]
			X2t <- X[-(1:n), ]
			Y1t <- oosMDS(dist(rbind(X1, Y1)), X1t)
			Y20t <- oosMDS(dist(rbind(X2, Y20)), X2t)
			Y2At <- oosMDS(dist(rbind(X2, Y2A)), X2t)
			T0 <- rowSums((Y1t - Y20t)^2)
			TA <- rowSums((Y1t - Y2At)^2)
		} else {
			M0 <- omnibusM(D10A, D20, W=(D10A+D20)/2)
			MA <- omnibusM(D10A, D2A, W=(D10A+D2A)/2)
			X0 <- cmdscale(M0, d)
			XA <- cmdscale(MA, d)
			T0 <- rowSums((X0[(n+1):(n+m), ] - X0[(n+m+n+1):(n+m+n+m), ])^2)
			TA <- rowSums((XA[(n+1):(n+m), ] - XA[(n+m+n+1):(n+m+n+m), ])^2)
		}
		power$jofc[mc, ] <- get_power(T0, TA, size)
		
		Test.stats.iden <- ThreewayMDS.Embed.Hyp.Test(D1,D2,X1,X2,Y1,Y20,Y2A,"identity",d)
		power$threeway.iden[mc, ] <- get_power(Test.stats.iden$T0, Test.stats.iden$TA, size)
		
		Test.stats.diag <- ThreewayMDS.Embed.Hyp.Test(D1,D2,X1,X2,Y1,Y20,Y2A,"diagonal",d)
		power$threeway.diag[mc, ] <- get_power(Test.stats.diag$T0, Test.stats.diag$TA, size)
		
		Test.stats.idio <- ThreewayMDS.Embed.Hyp.Test(D1,D2,X1,X2,Y1,Y20,Y2A,"idioscal",d)
		power$threeway.idio[mc, ] <- get_power(Test.stats.idio$T0, Test.stats.idio$TA, size)
	}    
	if (eigen.spectrum){
		plot.MC.evalues.with.CI(t(eigenvalue.sample.x1),"X1","blue",TRUE)
		plot.MC.evalues.with.CI(t(eigenvalue.sample.x2),"X2","red",TRUE)
	}
	return(list(power=power,grassmannian=grass.metric,geodesic=geo.dist,
					fc.pom = epsilon.pom,Haursdorf=Haursdorf))
}



gaussian_simulation_jofc_tradeoff <- function(p, r, q, c.val,
		d           = p-1,
		pprime1     = p+q,   # cca arguments
		pprime2     = p+q,   # cca arguments
		Wchoice     = "avg", 
		pre.scaling = TRUE,
		oos         = TRUE,
		alpha       = NULL,
		n = 100, m = 100, nmc = 100,
		sim.grass=FALSE,
		eigen.spectrum=FALSE,
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
		verbose) {
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
	config.dist<- array(0,dim=c(nmc,w.max.index,3))
	
	agg.cont.table <- matrix(0,2,2)
	empty.cont.tab<- list(matrix(0,2,2))
	cont.tables<-rep(empty.cont.tab,nmc)
	min.stress<- array(0,dim=c(nmc,w.max.index+1))
	
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
						sim.grass=sim.grass,
						eigen.spectrum=eigen.spectrum,
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
						verbose=verbose))
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
	
		
		
		
		
	}
	return (list(power=power,power.cmp=power.cmp, conting.table=agg.cont.table,conting.table.list=cont.tables,
					config.dist=config.dist,min.stress=min.stress))
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
		sim.grass=FALSE,
		eigen.spectrum=FALSE,
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
	
	F.to.C.ratio <-  c()
	wtF.to.C.ratio <- c()
	F.bar.to.C.bar.ratio <-c()
	
	seeds<-rep(list(),nmc)
	par.mc.result<- foreach(mc =1:nmc,.packages=c("MASS","MCMCpack","smacof") ) %dopar% { #,.export=c("power.comparison.test")
		
		source("./src/simulation_math_util_fn.R")
		
		source("./src/oosMDS.R")
		source("./src/smacofM.R")
		source("./src/oosIM.R")
		sink(file=file.path('log',paste("debug-G-",mc,".txt",collapse="")))
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
		sink(file=file.path('log',paste("traceback-debug-G-",mc,".txt",collapse="")))
		traceback()
		sink()
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
