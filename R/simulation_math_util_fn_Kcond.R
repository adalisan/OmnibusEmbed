## functions
##




#' Title
#'
#' @param model
#' @param p
#' @param r
#' @param q
#' @param c.val
#' @param K
#' @param d
#' @param p.prime.cond
#' @param Wchoice
#' @param pre.scaling
#' @param oos
#' @param alpha
#' @param n
#' @param m
#' @param hardest.alt
#' @param old.gauss.model.param
#' @param separability.entries.w
#' @param compare.pom.cca
#' @param oos.use.imputed
#' @param level.mcnemar
#' @param def.w
#' @param rival.w
#' @param proc.dilation
#' @param assume.matched.for.oos
#' @param w.vals
#' @param wt.equalize
#' @param verbose
#' @param power.comparison.test
#'
#' @return
#' @export
#'
#' @examples
run.mc.replicate.Kcond<-function(model,p, r, q, c.val,K,
		d           = p-1,
		p.prime.cond     = ifelse(model=="gaussian",rep(p+q-1,K),rep(p+q+1,K)),   # cca arguments , signal+noise dimension
		Wchoice     = "avg", #How to impute L
		pre.scaling = TRUE,  #Make the measurement spaces have the same scale
		oos         = TRUE,  #embed test observations by Out-of-sampling  ?
		alpha       = NULL,
		n = 100, m = 100,    #Number of training and test observations
		hardest.alt= FALSE,
		old.gauss.model.param=FALSE,
		separability.entries.w,
		compare.pom.cca=TRUE,  # Run PoM and CCA to compare with JOFC?
		oos.use.imputed,
		level.mcnemar=0.01,  #At what alpha, should unweighted(w=0.5) and optimal w^* be compared
		def.w=0.5,           #The null hypothesis is that power(def.w) >= power(rival.w) (by default ,def.w is the w for the unweighted case which equals 0.5)
		rival.w=NULL,
		proc.dilation=FALSE, #when investigating convergence of JOFC to PoM, should Procrustes analysis of configurations include the dilation component?
		assume.matched.for.oos,
		w.vals,				  #w values to use for JOFC
		wt.equalize,

		verbose=FALSE,
		power.comparison.test=TRUE){


	if (verbose) print("verbose logging")
#	try.flag <- try({
				print(paste("random ",runif(1)))
				print("run.mc.replicate")
				print("running parametres")
				print(paste(model,p, r, q, c.val,K,
		d          ,
		p.prime.cond    ,   # cca arguments , signal+noise dimension
		Wchoice    , #How to impute L
		pre.scaling,  #Make the measurement spaces have the same scale
		oos        ,  #embed test observations by Out-of-sampling  ?
		alpha      ,
		n , m ,    #Number of training and test observations
		hardest.alt,
		old.gauss.model.param,
		separability.entries.w,
		compare.pom.cca,  # Run PoM and CCA to compare with JOFC?
		oos.use.imputed,
		level.mcnemar,  #At what alpha, should unweighted(w=0.5) and optimal w^* be compared
		def.w,           #The null hypothesis is that power(def.w) >= power(rival.w) (by default ,def.w is the w for the unweighted case which equals 0.5)
		rival.w,
		proc.dilation, #when investigating convergence of JOFC to PoM, should Procrustes analysis of configurations include the dilation component?
		assume.matched.for.oos,
		w.vals,				  #w values to use for JOFC
		wt.equalize,

		verbose,
		power.comparison.test,sep=" ",collapse="    "))
				#
				# The followin if statement is Not really necessary, unless we change our mind about rival.w being the best in every MC replicate
				# and want to make rival.w a constant (preferably the best overall)

				if (is.null(rival.w)){
					if (0.95 %in% w.vals){
						rival.w=0.95
					}
					else if (0.9 %in% w.vals){
						rival.w=0.9
					}
					else if (0.99 %in% w.vals){
						rival.w=0.99
					}
				}
				w.max.index <- length(w.vals)
				size <- seq(0, 1, 0.01)
				len <- length(size)


				power.mc= array(0,dim=c(w.max.index,len))  #power values for JOFC in this MC replicate
				power.cca.mc = array(0,dim=c(len))         #power values for CCA in this MC replicate
				power.pom.mc = array(0,dim=c(len))         #power values for PoM in this MC replicate
				power.cca.reg.mc = array(0,dim=c(len))     #power values for reg CCA in this MC replicate


				config.mismatch <-  list(frob.norm=array(0,dim=c(w.max.index))) #Frob. norm of configuration difference
				#between PoM and JOFC with smallest w
				min.stress.for.w.val = array(0,dim=c(w.max.index))   #minimum stress value for  smacof algorithm
				pom.stress <- 0


				cont.table  <- matrix(0,2,2)


				sigma <- matrix(0,p,p)
				means <- array(0 , dim=c(w.max.index,2*d))

				if (verbose) print("generating alphas")
				if (verbose) print(alpha)
				if (verbose)  print(model)
				if (verbose)  print(p)
				if (verbose)  print(K)
				if (verbose)  print(Posdef)
				if (verbose)  print(mvrnorm)
				if (verbose)  print(old.gauss.model.param)


				if (verbose)  print(K)


				if (is.null(alpha)) {
					if (model=="gaussian"){
						sigma<- diag(p)
						if (old.gauss.model.param) sigma <-Posdef(p,r)
						alpha.mc <- mvrnorm(n+(K*m), rep(0,p),sigma)
					} else if (model=="dirichlet"){
						alpha.mc <- rdirichlet(n+(K*m), rep(1,p+1))
					} else stop("unknown model")


				} else {
					alpha.mc <- alpha
				}

				if (verbose) print("n pairs of matched points")

				## n pairs of matched points
				if (model=="gaussian"){
					xlist <- matched_rnorm_Kcond(n, p, q, c.val, r,K, alpha=alpha.mc[1:n, ],sigma.alpha=sigma,
							old.gauss.model.param=old.gauss.model.param,sigma.beta=NULL)
				} else{
					xlist <- matched_rdirichlet_Kcond(n, p, r, q, c.val,K, alpha.mc[1:n, ])
				}
				if (verbose) print("n pairs of K-matched  points  generated")
				x.config<- xlist$X

				if (verbose) print("Dissimilarities of matched points")
				D.cond.list<-list()
				for (k in 1:K){
					X.cond<-x.config[,,k]

					D.cond <- dist(X.cond)
					D.cond.list <- c(D.cond.list, list(D.cond))
				}

				if (model=="gaussian")
					sigma.mc<-xlist$sigma.beta


				if (verbose) print("random matched pairs generated\n")

				#prescaling
				sc.cond<-rep(1,K)
				if (pre.scaling) {
					for (cond.index in 2:K){
						sc.cond[cond.index] <- lm(as.vector(D.cond.list[[1]]) ~ as.vector(D.cond.list[[cond.index]]) + 0)$coefficients
						#D.cond.list[[cond.index]] <- D.cond.list[[cond.index]]*sc.cond[cond.index]
					}


				}



				Y.cond.A<-list()
				#m pairs of unmatched points
				if (model=="gaussian"){
					## test observations -- m pairs of matched and m pairs of unmatched
					ylist <- matched_rnorm_Kcond(m, p, q, c.val, r, K,alpha=alpha.mc[(n+1):(n+m), ],
							sigma.alpha=sigma,old.gauss.model.param=old.gauss.model.param, sigma.beta=sigma.mc)
					if (!hardest.alt) {
						for (cond.index in 1:K){
							Y.cond.A<-c(Y.cond.A,list(matched_rnorm_Kcond(m, p, q, c.val, r,K, alpha=alpha.mc[(n+(cond.index-1)*m+1):(n+(cond.index)*m), ],
													sigma.alpha=sigma,old.gauss.model.param=old.gauss.model.param, sigma.beta=sigma.mc)$X[,,cond.index]))
						}
					}


				} else{
					ylist <- matched_rdirichlet_Kcond(m, p, r, q, c.val,K, alpha.mc[(n+1):(n+m), ])
					if (!hardest.alt) {
						for (cond.index in 1:K){
							Y.cond.A<-c(Y.cond.A,list(matched_rdirichlet_Kcond(m, p, r, q, c.val,K,
													alpha=alpha.mc[(n+(cond.index-1)*m+1):(n+(cond.index)*m), ])$X[,,cond.index]))
						}
					}
				}

				if (hardest.alt) {
					alt.matched.tuple <- matched_rnorm_Kcond(m, p, q, c.val, r,K, alpha=alpha.mc[(n+1):(n+(m)), ],
							sigma.alpha=sigma,old.gauss.model.param=old.gauss.model.param, sigma.beta=sigma.mc)
					alt.unmatched.tuple.unmatched.source <- matched_rnorm_Kcond(m, p, q, c.val, r,K, alpha=alpha.mc[(n+m):(n+(2*m)), ],
							sigma.alpha=sigma,old.gauss.model.param=old.gauss.model.param, sigma.beta=sigma.mc)

					for (cond.index in 1:K){
						Y.cond.A<-c(Y.cond.A,list(alt.matched.tuple$X[,,cond.index]))
					}
					unmatched.cond <-sample(1:K,1)
					Y.cond.A[[unmatched.cond]] <- alt.unmatched.tuple.unmatched.source$X[,,unmatched.cond]
				}
				y.config<-ylist$X
				#
				# Dissimilarity matrices for in-sample + out-of-sample
				#
				# Dissimilarity matrices
				D.in.oos.list.0<-list()
				D.in.oos.list.A<-list()


				for (cond.idx in 1:K){
					diss.in.oos.cond.idx.null <- as.matrix(dist(rbind(x.config[,,cond.idx],y.config[,,cond.idx])))*sc.cond[cond.idx]
					D.in.oos.list.0<-c(D.in.oos.list.0,list(diss.in.oos.cond.idx.null))
				}
				for (cond.idx in 1:K){
					diss.in.oos.cond.idx.alt <- as.matrix(dist(rbind(x.config[,,cond.idx],Y.cond.A[[cond.idx]])))*sc.cond[cond.idx]
					D.in.oos.list.A<-c(D.in.oos.list.A,list(diss.in.oos.cond.idx.alt))
				}






				if (verbose) print("PoM and CCA embedding")
				if (compare.pom.cca) {

					# cca ====
					#embed in-sample measurements
					CCA.results <- run.cca.Kcond(D.cond.list,D.in.oos.list.0,D.in.oos.list.A,
							pprime.cond=p.prime.cond,d=d,hardest.alt = hardest.alt,size=size)
					T0.cca <- CCA.results$T0
					TA.cca <- CCA.results$TA
					power.cca.mc <- CCA.results$power


					if (verbose) print("CCA test statistic complete\n")

					## ==== pom ====
					Pom.results  <- run.pom.Kcond(D.cond.list,D.in.oos.list.0,D.in.oos.list.A,
							embed.dim=d,hardest.alt = hardest.alt,size=size)
					power.pom.mc <- Pom.results$power

				}
				JOFC.results <- run.jofc.Kcond(D.cond.list,w.vals,x.config,y.config,Y.cond.A,D.in.oos.list.0,D.in.oos.list.A,
						K,n,m,d,
						oos,separability.entries.w,init.conf,wt.equalize,Wchoice,assume.matched.for.oos,oos.use.imputed,
						verbose)
				if (verbose) print("JOFC test statistic complete \n")

				T0.w <- JOFC.results$T0
				TA.w <- JOFC.results$TA

				#if (verbose) print("T0.w")
				#if (verbose) print(T0.w)
				#if (verbose) print("TA.w")
				#if (verbose) print(TA.w)
				for (l in 1:w.max.index){
					power.mc[l, ] <- get_power(T0.w[l,], TA.w[l,], size)
				}


				# find which of the w values has the highest power at level.mcnemar
				# Use that w value for mcnemar's tes
				power.w.star <- 0



				for (l in 1:w.max.index) {

					power.mcnemar.l <- get_power(T0.w[l,],TA.w[l,],level.mcnemar)
					if (power.mcnemar.l>power.w.star){
						rival.w <- w.vals[l]
						power.w.star <- power.mcnemar.l
						w.val.rival.idx <- l
					}
				}








				################################################################
				################################################################
				# Power comparison test
				# In order to compare the best w^* vs w=0.5 in an unbiased way
				# re-run the simulation only for w= w^* and w=0.5
				# compute the contingency table using those results
				if (verbose) print("Power comparison test")

				if (power.comparison.test){

					if (verbose) print("Power comparison test starting")
					## n pairs of matched points
					D.cond.list<-list()

					w.vals.best.vs.equal.wt  <-c(w.vals[w.val.rival.idx],def.w)

					if (verbose) print("n matched points starting")

					## n pairs of matched points
					if (model=="gaussian"){
						xlist <- matched_rnorm_Kcond(n, p, q, c.val, r, K, alpha=alpha.mc[1:n, ],sigma.alpha=sigma,
								old.gauss.model.param=old.gauss.model.param,sigma.beta=NULL)
					} else{
						xlist <- matched_rdirichlet_Kcond(n, p, r, q, c.val,K, alpha.mc[1:n, ])
					}
					x.config<-xlist$X
					for (k in 1:K){
						X.cond<-x.config[,,k]
						D.cond <- dist(X.cond)
						D.cond.list <- c(D.cond.list, list(D.cond))
					}

					if (model=="gaussian")
						sigma.mc<-xlist$sigma.beta


					if (verbose) print("random matched pairs generated\n")

					#prescaling
					sc.cond<-rep(1,K)
					if (pre.scaling) {
						for (cond.index in 1:K)
							sc.cond[cond.index] <- lm(as.vector(D.cond.list[[1]]) ~ as.vector(D.cond.list[[cond.index]]) + 0)$coefficients

						#Apply the scaling coefficients
					}
					Y.cond.A<-list()





					#m pairs of unmatched points
					if (model=="gaussian"){
						## test observations -- m pairs of matched and m pairs of unmatched
						ylist <- matched_rnorm_Kcond(m, p, q, c.val, r,K, alpha=alpha.mc[(n+1):(n+m), ],
								sigma.alpha=sigma,old.gauss.model.param=old.gauss.model.param, sigma.beta=sigma.mc)
						for (cond.index in 1:K){
							Y.cond.A<-c(Y.cond.A,list(matched_rnorm_Kcond(m, p, q, c.val, r,K, alpha=alpha.mc[(n+(cond.index-1)*m+1):(n+(cond.index)*m), ],
													sigma.alpha=sigma,old.gauss.model.param=old.gauss.model.param, sigma.beta=sigma.mc)$X[,,cond.index]))
						}


					} else{
						ylist <- matched_rdirichlet_Kcond(m, p, r, q, c.val,K, alpha.mc[(n+1):(n+m), ])
						for (cond.index in 1:K){
							Y.cond.A<-c(Y.cond.A,list(matched_rdirichlet_Kcond(m, p, r,q, c.val, K,alpha=alpha.mc[(n+(cond.index-1)*m+1):(n+(cond.index)*m), ],
													sigma.alpha=sigma,old.gauss.model.param=old.gauss.model.param, sigma.beta=sigma.mc)$X[,,cond.index]))
						}
					}


					Y.cond.A<-list()
					#m pairs of unmatched points
					if (model=="gaussian"){
						## test observations -- m pairs of matched and m pairs of unmatched
						ylist <- matched_rnorm_Kcond(m, p, q, c.val, r, K,alpha=alpha.mc[(n+1):(n+m), ],
								sigma.alpha=sigma,old.gauss.model.param=old.gauss.model.param, sigma.beta=sigma.mc)
						if (!hardest.alt) {
							for (cond.index in 1:K){
								Y.cond.A<-c(Y.cond.A,list(matched_rnorm_Kcond(m, p, q, c.val, r,K, alpha=alpha.mc[(n+(cond.index-1)*m+1):(n+(cond.index)*m), ],
														sigma.alpha=sigma,old.gauss.model.param=old.gauss.model.param, sigma.beta=sigma.mc)$X[,,cond.index]))
							}
						}


					} else{
						ylist <- matched_rdirichlet_Kcond(m, p, r, q, c.val,K, alpha.mc[(n+1):(n+m), ])
						if (!hardest.alt) {
							for (cond.index in 1:K){
								Y.cond.A<-c(Y.cond.A,list(matched_rdirichlet_Kcond(m, p, r, q, c.val,K,
														alpha=alpha.mc[(n+(cond.index-1)*m+1):(n+(cond.index)*m), ])$X[,,cond.index]))
							}
						}
					}

					if (hardest.alt) {
						alt.matched.tuple <- matched_rnorm_Kcond(m, p, q, c.val, r,K, alpha=alpha.mc[(n+1):(n+(m)), ],
								sigma.alpha=sigma,old.gauss.model.param=old.gauss.model.param, sigma.beta=sigma.mc)
						alt.unmatched.tuple.unmatched.source <- matched_rnorm_Kcond(m, p, q, c.val, r,K, alpha=alpha.mc[(n+m):(n+(2*m)), ],
								sigma.alpha=sigma,old.gauss.model.param=old.gauss.model.param, sigma.beta=sigma.mc)

						for (cond.index in 1:K){
							Y.cond.A<-c(Y.cond.A,list(alt.matched.tuple$X[,,cond.index]))
						}
						unmatched.cond <-sample(1:K,1)
						Y.cond.A[[unmatched.cond]] <- alt.unmatched.tuple.unmatched.source$X[,,unmatched.cond]
					}
					y.config<-ylist$X







					#
					# Dissimilarity matrices for in-sample + out-of-sample
					#
					# Dissimilarity matrices
					D.in.oos.list.0<-list()
					D.in.oos.list.A<-list()


					for (cond.idx in 1:K){
						diss.in.oos.cond.idx.null <- as.matrix(dist(rbind(x.config[,,cond.idx],y.config[,,cond.idx])))*sc.cond[cond.idx]
						D.in.oos.list.0<-c(D.in.oos.list.0,list(diss.in.oos.cond.idx.null))
					}
					for (cond.idx in 1:K){
						diss.in.oos.cond.idx.alt  <- as.matrix(dist(rbind(x.config[,,cond.idx],Y.cond.A[[cond.idx]])))*sc.cond[cond.idx]
						D.in.oos.list.A<-c(D.in.oos.list.A,list(diss.in.oos.cond.idx.alt))
					}



					print("Running JOFC")

					JOFC.results.new <- run.jofc.Kcond(D.cond.list,w.vals.best.vs.equal.wt,
							x.config,y.config,Y.cond.A,D.in.oos.list.0,D.in.oos.list.A,
							K,n,m,d,
							oos,separability.entries.w,init.conf,wt.equalize,Wchoice,assume.matched.for.oos,oos.use.imputed,
							verbose)

					print("Ended JOFC")

					T0.best.w<- JOFC.results.new$T0
					TA.best.w<- JOFC.results.new$TA

					crit.value<-get_crit_val(T0.best.w[1,],level.mcnemar)
					crit.value.2<-get_crit_val(T0.best.w[2,],level.mcnemar)
					if (verbose){
						print("crit.values")
						print(crit.value)
						print(crit.value.2)
					}
					cont.table[1,1] <- sum(T0.best.w[1,]<=crit.value & T0.best.w[2,]<=crit.value.2) +
							sum(TA.best.w[1,]>crit.value & TA.best.w[2,]>crit.value.2)
					cont.table[1,2] <- sum(T0.best.w[1,]>crit.value & T0.best.w[2,]<=crit.value.2)  +
							sum(TA.best.w[1,]<=crit.value & TA.best.w[2,]>crit.value.2)
					cont.table[2,1] <- sum(T0.best.w[1,]<=crit.value & T0.best.w[2,]>crit.value.2)  +
							sum(TA.best.w[1,]>crit.value & TA.best.w[2,]<=crit.value.2)
					cont.table[2,2] <- sum(T0.best.w[1,]>crit.value & T0.best.w[2,]>crit.value.2)   +
							sum(TA.best.w[1,]<=crit.value & TA.best.w[2,]<=crit.value.2)

				}




				print("end run.mc.replicate")
				func.return<- list(power.mc=power.mc,power.cmp=list(cca = power.cca.mc,pom = power.pom.mc,cca.reg =power.cca.reg.mc),
						cont.tables=cont.table,
						config.dist= config.mismatch,
						min.stress=c(min.stress.for.w.val,pom.stress)
				#		,means=means,FidComm.Terms=FidComm.Terms,
				#		FidComm.Sum.Terms = FidComm.Sum.Terms,F.to.C.ratio = FC.ratio, wtF.to.C.ratio=FC.ratio.2,
				#		F.bar.to.C.bar.ratio= FC.ratio.3
				)
	#		}) #end try

	#if (inherits(try.flag,"try-error")) {
	#	print("run.mc.replicate error")
	#	sink(file.path('logs',"run.mc.rep-traceback.txt"))
	#	traceback()
	#	sink()
	#}
	return(func.return)
}


#' Title
#'
#' @param D.cond.list
#' @param w.vals
#' @param x.config
#' @param y.config
#' @param Y.cond.A
#' @param D.in.oos.list.0
#' @param D.in.oos.list.A
#' @param K
#' @param n
#' @param m
#' @param d
#' @param oos
#' @param separability.entries.w
#' @param init.conf
#' @param wt.equalize
#' @param Wchoice
#' @param assume.matched.for.oos
#' @param oos.use.imputed
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
run.jofc.Kcond<-function(D.cond.list,w.vals,x.config,y.config,Y.cond.A,D.in.oos.list.0,D.in.oos.list.A,
		K,n,m,d,
		oos,separability.entries.w,init.conf,wt.equalize,Wchoice,assume.matched.for.oos,oos.use.imputed,
		verbose){



	w.max.index<-length(w.vals)
	T0 <- matrix(0,w.max.index,m)   #Test statistics for JOFC under null
	TA <- matrix(0,w.max.index,m)    #Test statistics for JOFC under alternative

	#
	#  Forming omnibus dissimilarity matrix  and weight matrix for embedding via MDS
	#

# O/----------------M-----------------------------\              ideal.omnibus.0[1:(K*n),(K*n)+(1:(K*m))]
# O|                                              | _ ..----------------------------^----------------------.._
# OO----------------,-----------------------------:/         or  ideal.omnibus.A[1:(K*n),(K*n)+(1:(K*m))]     `-.
# OO----------------,-----------------------------|-------------------------------------------------------------:
# ||                |             ||              ||               ||            ||                             |
# ||                |             ||              ||               ||            ||                             |
# ||                |             ||              ||               .|            ||                             |
# ||                |intercondition.diss[1,,]     ||               ||            ||                             |
# || D_cond[[1]]    |             ||intercondition.diss[2,,]       ||            ||                             |
# ||                |             ||              ||               ||            ||                             |
# ||                |             ||              ||               ||            ||                             |
# OO----------------+-------------++--------------++------------------------------------------------------------|
# OO----------------+-------------++--------------++------------------------------------------------------------|
# ||                |             ||              ||               ||            ||                             |
# ||    .           |             ||              ||               ||            ||                             |
# ||    .           |             || intercondition.diss[3,,]      ||            ||                             |
# ||    .           |  D_cond[[2]]||              ||               ||            ||                             |
# ||                |             ||              ||               ||            ||                             |
# ||                |             ||              ||               ||            ||                             |
# ||                O-------------OO--------------OO              ```............/|-----------------------------|
# ||                O-------------OO--------------OO                               ---------------------------- |
# ||                              ||              ||                           M.oos.0                          |
# ||     .                 .      ||              ||                           M.oos.A                          |
# ||     .                 .      ||  D_cond[[3]] ||                            ^                               |
# ||     .                 .      ||              ||                         ___|___                            |
# ||                              ||              ||     __,,...------'''''''       `'''''''------.....___      |
# ++------------------------------OO--------------O.,--''-------------------------------------------------`'''---
# ++------------------------------OO--------------/O------------------------------------------------------------|
# ||                                             /||               ||          ||                               |
# ||                                             |||               ||          ||                               |
# ||                                             |||                intercondition.diss.tilde.OOS.null[]        |
# ||                   .                         |||D.oos.null[[1]]            ||                               |
# ||                   .                         |||D.oos.alt[[1]]|||          |intercondition.diss.tilde.OOS.nul
# ||                   .                         |||               ||          ||                               |
# ||                   .                         |||               ||          ||                               |
# ||                                             |||               ||          ||                               |
# ||                                             |OO-------------- ++----------++                               |
# ||                   .                         |OO-------------- ++----------++                               |
# ||                   .                         |||               ||          ||                               |
# ||                   .                         |||               ||         intercondition.diss.tilde.OOS.null|
# ||                                             |||               D.oos.null[[2]]                              |
# ||                                             |||               D.oos.alt[[2]]                               |
# ||                                             |||               ||          ||                               |
# ||                   .                         |||               ||          ||                               |
# ||                   .                         |||               OO----------OO-------------------------------|
# ||                   .                         |||                           ||                               |
# ||                                             |||                           ||                               |
# ||                                             |||                           ||D.oos.null[[3]]                |
# ||                                             \||                           || D.oos.alt[[3]]                |
# ++----------------------------------------------\+---------------------------OO-------------------------------/
#



	## ==== jofc ====
	intercondition.diss<- array(0,dim=c(K*(K-1)/2,n,n))
	# Impute "between-condition" dissimilarities from different objects
	for (l in 2:K)
		for (k in 1:(l-1))  {
			diss.matrix.index <- ((l-1)*(l-2)/2)+k
			#print(c(l,k,diss.matrix.index))
			if (Wchoice == "avg") {
				intercondition.diss[diss.matrix.index,,] <- as.matrix((D.cond.list[[l]] + D.cond.list[[k]])/2)
			} else if (Wchoice == "sqrt") {
				intercondition.diss[diss.matrix.index,,] <- as.matrix(sqrt((D.cond.list[[l]]^2 + D.cond.list[[k]]^2)/2))

			} else if (Wchoice == "NA+diag(0)") {
				intercondition.diss[diss.matrix.index,,] <- matrix(NA,n,n)
				diag(intercondition.diss[diss.matrix.index,,])<- 0
			}
		}


	if (oos == TRUE) {

		#In sample embedding
		# Form omnibus dissimilarity matrix
		M <- omnibusM.Kcond(D.cond.list, intercondition.diss)
		init.conf<-NULL

#		if (compare.pom.cca) init.conf<- pom.config
		if (verbose) print("startingJOFC in-sample embedding")
		# Embed in-sample using different weight matrices (different w values)
		X.in.JOFC.embeds <- in.embed.MDS.smacof.Kcond(n*K,M,d,K,w.vals,

				separability.entries.w,init.conf=init.conf,wt.equalize=wt.equalize)
		if (verbose) print("JOFC in-sample embedded")
#
#		Fid.Err.Terms <- X.in.JOFC.embeds[[w.max.index+2]]
#
#		Comm.Err.Terms  <- X.in.JOFC.embeds[[w.max.index+3]]
#
#		Fid.Err.Sum.Terms <- X.in.JOFC.embeds[[w.max.index+4]]
#
#		Comm.Err.Sum.Terms  <- X.in.JOFC.embeds[[w.max.index+5]]
#		FC.ratio  <- X.in.JOFC.embeds[[w.max.index+6]]
#		FC.ratio.2  <- X.in.JOFC.embeds[[w.max.index+7]]
#		FC.ratio.3  <- X.in.JOFC.embeds[[w.max.index+8]]
#
#		min.stress.for.w.val <- X.embeds[[w.max.index+1]]
		if (verbose) print("JOFC embeddings complete\n")


		#
		# OOS Dissimilarity matrices
		#

		D.oos.null <- list()
		D.oos.alt  <- list()

		for (k in 1:K) {
			D.oos.null.k<-dist(y.config[,,k])
			D.oos.null<-c(D.oos.null,list(D.oos.null.k))
		}
		for (k in 1:K) {
			D.oos.alt.k<-dist(Y.cond.A[[k]])
			D.oos.alt<-c(D.oos.alt,list(D.oos.alt.k))
		}



		intercondition.diss.tilde.OOS.null <- array(0,dim=c(K*(K-1)/2,m,m))
		intercondition.diss.tilde.OOS.alt <- array(0,dim=c(K*(K-1)/2,m,m))
		#Imputing dissimilarity  entries for OOS
		for (l in 2:K)
			for (k in 1:(l-1))  {
				diss.matrix.index <- (l-1)*(l-2)/2+k
				if (Wchoice == "avg") {
					intercondition.diss.tilde.OOS.null[diss.matrix.index,,] <- as.matrix((D.oos.null[[l]] + D.oos.null[[k]])/2)
					intercondition.diss.tilde.OOS.alt[diss.matrix.index,,]  <- as.matrix((D.oos.alt[[l]] + D.oos.alt[[k]])/2)
				} else if (Wchoice == "sqrt") {
					intercondition.diss.tilde.OOS.null[diss.matrix.index,,] <- as.matrix(sqrt((D.oos.null[[l]]^2 + D.oos.null[[k]]^2)/2))
					intercondition.diss.tilde.OOS.alt[diss.matrix.index,,]  <- as.matrix(sqrt((D.oos.alt[[l]]^2 + D.oos.alt[[k]]^2)/2))

				} else if (Wchoice == "NA+diag(0)") {
					NA.diag.0.mat <- matrix(NA,m,m)
					diag(NA.diag.0.mat)<- 0
					intercondition.diss.tilde.OOS.null[diss.matrix.index,,] <- NA.diag.0.mat
					intercondition.diss.tilde.OOS.alt[diss.matrix.index,,] <- NA.diag.0.mat

				}
			}

		#Form OOS omnibus matrices
		M.oos.0<- omnibusM.Kcond(D.oos.null,intercondition.diss.tilde.OOS.null)
		M.oos.A<- omnibusM.Kcond(D.oos.alt,intercondition.diss.tilde.OOS.alt)



		for (l in 1:w.max.index){
			if (verbose) print("OOS embedding for JOFC for w= \n")
			if (verbose) print(w.vals[l])

			w.val.l <- w.vals[l]
			X <- X.in.JOFC.embeds[[l]]


			oos.obs.flag<- c(rep(1,K*n),rep(0,K*m))

			#Forming omnibus weight matrix for oos embedding

# O/----------------oos.Weight.mat.1--------------\     |
# O|                                               |
# OO---------------------------------------------- |,..----m-----..
# OO---------------\,-----------------------------.)-----------------------------------------------..------.
# ||               |w  **2**   0   |w           0  |              |||             |                | \     .
# ||               |  w       0 0  |  w        0 0 |              |||             |                | |     |
# ||               |    w     0 0  |    w      0 0 |              |.| **1**       |                |  |    |
# ||               |      w    0   |      w     0  |   1-w        |||             |    **1**       |  |    |
# ||    1-w        |    0   w      |    0   w      |              |||             |                |  |->n |
# ||               |   0 0    w    |   0 0    w    |              |||             |                |  |    '.
# ||               |   0 0      w  |   0 0      w  |              |||             |                | |      | o
# OO---------------O    0         w|    0         w|              |||             |                | /      | o
# OO---------------O -------------Ow------------------------------++-------------------------------|.       | s
# ||               |              |  w         0  ||              ||              |                |        | .
# ||               |              |    w      0 0 ||              ||              |                |        | W
# ||      .        |              |      w    0 0       .         ||              |                |        | e
# ||      .        |     1-w      |     0  w   0  ||    .         ||     1-w      |    **1**       |        | i
# ||      .        |              |    0 0   w    ||    .         ||              |                |        | g
# ||               |              |    0 0     w  ||              ||              |                |        | h
# ||               O--------------|     0        wOO              ||              |                |        . t
# ||               O--------------+---------------OO-----------------------------------------------|       /  .
# ||                              ||              ||                              |                |       |  m
# ||     .              .         ||              ||     .                 .      |                |       |  a
# ||     .              .         ||    1-w       ||     .                 .      |      1-w       |       |  t
# ||     .              .         ||              ||     .                 .      |                |       |  .
# ||                              ||              ||                              |                |     _,   w
# +\--------------/---------------OO--------------OO-----------------------------------------------|.---'
# ++'.-----.----.'.---------------OO--------------/O-------------------oos.Weight.mat.2------------|
# ||       n                                     /||               .w              .w           0  |
# ||                                             ||||              |  w **3**      |  w        0 0 |
# ||                                             ||||              |    w          |    w      0 0 |
# ||                    .                        ||||              |      w        |   0  w     0  |
# ||                    .                        ||     1-w        |   0    w      |  0 0   w      |
# ||                    .                        ||                |  0 0     w    |  0 0     w    |
# ||                                             ||||              |  0 0       w  |   0        w  |
# ||                    .                        ||::              |   0  **2**   w|              w|
# ||                    .                        ||--------------------------------+-w-------------|
# ||                    .                        ||||              ||              |   w        0  |
# ||                                             ||||              ||              |     w     0 0 |
# ||                                             ||      .         ||              |   0   w   0 0 |
# ||                    .                        ||||    .         |    1-w        |  0 0    w  0  |
# ||                    .                        ||||    .         |               |  0 0      w   |
# ||                    .                        ||||              ||              |   0         w |
# ||                                             ||OO--------------||----------------------------- +.
# ||                    .                        ||OO  **3**If  assume.matched.for.oos is F        | \
# ||                    .                        ||||         w terms in oos.Weight.mat.2          | `.
# ||                    .                        ||||    .         |         are replaced by 0     |   |
# ||                                             ||||    .         |      .        |               |   |->m
# ||                                             ||||    .         |      .        |       1-w     |   |
# ++-----------------------------------------------\.-------------/--------------------------------|- '
#                                                   \....._,,,.../
#                                                         m
# **1** If oos.use.imputed is T, this part of    oos.Weight.mat.w is 1-w  else they are 0
# **2**  If sep.err.w is F, these sep terms are 0   otherwise they are 1-w like fid terms.
# **3**If  assume.matched.for.oos is F        w terms in oos.Weight.mat.2  are replaced by 0



			#Compute Weight matrix corresponding in-sample  entries
			oos.Weight.mat.1<-w.val.to.W.mat.Kcond(w.val.l,(K*n),K,separability.entries.w,wt.equalize)

			#Compute Weight matrix corresponding OOS  entries
			oos.Weight.mat.2<-w.val.to.W.mat.Kcond(w.val.l,(K*m),K,separability.entries.w,wt.equalize)

			# If assume.matched.for.oos is true, we assume OOS dissimilarities are matched(in reality,
			# they are matched for the matched pairs, but unmatched for the unmatched pairs)
			# Alternatively, If assume.matched.for.oos is false, we ignore the dissimilarities between matched/unmatched
			# pairs
			if (!assume.matched.for.oos){

				for (k in 1:(K))
					for (u in 1:(K)){
						if (k==u) next
						diag.entries.of.submatrices<-cbind(1:m+(k-1)*m,(1:m)+(u-1)*m)

						oos.Weight.mat.2[diag.entries.of.submatrices]<-0

					}
			}
			# if (oos.use.imputed is true) we treat the dissimiilarities between  in-sample and out-of-sample measurements
			# from different conditions like fidelity terms
			# otherwise they are ignored
			if (oos.use.imputed){
				oos.Weight.mat.w <- matrix(1-w.val.l,K*n,K*m)
			} else{
				block.diag.mat <- matrix(0,K*n,K*m)
				for (k in 0:(K-1)){
					block.diag.mat[1:n+(k*n),1:m+(k*m)] <- 1-w.val.l
				}
				oos.Weight.mat.w <- block.diag.mat

			}
			oos.Weight.mat<-omnibusM(oos.Weight.mat.1,oos.Weight.mat.2,oos.Weight.mat.w)
			# Since we are going to oos-embedding, set the weights  of in-sample embedding of stress
			# We are using previous in-sample embeddings, anyway
			oos.Weight.mat[1:(K*n),1:(K*n)]<-0
			if (verbose) print("dim(M.oos.0)")
			if (verbose) print(dim(M.oos.0))
			if (verbose) print("dim(M.oos.A)")
			if (verbose) print(dim(M.oos.A))
			if (verbose) print("dim(oos.Weight.mat)")
			if (verbose) print(dim(oos.Weight.mat))
			if (verbose) print("dim(X)")
			if (verbose) print(dim(X))

			if (verbose) print("ideal omnibus matrices(not used for embedding including in-oos diss. from orig.measurements")

			ideal.omnibus.0  <- as.matrix(diss.matrix(x.config,y.config))
			ideal.omnibus.A  <- as.matrix(diss.matrix(x.config,Y.cond.A))

			if (verbose) print("forming  omnibus matrices")
			omnibus.oos.D.0 <- omnibusM(M,M.oos.0,ideal.omnibus.0[1:(K*n),(K*n)+(1:(K*m))])
			omnibus.oos.D.A <- omnibusM(M,M.oos.A, ideal.omnibus.A[1:(K*n),(K*n)+(1:(K*m))])
			oos.Weight.mat[is.na(omnibus.oos.D.0)]<-0
			# This commented line is unnecessary, since it accomplishes the same task as the line above, setting weights for ignored dissimilarities(NAs) to 0
			#oos.Weight.mat[is.na(omnibus.oos.D.A)]<-0
			omnibus.oos.D.0[is.na(omnibus.oos.D.0)]<-1
			omnibus.oos.D.A[is.na(omnibus.oos.D.A)]<-1
			if (verbose) print("JOFC null omnibus OOS embedding \n")
			Y.0t<-oosIM(D=omnibus.oos.D.0,
					X=X,
					init     = "random",
					verbose  = FALSE,
					itmax    = 1000,
					eps      = 1e-8,
					W        = oos.Weight.mat,
					isWithin = oos.obs.flag,
					bwOos    = TRUE)

			if (verbose) print("JOFC alternative omnibus OOS embedding \n")
			Y.At<-oosIM(D=omnibus.oos.D.A,
					X=X,
					init     = "random",
					verbose  = FALSE,
					itmax    = 1000,
					eps      = 1e-8,
					W        = oos.Weight.mat,
					isWithin = oos.obs.flag,
					bwOos    = TRUE)
			#print(Y.0t)
			#print(Y.At)

			print("Computation of embedded distances between multiple conditions ")

			T0.mat<-array(0,dim=c(K,K,m))
			TA.mat<-array(0,dim=c(K,K,m))
			for  (g in (1:(K-1))) {
				for (h in (g+1):K){
					T0.mat[g,h,] <-rowSums((Y.0t[((g-1)*m)+(1:m),] -Y.0t[((h-1)*m)+(1:m),])^2)
					TA.mat[g,h,] <-rowSums((Y.At[((g-1)*m)+(1:m),] -Y.At[((h-1)*m)+(1:m),])^2)
					#The test stat with g and h being compared is same as h and g compared
					#Don't really need to symmetrize if we are choosing the max of test stats
					# But in case we choose a different alternative to test against
					T0.mat[h,g,]<-T0.mat[g,h,]
					TA.mat[h,g,]<-TA.mat[g,h,]
				}
			}

			print("Computation of test statistics")

			for (inst.it in 1:m){
				T0[l,inst.it] <- max(T0.mat[,,inst.it])
				TA[l,inst.it] <- max(TA.mat[,,inst.it])
			}


		}
	}

	print(traceback())
	return (list(T0=T0,TA=TA ))


}

#' Title
#'
#' @param D.cond.list
#' @param D.in.oos.list.0
#' @param D.in.oos.list.A
#' @param embed.dim
#' @param hardest.alt
#' @param size
#'
#' @return
#' @export
#'
#' @examples
run.pom.Kcond <- function(D.cond.list,D.in.oos.list.0,D.in.oos.list.A,embed.dim,hardest.alt,size) {


	in.n<-attr(D.cond.list[[1]],"Size")
	oos.n<-dim(D.in.oos.list.0[[1]])[1]-in.n
	d<-embed.dim
	## ==== pom = procrustes o mds ====
	if (oos == TRUE) {
		#Embed in-sample

		embed.joint.config<-list()
		for (k in 1:K){
			embed.joint.config<-c(embed.joint.config,list(smacofM(as.matrix(D.cond.list[[k]]), ndim=d,verbose=FALSE)))
		}
		Y.oos.config.null<-array(0,dim=c(oos.n,d,K))
		Y.oos.config.alt<-array(0,dim=c(oos.n,d,K))

#		# Compute Proc from in-sample embeddings
#		#proc <- procGPA(embed.joint.config,reflect=TRUE)
		Id.mat<-matrix(0,d,d)
		proc<- list(Id.mat)
		for (k in 2:K){
			new.proc <- procrustes(embed.joint.config[[k]],embed.joint.config[[1]],translation=TRUE,dilation=TRUE)
			proc<-c(proc,list(new.proc))
		}
		# Out-of sample embed and Proc Transform dissimilarities

		# TODO: Need to compute scale for each condition
		for (k in 2:K){
			Y.null.k<- oosMDS(D=D.in.oos.list.0[[k]],  X=embed.joint.config[[k]])
			Y.alt.k <- oosMDS(D=D.in.oos.list.A[[k]],  X=embed.joint.config[[k]])

			print("str(Y.null.k)")
			print(str(Y.null.k))


			Y.oos.config.null[,,k]  <- Y.null.k%*%proc[[k]]$R*proc[[k]]$s+(matrix(proc[[k]]$tt,nrow=oos.n,ncol=d,byrow=TRUE))
			Y.oos.config.alt[,,k]  <-  Y.alt.k%*%proc[[k]]$R*proc[[k]]$s+(matrix(proc[[k]]$tt,nrow=oos.n,ncol=d,byrow=TRUE))
		}

		d.btw.pairs.null<- array(0,dim=c(oos.n,K))
		d.btw.pairs.alt <- array(0,dim=c(oos.n,K))
		for (k in 2:K){
			d.btw.pairs.null[,k]<- rowSums((Y.oos.config.null[,,k] - Y.oos.config.null[,,1])^2)
			d.btw.pairs.alt[,k]<- rowSums((Y.oos.config.alt[,,k] - Y.oos.config.alt[,,1])^2)

		}
		T0.pom<-rep(0,oos.n)
		TA.pom<-rep(0,oos.n)

		if (hardest.alt){
			T0.pom <- apply(d.btw.pairs.null,1,max)
			TA.pom <- apply(d.btw.pairs.alt,1,max)
		} else{
			T0.pom <- apply(d.btw.pairs.null,1,sum)
			TA.pom <- apply(d.btw.pairs.alt,1,sum)

		}

		 print("PoM embedding complete\n")

	}
	power.pom.mc <- get_power(T0.pom, TA.pom, size)
	 print("PoM test statistic complete \n")

	return(list(T0=T0.pom,TA=TA.pom,power=power.pom.mc))
}


#' Title
#'
#' @param D.cond.list
#' @param D.in.oos.list.0
#' @param D.in.oos.list.A
#' @param pprime.cond
#' @param d
#' @param hardest.alt
#' @param size
#'
#' @return
#' @export
#'
#' @examples
run.cca.Kcond<-function(D.cond.list,D.in.oos.list.0,D.in.oos.list.A,pprime.cond,d,hardest.alt,size) {

	in.n<-attr(D.cond.list[[1]],"Size")
	oos.n<-dim(D.in.oos.list.0[[1]])[1]-in.n
	#info(logger,paste("values of m,n,d ",m,n,d))
	#info(logger,"str(D.cond.list)")
	#info(logger,str(D.cond.list))
	#info(logger,"D.in.oos.list.0")
	#info(logger,str(D.in.oos.list.0))
	print(str(D.cond.list))
	print(str(D.in.oos.list.0))
	print(str(D.in.oos.list.A))

	## ==== pom = procrustes o mds ====
	if (oos == TRUE) {
		#Embed in-sample
		embed.Config.list<-list()
		embed.joint.config<-list()
		for (k in 1:K){
			embed.joint.config<-c(embed.joint.config,list(smacofM(as.matrix(D.cond.list[[k]]), ndim=pprime.cond[k],verbose=FALSE)))

		}
		Y.oos.config.null<-array(0,dim=c(oos.n,d,K))
		Y.oos.config.alt<-array(0,dim=c(oos.n,d,K))
		print("in.n,oos.n,d")
		print(in.n)
		print(oos.n)
		print(d)

		canon.vecs<-regCCA(embed.joint.config)$eigvecs
		print("str(canon.vecs[[1]])")
		print(str(canon.vecs[[1]]))
#
#		# TODO: Need to compute scale for each condition
		for (k in 1:K){
			Y.null.k<- oosMDS(D=D.in.oos.list.0[[k]],  X=embed.joint.config[[k]])
			print("Y.null.k")
			print(str(Y.null.k))

			Y.alt.k <- oosMDS(D=D.in.oos.list.A[[k]],  X=embed.joint.config[[k]])
			Y.oos.config.null[,,k]  <- Y.null.k%*%canon.vecs[[k]][,1:d]
			Y.oos.config.alt[,,k]  <-  Y.alt.k%*%canon.vecs[[k]][,1:d]
		}
#
		T0.mat<-array(0,dim=c(K,K,oos.n))
		TA.mat<-array(0,dim=c(K,K,oos.n))

		# Computing KxK dist matrix for each K-tuple(matched or unmatched) of projected vectors
		# The transpose function t(.) is used because Y.oos.config...[m.i,,] is d \times K  instead of K \times d
		for (m.i in 1:oos.n){
			T0.mat[,,m.i]<- as.matrix(dist(t(Y.oos.config.null[m.i,,])))
			TA.mat[,,m.i]<- as.matrix(dist(t(Y.oos.config.alt[m.i,,])))
		}

		print("T0.mat")
		print(str(T0.mat))

		T0.cca<-rep(0,oos.n)
		TA.cca<-rep(0,oos.n)
		if (hardest.alt){
			for (inst.it in 1:oos.n){
				T0.cca[inst.it] <- max(T0.mat[,,inst.it])
				TA.cca[inst.it] <- max(TA.mat[,,inst.it])
			}
		} else {
			for (inst.it in 1:oos.n){
				T0.cca[inst.it] <- max(T0.mat[,,inst.it])
				TA.cca[inst.it] <- max(TA.mat[,,inst.it])
			}
		}
#
		print("T0.cca")
		#print(T0.cca)
		print("TA.cca")
		#print(TA.cca)
#		d.btw.pairs.null<- array(0,dim=c(m,K-1))
#		d.btw.pairs.alt <- array(0,dim=c(m,K-1))
#		for (k in 2:K){
#			d.btw.pairs.null[,k]<- (Y.oos.config.null[,,k] - Y.oos.config.null[,,1])^2
#			d.btw.pairs.alt[,k]<- (Y.oos.config.alt[,,k] - Y.oos.config.alt[,,1])^2
#
#		}
#		if (hardest.alt){
#			T0.pom <- apply(d.btw.pairs.null,max,1)
#			TA.pom <- apply(d.btw.pairs.alt,max,1)
#		} else{
#			T0.pom <- apply(d.btw.pairs.null,sum,1)
#			TA.pom <- apply(d.btw.pairs.alt,sum,1)
#
#		}
#
#		if (verbose) print("PoM embedding complete\n")

	}
	power.cca.mc <- get_power(T0.cca, TA.cca, size)
	print("CCA test statistic complete \n")

	return(list(T0=T0.cca,TA=TA.cca,power=power.cca.mc))

}




run.reg.cca.Kcond <-function()
{
	#TODO

	##low-dimensional (regularized) CCA
	## ==== cca ====
	#embed in-sample measurements




}





#Input parameters
# n : numbers of matched K-tuples
# p : signal dimension of matched K-tuples
# q, c, r,K, alpha,sigma.alpha,old.gauss.model.param,sigma.beta=NULL)


#' Title
#'
#' @param n
#' @param p
#' @param q
#' @param c
#' @param r
#' @param K
#' @param alpha
#' @param sigma.alpha
#' @param old.gauss.model.param
#' @param sigma.beta
#'
#' @return
#' @export
#'
#' @examples
matched_rnorm_Kcond<- function(n, p,  q, c, r,K, alpha,sigma.alpha,old.gauss.model.param,sigma.beta=NULL) {

	## Return n K-tuples of matched Normal distributed random vectors, given by
	## X_{ik} ~ (1-c) Norm(alpha[i] ,I/r) + c Norm(0,(1+1/r)I), i = 1, ..., n; k = 1, 2,...,K,
	## where alpha[i]  gaussian distributed common means, Norm(I(1+1/r)) is gaussian
	## noise on R^q

	signals <- array(0, dim=c( n, p,K))
	noise <- array(0, dim=c( n, q,K))


	if (is.null(sigma.beta))
		sigma.beta<- Posdef(p,1/r)
	sigma.eta <- sigma.beta

	print("dims of  alpha,sigma")
	print(dim(alpha))
	print(dim(sigma.beta))
	for (i in 1:n) {

		signals[i, ,] <- t(mvrnorm(K, alpha[i,],sigma.beta))

	}

	noise.sigma.1<- (1+1/r)*diag(q)
	noise.sigma.2<- noise.sigma.1

	for (k in 1:K){
		noise[,,k] <-mvrnorm(n, rep(0,q), noise.sigma.1)
	}

	print("noise dimensions generated")
	signals.wt.noise <-list()
	if (c == 0) {
		return(list(X=signals	,sigma.beta=sigma.beta	)	)

	} else {

		print(abind)
		signals.wt.noise <-abind((1-c)*signals, c*noise,along=2)
		print("combined signals.with.noise")

		return(list(X=signals.wt.noise	,sigma.beta=sigma.beta	)	)
	}

}





matched_rdirichlet_Kcond <- function(n, p, r, q, c,K, alpha) {
	## Return n pairs of matched Dirichlet distribtued random vectors, given by
	## X_{ik} ~ (1-c) Dir(r alpha[i] + 1) + c Dir(1), i = 1, ..., n; k = 1, 2,
	## where alpha are given, Dir(r alpha[i] + 1) is on Delta^p, Dir(1) is uniform
	## noise on Delta^q

	signals <- array(0,dim=c( n, p+1,K))
	noise <- array(0, dim=c( n, q+1,K))

	for (i in 1:n) {
		signals[i, ,] <- t(rdirichlet(K, r*alpha[i,]+1))

	}

	for (k in 1:K){
		noise[,,k] <- rdirichlet(n, rep(1, q+1))
	}


	if (c == 0) {
		return(list(X=signals,sigma.beta=sigma.beta))
	} else {
		signals.wt.noise <-list()
		for (k in 1:K){
			signals.wt.noise <-c(signals.wt.noise, list(cbind((1-c)*signals[,,k], c*noise[,,k])))
		}
		return(list(X=signals.wt.noise	,sigma.beta=sigma.beta	)	)
	}


}







w.val.to.W.mat<-function(w,n,sep.err.w,wt.equalize){
	Weight.Mat<-matrix(1-w,n,n)

	num.pt.pairs<- n/2
	commens.entries <- cbind(1:num.pt.pairs,num.pt.pairs+(1:num.pt.pairs))
	commens.entries <- rbind(commens.entries,cbind(num.pt.pairs+1:num.pt.pairs,1:num.pt.pairs))
	correction.factor <- 1
	if (sep.err.w==FALSE){
		Weight.Mat[1:num.pt.pairs,][,num.pt.pairs+(1:num.pt.pairs)]<- 0
		Weight.Mat[num.pt.pairs+(1:num.pt.pairs),][,(1:num.pt.pairs)]<- 0
		correction.factor<-(1/2)*(n-2)
	}
	else {correction.factor<-(n-2)
	}
	if (wt.equalize==FALSE)
		correction.factor <- 1

	diag(Weight.Mat)<-0
	normalized.w<- w*correction.factor
	weights.sum <- normalized.w+(1-w)

	Weight.Mat[commens.entries]<-normalized.w
	Weight.Mat <- Weight.Mat / weights.sum
	return(Weight.Mat)

}


#JOFC.data.test.stats()

w.val.to.W.mat.Kcond<-function(w,n,K,sep.err.w,wt.equalize){
	Weight.Mat<-matrix(1-w,n,n)

	num.matched.pts<- n/K
	commens.entries <- cbind(1:num.matched.pts,num.matched.pts+(1:num.matched.pts))
	for (k in 1:(K-1))
		for (l in (k+1):(K)){
			commens.entries <-  rbind(commens.entries,cbind(1:num.matched.pts+(k-1)*num.matched.pts,
							1:num.matched.pts+(l-1)*num.matched.pts))
			commens.entries <- rbind(commens.entries,cbind(1:num.matched.pts+(l-1)*num.matched.pts,
							1:num.matched.pts+(k-1)*num.matched.pts))
		}
	correction.factor <- 1

	if (sep.err.w==FALSE){
		for (k in 1:(K-1))
			for (l in (k+1):(K)){
				Weight.Mat[1:num.matched.pts+(k-1)*num.matched.pts,][,1:num.matched.pts+(l-1)*num.matched.pts]<- 0
				Weight.Mat[1:num.matched.pts+(l-1)*num.matched.pts,][,1:num.matched.pts+(k-1)*num.matched.pts]<- 0
			}
		correction.factor<-(n-1)
	}
	else {correction.factor<-(n-1)*K/(K-1)
	}
	if (wt.equalize==FALSE)
		correction.factor <- 1

	diag(Weight.Mat)<-0
	normalized.w<- w*correction.factor
	weights.sum <- normalized.w+(1-w)

	Weight.Mat[commens.entries]<-normalized.w
	Weight.Mat <- Weight.Mat / weights.sum
	return(Weight.Mat)

}




in.embed.MDS.smacof.Kcond <-function(dim.Diss,D,ndimens,K,w.vals,sep.err.w,init.conf,wt.equalize){
	n<-dim.Diss
	smacof.embed<-list()
#		stress.vec<-c()
#		comm.sum.vec<-c()
#		fid1.sum.vec<-c()
#		fid2.sum.vec<-c()
#
#		comm.vec<-c()
#		fid1.vec<-c()
#		fid2.vec<-c()
	#print("Class and Str of Diss matrix D")
	#	print(class(D))
	#	print(str(D))
	#	print(D)

	num.matched.pts<- n/K
	print("in-sample embeddings starting")
	for (w in w.vals){
		Weight.Mat<-w.val.to.W.mat.Kcond(w,n,K,sep.err.w,wt.equalize)
		Weight.Mat[is.na(D)]<-0
		D[is.na(D)] <-1

		new.embed <- smacofM(D,ndimens    ,	W=Weight.Mat        ,
				init    = init.conf,
				verbose = FALSE,
				itmax   = 1000,
				eps     = 1e-6)
		smacof.embed<-c(smacof.embed,list(new.embed ))
#			stress.mat <- (as.dist(D) - dist(new.embed))^2


#			comm.term  <- 0
#			fid.term.1 <-0
#			fid.term.2 <-0
#			for (i in 1:(half.n-1)) {
#				comm.term <- comm.term + (stress.mat [n*(i-1) - i*(i-1)/2 + half.n])
#
#				for (j in (i+1):half.n) {
#					fid.term.1 <- fid.term.1  + (stress.mat [n*(i-1) - i*(i-1)/2 + j-i])
#				}
#			}
#			i <- half.n
#			comm.term <- comm.term + (stress.mat [n*(i-1) - i*(i-1)/2 + half.n])
#			for (i in (half.n+1):(n-1)) {
#				for (j in (i+1):n) {
#					fid.term.2 <- fid.term.2  + (stress.mat [n*(i-1) - i*(i-1)/2 + j-i])
#				}
#			}
#
#			fid1.sum.vec <- c(fid1.sum.vec,fid.term.1)
#			fid2.sum.vec <- c(fid2.sum.vec,fid.term.2)
#			comm.sum.vec <- c(comm.sum.vec,comm.term)
#
		#	stress.mat<-as.dist(Weight.Mat) * stress.mat
		#		stress <- sum(stress.mat)
#			stress.vec<-c(stress.vec,stress)
#			num.fid.terms<-half.n*(half.n-1)/2
#			fid1.vec <- c(fid1.vec,fid.term.1/num.fid.terms)
#			fid2.vec <- c(fid2.vec,fid.term.2/num.fid.terms)
#			comm.vec <- c(comm.vec,comm.term/half.n)

	}
#		FC.ratio   <- (fid1.sum.vec+fid2.sum.vec)/comm.sum.vec
#		FC.ratio.2 <- ((1-w.vals)/w.vals)*(fid1.sum.vec+fid2.sum.vec)/comm.sum.vec
#		FC.ratio.3 <- (fid1.vec + fid2.vec) / comm.vec
	#smacof.embed<-c(smacof.embed),list(stress.vec),list(fid1.vec),list(fid2.vec),list(comm.vec) ,
	#list(fid1.sum.vec),list(fid2.sum.vec),list(comm.sum.vec),list(FC.ratio),list(FC.ratio.2),list(FC.ratio.3))
	#print("length(smacof.embed)")
	#print(length(smacof.embed))
	#print(sum(smacof.embed[[1]]-smacof.embed[[2]]))
	#print(sum(smacof.embed[[1]]-smacof.embed[[length(w.vals)]]))
	print("in-sample embeddings ended")
	return(smacof.embed)
}


diss.matrix<-function(x.config,y.config){
	n<-dim(x.config)[1]
	d<-dim(x.config)[2]
	if (is.list(y.config)){
		m<-dim(y.config[[1]])[1]
	} else{

		m<-dim(y.config)[1]
	}
	K<-dim(x.config)[3]
	x.config.agg<- matrix(0,K*n,d)
	y.config.agg<- matrix(0,K*m,d)


	for (k in 1:K){
		x.config.agg[((k-1)*n)+(1:n),]<-x.config[,,k]
		if (is.list(y.config)){
			y.config.agg[((k-1)*m)+(1:m),]<-y.config[[k]]
		} else{

			y.config.agg[((k-1)*m)+(1:m),]<-y.config[,,k]

		}
		#diss.matrix.list<-c(diss.matrix.list,list(dist(rbind(x.config[,,k],y.config.cond) ) ) )
	}
	return(dist(rbind(x.config.agg,y.config.agg)))
}





get_crit_val<- function(T0,size)
{
	n <- length(T0)
	T0 <- sort(T0)
	return(T0[ceiling(n*(1-size))])
}

get_power <- function(T0, TA, size)
## T0: values of test statistic under H0
## TA: values of test statistic under HA
{
	n <- length(T0)
	m <- length(size)
	T0 <- sort(T0)
	power <- rep(0, m)
	for (i in 1:m) {
		if (size[i] == 0) {
			power[i] <- 0
		} else if(size[i] == 1) {
			power[i] <- 1
		} else {
			power[i] <- sum(TA > T0[ceiling(n*(1-size[i]))]) / n
		}
	}
	power
}

omnibusM <- function(D1, D2, W)
{
	D1 <- as.matrix(D1)
	D2 <- as.matrix(D2)
	W <- as.matrix(W)
	rbind(cbind(D1, W), cbind(t(W), D2))
}

# Returns an omnibus dissimilarity matrix with diagonal block matrices from D.cond.list
# and  off-diagonal block matrices intercondition.diss
omnibusM.Kcond<-function(D.cond.list, intercondition.diss)
{
	K<-length(D.cond.list)
	n<- dim(intercondition.diss)[2]
	omnibusM<- matrix(0,K*n,K*n)

	for (l in 2:(K))
		for (k in (1:(l-1)))  {
			diss.matrix.index <- ((l-1)*(l-2)/2)+k
			omnibusM[(l-1)*n+(1:n),(k-1)*n+(1:n)]<- intercondition.diss[diss.matrix.index,,]
			omnibusM[(k-1)*n+(1:n),(l-1)*n+(1:n)]<- intercondition.diss[diss.matrix.index,,]
		}
	for (l in 1:K) {
		omnibusM[(l-1)*n+(1:n),(l-1)*n+(1:n)]<-as.matrix(D.cond.list[[l]])
	}
	return(omnibusM)

}
#
#	plot.MC.evalues.with.CI<-function(evalues.mc,plot.title,plot.col,conf.int=TRUE,add=FALSE){
#
#
#		num.sims<-dim(evalues.mc)[1]
#		evalue.count <- dim(evalues.mc)[2]
#
#		fp.points <- 1:evalue.count
#		num.plot.points <-evalue.count
#		y.points<-colMeans(evalues.mc,na.rm=TRUE)
#		var.y.points <-rep (0,num.plot.points)
#		valid.sample.count <- rep (0,num.plot.points)
#		for (i in 1:num.sims){
#
#			err.points <- evalues.mc[i,]-y.points
#			err.points <- err.points^2
#			for (j in 1:num.plot.points){
#				if (!is.na(evalues.mc[i,j])){
#					var.y.points[j] <- var.y.points[j] + err.points[j]
#					valid.sample.count[j] <- valid.sample.count[j] + 1
#				}
#			}
#
#		}
#		var.y.points <- var.y.points/valid.sample.count
#		std.y.points <- 2*sqrt(var.y.points)
#		ucl <- y.points+std.y.points
#		lcl <- y.points-std.y.points
#		if (add){
#			lines(x=fp.points,y= y.points,main=plot.title,
#					xlab="Eigenvalues",col=plot.col,xlim=c(0,1),ylim=c(0,1),lwd=2.5)
#
#		}
#		else{
#			plot(x=fp.points,y= y.points,main=plot.title,
#					xlab="Eigenvalues",ylab="",col=plot.col,xlim=c(0,1),ylim=c(0,1),type='l',lwd=2.5)
#		}
#		if (conf.int){
#			arrows(fp.points,ucl,fp.points,lcl,length=.05,angle=90,code=3, lty=3,col=plot.col)
#		}
#		par(lty=3)
#		abline(0,1,col="blue")
#		par(lty=1)
#	}
#
#
#
#	plot.ROC.with.CI<-function(plot.roc.points,plot.title,plot.col,conf.int=TRUE,add=FALSE,fp.points=seq(0,1,0.01),
#			ispowercurve=TRUE,linewd=2.5,xlim=1,ylim=1){
#
#		num.sims<-dim(plot.roc.points)[1]
#		fp.points <- fp.points[fp.points<=xlim]
#		num.x.pts <- length(fp.points)
#		plot.roc.points <- plot.roc.points[,1:num.x.pts]
#		y.points<-colMeans(plot.roc.points,na.rm=TRUE)
#		var.y.points <-rep (0,length(fp.points))
#
#		var.y.points <- colVars(plot.roc.points,na.rm=TRUE)
#		std.y.points <- 2*sqrt(var.y.points)
#		ucl <- y.points+std.y.points
#		lcl <- y.points-std.y.points
#		#if (is.finite(max(y.points))){
#		#	if (max(ucl) < ylim)
#		#		ylim <-  max(y.points)
#		#}
#
#		if (add){
#			lines(x=fp.points,y= y.points,main=plot.title,
#					xlab=expression(alpha),ylab=expression(beta),col=plot.col,xlim=c(0,xlim),ylim=c(0,1),lwd=linewd)
#
#		}
#		else{
#			plot(x=fp.points,y= y.points,main=plot.title,
#					xlab=expression(alpha),ylab=expression(beta),col=plot.col,xlim=c(0,xlim),ylim=c(0,1),type='l',lwd=linewd)
#		}
#		if (conf.int){
#			arrows(fp.points,ucl,fp.points,lcl,length=.05,angle=90,code=3, lty=3,col=plot.col)
#		}
#		par(lty=3)
#		abline(0,1,col="blue")
#		par(lty=1)
#	}
#
#	plot.graph.with.CI<-function(plot.roc.points,plot.title,plot.col,conf.int=TRUE,add=FALSE,fp.points=seq(0,1,0.01),customx.labels=NULL,customy.labels=NULL,ispowercurve=TRUE){
#
#		standardx.axis <- FALSE
#		standardy.axis <- FALSE
#		if (is.null(customx.labels))
#			standardx.axis<-TRUE
#		if (is.null(customy.labels))
#			standardy.axis<-TRUE
#
#		num.sims<-dim(plot.roc.points)[1]
#
#		y.points<-colMeans(plot.roc.points,na.rm=TRUE)
#		var.y.points <-rep (0,length(fp.points))
#		var.y.points <- colVars(plot.roc.points,na.rm=TRUE)
#		std.y.points <- 2*sqrt(var.y.points)
#		ucl <- y.points+std.y.points
#		lcl <- y.points-std.y.points
#
#		if (add){
#			lines(x=fp.points,y= y.points,main=plot.title,
#					col=plot.col,xaxt=ifelse(standardx.axis,"s","n"),
#					yaxt=ifelse(standardy.axis,"s","n"), lwd=2.5,xlab="",ylab="")
#
#		}
#		else{
#			plot(x=fp.points,y= y.points,main=plot.title,xaxt=ifelse(standardx.axis,"s","n"),
#					yaxt=ifelse(standardy.axis,"s","n"), col=plot.col,type='l',lwd=2.5,xlab="",ylab="")
#		}
#
#		if (!standardx.axis)
#			axis(1, at=fp.points,labels=customx.labels)
#		if (!standardy.axis)
#			axis(2, at=y.points,labels=customy.labels)
#
#
#
#		if (conf.int){
#			arrows(fp.points,ucl,fp.points,lcl,length=.05,angle=90,code=3, lty=3,col=plot.col)
#		}
#
#		par(lty=1)
#	}
#
#
#
#
#	get_epsilon_c <- function(X, Y)
#	## Return commensurability error
#	{
#		sum((X - Y)^2) / nrow(X)
#	}
#
#	get_epsilon_f <- function(D, DX)
#	## Return fidelity error
#	{
#		mean((as.dist(D) - as.dist(DX))^2)
#	}
#
#	get_epsilon <- function(D1, D1X, D2, D2X, X1t, X2t)
#	{
#		c(get_epsilon_f(D1, D1X),
#				get_epsilon_f(D2, D2X),
#				get_epsilon_c(X1t, X2t))
#	}
#
#	get_power <- function(T0, TA, size)
#	## T0: values of test statistic under H0
#	## TA: values of test statistic under HA
#	{
#		n <- length(T0)
#		m <- length(size)
#		T0 <- sort(T0)
#		power <- rep(0, m)
#		for (i in 1:m) {
#			if (size[i] == 0) {
#				power[i] <- 0
#			} else if(size[i] == 1) {
#				power[i] <- 1
#			} else {
#				power[i] <- sum(TA > T0[round(n*(1-size[i]))]) / n
#			}
#		}
#		power
#	}
#
#	weight <- function(n, c = 1)
#	## Create the weight matrix for W=diag(0)+NA
#	{
#		rbind(cbind(matrix(1,n,n), c*diag(n)), cbind(c*diag(n), matrix(0,n,n)))
#	}
#
#
#	grassmannian <- function(Q1, Q2) {
#		## Q1 and Q2 are two pxd projection matrices
#		svd(Q1 - Q2)$d[1]
#	}
#
#	## theta <- acos(svd(t(Q1)%*%Q2)$d)
#	##
#	## then geodesic distance is sqrt(sum(theta^2)) (and there are a
#	## boatload of other distances computable from theta).
#	geo_dist <- function(Q1, Q2) {
#		theta <- acos(svd(t(Q1) %*% Q2)$d)
#		sqrt(sum(theta^2))
#		## sum(theta^2)
#	}
#
#
#	Haursdorf_dist <- function(Q1,Q2)
#	{
#		sin(geo_dist(Q1,Q2)/2)
#	}
#
#
	Posdef <- function (dim, maxev = 1)
	## Generating a random positive-definite matrix
	## Eigenvalues are generated from uniform(0, maxev)
	{
		ev = runif(dim-1, 0, maxev)
		ev <- c(ev,maxev)
		Z <- matrix(ncol=dim, rnorm(dim^2))
		decomp <- qr(Z)
		Q <- qr.Q(decomp)
		R <- qr.R(decomp)
		d <- diag(R)
		ph <- d / abs(d)
		O <- Q %*% diag(ph)
		Z <- t(O) %*% diag(ev) %*% O
		return(Z)
	}
#
#	## polarity <- function(X, Xstar)
#	## ## Change the signs of each column of X to best match Xstar
#	## ## in the sum of squared difference sense
#	## ##
#	## ## Return a ncol(X) by ncol(X) diagonal matrix Q, with {1, -1}
#	## ## entries, such that ||XQ - Xstar|| is minimized
#	## {
#	##   d <- ncol(X)
#	##   ss <- rep(0, 2^d)
#	##   diag.entries <- as.matrix(expand.grid(lapply(1:d, function(i) c(1,-1))))
#
#	##   for (i in 1:2^d) {
#	##     ss[i] <- sum((X %*% diag(diag.entries[i, ]) - Xstar)^2)
#	##   }
#
#	##   diag(diag.entries[which.min(ss), ])
#	## }
#
#	polarity <- function(X, Xstar)
#	## Change the signs of each column of X to best match Xstar
#	## in the sum of squared difference sense
#	##
#	## Return a ncol(X) by ncol(X) diagonal matrix Q, with {1, -1}
#	## entries, such that ||XQ - Xstar|| is minimized
#	{
#		d <- ncol(X)
#		diagv <- rep(1L, d)
#		for (i in 1:d) {
#			if (sum((X[, i] - Xstar[, i])^2) > sum((-X[, i] - Xstar[, i])^2))
#				diagv[i] <- -1L
#		}
#		diag(diagv)
#	}
#
#
#	impute_dYX <- function(D, W, dYX, k=3) {
#		## get imputed distances between oos objects Y and within-sample
#		## objects X of different conditions. That is, d(Y2, X1) or d(Y1, X2)
#		##
#		## D = dist(X1) or dist(X2)
#		## W = imputed distance between X1 and X2
#		## dYX = dist(Y1, X1) or dist(Y2, X2)
#		D <- as.matrix(D)
#		W <- as.matrix(W)
#		n <- ncol(D)
#		ord <- t(apply(dYX, 1, function(dyx) order(dyx, rnorm(n))))[, 1:k]
#		imputed.dYX <- t(apply(ord, 1, function(ii) colMeans(W[ii, ])))
#		imputed.dYX
#	}
#
#
#	generateX <- function(alpha, r, q, c) {
#		## Dir(r alpha_i + 1) is on Delta^p, p=ncol(alpha)
#		## Consider uniform noise Dir(1) on Delta^q
#		## X_{ik} ~ (1-c) Dir(r alpha_i + 1) + c Dir(1)
#		n <- nrow(alpha)
#		p <- ncol(alpha)
#		Signal1 <- matrix(0, n, p)
#		Signal2 <- matrix(0, n, p)
#		for (i in 1:n) {
#			Signal1[i, ] <- rdirichlet(1, r*alpha[i, ]+1)
#			Signal2[i, ] <- rdirichlet(1, r*alpha[i, ]+1)
#		}
#		Noise1 <- rdirichlet(n, rep(1, q+1))
#		Noise2 <- rdirichlet(n, rep(1, q+1))
#		if (c == 0) {
#			return(list(X1=Signal2, X2=Signal2))
#		} else {
#			return(list(X1=cbind((1-c)*Signal1, c*Noise1),
#							X2=cbind((1-c)*Signal2, c*Noise2)))
#		}
#	}
#
#
#	colVars <- function(x, na.rm=FALSE, dims=1, unbiased=TRUE, SumSquares=FALSE,
#			twopass=FALSE) {
#		if (SumSquares) return(colSums(x^2, na.rm, dims))
#		N <- colSums(!is.na(x), FALSE, dims)
#		Nm1 <- if (unbiased) N-1 else N
#		if (twopass) {x <- if (dims==length(dim(x))) x - mean(x, na.rm=na.rm) else
#						sweep(x, (dims+1):length(dim(x)), colMeans(x,na.rm,dims))}
#		(colSums(x^2, na.rm, dims) - colSums(x, na.rm, dims)^2/N) / Nm1
#	}
#
#	ThreewayMDS.Embed.Hyp.Test <- function(D1,D2,X1,X2,Y1,Y20,Y2A,model,ndim){
#
#
#		Threeway.Embed<- smacofIndDiff(delta=list(D1,D2), ndim = d, weightmat = NULL, init = NULL, metric = TRUE,
#				ties = "primary", constraint = model, verbose = FALSE, modulus = 1,
#				itmax = 1000, eps = 1e-6)
#
#		X1t <- Threeway.Embed$conf[[1]]
#		X2t <- Threeway.Embed$conf[[2]]
#		Y1t <- oosMDS(dist(rbind(X1, Y1)), X1t)%*% (solve(Threeway.Embed$cweights[[1]]))  # Check row column matchup
#		Y20t <- oosMDS(dist(rbind(X2, Y20)), X2t) %*% (solve(Threeway.Embed$cweights[[2]]))
#		Y2At <- oosMDS(dist(rbind(X2, Y2A)), X2t) %*% (solve(Threeway.Embed$cweights[[2]]))
#		T0 <- rowSums((Y1t - Y20t)^2)
#		TA <- rowSums((Y1t - Y2At)^2)
#		return (list(T0=T0,TA=TA))
#
#	}
#
#	sign.test.cont.table <-function(cont.table.list){
#		suc<-0
#		num.trials <-0
#		for (tab in cont.table.list){
#			if (tab[1,2]>tab[2,1]){
#				suc <- suc + 1
#			}
#			num.trials <- num.trials+1
#		}
#		fail<-length(cont.table.list)-suc
#		return(binom.test(suc,num.trials,p=0.5,alternative="greater"))
#	}
#	sign.rank.sum.test.cont.table <-function(cont.table.list){
#
#		num.trials <- length(cont.table.list)
#		x<-rep(0,num.trials)
#		y<-rep(0,num.trials)
#		i<-1
#		for (tab in cont.table.list){
#			x[i]<-tab[1,2]
#			y[i]<-tab[2,1]
#			i <- i + 1
#
#		}
#
#		return(wilcox.test(x,y,paired=TRUE,alternative="greater"))
#	}
#	binom.out <-function(cont.table.list){
#		num.trials <- length(cont.table.list)
#		binomial.v<-rep(0,num.trials)
#		i<-1
#		for (tab in cont.table.list){
#			if (tab[1,2]>tab[2,1]){
#				binomial.v[i] <- 1
#
#			}
#			i<- i+1
#		}
#		return (binomial.v)
#	}
#
#omnibusM.inoos <- function(D1, D2, W)
#{
#	D1 <- as.matrix(D1)
#	D2 <- as.matrix(D2)
#	W <- as.matrix(W)
#	rbind(cbind(D1, W), cbind(W, D2))
#}
#
