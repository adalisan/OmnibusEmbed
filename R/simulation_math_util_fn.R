## functions
run.mc.replicate<-function(model,p, r, q, c.val,
                           d            =  p-1,
                           pprime1      =  ifelse(model == "gaussian",p+q,p+q+2),   # cca arguments , signal+noise dimension
                           pprime2      =  ifelse(model == "gaussian",p+q,p+q+2),   # cca arguments, signal+noise dimension
                           Wchoice      =  "avg", #How to impute L
                           pre.scaling  =  TRUE,  #Make the measurement spaces have the same scale
                           oos          =  TRUE,  #embed test observations by Out-of-sampling  ?
                           alpha        =  NULL,
                           n  =  100, m  =  100,    #Number of training and test observations

                           old.gauss.model.param = FALSE,
                           separability.entries.w,
                           compare.pom.cca = TRUE,  # Run PoM and CCA to compare with JOFC?
                           oos.use.imputed,
                           level.mcnemar = 0.01,  #At what alpha, should unweighted(w = 0.5) and optimal w^* be compared
                           def.w = 0.5,           #The null hypothesis is that power(def.w)  >= power(rival.w) (by default ,def.w is the w for the unweighted case which equals 0.5)
                           rival.w = NULL,
                           proc.dilation = FALSE, #when investigating convergence of JOFC to PoM, should Procrustes analysis of configurations include the dilation component?
                           assume.matched.for.oos,
                           w.vals,				  #w values to use for JOFC
                           wt.equalize,

                           verbose = FALSE,
                           power.comparison.test = TRUE,
                           cca.reg = FALSE,
                           threewaymds.test=FALSE){

  print(paste("random ",runif(1)))
  print("run.mc.replicate")
  #
  # The followin if statement is Not really necessary, unless we change our mind about rival.w being the best in average in every MC replicate
  # and want to make rival.w a constant (preferably the best overall)

  if (is.null(rival.w)){
    if (0.95 %in% w.vals){
      rival.w = 0.95
    }
    else if (0.9 %in% w.vals){
      rival.w = 0.9
    }
    else if (0.99 %in% w.vals){
      rival.w = 0.99
    }
  }
  w.max.index <- length(w.vals)
  size <- seq(0, 1, 0.01)
  len <- length(size)

  power.w.star <- 0


  power.mc =  array(0,dim = c(w.max.index,len))  #power values for JOFC in this MC replicate
  power.cca.mc  =  array(0,dim = c(len))         #power values for CCA in this MC replicate
  power.pom.mc  =  array(0,dim = c(len))         #power values for PoM in this MC replicate
  power.cca.reg.mc  =  array(0,dim = c(len))     #power values for reg CCA in this MC replicate
  power.indscal.mc  =  array(0,dim = c(len))     #power values for reg CCA in this MC replicate
  power.idioscal.mc  =  array(0,dim = c(len))     #power values for reg CCA in this MC replicate


  config.mismatch <-  list(frob.norm = array(0,dim = c(w.max.index))) #Frob. norm of configuration difference
  #between PoM and JOFC with smallest w
  min.stress.for.w.val  =  array(0,dim = c(w.max.index))   #minimum stress value for  smacof algorithm
  pom.stress <- 0
  auc.meas <-rep(0,w.max.index)

  T0.best.w <- matrix(0,2,m)    #Test statistics for JOFC (comparison of w = 0.5 with optimal w*
  TA.best.w <- matrix(0,2,m)

  cont.table  <- matrix(0,2,2)

  Fid.Err.Term.1 <- array(0,dim = c(w.max.index))
  Fid.Err.Term.2 <- array(0,dim = c(w.max.index))
  Comm.Err.Term <- array(0,dim = c(w.max.index))

  means <- array(0 , dim = c(w.max.index,2*d))

  dissim.list <- generate.dissim(p  =  p,  d  =  d, alpha  =  alpha,
                                 model  =  model, old.gauss.model.param  =  old.gauss.model.param,
                                 r  =  r, n  =  n, m  =  m, mc  =  mc, size  =  size, q  =  q, c.val  =  c.val,
                                 verbose  =  verbose, pre.scaling  =  pre.scaling)

  if (verbose) print("PoM and CCA embedding\n")


  pom.config <-NULL
  regCCA.teststats <- list()
  PoM.teststats<- list()
  CCA.teststats <- list()
  if (compare.pom.cca) {




    PoM.teststats <- with(dissim.list,
                          run.pom(D1, D2, D10A,D20,D2A,

                                  p,q,d,c.val,
                                  n,m,

                                  model,oos,proc.dilation,
                                  size,
                                  verbose)
    )



    pom.config<-PoM.teststats$pom.config
    if (verbose) print(PoM.teststats$T0[1:3])
    if (verbose) print(PoM.teststats$TA[1:3])



    if (verbose) print("PoM test statistic complete\n")
    CCA.teststats <- with(dissim.list,
                          run.cca(D1, D2, D10A,D20,D2A,
                                  p,q,d,c.val,
                                  pprime1,pprime2,

                                  n,m,

                                  model,oos,
                                  size,
                                  verbose)
    )
    if (verbose) print("CCA test statistic complete\n")
    if (cca.reg){
      regCCA.teststats <- with(dissim.list,
                               run.reg.cca(D1, D2, D10A,D20,D2A,
                                           p,q,d,c.val,
                                           pprime1,pprime2,
                                           d.super = floor((d /4 + (p+q*as.numeric(c.val>0)*3/4))),
                                           n,m,

                                           model,oos,
                                           size,
                                           verbose)
      )

    }


  }
  if (threewaymds.test){
    indscal.teststats <- with(dissim.list,
                          run.indscal(D1, D2, D10A,D20,D2A,
                                  d,
                                  m,

                                  "indscal",oos,
                                  size,
                                  verbose)
    )
    idioscal.teststats <- with(dissim.list,
                              run.indscal(D1, D2, D10A,D20,D2A,
                                          d,
                                          m,

                                          "idioscal",oos,
                                          size,
                                          verbose)
    )

  }









  JOFC.results <- with(dissim.list,
                       run.jofc(D1, D2, D10A,D20,D2A,
                                D.oos.1,
                                D.oos.2.null ,
                                D.oos.2.alt ,
                                L.in.oos.0  ,
                                L.in.oos.A ,

                                n,m,
                                d,
                                model,oos,Wchoice,separability.entries.w,wt.equalize,assume.matched.for.oos,oos.use.imputed,
                                pom.config = pom.config,
                                w.vals = w.vals,
                                size,
                                verbose)
  )

  T0<- JOFC.results$T0
  TA<- JOFC.results$TA

  if (verbose) print("JOFC test statistic complete \n")
  for (l in 1:w.max.index){

    power.mcnemar.l <- get_power(T0[l,],TA[l,],level.mcnemar)
    if (power.mcnemar.l>power.w.star){
      rival.w <- w.vals[l]
      power.w.star <- power.mcnemar.l
      w.val.rival.idx <- l
    }
  }


  FidComm.Terms<- list(F1 = JOFC.results$Fid.Err.Term.1,F2 = JOFC.results$Fid.Err.Term.2,C = JOFC.results$Comm.Err.Term)
  FidComm.Sum.Terms <- list(F1 = JOFC.results$Fid.Err.Sum.Term.1,F2 = JOFC.results$Fid.Err.Sum.Term.2,C = JOFC.results$Comm.Err.Sum.Term)


  if(verbose) print(str(FidComm.Terms))
  if(verbose) print("FC.ratio")
  if(verbose) print(str(JOFC.results$FC.ratio))
  if(verbose) print("FC.ratio.2")
  if(verbose) print(str(JOFC.results$FC.ratio.2))
  if(verbose) print("FC.ratio.3")
  if(verbose) print(str(JOFC.results$FC.ratio.3))

  F.to.C.ratio  =  JOFC.results$FC.ratio
  wtF.to.C.ratio = JOFC.results$FC.ratio.2
  F.bar.to.C.bar.ratio =  JOFC.results$FC.ratio.3


  # Power comparison test
  # In order to compare the best w^* vs w = 0.5 in an unbiased way
  # re-run the simulation only for w =  w^* and w = 0.5
  # compute the contingency table using those results
  if (power.comparison.test){
    ## n pairs of matched points
    dissim.list.2<- generate.dissim(p  =  p, d  =  d, alpha  =  alpha,
                                    model  =  model, old.gauss.model.param  =  old.gauss.model.param,
                                    r  =  r, n  =  n, m  =  m, mc  =  mc, size  =  size, q  =  q, c.val  =  c.val,
                                    verbose  =  verbose, pre.scaling  =  pre.scaling)

    JOFC.results <-with(dissim.list.2,run.jofc(D1, D2, D10A,D20,D2A,
                                               D.oos.1,
                                               D.oos.2.null ,
                                               D.oos.2.alt ,

                                               L.in.oos.0 ,
                                               L.in.oos.A ,

                                               n,m,
                                               d,c.val,
                                               model,oos,Wchoice,separability.entries.w,wt.equalize,assume.matched.for.oos,oos.use.imputed,
                                               compare.pom.cca  = TRUE,
                                               w.vals = c(rival.w,def.w),
                                               verbose)
    )


    T0.best.w <- JOFC.results$T0
    TA.best.w <- JOFC.results$TA
    if (verbose) print("JOFC test statistic complete \n")




    w.val.def.idx <- which(w.vals == def.w)
    w.val.rival.idx<- which(w.vals == rival.w)
    crit.value<-get_crit_val(T0.best.w[1,],level.mcnemar)
    crit.value.2<-get_crit_val(T0.best.w[2,],level.mcnemar)
    if (verbose){
      print("crit.values")
      print(crit.value)
      print(crit.value.2)
    }
    cont.table[1,1] <- sum(T0.best.w[1,] <=crit.value & T0.best.w[2,] <=crit.value.2) +
      sum(TA.best.w[1,]>crit.value & TA.best.w[2,]>crit.value.2)
    cont.table[1,2] <- sum(T0.best.w[1,]>crit.value & T0.best.w[2,] <=crit.value.2)  +
      sum(TA.best.w[1,] <=crit.value & TA.best.w[2,]>crit.value.2)
    cont.table[2,1] <- sum(T0.best.w[1,] <=crit.value & T0.best.w[2,]>crit.value.2)  +
      sum(TA.best.w[1,]>crit.value & TA.best.w[2,] <=crit.value.2)
    cont.table[2,2] <- sum(T0.best.w[1,]>crit.value & T0.best.w[2,]>crit.value.2)   +
      sum(TA.best.w[1,] <=crit.value & TA.best.w[2,] <=crit.value.2)
    if (verbose) print("Cont table computed \n")
    if (verbose) print(cont.table)

  }
  for (l in 1:w.max.index){
    power.mc[l, ] <- get_power(T0[l,], TA[l,], size)
	auc.meas[l] <- get_AUC(T0[l,], TA[l,])
  }
  sink(paste("auc.meas",as.integer(round(10000*runif(1,0,1))), ".txt"))
       print(auc.meas)
       n1<-length(Tnull<-T0[1,])
       n2<-length(Talt<-TA[1,])

       k<-outer(Talt,Tnull,'-')
       auc<- (sum(k>0)+0.5*sum(k==0))/(n1*n2)
       print(auc)
       sink()




  print("end run.mc.replicate")
  list(power.mc = power.mc,
       power.cmp = list(cca  =  CCA.teststats$power,
                        pom  =  PoM.teststats$power,
                        regCCA   =  regCCA.teststats$power,
                        indscal  =  indscal.teststats$power,
                        idioscal =  idioscal.teststats$power),
       cont.tables = cont.table,
       config.dist =  config.mismatch, min.stress = c(min.stress.for.w.val,pom.stress),means = means,
       FidComm.Terms = FidComm.Terms,	FidComm.Sum.Terms  =  FidComm.Sum.Terms,
       F.to.C.ratio  =  F.to.C.ratio, wtF.to.C.ratio = wtF.to.C.ratio,
       F.bar.to.C.bar.ratio =  F.bar.to.C.bar.ratio,
       optim.power = dissim.list$optim.power,
	   best.w  = rival.w,
	   auc.meas=auc.meas
  )

}


run.bootstrapped.JOFC<-function(model,p, r, q, c.val,
                                d            =  p-1,
                                pprime1      =  ifelse(model == "gaussian",p+q,p+q+2),   # cca arguments , signal+noise dimension
                                pprime2      =  ifelse(model == "gaussian",p+q,p+q+2),   # cca arguments, signal+noise dimension
                                Wchoice      =  "avg", #How to impute L
                                pre.scaling  =  TRUE,  #Make the measurement spaces have the same scale
                                oos          =  TRUE,  #embed test observations by Out-of-sampling  ?
                                alpha        =  NULL,
                                n  =  300, m  =  100,    #Number of training and test observations
                                nmc = 150,  #resample this many times to estimate w^* (Placeholder, not yet implemented)
                                old.gauss.model.param = FALSE,
                                separability.entries.w,
                                compare.pom.cca = TRUE,  # Run PoM and CCA to compare with JOFC?
                                oos.use.imputed,
                                level.mcnemar = 0.01,  #At what alpha, should unweighted(w = 0.5) and optimal w^* be compared
                                def.w = 0.5,           #The null hypothesis is that power(def.w)  >= power(rival.w) (by default ,def.w is the w for the unweighted case which equals 0.5)
                                rival.w = NULL,
                                proc.dilation = FALSE, #when investigating convergence of JOFC to PoM, should Procrustes analysis of configurations include the dilation component?
                                assume.matched.for.oos,
                                w.vals,				  #w values to use for JOFC
                                wt.equalize,
                                verbose = FALSE,
                                power.comparison.test =  TRUE,cca.reg = FALSE,
                                ...) {
  print("run.bootstrapped.JOFC start")
  w.max.index <- length(w.vals)
  size <- seq(0, 1, 0.01)
  len <- length(size)

  power.w.star <- 0
  repl.n <- 5

  power.mc =  array(0,dim = c(w.max.index,len))  #power values for JOFC in this MC replicate


  T0.best.w <- matrix(0,2,m)    #Test statistics for JOFC (comparison of w = 0.5 with optimal w*
  TA.best.w <- matrix(0,2,m)

  cont.table  <- matrix(0,2,2)

  Fid.Err.Term.1 <- array(0,dim = c(w.max.index))
  Fid.Err.Term.2 <- array(0,dim = c(w.max.index))
  Comm.Err.Term <- array(0,dim = c(w.max.index))

  dissim.list <- generate.dissim(p = p,  d = d, alpha = alpha,
                                 model = model, old.gauss.model.param = old.gauss.model.param,
                                 r = r, n = n, m = m, mc = mc, size = size, q = q, c.val = c.val,
                                 verbose = verbose, pre.scaling = pre.scaling)
  JOFC.results.list<-list()
  for (repl.i in 1:repl.n ){
  bootstrap.parts<- resample.for.param.estim(dissim.list)

  JOFC.results <-  run.jofc(bootstrap.parts$dissim.list.param.est$D1,
								bootstrap.parts$dissim.list.param.est$D2,
								bootstrap.parts$dissim.list.param.est$D10A,
								bootstrap.parts$dissim.list.param.est$D20,
								bootstrap.parts$dissim.list.param.est$D2A,
								bootstrap.parts$dissim.list.param.est$D.oos.1,
								bootstrap.parts$dissim.list.param.est$D.oos.2.null ,
								bootstrap.parts$dissim.list.param.est$D.oos.2.alt ,
								bootstrap.parts$dissim.list.param.est$L.in.oos.0  ,
								bootstrap.parts$dissim.list.param.est$L.in.oos.A ,

								n=bootstrap.parts$dissim.list.param.est$estim.sample.n,
								m=bootstrap.parts$dissim.list.param.est$estim.sample.m,
                                 d,
                                 model,oos,Wchoice,separability.entries.w,wt.equalize,
                                 assume.matched.for.oos,oos.use.imputed,
                                 pom.config = NULL,
                                 w.vals = w.vals,
                                 size,
                                 verbose)
						 JOFC.results.list<- c(   JOFC.results.list ,list( JOFC.results ))
					 }


  best.w.est <- find.best.w.for.power(JOFC.results.list,w.vals,level.mcnemar)

  print("Starting the JOFC test using optimal w^*")
  n.test<- n- (bootstrap.parts$dissim.list.param.est$estim.sample.size)
  bootstrap.parts$dissim.list.test
  rival.w.l<- 1-sqrt(1-best.w.est)
  rival.w.u <- 1- (1-best.w.est)^2
  if (best.w.est>0.925) rival.w <- rival.w.u
  else rival.w <- rival.w.l
  JOFC.results.test <-   run.jofc(bootstrap.parts$dissim.list.test$D1,
		  bootstrap.parts$dissim.list.test$D2,
		  bootstrap.parts$dissim.list.test$D10A,
		  bootstrap.parts$dissim.list.test$D20,
		  bootstrap.parts$dissim.list.test$D2A,
		  bootstrap.parts$dissim.list.test$D.oos.1,
		  bootstrap.parts$dissim.list.test$D.oos.2.null ,
		  bootstrap.parts$dissim.list.test$D.oos.2.alt ,
		  bootstrap.parts$dissim.list.test$L.in.oos.0  ,
		  bootstrap.parts$dissim.list.test$L.in.oos.A ,


                                      n=n,m=m,
                                      d,
                                      model,oos,Wchoice,separability.entries.w,wt.equalize,
                                      assume.matched.for.oos,oos.use.imputed,
                                      pom.config = NULL,
                                      w.vals =  c(best.w.est,rival.w),  #Use estimate from bootstrapping
                                      size,
                                      verbose)



  T0.best.w.test <- JOFC.results.test$T0
  TA.best.w.test <- JOFC.results.test$TA



  crit.value<-get_crit_val(T0.best.w.test[1,],level.mcnemar)
  crit.value.2<-get_crit_val(T0.best.w.test[2,],level.mcnemar)
  if (verbose){
    print("crit.values")
    print(crit.value)
    print(crit.value.2)
  }
  cont.table[1,1] <- sum(T0.best.w.test[1,] <=crit.value & T0.best.w.test[2,] <=crit.value.2) +
    sum(TA.best.w[1,]>crit.value & TA.best.w[2,]>crit.value.2)
  cont.table[1,2] <- sum(T0.best.w.test[1,]>crit.value & T0.best.w.test[2,] <=crit.value.2)  +
    sum(TA.best.w[1,] <=crit.value & TA.best.w[2,]>crit.value.2)
  cont.table[2,1] <- sum(T0.best.w.test[1,] <=crit.value & T0.best.w.test[2,]>crit.value.2)  +
    sum(TA.best.w[1,]>crit.value & TA.best.w[2,] <=crit.value.2)
  cont.table[2,2] <- sum(T0.best.w.test[1,]>crit.value & T0.best.w.test[2,]>crit.value.2)   +
    sum(TA.best.w[1,] <=crit.value & TA.best.w[2,] <=crit.value.2)
  if (verbose) print("Cont table computed \n")
  if (verbose) print(cont.table)


  power.mc <- get_power(T0.best.w.test[1,], TA.best.w.test[1,], size)
  print("run.bootstrapped.JOFC end")
  return(list(power.mc = power.mc, cont.table = cont.table, 	optim.power = dissim.list$optim.power,best.w  = best.w.est ))

}


resample.for.param.estim <- function(dissim.list) {

	#Size of given training set
  n <- nrow(dissim.list$D1)
  total_n <- nrow(dissim.list$D10A)
  # m: number of test pairs
  m <- total_n - n

  estim.sample.size <- floor(n/3)
  estim.sample.n <-  estim.sample.size
  estim.sample.m <- n-estim.sample.n
  sample.for.param.est	<- 1:n
  #estim.sample.n <- floor(2*n/9)
  #estim.sample.m <- estim.sample.size-estim.sample.n
  #separate 1/3 of matched insample for estimation of w^*
	sample.for.param.est.insample <-sample(x=sample.for.param.est, size = estim.sample.n)
  #separate 2/9 of matched insample , (2/3) of sample.for.param.est as insample,
  #the remaining will be used for oos embedding
  #This can be done multiple times, and w^* can be estimated as the one that gives the best power
  #sample.for.param.est.insample <- sample (x=sample.for.param.est, size =  estim.sample.n)

  sample.log <- rep(FALSE,n)
  sample.log[sample.for.param.est.insample]<- TRUE
  # OOS embedding sample indices that will be matched pair set
  oos.sample.0 <- which(!sample.log)
  shuff.sample <- sample( oos.sample.0)
  # OOS embedding sample indices that will be unmatched pair set
  sample.alt <- c(sample.for.param.est.insample,shuff.sample)

  L.in.oos.0.imp <- (dissim.list$D1[sample.for.param.est.insample,oos.sample.0]+
			         dissim.list$D2[sample.for.param.est.insample,oos.sample.0])/2
  L.in.oos.A.imp <- (dissim.list$D1[sample.for.param.est.insample,shuff.sample]+
			         dissim.list$D2[sample.for.param.est.insample,shuff.sample])/2

  L.in.oos.0.param.est <- omnibusM.inoos( L.in.oos.0.imp , L.in.oos.0.imp, matrix(NA,estim.sample.n,estim.sample.m))
  L.in.oos.A.param.est  <- omnibusM.inoos( L.in.oos.A.imp , L.in.oos.A.imp, matrix(NA,estim.sample.n,estim.sample.m))

  #The data for estimating w^*
  data.for.param.est  <- list(D1 = dissim.list$D1[sample.for.param.est.insample,sample.for.param.est.insample],
                              D2 =  dissim.list$D2 [sample.for.param.est.insample,sample.for.param.est.insample] ,
                              D10A = dissim.list$D1[sample.for.param.est,sample.for.param.est] ,
                              D20  =  dissim.list$D2[sample.for.param.est,sample.for.param.est],
                              D2A  =  dissim.list$D2[sample.alt,sample.alt],
                              D.oos.1 =  dissim.list$D1[oos.sample.0,oos.sample.0] ,
                              D.oos.2.null  =  dissim.list$D2[oos.sample.0,oos.sample.0],
                              D.oos.2.alt =   dissim.list$D2[shuff.sample,shuff.sample],
                              L.in.oos.0  =  L.in.oos.0.param.est,
                              L.in.oos.A  =  L.in.oos.A.param.est,
                              estim.sample.size  =   estim.sample.size,
                              estim.sample.n  =  estim.sample.n,
                              estim.sample.m  =  estim.sample.m
  )

  #The data for testing predictor with w*
  data.for.test  <- list(D1 = dissim.list$D1,#[-sample.for.param.est,-sample.for.param.est],
                         D2 =  dissim.list$D2,#[-sample.for.param.est,-sample.for.param.est] ,
                         D10A =  dissim.list$D10A,#[-sample.for.param.est,-sample.for.param.est] ,
                         D20  =  dissim.list$D20,#[-sample.for.param.est,-sample.for.param.est],
                         D2A  =  dissim.list$D2A,#[-sample.for.param.est,-sample.for.param.est],
                         D.oos.1 =  dissim.list$D.oos.1 ,	D.oos.2.null  =  dissim.list$D.oos.2.null,
                         D.oos.2.alt =   dissim.list$D.oos.2.alt,
                         L.in.oos.0  =  dissim.list$L.in.oos.0, L.in.oos.A  =  dissim.list$L.in.oos.A
  )


  diss.list<-list(dissim.list.param.est = data.for.param.est,dissim.list.test = data.for.test)

   #end with

  return(diss.list) #end return


}




find.best.w.for.power <- function(JOFC.res.reps,w.vals.vec,level.for.w.compare = 0.05) {
	#w
	w.vals.best.count <- rep(0,length(w.vals.vec))
	for (JOFC.res  in JOFC.res.reps){
		w.val.len <- length(w.vals.vec)
		beta.w<- rep(0,w.val.len)
		for (l in 1:w.val.len){
			beta.w[l]<-get_power(JOFC.res$T0[l,], JOFC.res$TA[l,], level.for.w.compare)

		}
		best.w.list <- which.max((beta.w))
		w.vals.best.count[best.w.list]<- w.vals.best.count[best.w.list] + 1
	}


	max.w <- which.max(w.vals.best.count)

	if   (sum(w.vals.best.count==w.vals.best.count[max.w])>1){
		auc.sum <-  rep(0,length(w.vals.vec))
		for (JOFC.res  in JOFC.res.reps) {
			w.val.len <- length(w.vals.vec)
			beta.w<- rep(0,w.val.len)
			for (l in 1:w.val.len){

				auc[l]<-sum(outer(JOFC.res$T0[l,], JOFC.res$TA[l,],function(x,y){x<y}))

			}
			auc.sum <- auc.sum+auc


		}

		max.w <- which.max(auc.sim)
	}



	best.w.value <- w.vals.vec [max.w]
	if (verbose) print("Estimate of wstar for average power curve")
	if (verbose)  print(best.w.value)
	return(best.w.value)


}

run.pom <- function(D1, D2, D10A,D20,D2A,
                    p,q,d,c.val,
                    n,m,
                    model,oos,proc.dilation,
                    size,
                    verbose){
  # if (verbose) print (D1[1:5,1:5])
  #	if (verbose) print (D2[1:5,1:5])


  T0.pom <- array(0,dim = c(m))    #Test statistics for PoM under null
  TA.pom <- array(0,dim = c(m))    #Test statistics for JOFC under alternative

  ##  ==  ==  pom  =  procrustes o mds  ==  ==
  if (oos  ==  TRUE) {
    #Embed in-sample
    X1t <- smacofM(D1, ndim = d,verbose = FALSE)
    X2t <- smacofM(D2, ndim = d,verbose = FALSE)
    #      if (verbose) print ("Configs")
    #  		if (verbose) print (X1t)
    #			if (verbose) print (X2t)


    if (verbose) print ("Config. Means")
    #		if (verbose) print (colMeans(X1t))
    #		if (verbose) print (colMeans(X2t))
    # Compute Proc from in-sample embeddings
    proc <- MCMCpack::procrustes(X2t, X1t, dilation = proc.dilation)
    # Out-of sample embed and Proc Transform dissimilarities
    #if (profile.mode)			Rprof("profile-oosMDS.out",append = TRUE)
    Y1t  <- oosMDS(D10A, X1t)
    Y20t <- oosMDS(D20, X2t) %*% proc$R * proc$s
    Y2At <- oosMDS(D2A, X2t) %*% proc$R * proc$s
    #if (profile.mode)			Rprof(NULL)
    X2tp<-X2t %*% proc$R * proc$s
    pom.config<-rbind(X1t,X2tp)
    pom.stress<- sum((as.dist(D1) - dist(X1t))^2)
    pom.stress<- pom.stress+ sum((as.dist(D2) - dist(X2tp))^2)
    if (verbose) print("PoM embedding complete\n")
  } else {
    X1t <- smacofM(D10A,ndim =  d,verbose = FALSE,init = cmdscale(D10A,d))
    D20A <-dist(rbind(X2, Y20, Y2A))
    X2t <- smacofM(D20A,ndim =  d,verbose = FALSE,init = cmdscale(D20A,d))
    center1 <- colMeans(X1t[1:n, ])
    center2 <- colMeans(X2t[1:n, ])
    X1t <- X1t - matrix(center1, n+m, d, byrow = TRUE) # column-center training only
    X2t <- X2t - matrix(center2, n+2*m, d, byrow = TRUE)
    proc <- MCMCpack::procrustes(X2t[1:n, ], X1t[1:n, ], dilation = proc.dilation)
    Y1t <- X1t[(n+1):(n+m), ]
    Y20t <- X2t[(n+1):(n+m), ] %*% proc$R * proc$s
    Y2At <- X2t[(n+m+1):(n+2*m), ] %*% proc$R * proc$s
  }

  T0.pom <- rowSums((Y1t - Y20t)^2)
  TA.pom <- rowSums((Y1t - Y2At)^2)
  power.pom <- get_power(T0.pom, TA.pom, size)
  if (verbose) print("PoM test statistic complete \n")
  return(list(T0 = T0.pom,TA = TA.pom,power = power.pom,pom.config = pom.config))



}

run.cca<-function(D1, D2, D10A,D20,D2A,
                  p,q,d,c.val,
                  pprime1,pprime2,
                  n,m,
                  model,oos,
                  size,
                  verbose){

  T0.cca <- array(0,dim = c(m))     #Test statistics for CCA under null
  TA.cca <- array(0,dim = c(m))		#Test statistics for CCA under alternative
  if (model == "gaussian"){
    pprime1 <- p+q
    pprime2 <- p+q
  }
  else{
    pprime1 <- p+q+2
    pprime2 <- p+q+2

  }
  ##  ==  ==  cca  ==  ==
  #embed in-sample measurements
  if (oos  ==  TRUE) {
    if (c.val == 0){
      if (model == "gaussian"){

        X1t <- smacofM(D1,ndim  =  p,verbose = FALSE)
        X2t <- smacofM(D2,ndim  =  p,verbose = FALSE)
      } else{
        X1t <- smacofM(D1,ndim  =  p+1,verbose = FALSE)
        X2t <- smacofM(D2,ndim  =  p+1,verbose = FALSE)
      }
    } else{
      X1t <- smacofM(D = D1,ndim =  pprime1,verbose = FALSE)
      X2t <- smacofM(D = D2,ndim =  pprime2,verbose = FALSE)
    }

    xcca <- cancor(X1t, X2t)

    #project using projection vectors computed by CCA
    #if (profile.mode)			Rprof("profile-oosMDS.out",append = TRUE)
    feat.dim <- dim(X1t)[2]
    feat.dim.2 <- dim(X2t)[2]

    print (m)
    print (feat.dim)
    Y1t  <- matrix(matrix(oosMDS(D10A, X1t),m,feat.dim) %*% xcca$xcoef,m,feat.dim)[, 1:d]
    Y20t <- matrix(matrix(oosMDS(D20, X2t),m,feat.dim.2) %*% xcca$ycoef,m,feat.dim.2)[, 1:d]
    Y2At <- matrix(matrix(oosMDS(D2A, X2t),m,feat.dim.2) %*% xcca$ycoef,m,feat.dim.2)[, 1:d]
    #if (profile.mode)			Rprof(NULL)
    #cca.config<-rbind(X1t,X2t)

  } else {
    if (c.val == 0){
      if (model == "gaussian"){
        X1t <- smacofM(D10A, ndim = p,verbose = FALSE)
        D20A <-dist(rbind(X2, Y20, Y2A))
        X2t <- smacofM(D20A, ndim = p,verbose = FALSE)
      }
      else{
        X1t <- smacofM(D10A, ndim = p+1,verbose = FALSE)
        D20A <-dist(rbind(X2, Y20, Y2A))
        X2t <- smacofM(D20A, ndim = p+1,verbose = FALSE)
      }
    } else{

      X1t <- smacofM(D10A, ndim = pprime1,verbose = FALSE,init = cmdscale(D10A,pprime1))
      D20A <-dist(rbind(X2, Y20, Y2A))
      X2t <- smacofM(D20A, ndim = pprime2,verbose = FALSE,init = cmdscale(D20A,pprime2))

    }


    if (verbose) print("CCA embedding complete\n")
    center1 <- colMeans(X1t[1:n, ])   # column means of training obs
    center2 <- colMeans(X2t[1:n, ])
    X1t <- X1t - matrix(center1, n+m, pprime1, byrow = TRUE) # column-center training only
    X2t <- X2t - matrix(center2, n+2*m, pprime2, byrow = TRUE)
    cca <- cancor(X1t[1:n, ], X2t[1:n, ])
    Y1t <-  (X1t[(n+1):(n+m), ] %*% cca$xcoef )[, 1:d]
    Y20t <- (X2t[(n+1):(n+m), ] %*% cca$ycoef)[, 1:d]
    Y2At <- (X2t[(n+m+1):(n+2*m), ] %*% cca$ycoef)[, 1:d]
  }
  if (m==1){
  T0.cca <- sum((Y1t - Y20t)^2)
  TA.cca <- sum((Y1t - Y2At)^2)
  } else{

  T0.cca <- rowSums((Y1t - Y20t)^2)
  TA.cca <- rowSums((Y1t - Y2At)^2)
  }
  power.cca.mc <- get_power(T0.cca, TA.cca, size)
  return(list(power = power.cca.mc, T0 = T0.cca, TA = TA.cca	))


}


run.indscal<-function(D1, D2, D10A,D20,D2A,
                  d,
                  m,
                  model,
                  oos,
                  size,
                  verbose){
  require(smacof)

  T0.indscal <- array(0,dim = c(m))     #Test statistics for indscal under null
  TA.indscal <- array(0,dim = c(m))		#Test statistics for indscal under alternative

  ##  ==  ==  indscal  ==  ==
  #embed in-sample measurements
  if (oos  ==  TRUE) {
    indscal.model <- smacof::smacofIndDiff(list(D1,D2),ndim=d,constraint = model)
    X1t= indscal.model$conf[[1]]
    X2t= indscal.model$conf[[2]]
    group2cond_tform_1=indscal.model$cweights[[1]]
    group2cond_tform_2=indscal.model$cweights[[2]]
    #project using projection vectors computed by indscal
    #if (profile.mode)			Rprof("profile-oosMDS.out",append = TRUE)
    Y1t  <- (oosMDS(D10A, X1t) %*% solve(indscal.model$cweights[[1]]))
    Y20t <- (oosMDS(D20, X2t) %*% solve(indscal.model$cweights[[2]]))
    Y2At <- (oosMDS(D2A, X2t) %*% solve(indscal.model$cweights[[2]]))
    #if (profile.mode)			Rprof(NULL)
    #indscal.config<-rbind(X1t,X2t)
    }
  if (m==1){
    T0.indscal <- sum((Y1t - Y20t)^2)
    TA.indscal <- sum((Y1t - Y2At)^2)
  } else{

    T0.indscal <- rowSums((Y1t - Y20t)^2)
    TA.indscal <- rowSums((Y1t - Y2At)^2)
  }
  power.indscal.mc <- get_power(T0.indscal, TA.indscal, size)
  return(list(power = power.indscal.mc, T0 = T0.indscal, TA = TA.indscal	))

}




run.reg.cca<-function(D1, D2, D10A,D20,D2A,
                      p,q,d,c.val,
                      pprime1,pprime2,
                      d.super = floor((d+p+q*as.numeric(c.val>0))/2),
                      n,m,
                      model,oos,
                      size,
                      verbose){



  T0.cca.reg <- array(0,dim = c(m))     #Test statistics for regularized CCA under null
  TA.cca.reg <- array(0,dim = c(m))		#Test statistics for regularized CCA under alternative


  ##low-dimensional (regularized) CCA
  ##  ==  ==  cca  ==  ==
  #embed in-sample measurements
  if (oos  ==  TRUE) {

    if (model == "gaussian"){

      X1t <- smacofM(D1,ndim  = d.super,verbose = FALSE)
      X2t <- smacofM(D2,ndim  = d.super,verbose = FALSE)
    } else{
      X1t <- smacofM(D1,ndim  =  d.super+1,verbose = FALSE)
      X2t <- smacofM(D2,ndim  =  d.super+1,verbose = FALSE)
    }


    xcca <- cancor(X1t, X2t)

    #project using projection vectors computed by CCA
    #if (profile.mode)			Rprof("profile-oosMDS.out",append = TRUE)
    Y1t  <- (oosMDS(D10A, X1t) %*% xcca$xcoef)[, 1:d]
    Y20t <- (oosMDS(D20, X2t) %*% xcca$ycoef)[, 1:d]
    Y2At <- (oosMDS(D2A, X2t) %*% xcca$ycoef)[, 1:d]
    #if (profile.mode)			Rprof(NULL)
    #cca.config<-rbind(X1t,X2t)

  } else {
    if (c.val == 0){
      if (model == "gaussian"){
        X1t <- smacofM(D10A, ndim = d.super,verbose = FALSE)
        D20A <-dist(rbind(X2, Y20, Y2A))
        X2t <- smacofM(D20A, ndim = d.super,verbose = FALSE)
      }
      else{
        X1t <- smacofM(D10A, ndim =  d.super+1,verbose = FALSE)
        D20A <-dist(rbind(X2, Y20, Y2A))
        X2t <- smacofM(D20A, ndim =  d.super+1,verbose = FALSE)


      }
    } else{
      if (model == "gaussian"){
        pprime1 <- d.super
        pprime2 <- d.super
      }
      else{
        pprime1 <- d.super+2
        pprime2 <- d.super+2

      }
      X1t <- smacofM(D10A, ndim = pprime1,verbose = FALSE,init = cmdscale(D10A,pprime1))
      D20A <-dist(rbind(X2, Y20, Y2A))
      X2t <- smacofM(D20A, ndim = pprime2,verbose = FALSE,init = cmdscale(D20A,pprime2))
    }
    if (verbose) print("CCA embedding complete\n")
    center1 <- colMeans(X1t[1:n, ])   # column means of training obs
    center2 <- colMeans(X2t[1:n, ])
    X1t <- X1t - matrix(center1, n+m, pprime1, byrow = TRUE) # column-center training only
    X2t <- X2t - matrix(center2, n+2*m, pprime2, byrow = TRUE)
    cca <- cancor(X1t[1:n, ], X2t[1:n, ])
    Y1t <-  (X1t[(n+1):(n+m), ] %*% cca$xcoef )[, 1:d]
    Y20t <- (X2t[(n+1):(n+m), ] %*% cca$ycoef)[, 1:d]
    Y2At <- (X2t[(n+m+1):(n+2*m), ] %*% cca$ycoef)[, 1:d]
  }
  T0.cca.reg <- rowSums((Y1t - Y20t)^2)
  TA.cca.reg <- rowSums((Y1t - Y2At)^2)
  power.cca.reg.mc <- get_power(T0.cca.reg, TA.cca.reg, size)

  return(list(power = power.cca.reg.mc,T0 = T0.cca.reg,TA = TA.cca.reg	))
}


run.jofc <- function(D1, D2,
		             D10A, D20, D2A,
                     D.oos.1,
                     D.oos.2.null ,
                     D.oos.2.alt ,

                     L.in.oos.0  ,
                     L.in.oos.A ,

                     n,m,
                     d,
                     model,oos,Wchoice="avg",separability.entries.w,wt.equalize,assume.matched.for.oos,oos.use.imputed,

                     pom.config = NULL,
                     w.vals,
                     size,
                     verbose = FALSE)   {

 # if (dim(D1)[1]  !=   n) {print(paste(dim(D1)[1],"is dim of D1:dimension should be",n))}
  #if (dim(D2)[1]  !=   n) {print(paste(dim(D2)[1],"is dim of D2:dimension should be",n))}
 # if (dim(D10A)[1]  !=   (n+m)) {print(paste(dim(D10A)[1],"is dim of D10A:dimension should be",n+m))}
#  if (dim(D20)[1]  !=   (n+m)) {print(paste(dim(D20)[1],"is dim of D20:dimension should be",n+m))}
 # if (dim(D2A)[1]  !=   (n+m)) {print(paste(dim(D2A)[1],"is dim of D2A:dimension should be",n+m))}
  #if (attr(D.oos.1,"Size")  !=   (m))   {print(paste(attr(D.oos.1,"Size"),"is dim of D.oos.1:dimension should be",m))}
  #if (attr(D.oos.2.null,"Size")  !=   (m))  {print(paste(attr(D.oos.2.null,"Size") ,"is dim of D.oos.2.null:dimension should be",m))}
#  if (attr(D.oos.2.alt,"Size") !=   (m)) {print(paste(attr(D.oos.2.alt,"Size"),"is dim of D.oos.2.alt:dimension should be",m))}
#  if ((dim(  L.in.oos.0)[1]  !=   (n)) |  (dim(  L.in.oos.0)[2]  !=   (m))) {
#    print(paste(dim(L.in.oos.0)[1] , dim(  L.in.oos.0)[2],"is dim of D.oos.2.alt:dimension should be",n,"x",m))
#    }

 #   if ((dim(  L.in.oos.A)[1]  !=   (n)) |  (dim(  L.in.oos.A)[2]  !=   (m))) {
 #     print(paste(dim(L.in.oos.A)[1] , dim(  L.in.oos.A)[2],"is dim of D.oos.2.alt:dimension should be",n,"x",m))}
	if (verbose) sink("JOFC.log.txt")
	if (verbose) print("run.jofc run start")
      w.max.index <- length(w.vals)
  T0<-matrix(0,w.max.index,m)
  TA<-matrix(0,w.max.index,m)

  Fid.Err.Term.1 <- c()
  Fid.Err.Term.2 <- c()
  Comm.Err.Term  <- c()

  Fid.Err.Sum.Term.1 <-c()
  Fid.Err.Sum.Term.2 <- c()
  Comm.Err.Sum.Term  <- c()
  FC.ratio    <- c()
  FC.ratio.2  <- c()
  FC.ratio.3  <- c()



  ##  ==  ==  jofc  ==  ==
  L<-matrix(NA,n,n)
  # Impute "between-condition" dissimilarities from different objects
  if (Wchoice  ==  "avg") {
    L <- (D1 + D2)/2
  } else if (Wchoice  ==  "sqrt") {
    L <- sqrt((D1^2 + D2^2)/2)
  } else if (Wchoice  ==  "NA+diag(0)") {
    L <- matrix(NA,n,n)
    diag(L)<- 0
  }


  if (oos  ==  TRUE) {

    #In sample embedding
    # Form omnibus dissimilarity matrix
    M <- omnibusM(D1, D2, L)

    init.conf<-pom.config

    if (verbose) print("M and init.conf")
    if (verbose) print(str(M))
    if (verbose) print(str(init.conf))
    # Embed in-sample using different weight matrices (differentw values)
    X.embeds<-JOFC.Insample.Embed(M,d,w.vals,separability.entries.w,
                                  init.conf = init.conf,wt.equalize = wt.equalize)

    Fid.Err.Term.1 <- X.embeds[[w.max.index+2]]
    Fid.Err.Term.2 <- X.embeds[[w.max.index+3]]
    Comm.Err.Term  <- X.embeds[[w.max.index+4]]

    Fid.Err.Sum.Term.1 <- X.embeds[[w.max.index+5]]
    Fid.Err.Sum.Term.2 <- X.embeds[[w.max.index+6]]
    Comm.Err.Sum.Term  <- X.embeds[[w.max.index+7]]
    FC.ratio    <- X.embeds[[w.max.index+8]]
    FC.ratio.2  <- X.embeds[[w.max.index+9]]
    FC.ratio.3  <- X.embeds[[w.max.index+10]]

    print("Fid.Err.Term.1" )
    print(Fid.Err.Term.1 )
    print("Comm.Err.Term ")
    print(Comm.Err.Term )
    min.stress.for.w.val <- X.embeds[[w.max.index+1]]
    if (verbose) print("JOFC embeddings complete\n")


    #
    # OOS Dissimilarity matrices
    #



    #Imputing dissimilarity  entries for OOS
    if (Wchoice  ==  "avg") {
      L.tilde.null <- (D.oos.1 + D.oos.2.null)/2
      L.tilde.alt  <- (D.oos.1 + D.oos.2.alt)/2
    } else if (Wchoice  ==  "sqrt") {
      L.tilde.null <- sqrt((D.oos.1^2 + D.oos.2.null^2)/2)
      L.tilde.alt  <- sqrt((D.oos.1^2 + D.oos.2.alt^2)/2)

    } else if (Wchoice  ==  "NA+diag(0)") {
      L.tilde.null <- matrix(NA,m,m)
      L.tilde.alt <- matrix(NA,m,m)
      diag(L.tilde.null)<- 0
      diag(L.tilde.alt)<- 0
    }
    #Form OOS omnibus matrices
    M.oos.0 <- omnibusM(D.oos.1,D.oos.2.null, L.tilde.null)
    M.oos.A <- omnibusM(D.oos.1,D.oos.2.alt,  L.tilde.alt)



    for (l in 1:w.max.index){
      if (verbose) print("OOS embedding for JOFC for w =  \n")
      if (verbose) print(w.vals[l])

      w.val.l <- w.vals[l]
      X <- X.embeds[[l]]


      oos.obs.flag<- c(rep(1,2*n),rep(0,2*m))

      #Compute Weight matrix corresponding in-sample  entries
      oos.Weight.mat.1<-w.val.to.W.mat(w.val.l,(2*n),separability.entries.w,wt.equalize)

      #Compute Weight matrix corresponding OOS  entries
      oos.Weight.mat.2<-w.val.to.W.mat(w.val.l,(2*m),separability.entries.w,wt.equalize)

      # If assume.matched.for.oos is true, we assume OOS dissimilarities are matched(in reality,
      # they are matched for the matched pairs, but unmatched for the unmatched pairs)
      # If assume.matched.for.oos is true, we ignore the dissimilarities between matched/unmatched
      # pairs
      if (!assume.matched.for.oos){
        oos.Weight.mat.2[1:m,m+(1:m)]<-0
        oos.Weight.mat.2[m+(1:m),(1:m)]<-0
      }
      # if (oos.use.imputed is true) we treat the dissimiilarities between  in-sample and out-of-sample measurements
      # from different conditions like fidelity terms
      # otherwise they are ignored
      if (oos.use.imputed){
        oos.Weight.mat.w <- matrix(1-w.val.l,2*n,2*m)
      } else{
        oos.Weight.mat.w <- rbind(cbind(matrix(1-w.val.l,n,m), matrix(0,n,m) ),
                                  cbind(matrix(0,n,m),matrix(1-w.val.l,n,m))
        )
      }
      oos.Weight.mat<-omnibusM(oos.Weight.mat.1,oos.Weight.mat.2,oos.Weight.mat.w)
      # Since we are going to oos-embedding, set the weights  of in-sample embedding of stress
      # We are using previous in-sample embeddings, anyway
      oos.Weight.mat[1:(2*n),1:(2*n)]<-0
      if (verbose) print("dim(M.oos.0)")
      if (verbose) print(dim(M.oos.0))
      if (verbose) print("dim(M.oos.A)")
      if (verbose) print(dim(M.oos.A))
      if (verbose) print("dim(oos.Weight.mat)")
      if (verbose) print(dim(oos.Weight.mat))
      if (verbose) print("dim(X)")
      if (verbose) print(dim(X))
      #if (verbose) {print("oos.obs.flag")
      if (verbose)  print(head(colSums(oos.Weight.mat)))



      omnibus.oos.D.0 <- omnibusM(M,M.oos.0,L.in.oos.0 )
      omnibus.oos.D.A <- omnibusM(M,M.oos.A, L.in.oos.A)
      oos.Weight.mat[is.na(omnibus.oos.D.0)]<-0
      omnibus.oos.D.0[is.na(omnibus.oos.D.0)]<-1
      omnibus.oos.D.A[is.na(omnibus.oos.D.A)]<-1
      if (verbose) print("JOFC null omnibus OOS embedding \n")
       print(paste("w is ",w.val.l))
      #if (profile.mode)			Rprof("profile-oosIM.out",append = TRUE)
      Y.0t<-oosIM(D = omnibus.oos.D.0,
                  X = X,
                  init      =  "random",
                  verbose   =  FALSE,
                  itmax     =  1000,
                  eps       =  1e-8,
                  W         =  oos.Weight.mat,
                  isWithin  =  oos.obs.flag,
                  bwOos     =  TRUE)

      if (verbose) print("JOFC alternative omnibus OOS embedding \n")
      Y.At<-oosIM(D = omnibus.oos.D.A,
                  X = X,
                  init      =  "random",
                  verbose   =  FALSE,
                  itmax     =  1000,
                  eps       =  1e-8,
                  W         =  oos.Weight.mat,
                  isWithin  =  oos.obs.flag,
                  bwOos     =  TRUE)
      #if (profile.mode)				Rprof(NULL)
      Y1t<-Y.0t[1:m,]
      Y2t<-Y.0t[m+(1:m),]
      Y1t.A<-Y.At[1:m,]
      Y2At<-Y.At[m+(1:m),]

      if (m==1){
        print(str(Y1t))
        T0[l,] <- sum((Y1t - Y2t)^2)
        TA[l,] <- sum((Y1t.A - Y2At)^2)
      } else {
      T0[l,] <- rowSums((Y1t - Y2t)^2)
      TA[l,] <- rowSums((Y1t.A - Y2At)^2) }
    }
  }

  if (verbose) print("run.jofc  ended")
   if (verbose) sink()
  return(list(T0 = T0,TA = TA,
              Fid.Err.Term.1 = Fid.Err.Term.1 ,		Fid.Err.Term.2 = Fid.Err.Term.2,  Comm.Err.Term = Comm.Err.Term  ,
              Fid.Err.Sum.Term.1  =  Fid.Err.Sum.Term.1 ,		Fid.Err.Sum.Term.2  =  Fid.Err.Sum.Term.2 ,		Comm.Err.Sum.Term   = Comm.Err.Sum.Term  ,
              FC.ratio  =  FC.ratio ,		FC.ratio.2  = FC.ratio.2 ,		FC.ratio.3  =  FC.ratio.3
  )
  )

}



generate.dissim <- function(p,  d, alpha, model, old.gauss.model.param, r, n, m, mc, size, q, c.val, verbose, pre.scaling) {
  sigma <- matrix(0,p,p)


  if (is.null(alpha)) {
    if (model == "gaussian"){
      sigma<- diag(p)
      if (old.gauss.model.param) sigma <-Posdef(p,r)
      alpha.mc <- mvrnorm(n+(2*m), rep(0,p),sigma)
    } else if (model == "dirichlet"){
      alpha.mc <- rdirichlet(n+2*m, rep(1,p+1))
    } else stop("unknown model")


  } else {
    alpha.mc <- alpha[[mc]]
  }
  ## optimal power
  optim.power<- c()
  if (model == "gaussian"){
    for  (aleph in size){
      crit.val.1<-qgamma(aleph,(p)/2,scale = 2/r,lower.tail = FALSE)
      crit.val.2<-crit.val.1
      type.2.err<-pgamma(crit.val.2,shape = (p)/2,scale = 2*(1+1/r))
      beta<- 1-type.2.err
      optim.power<- c(optim.power,beta)
    }
  }


  ## n pairs of matched points
  if (model == "gaussian"){
    xlist <- matched_rnorm(n, p, q, c.val, r, alpha = alpha.mc[1:n, ],
                           sigma.alpha = sigma,old.gauss.model.param = old.gauss.model.param)
  } else{
    xlist <- matched_rdirichlet(n, p, r, q, c.val, alpha.mc[1:n, ])
  }
  X1 <- xlist$X1
  X2 <- xlist$X2
  if (model == "gaussian")
    sigma.mc<-xlist$sigma.beta

  D1 <- dist(X1)
  D2 <- dist(X2)

  if (verbose) print("random matched pairs generated\n")

  #prescaling
  if (pre.scaling) {
    s <- lm(as.vector(D1) ~ as.vector(D2) + 0)$coefficients
  } else {
    s <- 1
  }

  #m pairs of unmatched points
  if (model == "gaussian"){
    ## test observations -- m pairs of matched and m pairs of unmatched
    ylist <- matched_rnorm(m, p, q, c.val, r, alpha = alpha.mc[(n+1):(n+m), ],
                           sigma.alpha = sigma,old.gauss.model.param = old.gauss.model.param, sigma.beta = sigma.mc)
    Y2A <- matched_rnorm(m, p, q, c.val, r, alpha = alpha.mc[(n+m+1):(n+m+m), ],
                         sigma.alpha = sigma,old.gauss.model.param = old.gauss.model.param, sigma.beta = sigma.mc)$X2
  } else{
    ylist <- matched_rdirichlet(m, p, r, q, c.val, alpha.mc[(n+1):(n+m), ])
    Y2A <- matched_rdirichlet(m, p, r, q, c.val, alpha.mc[(n+m+1):(n+m+m), ])$X2
  }
  Y1 <- ylist$X1
  Y20 <- ylist$X2

  # Dissimilarity matrices for in-sample +out-of-sample
  D10A <- as.matrix(dist(rbind(X1, Y1)))

  D20 <- as.matrix(dist(rbind(X2, Y20))) * s
  D2A <- as.matrix(dist(rbind(X2, Y2A))) * s
  D1<-as.matrix(D1)
  D2<-as.matrix(D2)
  pom.config<-c()
  cca.config<-c()
  if (verbose) print(D2[,1:3])
  if (verbose) print("s")
  if (verbose) print(s)

  D2<-D2*s


  D.oos.1<-as.matrix(dist(Y1))
  D.oos.2.null <- as.matrix(dist(Y20))
  D.oos.2.alt <- as.matrix(dist(Y2A))

  ideal.omnibus.0  <- as.matrix(dist(rbind(X1,X2,Y1,Y20)))
  ideal.omnibus.A  <- as.matrix(dist(rbind(X1,X2,Y1,Y2A)))

  L.in.oos.0 <- ideal.omnibus.0[1:(2*n),(2*n)+(1:(2*m))]
  L.in.oos.A <- ideal.omnibus.A[1:(2*n),(2*n)+(1:(2*m))]


  return(list(D1 =  D1,D2 =  D2 ,  D10A =  D10A,D20  =  D20,D2A = D2A,
              D.oos.1 =  D.oos.1, D.oos.2.null =  D.oos.2.null,D.oos.2.alt =  D.oos.2.alt,
              L.in.oos.0  =  L.in.oos.0, L.in.oos.A  =  L.in.oos.A,
              optim.power = optim.power))

}



w.val.to.W.mat<-function(w,n,sep.err.w,wt.equalize){
  Weight.Mat<-matrix(1-w,n,n)

  num.pt.pairs<- n/2
  commens.entries <- cbind(1:num.pt.pairs,num.pt.pairs+(1:num.pt.pairs))
  commens.entries <- rbind(commens.entries,cbind(num.pt.pairs+1:num.pt.pairs,1:num.pt.pairs))
  correction.factor <- 1
  if (sep.err.w == FALSE){
    Weight.Mat[1:num.pt.pairs,num.pt.pairs+(1:num.pt.pairs)]<- 0
    Weight.Mat[num.pt.pairs+(1:num.pt.pairs),(1:num.pt.pairs)]<- 0
    correction.factor<-(1/2)*(n-1)
  }	else {
    correction.factor<-(n-1)
  }
  if (wt.equalize == FALSE)
    correction.factor <- 1

  diag(Weight.Mat)<-0
  normalized.w<- w*correction.factor
  weights.sum <- normalized.w+(1-w)

  Weight.Mat[commens.entries]<-normalized.w
  Weight.Mat <- Weight.Mat / weights.sum
  return(Weight.Mat)

}


#JOFC.data.test.stats()




JOFC.Insample.Embed <-function(D,ndimens,w.vals,sep.err.w,init.conf,wt.equalize){
  #	if (profile.mode) Rprof("JOFC.FC.out",append = TRUE)
  n<- nrow(D)
  smacof.embed<-list()
  stress.vec<-c()
  comm.sum.vec<-c()
  fid1.sum.vec<-c()
  fid2.sum.vec<-c()

  #avg comm and fid errors
  comm.avg.vec<-c()
  fid1.avg.vec<-c()
  fid2.avg.vec<-c()


  half.n<- n/2
  for (w in w.vals){
  sink(paste(w,"insample_embed.txt"))
    Weight.Mat<-w.val.to.W.mat(w,n,sep.err.w,wt.equalize)
    Weight.Mat[is.na(D)]<-0
    D[is.na(D)] <-1
    print(w)
    print(head(Weight.Mat))
    print(head(D))
    new.embed <- smacofM(D,ndim = ndimens    ,	W = Weight.Mat        ,
                         init     =  init.conf,
                         verbose  =  FALSE,
                         itmax    =  1000,
                         eps      =  1e-6)

    sink()
    smacof.embed<-c(smacof.embed,list(new.embed ))
    stress.mat <- (as.dist(D) - dist(new.embed))^2


    comm.term  <- 0
    fid.term.1 <-0
    fid.term.2 <-0
    for (i in 1:(half.n-1)) {
      comm.term <- comm.term + (stress.mat [n*(i-1) - i*(i-1)/2 + half.n])

      for (j in (i+1):half.n) {
        fid.term.1 <- fid.term.1  + (stress.mat [n*(i-1) - i*(i-1)/2 + j-i])
      }
    }
    i <- half.n
    comm.term <- comm.term + (stress.mat [n*(i-1) - i*(i-1)/2 + half.n])
    for (i in (half.n+1):(n-1)) {
      for (j in (i+1):n) {
        fid.term.2 <- fid.term.2  + (stress.mat [n*(i-1) - i*(i-1)/2 + j-i])
      }
    }

    num.fid.terms<-half.n*(half.n-1)/2
    fid1.sum.vec <- c(fid1.sum.vec,fid.term.1)
    fid2.sum.vec <- c(fid2.sum.vec,fid.term.2)
    comm.sum.vec <- c(comm.sum.vec,comm.term)

    stress.mat<-as.dist(Weight.Mat) * stress.mat
    stress <- sum(stress.mat)
    stress.vec<-c(stress.vec,stress)

    fid1.avg.vec <- c(fid1.avg.vec,fid.term.1/num.fid.terms)
    fid2.avg.vec <- c(fid2.avg.vec,fid.term.2/num.fid.terms)
    comm.avg.vec <- c(comm.avg.vec,comm.term/half.n)

  }
  FC.ratio   <- (fid1.sum.vec+fid2.sum.vec)/comm.sum.vec
  FC.ratio.2 <- ((1-w.vals)/w.vals)*(fid1.sum.vec+fid2.sum.vec)/comm.sum.vec
  FC.ratio.3 <- (fid1.avg.vec + fid2.avg.vec) / comm.avg.vec
  smacof.embed<-c(smacof.embed,list(stress.vec),
                  list(fid1.avg.vec),list(fid2.avg.vec),list(comm.avg.vec) ,
                  list(fid1.sum.vec),list(fid2.sum.vec),list(comm.sum.vec),
                  list(FC.ratio),list(FC.ratio.2),list(FC.ratio.3))

  return(smacof.embed)
}


JOFC.Insample.tSNE.Embed <-function(D,ndimens,w.vals,sep.err.w,init.conf,wt.equalize){
  #  if (profile.mode) Rprof("JOFC.FC.out",append = TRUE)
  n<- nrow(D)
  smacof.embed<-list()
  stress.vec<-c()
  comm.sum.vec<-c()
  fid1.sum.vec<-c()
  fid2.sum.vec<-c()

  #avg comm and fid errors
  comm.avg.vec<-c()
  fid1.avg.vec<-c()
  fid2.avg.vec<-c()


  half.n<- n/2
  for (w in w.vals){
    sink(paste("debug_log_",w,"insample_embed.txt"))
    Weight.Mat<-w.val.to.W.mat(w,n,sep.err.w,wt.equalize)
    Weight.Mat[is.na(D)]<-0
    D[is.na(D)] <-2*max(D,na.rm=T)
    print(w)
    print(head(Weight.Mat))
    print(head(D))
    new.embed <- tsne(as.dist(D),k = ndimens ,	W = Weight.Mat        ,
                         initial_config     =  init.conf
                         )

    sink()
    smacof.embed<-c(smacof.embed,list(new.embed ))
    stress.mat <- (as.dist(D) - dist(new.embed))^2


    comm.term  <- 0
    fid.term.1 <-0
    fid.term.2 <-0
    for (i in 1:(half.n-1)) {
      comm.term <- comm.term + (stress.mat [n*(i-1) - i*(i-1)/2 + half.n])

      for (j in (i+1):half.n) {
        fid.term.1 <- fid.term.1  + (stress.mat [n*(i-1) - i*(i-1)/2 + j-i])
      }
    }
    i <- half.n
    comm.term <- comm.term + (stress.mat [n*(i-1) - i*(i-1)/2 + half.n])
    for (i in (half.n+1):(n-1)) {
      for (j in (i+1):n) {
        fid.term.2 <- fid.term.2  + (stress.mat [n*(i-1) - i*(i-1)/2 + j-i])
      }
    }

    num.fid.terms<-half.n*(half.n-1)/2
    fid1.sum.vec <- c(fid1.sum.vec,fid.term.1)
    fid2.sum.vec <- c(fid2.sum.vec,fid.term.2)
    comm.sum.vec <- c(comm.sum.vec,comm.term)

    stress.mat<-as.dist(Weight.Mat) * stress.mat
    stress <- sum(stress.mat)
    stress.vec<-c(stress.vec,stress)

    fid1.avg.vec <- c(fid1.avg.vec,fid.term.1/num.fid.terms)
    fid2.avg.vec <- c(fid2.avg.vec,fid.term.2/num.fid.terms)
    comm.avg.vec <- c(comm.avg.vec,comm.term/half.n)

  }
  FC.ratio   <- (fid1.sum.vec+fid2.sum.vec)/comm.sum.vec
  FC.ratio.2 <- ((1-w.vals)/w.vals)*(fid1.sum.vec+fid2.sum.vec)/comm.sum.vec
  FC.ratio.3 <- (fid1.avg.vec + fid2.avg.vec) / comm.avg.vec
  smacof.embed<-c(smacof.embed,list(stress.vec),
                  list(fid1.avg.vec),list(fid2.avg.vec),list(comm.avg.vec) ,
                  list(fid1.sum.vec),list(fid2.sum.vec),list(comm.sum.vec),
                  list(FC.ratio),list(FC.ratio.2),list(FC.ratio.3))

  return(smacof.embed)
}





matched_rnorm_old_form<- function(n, p,  q, c, r, alpha,sigma.alpha) {
  ## Return n pairs of matched Normal distributed random vectors, given by
  ## X_{ik} ~ (1-c) Norm(alpha[i] ,I) + c Norm(0,SIGMA+I), i  =  1, ..., n; k  =  1, 2,
  ## where alpha[i]  gaussian distributed common means, Norm(0,SIGMA+I)is gaussian
  ## noise on R^q

  signal1 <- matrix(0, n, p)
  signal2 <- matrix(0, n, p)

  sigma.beta<- diag(p)
  sigma.eta <- diag(p)

  for (i in 1:n) {

    signal1[i, ] <- mvrnorm(1, alpha[i,],sigma.beta)
    signal2[i, ] <- mvrnorm(1, alpha[i,],sigma.eta)
  }
  sigma.X <-sigma.alpha+sigma.beta
  sigma.Y <-sigma.alpha+sigma.eta
  #  eig.1 <- eigen(sigma.X)
  #  eig.2 <- eigen(sigma.Y)

  noise.sigma.1<- Posdef(q,
                         max(eigen(sigma.alpha, symmetric = TRUE, only.values  =  TRUE)$values))
  noise.sigma.2<- noise.sigma.1

  noise1 <- mvrnorm(n, rep(0,q), noise.sigma.1)
  noise2 <- mvrnorm(n, rep(0,q), noise.sigma.2)
  if (c  ==  0) {
    return(list(X1 = signal1, X2 = signal2))
  } else {
    return(list(X1 = cbind((1-c)*signal1, c*noise1),
                X2 = cbind((1-c)*signal2, c*noise2)))
  }
}


# Leave one row/column out from each dissimilarity matrix
# compute test statistic for matching.
#
#jofc.simple.leaveoneout<-function(D1,D2,D.2.all,in.sample.ind,w.vals = 0.95){
#
#  D.oos.1 <-D.1.all[!in.sample.ind,!in.sample.ind]
#  D.oos.2 <-D.2.all[!in.sample.ind,!in.sample.ind]
#
#  Wchoice = "NA+diag(0)"
#  separability.entries.w <- FALSE
#  wt.equalize = FALSE
#
#  w.max.index<-length(w.vals)
#  assume.matched.for.oos<-FALSE
#  oos.use.imputed<-FALSE
#  bwOOS<-FALSE
#
#
#  # Impute "between-condition" dissimilarities from different objects
#  if (Wchoice  ==  "avg") {
#    L <- (D1 + D2)/2
#  } else if (Wchoice  ==  "sqrt") {
#    L <- sqrt((D1^2 + D2^2)/2)
#  } else if (Wchoice  ==  "NA+diag(0)") {
#    L <- matrix(NA,n,n)
#    diag(L)<- 0
#  }
#
#
#
#
#  #In sample embedding
#  # Form omnibus dissimilarity matrix
#  M <- omnibusM(D1, D2, L)
#
#  init.conf<-pom.config
#
#  if (verbose) print("M and init.conf")
#  if (verbose) print(str(M))
#  if (verbose) print(str(init.conf))
#  # Embed in-sample using different weight matrices (differentw values)
#  X.embeds<-JOFC.Insample.Embed(M,d,w.vals,separability.entries.w,init.conf = init.conf,wt.equalize = wt.equalize)
#
#
#
#  #
#  # OOS Dissimilarity matrices
#  #
#
#
#
#  #Imputing dissimilarity  entries for OOS
#  if (Wchoice  ==  "avg") {
#    L.tilde.null <- (D.oos.1 + D.oos.2.null)/2
#    L.tilde.alt  <- (D.oos.1 + D.oos.2.alt)/2
#  } else if (Wchoice  ==  "sqrt") {
#    L.tilde.null <- sqrt((D.oos.1^2 + D.oos.2.null^2)/2)
#    L.tilde.alt  <- sqrt((D.oos.1^2 + D.oos.2.alt^2)/2)
#
#  } else if (Wchoice  ==  "NA+diag(0)") {
#    L.tilde.null <- matrix(NA,m,m)
#    L.tilde.alt <- matrix(NA,m,m)
#    diag(L.tilde.null)<- 0
#    diag(L.tilde.alt)<- 0
#  }
#  #Form OOS omnibus matrices
#  M.oos.0 <- omnibusM(D.oos.1,D.oos.2.null, L.tilde.null)
#  M.oos.A <- omnibusM(D.oos.1,D.oos.2.alt,  L.tilde.alt)
#
#
#
#  for (l in 1:w.max.index){
#    if (verbose) print("OOS embedding for JOFC for w =  \n")
#    if (verbose) print(w.vals[l])
#
#    w.val.l <- w.vals[l]
#    X <- X.embeds[[l]]
#
#
#    oos.obs.flag<- c(rep(1,2*n),rep(0,2*m))
#
#    #Compute Weight matrix corresponding in-sample  entries
#    oos.Weight.mat.1<-w.val.to.W.mat(w.val.l,(2*n),separability.entries.w,wt.equalize)
#
#    #Compute Weight matrix corresponding OOS  entries
#    oos.Weight.mat.2<-w.val.to.W.mat(w.val.l,(2*m),separability.entries.w,wt.equalize)
#
#    # If assume.matched.for.oos is true, we assume OOS dissimilarities are matched(in reality,
#    # they are matched for the matched pairs, but unmatched for the unmatched pairs)
#    # If assume.matched.for.oos is true, we ignore the dissimilarities between matched/unmatched
#    # pairs
#    if (!assume.matched.for.oos){
#      oos.Weight.mat.2[1:m,m+(1:m)]<-0
#      oos.Weight.mat.2[m+(1:m),(1:m)]<-0
#    }
#    # if (oos.use.imputed is true) we treat the dissimiilarities between  in-sample and out-of-sample measurements
#    # from different conditions like fidelity terms
#    # otherwise they are ignored
#    if (oos.use.imputed){
#      oos.Weight.mat.w <- matrix(1-w.val.l,2*n,2*m)
#    } else{
#      oos.Weight.mat.w <- rbind(cbind(matrix(1-w.val.l,n,m), matrix(0,n,m) ),
#                                cbind(matrix(0,n,m),matrix(1-w.val.l,n,m))
#      )
#    }
#    oos.Weight.mat<-omnibusM(oos.Weight.mat.1,oos.Weight.mat.2,oos.Weight.mat.w)
#    # Since we are going to oos-embedding, set the weights  of in-sample embedding of stress
#    # We are using previous in-sample embeddings, anyway
#    oos.Weight.mat[1:(2*n),1:(2*n)]<-0
#    if (verbose) print("dim(M.oos.0)")
#    if (verbose) print(dim(M.oos.0))
#    if (verbose) print("dim(M.oos.A)")
#    if (verbose) print(dim(M.oos.A))
#    if (verbose) print("dim(oos.Weight.mat)")
#    if (verbose) print(dim(oos.Weight.mat))
#    if (verbose) print("dim(X)")
#    if (verbose) print(dim(X))
#    #if (verbose) {print("oos.obs.flag")
#
#
#
#    omnibus.oos.D.0 <- omnibusM(M,M.oos.0,ideal.omnibus.0[1:(2*n),(2*n)+(1:(2*m))])
#    omnibus.oos.D.A <- omnibusM(M,M.oos.A, ideal.omnibus.A[1:(2*n),(2*n)+(1:(2*m))])
#    oos.Weight.mat[is.na(omnibus.oos.D.0)]<-0
#    omnibus.oos.D.0[is.na(omnibus.oos.D.0)]<-1
#    omnibus.oos.D.A[is.na(omnibus.oos.D.A)]<-1
#    if (verbose) print("JOFC null omnibus OOS embedding \n")
#    #if (profile.mode)			Rprof("profile-oosIM.out",append = TRUE)
#    Y.0t<-oosIM(D = omnibus.oos.D.0,
#                X = X,
#                init      =  "random",
#                verbose   =  FALSE,
#                itmax     =  1000,
#                eps       =  1e-8,
#                W         =  oos.Weight.mat,
#                isWithin  =  oos.obs.flag,
#                bwOos     =  bwOOS)
#    Y.At<-oosIM(D = omnibus.oos.D.A,
#                X = X,
#                init      =  "random",
#                verbose   =  FALSE,
#                itmax     =  1000,
#                eps       =  1e-8,
#                W         =  oos.Weight.mat,
#                isWithin  =  oos.obs.flag,
#                bwOos     =  bwOOS)
#    #if (profile.mode)				Rprof(NULL)
#    Y1t<-Y.0t[1:m,]
#    Y2t<-Y.0t[m+(1:m),]
#    Y1t.A<-Y.At[1:m,]
#    Y2At<-Y.At[m+(1:m),]
#
#    T0[l,] <- rowSums((Y1t - Y2t)^2)
#    TA[l,] <- rowSums((Y1t.A - Y2At)^2)
#  }
#
#  return(list(T0 = T0,TA = TA))
#
#}



matched_rnorm<- function(n, p,  q, c, r, alpha,sigma.alpha,old.gauss.model.param,sigma.beta = NULL) {
  if (old.gauss.model.param)  return (matched_rnorm_old_form(n,p,q,c,r,alpha,sigma.alpha))
  ## Return n pairs of matched Normal distributed random vectors, given by
  ## X_{ik} ~ (1-c) Norm(alpha[i] ,I/r) + c Norm(0,(1+1/r)I), i  =  1, ..., n; k  =  1, 2,
  ## where alpha[i]  gaussian distributed common means, Norm(I(1+1/r)) is gaussian
  ## noise on R^q

  signal1 <- matrix(0, n, p)
  signal2 <- matrix(0, n, p)

  if (is.null(sigma.beta))
    sigma.beta<- Posdef(p,1/r)
  sigma.eta <- sigma.beta

  for (i in 1:n) {

    signal1[i, ] <- mvrnorm(1, alpha[i,],sigma.beta)
    signal2[i, ] <- mvrnorm(1, alpha[i,],sigma.eta)
  }

  noise.sigma.1<- (1+1/r)*diag(q)
  noise.sigma.2<- noise.sigma.1

  noise1 <- matrix(mvrnorm(n, rep(0,q), noise.sigma.1),n,q)
  noise2 <- matrix(mvrnorm(n, rep(0,q), noise.sigma.2),n,q)
  if (c  ==  0) {
    return(list(X1 = signal1, X2 = signal2,sigma.beta = sigma.beta))
  } else {
    return(list(X1 = cbind((1-c)*signal1, c*noise1),
                X2 = cbind((1-c)*signal2, c*noise2)
                ,sigma.beta = sigma.beta
    )
    )
  }
}





matched_rdirichlet <- function(n, p, r, q, c, alpha) {
  ## Return n pairs of matched Dirichlet distribtued random vectors, given by
  ## X_{ik} ~ (1-c) Dir(r alpha[i] + 1) + c Dir(1), i  =  1, ..., n; k  =  1, 2,
  ## where alpha are given, Dir(r alpha[i] + 1) is on Delta^p, Dir(1) is uniform
  ## noise on Delta^q

  signal1 <- matrix(0, n, p+1)
  signal2 <- matrix(0, n, p+1)
  for (i in 1:n) {
    signal1[i, ] <- rdirichlet(1, r*alpha[i, ]+1)
    signal2[i, ] <- rdirichlet(1, r*alpha[i, ]+1)
  }
  noise1 <- rdirichlet(n, rep(1, q+1))
  noise2 <- rdirichlet(n, rep(1, q+1))
  if (c  ==  0) {
    return(list(X1 = signal1, X2 = signal2))
  } else {
    return(list(X1 = cbind((1-c)*signal1, c*noise1),
                X2 = cbind((1-c)*signal2, c*noise2)))
  }
}


get_crit_val<- function(T0,size)
{
  n <- length(T0)
  T0 <- sort(T0)
  return(T0[ceiling(n*(1-size))])
}

get_reject_rate <- function(T0, crit.val ){

  return(sum(T0>=crit.val)/length(T0)	)
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
    if (size[i]  ==  0) {
      power[i] <- 0
    } else if(size[i]  ==  1) {
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


omnibusM.inoos <- function(D1, D2, W)
{
  D1 <- as.matrix(D1)
  D2 <- as.matrix(D2)
  W <- as.matrix(W)
  rbind(cbind(D1, W), cbind(W, D2))
}

plot.MC.evalues.with.CI<-function(evalues.mc,plot.title,plot.col,conf.int = TRUE,add = FALSE){


  num.sims<-dim(evalues.mc)[1]
  evalue.count <- dim(evalues.mc)[2]

  fp.points <- 1:evalue.count
  num.plot.points <-evalue.count
  y.points<-colMeans(evalues.mc,na.rm = TRUE)
  var.y.points <-rep (0,num.plot.points)
  valid.sample.count <- rep (0,num.plot.points)
  for (i in 1:num.sims){

    err.points <- evalues.mc[i,]-y.points
    err.points <- err.points^2
    for (j in 1:num.plot.points){
      if (!is.na(evalues.mc[i,j])){
        var.y.points[j] <- var.y.points[j] + err.points[j]
        valid.sample.count[j] <- valid.sample.count[j] + 1
      }
    }

  }
  var.y.points <- var.y.points/valid.sample.count
  std.y.points <- 2*sqrt(var.y.points)
  ucl <- y.points+std.y.points
  lcl <- y.points-std.y.points
  if (add){
    lines(x = fp.points,y =  y.points,main = plot.title,
          xlab = "Eigenvalues",col = plot.col,xlim = c(0,1),ylim = c(0,1),lwd = 2.5)

  }
  else{
    plot(x = fp.points,y =  y.points,main = plot.title,
         xlab = "Eigenvalues",ylab = "",col = plot.col,xlim = c(0,1),ylim = c(0,1),type = 'l',lwd = 2.5)
  }
  if (conf.int){
    arrows(fp.points,ucl,fp.points,lcl,length = .05,angle = 90,code = 3, lty = 3,col = plot.col)
  }
  par(lty = 3)
  abline(0,1,col = "blue")
  par(lty = 1)
}



plot.ROC.with.CI<-function(plot.roc.points,plot.title,plot.col,conf.int = TRUE,add = FALSE,fp.points = seq(0,1,0.01),
                           ispowercurve = TRUE,linewd = 2.5,xlim = 1,ylim = 1){

  num.sims<-dim(plot.roc.points)[1]
  fp.points <- fp.points[fp.points <=xlim]
  num.x.pts <- length(fp.points)
  plot.roc.points <- plot.roc.points[,1:num.x.pts]
  y.points<-colMeans(plot.roc.points,na.rm = TRUE)
  var.y.points <-rep (0,length(fp.points))

  var.y.points <- colVars(plot.roc.points,na.rm = TRUE)
  std.y.points <- 2*sqrt(var.y.points)
  ucl <- y.points+std.y.points
  lcl <- y.points-std.y.points
  #if (is.finite(max(y.points))){
  #	if (max(ucl) < ylim)
  #		ylim <-  max(y.points)
  #}

  if (add){
    lines(x = fp.points,y =  y.points,main = plot.title,
          xlab = expression(alpha),ylab = expression(beta),col = plot.col,xlim = c(0,xlim),ylim = c(0,1),lwd = linewd)

  }
  else{
    plot(x = fp.points,y =  y.points,main = plot.title,
         xlab = expression(alpha),ylab = expression(beta),col = plot.col,xlim = c(0,xlim),ylim = c(0,1),type = 'l',lwd = linewd)
  }
  if (conf.int){
    arrows(fp.points,ucl,fp.points,lcl,length = .05,angle = 90,code = 3, lty = 3,col = plot.col)
  }
  par(lty = 3)
  abline(0,1,col = "blue")
  par(lty = 1)
}

plot.graph.with.CI<-function(plot.roc.points,plot.title,plot.col,conf.int = TRUE,add = FALSE,fp.points = seq(0,1,0.01),
                             customx.labels = NULL,customy.labels = NULL,ispowercurve = TRUE,...){

  standardx.axis <- FALSE
  standardy.axis <- FALSE
  if (is.null(customx.labels))
    standardx.axis<-TRUE
  if (is.null(customy.labels))
    standardy.axis<-TRUE

  num.sims<-dim(plot.roc.points)[1]

  y.points<-colMeans(plot.roc.points,na.rm = TRUE)
  var.y.points <-rep (0,length(fp.points))
  var.y.points <- colVars(plot.roc.points,na.rm = TRUE)
  std.y.points <- 2*sqrt(var.y.points)
  ucl <- y.points+std.y.points
  lcl <- y.points-std.y.points

  if (add){
    lines(x = fp.points,y =  y.points,main = plot.title,
          col = plot.col,xaxt = ifelse(standardx.axis,"s","n"),
          yaxt = ifelse(standardy.axis,"s","n"), lwd = 2.5)

  }
  else{
    plot(x = fp.points,y =  y.points,main = plot.title,xaxt = ifelse(standardx.axis,"s","n"),
         yaxt = ifelse(standardy.axis,"s","n"), col = plot.col,type = 'l',lwd = 2.5,...)
  }

  if (!standardx.axis)
    axis(1, at = fp.points,labels = customx.labels)
  if (!standardy.axis)
    axis(2, at = y.points,labels = customy.labels)



  if (conf.int){
    arrows(fp.points,ucl,fp.points,lcl,length = .05,angle = 90,code = 3, lty = 3,col = plot.col)
  }

  par(lty = 1)
}




get_epsilon_c <- function(X, Y)
  ## Return commensurability error
{
  sum((X - Y)^2) / nrow(X)
}

get_epsilon_f <- function(D, DX)
  ## Return fidelity error
{
  mean((as.dist(D) - as.dist(DX))^2)
}

get_epsilon <- function(D1, D1X, D2, D2X, X1t, X2t)
{
  c(get_epsilon_f(D1, D1X),
    get_epsilon_f(D2, D2X),
    get_epsilon_c(X1t, X2t))
}


compute.crit.val <- function(T0.stats,test.level){
  stats.len <- length(T0.stats)
  upper.quant.est<-floor(test.level*stats.len)
  T.crit.val <- sort(T0.stats)[stats.len-upper.quant.est+1]
  return(T.crit.val)

}



weight <- function(n, c  =  1)
  ## Create the weight matrix for W = diag(0)+NA
{
  rbind(cbind(matrix(1,n,n), c*diag(n)), cbind(c*diag(n), matrix(0,n,n)))
}


grassmannian <- function(Q1, Q2) {
  ## Q1 and Q2 are two pxd projection matrices
  svd(Q1 - Q2)$d[1]
}

## theta <- acos(svd(t(Q1)%*%Q2)$d)
##
## then geodesic distance is sqrt(sum(theta^2)) (and there are a
## boatload of other distances computable from theta).
geo_dist <- function(Q1, Q2) {
  theta <- acos(svd(t(Q1) %*% Q2)$d)
  sqrt(sum(theta^2))
  ## sum(theta^2)
}


Haursdorf_dist <- function(Q1,Q2)
{
  sin(geo_dist(Q1,Q2)/2)
}


Posdef <- function (dim, maxev  =  1)
  ## Generating a random positive-definite matrix
  ## Eigenvalues are generated from uniform(0, maxev)
{
  ev  =  runif(dim-1, 0, maxev)
  ev <- c(ev,maxev)
  Z <- matrix(ncol = dim, rnorm(dim^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp)
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

## polarity <- function(X, Xstar)
## ## Change the signs of each column of X to best match Xstar
## ## in the sum of squared difference sense
## ##
## ## Return a ncol(X) by ncol(X) diagonal matrix Q, with {1, -1}
## ## entries, such that ||XQ - Xstar|| is minimized
## {
##   d <- ncol(X)
##   ss <- rep(0, 2^d)
##   diag.entries <- as.matrix(expand.grid(lapply(1:d, function(i) c(1,-1))))

##   for (i in 1:2^d) {
##     ss[i] <- sum((X %*% diag(diag.entries[i, ]) - Xstar)^2)
##   }

##   diag(diag.entries[which.min(ss), ])
## }

polarity <- function(X, Xstar)
  ## Change the signs of each column of X to best match Xstar
  ## in the sum of squared difference sense
  ##
  ## Return a ncol(X) by ncol(X) diagonal matrix Q, with {1, -1}
  ## entries, such that ||XQ - Xstar|| is minimized
{
  d <- ncol(X)
  diagv <- rep(1L, d)
  for (i in 1:d) {
    if (sum((X[, i] - Xstar[, i])^2) > sum((-X[, i] - Xstar[, i])^2))
      diagv[i] <- -1L
  }
  diag(diagv)
}


impute_dYX <- function(D, W, dYX, k = 3) {
  ## get imputed distances between oos objects Y and within-sample
  ## objects X of different conditions. That is, d(Y2, X1) or d(Y1, X2)
  ##
  ## D  =  dist(X1) or dist(X2)
  ## W  =  imputed distance between X1 and X2
  ## dYX  =  dist(Y1, X1) or dist(Y2, X2)
  D <- as.matrix(D)
  W <- as.matrix(W)
  n <- ncol(D)
  ord <- t(apply(dYX, 1, function(dyx) order(dyx, rnorm(n))))[, 1:k]
  imputed.dYX <- t(apply(ord, 1, function(ii) colMeans(W[ii, ])))
  imputed.dYX
}


generateX <- function(alpha, r, q, c) {
  ## Dir(r alpha_i + 1) is on Delta^p, p = ncol(alpha)
  ## Consider uniform noise Dir(1) on Delta^q
  ## X_{ik} ~ (1-c) Dir(r alpha_i + 1) + c Dir(1)
  n <- nrow(alpha)
  p <- ncol(alpha)
  Signal1 <- matrix(0, n, p)
  Signal2 <- matrix(0, n, p)
  for (i in 1:n) {
    Signal1[i, ] <- rdirichlet(1, r*alpha[i, ]+1)
    Signal2[i, ] <- rdirichlet(1, r*alpha[i, ]+1)
  }
  Noise1 <- rdirichlet(n, rep(1, q+1))
  Noise2 <- rdirichlet(n, rep(1, q+1))
  if (c  ==  0) {
    return(list(X1 = Signal2, X2 = Signal2))
  } else {
    return(list(X1 = cbind((1-c)*Signal1, c*Noise1),
                X2 = cbind((1-c)*Signal2, c*Noise2)))
  }
}


colVars <- function(x, na.rm = FALSE, dims = 1, unbiased = TRUE, SumSquares = FALSE,
                    twopass = FALSE) {
  if (SumSquares) return(colSums(x^2, na.rm, dims))
  N <- colSums(!is.na(x), FALSE, dims)
  Nm1 <- if (unbiased) N-1 else N
  if (twopass) {x <- if (dims == length(dim(x))) x - mean(x, na.rm = na.rm) else
    sweep(x, (dims+1):length(dim(x)), colMeans(x,na.rm,dims))}
  (colSums(x^2, na.rm, dims) - colSums(x, na.rm, dims)^2/N) / Nm1
}

ThreewayMDS.Embed.Hyp.Test <- function(D1,D2,X1,X2,Y1,Y20,Y2A,model,ndim){


  Threeway.Embed<- smacofIndDiff(delta = list(D1,D2), ndim  =  d, weightmat  =  NULL, init  =  NULL, metric  =  TRUE,
                                 ties  =  "primary", constraint  =  model, verbose  =  FALSE, modulus  =  1,
                                 itmax  =  1000, eps  =  1e-6)

  X1t <- Threeway.Embed$conf[[1]]
  X2t <- Threeway.Embed$conf[[2]]
  Y1t <- oosMDS(dist(rbind(X1, Y1)), X1t)%*% (solve(Threeway.Embed$cweights[[1]]))  # Check row column matchup
  Y20t <- oosMDS(dist(rbind(X2, Y20)), X2t) %*% (solve(Threeway.Embed$cweights[[2]]))
  Y2At <- oosMDS(dist(rbind(X2, Y2A)), X2t) %*% (solve(Threeway.Embed$cweights[[2]]))
  T0 <- rowSums((Y1t - Y20t)^2)
  TA <- rowSums((Y1t - Y2At)^2)
  return (list(T0 = T0,TA = TA))

}

sign.test.cont.table <-function(cont.table.list){
  suc<-0
  num.trials <-0
  for (tab in cont.table.list){
    if (tab[1,2]>tab[2,1]){
      suc <- suc + 1
    }
    num.trials <- num.trials+1
  }
  fail<-length(cont.table.list)-suc
  return(binom.test(suc,num.trials,p = 0.5,alternative = "greater"))
}
sign.rank.sum.test.cont.table <-function(cont.table.list){

  num.trials <- length(cont.table.list)
  x<-rep(0,num.trials)
  y<-rep(0,num.trials)
  i<-1
  for (tab in cont.table.list){
    x[i]<-tab[1,2]
    y[i]<-tab[2,1]
    i <- i + 1

  }

  return(wilcox.test(x,y,paired = TRUE,alternative = "greater"))
}
binom.out <-function(cont.table.list){
  num.trials <- length(cont.table.list)
  binomial.v<-rep(0,num.trials)
  i<-1
  for (tab in cont.table.list){
    if (tab[1,2]>tab[2,1]){
      binomial.v[i] <- 1

    }
    i<- i+1
  }
  return (binomial.v)
}

The.mode <- function(x, show_all_modes  =  F)
{
  x_freq <- table(x)
  mode_location <- which.max(x_freq)
  The_mode <- names(x_freq)[mode_location]
  Number_of_modes <- length(mode_location)
  #
  if(show_all_modes) {
    if(Number_of_modes >1) {
      warning(paste("Multiple modes exist - returning all",Number_of_modes,"of
									them"))}
		return(The_mode)
    } else {
      if(Number_of_modes >1) {
        warning(paste("Multiple modes exist - returning only the first one out
									of", Number_of_modes))}
		return(The_mode[1])
      }
    }

	get_AUC <- function(T0_stats,TA_stats) {
		n1<-length(T0_stats)
		n2<-length(TA_stats)

		#diff_stats<-outer(TA_stats,T0_stats,'-')
		#auc.estim<- (sum(diff_stats>0)+0.5*sum(diff_stats==0))/(n1*n2)

		library(caTools)
		auc.estim = colAUC(c(T0_stats,TA_stats),rep(c(0,1),each=c(n1,n2)), plotROC=FALSE, alg=c("ROC"))

		return(auc.estim)

	}

get_soft_AUC <- function(T0_stats,TA_stats) {
  n1<-length(T0_stats)
  n2<-length(TA_stats)

  diff_stats<-outer(TA_stats,T0_stats,'-')
  auc.estim<- sum(plogis(diff_stats,scale=2*min(abs(diff_stats))) )/(n1*n2)
  return(auc.estim)

}



## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr,quietly=TRUE)

  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }

  # This is does the summary; it's not easy to understand...
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun= function(xx, col, na.rm) {
                   c( N    = length2(xx[,col], na.rm=na.rm),
                      mean = mean   (xx[,col], na.rm=na.rm),
                      sd   = sd     (xx[,col], na.rm=na.rm)
                   )
                 },
                 measurevar,
                 na.rm
  )

  # Rename the "mean" column
  datac <- rename(datac, c("mean"=measurevar))

  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult

  return(datac)
}




