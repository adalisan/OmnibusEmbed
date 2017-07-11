#' Title
#'
#' @param Diss.E
#' @param Diss.F
#' @param m
#' @param d
#' @param oos
#' @param separability.entries.w
#' @param wt.equalize
#' @param assume.matched.for.oos
#' @param oos.use.imputed
#' @param pom.config
#' @param w.vals
#' @param size
#' @param jack.rep.count
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
run.jacknife<- function( Diss.E, Diss.F,
                         m =2 ,	d,
                         oos ,separability.entries.w,wt.equalize,
                         assume.matched.for.oos,oos.use.imputed,
                         pom.config=NULL,
                         w.vals, 	size,	jack.rep.count=100,verbose=verbose){



  N<-nrow(Diss.E)
  test.samp.size<-m
  test.stats.jk <-list()
  test.stats.jk.CCA <-list()
  test.stats.jk.POM <-list()
  test.stats.jk.INDSCAL <-list()
  test.stats.jk.IDIOSCAL <-list()
  w.val.len<- length(w.vals)
  n= N- (2*m)

  for (jack.i in 1:jack.rep.count){
    left.out.samp<- sample(1:N,2*test.samp.size)
    test.matched<- left.out.samp[1:test.samp.size]
    test.unmatched<- left.out.samp[test.samp.size+(1:test.samp.size)]
    test.matched   <- sort(test.matched)
    test.unmatched <- sort(test.unmatched)

    #sample.for.unmatched<- sample(1:n,test.samp.size)
    orig.indices<-1:N
    #sample.for.unmatched.orig.index<-orig.indices[-left.out.samp][sample.for.unmatched]

    T0 <- matrix(0,w.val.len,m)   #Test statistics for JOFC under null
    TA <- matrix(0,w.val.len,m)    #Test statistics for JOFC under alternative

    D1<-Diss.E[-left.out.samp,-left.out.samp]
    D2<-Diss.F[-left.out.samp,-left.out.samp]
    train.test.0<-orig.indices[-left.out.samp]
    train.test.0<-c(train.test.0,test.matched)
    train.test.A<-orig.indices[-left.out.samp]
    train.test.A<-c(train.test.A,test.unmatched)

    D10A<- Diss.E[train.test.0,train.test.0]
    D20<-  Diss.F[train.test.0,train.test.0]

    D2A<-  Diss.F[train.test.A, train.test.A]

    D.oos.1      <- as.matrix(Diss.E[test.matched,test.matched])
    D.oos.2.null <- as.matrix(Diss.F[test.matched,test.matched])
    D.oos.2.alt  <- as.matrix(Diss.F[test.unmatched,test.unmatched])


    L.in.oos.0 <- omnibusM.inoos(Diss.E[-left.out.samp,][,test.matched],Diss.F[-left.out.samp,][,test.matched],matrix(0,n,m))
    L.in.oos.A <- omnibusM.inoos(Diss.E[-left.out.samp,][,test.matched],Diss.F[-left.out.samp,][,test.unmatched],matrix(0,n,m))



    print(str(T0))
    print(str(TA))

    print(str(D1))
    print(str(D2))
    print(str(D10A))
    print(str(D20))


    power.w.star<- 0
    print("starting JOFC embedding ")


    test.stats.null.alt <- run.jofc(  D1, D2, D10A,D20,D2A,
                                      D.oos.1, D.oos.2.null ,    D.oos.2.alt ,

                                      L.in.oos.0 ,	L.in.oos.A,n, m,	d,
                                      model=NULL, oos, Wchoice="avg", separability.entries.w, wt.equalize,
                                      assume.matched.for.oos,oos.use.imputed,
                                      pom.config=NULL,
                                      w.vals, 	size,	verbose=verbose)

    p.dim  <- 10
    test.stats.null.alt.CCA<- run.cca( D1 = D1, D2 = D2 ,
                                       D10A = D10A, D20= D20, D2A=D2A,
                                       p = p.dim, q = 0, d = d, c.val = 0,
                                       m = m, n= n                          , pprime1 = p.dim, pprime2 = p.dim,
                                       model = "", oos=TRUE, size,verbose = verbose)

    test.stats.null.alt.POM<- run.pom( D1 = D1, D2 = D2 ,
                                       D10A = D10A, D20= D20,
                                       D2A=D2A,
                                       p = p.dim, q = 0,d = d,
                                       c.val = 0, proc.dilation = TRUE,
                                       m = m,
                                       n= n,

                                       model = "", oos=TRUE,
                                       size=size, verbose = verbose)

    test.stats.null.alt.INDSCAL<- run.indscal( D1 = D1, D2 = D2 ,
                                    D10A = D10A, D20= D20,
                                    D2A=D2A,
                                    d=d,
                                    m=m,
                                    model = "indscal", oos=TRUE,
                                    size=size, verbose = verbose)

    test.stats.null.alt.IDIOSCAL<- run.indscal( D1 = D1, D2 = D2 ,
                                                D10A = D10A, D20= D20,
                                                D2A=D2A,
                                                d=d,
                                                m=m,
                                                model = "idioscal", oos=TRUE,
                                                size=size, verbose = verbose)


    if (w.val.len==1){
    test.stats.jk$T0 <- c(test.stats.jk$T0, test.stats.null.alt$T0)
    test.stats.jk$TA <- c(test.stats.jk$TA, test.stats.null.alt$TA)

    test.stats.jk.CCA$T0 <- c(test.stats.jk.CCA$T0, test.stats.null.alt.CCA$T0)
    test.stats.jk.CCA$TA <- c(test.stats.jk.CCA$TA, test.stats.null.alt.CCA$TA)
    test.stats.jk.POM$T0 <- c(test.stats.jk.POM$T0, test.stats.null.alt.POM$T0)
    test.stats.jk.POM$TA <- c(test.stats.jk.POM$TA, test.stats.null.alt.POM$TA)

    test.stats.jk.IDIOSCAL$T0 <- c(test.stats.jk.IDIOSCAL$T0, test.stats.null.alt.IDIOSCAL$T0)
    test.stats.jk.IDIOSCAL$TA <- c(test.stats.jk.IDIOSCAL$TA, test.stats.null.alt.IDIOSCAL$TA)
    test.stats.jk.INDSCAL$T0 <- c(test.stats.jk.INDSCAL$T0, test.stats.null.alt.INDSCAL$T0)
    test.stats.jk.INDSCAL$TA <- c(test.stats.jk.INDSCAL$TA, test.stats.null.alt.INDSCAL$TA)

    } else {
      test.stats.jk$T0 <- cbind(test.stats.jk$T0, test.stats.null.alt$T0)
      test.stats.jk$TA <- cbind(test.stats.jk$TA, test.stats.null.alt$TA)
      test.stats.jk.CCA$T0 <- cbind(test.stats.jk.CCA$T0, test.stats.null.alt.CCA$T0)
      test.stats.jk.CCA$TA <- cbind(test.stats.jk.CCA$TA, test.stats.null.alt.CCA$TA)
      test.stats.jk.POM$T0 <- cbind(test.stats.jk.POM$T0, test.stats.null.alt.POM$T0)
      test.stats.jk.POM$TA <- cbind(test.stats.jk.POM$TA, test.stats.null.alt.POM$TA)
      test.stats.jk.IDIOSCAL$T0 <- cbind(test.stats.jk.IDIOSCAL$T0, test.stats.null.alt.IDIOSCAL$T0)
      test.stats.jk.IDIOSCAL$TA <- cbind(test.stats.jk.IDIOSCAL$TA, test.stats.null.alt.IDIOSCAL$TA)
      test.stats.jk.INDSCAL$T0 <- cbind(test.stats.jk.INDSCAL$T0, test.stats.null.alt.INDSCAL$T0)
      test.stats.jk.INDSCAL$TA <- cbind(test.stats.jk.INDSCAL$TA, test.stats.null.alt.INDSCAL$TA)


    }
  }
  return(list(jofc = test.stats.jk, cca = test.stats.jk.CCA,pom=test.stats.jk.POM,
              indscal=test.stats.jk.INDSCAL, idioscal=test.stats.jk.IDIOSCAL ))

}

#' Title
#'
#' @param m.i
#' @param N
#' @param test.samp.size
#' @param w.val.len
#' @param Diss.E
#' @param Diss.F
#' @param d
#' @param oos
#' @param separability.entries.w
#' @param wt.equalize
#' @param assume.matched.for.oos
#' @param oos.use.imputed
#' @param w.vals
#' @param size
#' @param verbose
#' @param level.critical.val
#'
#' @return
#' @export
#'
#' @examples
run.JOFC.match.jacknife.replicate<- function(m.i, N, test.samp.size, w.val.len, Diss.E, Diss.F
                                             , d,   oos,   separability.entries.w
                                             , wt.equalize, assume.matched.for.oos
                                             , oos.use.imputed, w.vals, size, verbose
                                             , level.critical.val) {




  m<-test.samp.size
  n<- N-2*test.samp.size
  power.mc <-array(0,dim=c(w.val.len,length(size)))
  left.out.samp<- sample(1:N,2*test.samp.size)
  test.matched<- left.out.samp[1:test.samp.size]
  test.unmatched<- left.out.samp[test.samp.size+(1:test.samp.size)]
  test.matched   <- sort(test.matched)
  test.unmatched <- sort(test.unmatched)

  #sample.for.unmatched<- sample(1:n,test.samp.size)
  orig.indices<-1:N
  #sample.for.unmatched.orig.index<-orig.indices[-left.out.samp][sample.for.unmatched]

  T0 <- matrix(0,w.val.len,m)   #Test statistics for JOFC under null
  TA <- matrix(0,w.val.len,m)    #Test statistics for JOFC under alternative

  D1<-Diss.E[-left.out.samp,-left.out.samp]
  D2<-Diss.F[-left.out.samp,-left.out.samp]
  train.test.0<-orig.indices[-left.out.samp]
  train.test.0<-c(train.test.0,test.matched)
  train.test.A<-orig.indices[-left.out.samp]
  train.test.A<-c(train.test.A,test.unmatched)

  D10A<- Diss.E[train.test.0,train.test.0]
  D20<-  Diss.F[train.test.0,train.test.0]

  D2A<-  Diss.F[train.test.A,train.test.A]

  D.oos.1      <- Diss.E[test.matched,test.matched]
  D.oos.2.null <- Diss.F[test.matched,test.matched]
  D.oos.2.alt  <- Diss.F[test.unmatched,test.unmatched]


  L.in.oos.0 <- omnibusM.inoos(Diss.E[-left.out.samp,][,test.matched],Diss.F[-left.out.samp,][,test.matched],matrix(0,n,m))
  L.in.oos.A <- omnibusM.inoos(Diss.E[-left.out.samp,][,test.matched],Diss.F[-left.out.samp,][,test.unmatched],matrix(0,n,m))

  ideal.omnibus.0 <- omnibusM(omnibusM (D1,D2,matrix(0,n,n)),omnibusM(D.oos.1,D.oos.2.null,matrix(0,m,m)),L.in.oos.0)
  ideal.omnibus.A <- omnibusM(omnibusM (D1,D2,matrix(0,n,n)),omnibusM(D.oos.1,D.oos.2.alt,matrix(0,m,m)),L.in.oos.A)


  if (verbose) print(str(T0))
  if (verbose) print(str(TA))

  if (verbose) print(str(D1))
  if (verbose) print(str(D2))
  if (verbose) print(str(D10A))
  if (verbose) print(str(D20))


  power.w.star<- 0
  print("starting JOFC embedding ")

  jacknife.res <- run.jacknife( D1, D2,  m =1 ,	d,
                oos ,separability.entries.w,wt.equalize,
                assume.matched.for.oos,oos.use.imputed,
                pom.config=NULL,
                w.vals, 	size,	jack.rep.count=100,verbose=verbose)
    T.crit.val <- compute.crit.val (jacknife.res$jofc$T0,level.critical.val)
  T.crit.val.cca <- compute.crit.val (jacknife.res$cca$T0,level.critical.val)
  T.crit.val.pom <- compute.crit.val (jacknife.res$pom$T0,level.critical.val)
  T.crit.val.indscal <- compute.crit.val (jacknife.res$indscal$T0,level.critical.val)
  T.crit.val.idioscal <- compute.crit.val (jacknife.res$idioscal$T0,level.critical.val)


    JOFC.results <- run.jofc(  D1, D2, D10A,D20,D2A,
                               D.oos.1, D.oos.2.null ,		D.oos.2.alt ,

                               L.in.oos.0 ,	L.in.oos.A, 	n,m,	d,
                               model=NULL,oos,Wchoice="avg" ,separability.entries.w,wt.equalize,
                               assume.matched.for.oos,oos.use.imputed,
                               pom.config=NULL,
                               w.vals=w.vals, 	size,	verbose=verbose)
  if (verbose) print("JOFC test statistic complete \n")

  CCA.results<- run.cca( D1 = D1, D2 = D2 ,
                                     D10A = D10A, D20= D20, D2A=D2A,
                                     p = 5, q = 0,d = d,c.val = 0,
                                     m = m, n= n , pprime1 = 5, pprime2 = 5,
                                     model = "", oos=TRUE, size=size,verbose = verbose)

  POM.results<- run.pom( D1 = D1, D2 = D2 ,
                                     D10A = D10A, D20= D20,
                                     D2A=D2A,
                                     p = d, q = 0,d = d,c.val = 0, proc.dilation = TRUE,
                                     m = m,
                                     n= n,

                                     model = "", oos=TRUE,
                                     size=size, verbose = verbose)
  INDSCAL.results <- run.indscal( D1 = D1, D2 = D2 ,
                                  D10A = D10A, D20= D20,
                                     D2A=D2A,
                                  d=d,
                                  m=m,
                                  model = "indscal", oos=TRUE,
                                  size=size, verbose = verbose)

  IDIOSCAL.results <- run.indscal( D1 = D1, D2 = D2 ,
                                   D10A = D10A, D20= D20,
                                   D2A=D2A,
                                   d=d,
                                   m=m,
                                   model = "idioscal", oos=TRUE,
                                   size=size, verbose = verbose)


    T0 <- JOFC.results$T0
    TA <- JOFC.results$TA

    power.mc <-array(0,dim=c(w.val.len,length(size)))
    real.size.mc <-array(0,dim=w.val.len)
    power.at.real.size.mc <-array(0,dim=w.val.len)
    for (l in 1:w.val.len){
      w.val.l <- w.vals[l]
      real.size.mc <- get_reject_rate(T0[l,],T.crit.val)
      power.at.real.size.mc <- get_reject_rate(TA[l,],T.crit.val)
      power.l <- get_power(T0[l,],TA[l,],size)
      power.mc[l,]<-power.l
    }


  T0.cca <- CCA.results$T0
  TA.cca <- CCA.results$TA

  power.mc.CCA <-array(0,dim=c(1,length(size)))
  real.size.mc.CCA <-array(0,dim=c(1))
  power.at.real.size.mc.CCA <-array(0,dim=1)


    real.size.mc.CCA <- get_reject_rate(T0.cca,T.crit.val.cca)
    power.at.real.size.mc.CCA <- get_reject_rate(TA.cca,T.crit.val.cca)

    power.mc.CCA<- get_power(T0.cca,TA.cca,size)


  T0.pom<- POM.results$T0
  TA.pom <- POM.results$TA

  power.mc.POM <-array(0,dim=c(1,length(size)))
  real.size.mc.POM <-array(0,dim=c(1))
  power.at.real.size.mc.POM <-array(0,dim=c(1))


    real.size.mc.POM <- get_reject_rate(T0.pom,T.crit.val.pom)
    power.at.real.size.mc.POM <- get_reject_rate(TA.pom,T.crit.val.pom)

    power.mc.POM <- get_power(T0.pom,TA.pom,size)

    T0.indscal <- INDSCAL.results$T0
    TA.indscal <- INDSCAL.results$TA

    power.mc.INDSCAL <-array(0,dim=c(1,length(size)))
    real.size.mc.INDSCAL <-array(0,dim=c(1))
    power.at.real.size.mc.INDSCAL <-array(0,dim=1)


    real.size.mc.INDSCAL <- get_reject_rate(T0.indscal,T.crit.val.indscal)
    power.at.real.size.mc.INDSCAL <- get_reject_rate(TA.indscal,T.crit.val.indscal)

    power.mc.INDSCAL<- get_power(T0.indscal,TA.indscal,size)


    T0.idioscal <- IDIOSCAL.results$T0
    TA.idioscal <- IDIOSCAL.results$TA

    power.mc.IDIOSCAL <-array(0,dim=c(1,length(size)))
    real.size.mc.IDIOSCAL <-array(0,dim=c(1))
    power.at.real.size.mc.IDIOSCAL <-array(0,dim=1)


    real.size.mc.IDIOSCAL <- get_reject_rate(T0.idioscal,T.crit.val.idioscal)
    power.at.real.size.mc.IDIOSCAL <- get_reject_rate(TA.idioscal,T.crit.val.idioscal)

    power.mc.IDIOSCAL<- get_power(T0.idioscal,TA.idioscal,size)




    return(list(T0=T0,TA=TA,power.mc=power.mc,
                real.size.mc=real.size.mc,real.power.mc=power.at.real.size.mc,
                cca=list(T0=T0.cca,TA=TA.cca,power.mc=power.mc.CCA,
                         real.size.mc=real.size.mc.CCA,real.power.mc=power.at.real.size.mc.CCA
                ),
                pom = list(T0=T0.pom,TA=TA.pom,power.mc=power.mc.POM,
                           real.size.mc=real.size.mc.POM,real.power.mc=power.at.real.size.mc.POM
                ),
                idioscal=list(T0=T0.idioscal,TA=TA.idioscal,power.mc=power.mc.IDIOSCAL,
                         real.size.mc=real.size.mc.IDIOSCAL,real.power.mc=power.at.real.size.mc.IDIOSCAL
                ),
                indscal=list(T0=T0.indscal,TA=TA.indscal,power.mc=power.mc.INDSCAL,
                         real.size.mc=real.size.mc.INDSCAL,real.power.mc=power.at.real.size.mc.INDSCAL
                )


                )
           )
  }



