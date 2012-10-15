
genG  <-   function(s  =  seed,n  =  1000,np  =  50,p  =  0.1,q  =  0.5) {
  x  <-   matrix(p, nrow  =  n, ncol  =  n)    ## null: (100-np)x98 matrix
  if (np>0) {
    egg  <-   matrix(q, nrow  =  np, ncol  =  np) ## alternative: npx98 matrix
    x[(n-np+1):n,(n-np+1):n]  <-   egg
  }
  ordering <-  sample.int(n)
  diag(x)  <-   0
  
  A  <-   matrix(0,nrow  =  n,ncol  =  n)
  for(i in 1:(n-1))
    for(j in i:n)
      A[i,j]  <-   A[j,i]  <-   sample(c(0,1),1,prob  =  c(1-x[i,j],x[i,j]))
  
  return(A)
}


perturbG <-  function(G,q){
  n <-  nrow(G)
  Flip.mat <-  matrix(0,nrow  =  n,ncol  =  n)
  for(i in 1:(n-1))
    for(j in i:n){
      Flip.mat[i,j] <-  sample(c(0,1),1,prob  =  c(1-q,q))
      
    }
  for(i in 1:(n-1))
    for(j in i:n){
      if (Flip.mat[i,j] == 1){
        G[i,j] <-  1-G[i,j]
        G[j,i] <-  1-G[j,i]
      }
      
    }
  return(G)
  
}


graph2dissimilarity <- function (G,Gp,
                                      in.sample.ind,
                                      d.dim,
                                      w.vals.vec,
                                 graph.mode,
                                 vert_diss_measure,
                                 T.param,
                                 num_v_to_embed_at_a_time  ,
                                 weighted.g) {
  n <-  nrow(G)
  Graph.1 <-  graph.adjacency(G, mode =  graph.mode,weighted=weighted.g)
  Graph.2 <-  graph.adjacency(Gp,mode  =  graph.mode,weighted=weighted.g)
  #A.M <-   diag(n)
  #A.M[!in.sample.ind[1:n],!in.sample.ind[1:n]] <-   0 # Make sure vertices that are NOT known to be matched
  #are not connected
  
  
  #Now that graphs are generated from adjacency or weight matrices,
  # compute dissimilarities using shortest.paths	
  
  
  #compute dissimilarities in separate graphs, then impute dissimilarity between different condition
  if (vert_diss_measure == 'default'){
    D.1 <-  shortest.paths(Graph.1)
    D.2 <-  shortest.paths(Graph.2)
  } else if (vert_diss_measure == "diffusion"){
    D.1 <-  diff.dist.fun(G,T.param)
    D.2 <-  diff.dist.fun(Gp,T.param)	
  } else if (vert_diss_measure == 'ECT'){
    if (any(!is.finite(G))) {
      print(G(!is.finite(G)))
    }
    D.1 <-   ectime(G)
    if (any(!is.finite(Gp)))
      print(Gp(!is.finite(G)))
    D.2 <-   ectime(Gp)
    
  } else if (vert_diss_measure == 'ell1'){
    D.1 <-   as.matrix(dist(G,'manhattan'))
    D.2 <-   as.matrix(dist(Gp,'manhattan'))
  } else if (vert_diss_measure == 'jaccard'){
    D.1  <-   similarity.jaccard( Graph.1)
    D.2  <-   similarity.jaccard( Graph.2)
    
    D.1  <-   2- 2*D.1
    D.2  <-   2- 2*D.2
  } else if (vert_diss_measure == 'dice'){
    D.1  <-   similarity.dice( Graph.1)
    D.2  <-   similarity.dice( Graph.2)
    D.1  <-   2- 2*D.1
    D.2  <-   2- 2*D.2
  }  else if (vert_diss_measure == 'invlogweighted'){
    D.1  <-   similarity.invlogweighted( Graph.1)
    D.2  <-   similarity.invlogweighted( Graph.2)
    D.1  <-   2- 2*D.1
    D.2  <-   2- 2*D.2
  } else if (vert_diss_measure == 'exp_minus'){
    D.1 <- exp(-1*G/2)
    D.2 <- exp(-1*Gp/2)    
  } else if (vert_diss_measure == 'C_dice_weighted'){
       
    D.1 <- C_dice_weighted(G)
    D.2 <- C_dice_weighted(Gp)
  }
  # In case dissimilarities blow up for diss. measure, put a limit on max 
  # value for diss.
  
  
  
  if (vert_diss_measure == "diffusion"){
  
    upper_diss_limit<-min(quantile(c(as.vector(D.1),as.vector(D.2)),
                                   0.99,na.rm=TRUE),1E3)
    
    D.1[D.1>upper_diss_limit] <- upper_diss_limit
    D.2[D.2>upper_diss_limit] <- upper_diss_limit
  }
  #Should we replace   infinite values of  D.1 and D.2  to large number
  # or ignore them in embeddin?
  #D.1.temp<-D.1
  #D.2.temp<-D.2
  D.1[is.infinite(D.1)] <-  2*max(D.1[is.finite(D.1)])
  D.2[is.infinite(D.2)] <-  2*max(D.2[is.finite(D.2)])
  D.imp <- (D.1+D.2)/2
  
  
  in.sample.ind.half <- in.sample.ind[1:n]
  btw.cond.matched.diss <- rep(0,n)
  btw.cond.matched.diss[!in.sample.ind.half] <- NA
  #in.sample.indices.half<-which(in.sample.ind.half)
  
  D.w<-matrix(NA,n,n)
  D.w.p<-matrix(NA,n,n)
  
  
  #diag(D.w) <- btw.cond.matched.diss
  #diag(D.w.p) <- btw.cond.matched.diss
  
  D.w   [in.sample.ind.half,]   <- D.2[in.sample.ind.half,]
  D.w.p [in.sample.ind.half,]   <- D.1[in.sample.ind.half,]
  
  
  
  D.M <-   rbind(cbind(D.1,D.w),
                 cbind(D.w.p,D.2)
                )
  return(D.M)
}



JOFC.graph.custom.dist  <-   function(G,Gp,
                                      in.sample.ind,
                                      d.dim,
                                      w.vals.vec,
                                      graph.is.directed  =  FALSE,
                                      vert_diss_measure  =  "default",
                                      T.param  =  NULL,
                                      num_v_to_embed_at_a_time   =   sum(!in.sample.ind)/2,
                                      graph.is.weighted=FALSE   )
{
  
  
  
  n <-  nrow(G)
  graph.mode <-   ifelse(graph.is.directed,"directed","undirected")
  weighted.g <- graph.is.weighted
  if (!weighted.g)
    weighted.g<-NULL
 # print("T.param")
 # print(T.param)
  
  
  
  D.M<-graph2dissimilarity (G,Gp,
                            in.sample.ind,
                            d.dim,
                            w.vals.vec,
                            graph.mode,
                            vert_diss_measure,
                            T.param,
                            num_v_to_embed_at_a_time  ,
                            weighted.g)
  
  
  #Given adjacency matrix, generate unweighted graph
  if (!graph.is.weighted)   print("Using adjacency for computing dissimilarities")
  else print("Using weighted for computing dissimilarities") 
  if (isSymmetric(G)) {print("G is symmetric")}
  else   {print("G is unsymmetric")}
  
  if (isSymmetric(Gp)) {print("Gp is symmetric")}
  else   {print("G is unsymmetric")}
  
  #num_v_to_embed_at_a_time   <-   min(floor(1.2*sum(in.sample.ind)),sum(!in.sample.ind))
  #num_v_to_embed_at_a_time   <-   sum(!in.sample.ind)
 # print("num_v_to_embed_at_a_time")
 # print(num_v_to_embed_at_a_time)
#   Embed.List <-  Embed.Nodes(D.M,  in.sample.ind ,oos  =  TRUE ,
#                              d.start  =  d.dim,
#                              wt.equalize =  FALSE,
#                              separability.entries.w  =  FALSE,
#                              assume.matched.for.oos   =   FALSE,
#                              w.vals  =  w.vals.vec,
#                              oos.embed.n.at.a.time  =  num_v_to_embed_at_a_time                             
#   )	
#   
  
  Embed.List <- Embed.Nodes.one.atat (D.M,
                                     in.sample.ind,
                                     oos=TRUE, 
                             d.start  =  d.dim,
                             wt.equalize =  FALSE,
                             separability.entries.w  =  FALSE,
                             assume.matched.for.oos   =   FALSE,
                                     w.vals  =  0.8,
                                     oos.embed.n.at.a.time   =  1,
                                     mds.init.method="gower")
  
  
  #print(str(Embed.List))
  J <-  list()
  for (Y.embed in Embed.List$Y.embeds){
    
    test.samp.size <-  nrow(Y.embed)/2
    Dist  <-  as.matrix(dist(Y.embed))[1:test.samp.size,(1:test.samp.size)+test.samp.size]
    sink("debug.txt")
    print(Y.embed)
    print(Dist)
    sink()
    J <-  c(J,list(Dist))
    
  }
  return(J)
  
  
}


JOFC.graph <-  function(G,Gp,
                        in.sample.ind,
                        d.dim,
                        w.vals.vec,
                        graph.is.directed  =  FALSE){
  
  return(JOFC.graph.custom.dist(G  =  G,Gp  =  Gp,in.sample.ind  =  in.sample.ind,d.dim  =  d.dim,w.vals.vec  =  w.vals.vec,
                                graph.is.directed  =  graph.is.directed,vert_diss_measure  =  "default"))
  
}

JOFC.graph.diff <-  function(G,Gp,
                             in.sample.ind,
                             d.dim,
                             w.vals.vec,
                             graph.is.directed  =  FALSE,
                             T.param  =  2){
  return(JOFC.graph.custom.dist(G  =  G,Gp  =  Gp,in.sample.ind  =  in.sample.ind,d.dim  =  d.dim,w.vals.vec  =  w.vals.vec,
                                graph.is.directed  =  graph.is.directed,vert_diss_measure  =  "diffusion",T.param  =  T.param))
  
}



JOFC.graph.with.opts <-  function(G,Gp,
                  in.sample.ind,
                  d.dim,
                  w.vals.vec,
                  graph.is.directed  =  FALSE,
                  oos  =  TRUE,
                  notconnect.wt  =  10,
                  use.weighted.graph  =  TRUE,
                  wt.matrix.1  =  NULL,
                  wt.matrix.2  =  NULL,
                  sep.graphs  =  TRUE, # if TRUE, treat two graphs separately to compute dissimilarities
                  #and impute W (off-diagonalblock matrix)
                  # if FALSE, join the graphs and compute dissimilarities from joint graph
                  matched.cost  =  0.01
){
  #obsolete
  print("use other JOFC.graph functions")
  return()
  n <-  nrow(G)
  graph.mode <-   ifelse(graph.is.directed,"directed","undirected")
  
  Graph.1 <-  graph.empty(directed  =  graph.is.directed)
  Graph.2 <-  graph.empty(directed  =  graph.is.directed)
  Graph.M <-  graph.empty(directed =  graph.is.directed)
  
  if (!use.weighted.graph) {
    #Given adjacency matrix, generate unweighted graph
    print("Using adjacency for computing dissimilarities")
    Graph.1 <-  graph.adjacency(G, mode  =  graph.mode)
    Graph.2 <-  graph.adjacency(Gp,mode  =  graph.mode)
    A.M =   diag(n)
    A.M[!in.sample.ind[1:n],!in.sample.ind[1:n]] <-   0 # Make sure vertices that are NOT known to be matched
    #are not connected
    G.comb <-  omnibusM(G,Gp,A.M)
    Graph.M  <-   graph.adjacency(G.comb,
                                  weighted  =   NULL ,mode  =  graph.mode)
  } else{
    print("Using weighted graph for computing dissimilarities")
    # If creating  a weighted graph 
    # make the weight matrix from adjacency matrix. Those with same-condition edges have weights of wt.connect/10
    #those  with no edges have weights of wt.connect. Those with "matched edges has weights of matched.cost
    A.M <-   matrix(1E8,n,n)
    
    diag(A.M)  <-   matched.cost
    
    A.M[!in.sample.ind[1:n],!in.sample.ind[1:n]] <-   0
    if (is.null(wt.matrix.1)){
      #Given adjacency matrix, generate weighted graph
      wt.matrix.1  <-   G
      wt.matrix.2  <-   Gp
      wt.matrix.1[G == 0]  <-   notconnect.wt
      wt.matrix.2[Gp == 0]  <-   notconnect.wt						
      wt.matrix.1[G == 1]  <-   notconnect.wt/10			
      wt.matrix.2[Gp == 1]  <-   notconnect.wt/10
      G.comb.w <-  omnibusM(wt.matrix.1,wt.matrix.2,A.M)
      
      if (sep.graphs){
        Graph.1 <-  graph.adjacency(wt.matrix.1, weighted  =   TRUE , mode  =  graph.mode)
        Graph.2 <-  graph.adjacency(wt.matrix.2, weighted  =   TRUE , mode  =  graph.mode)
      }
      else{
        Graph.1 <-  graph.adjacency(wt.matrix.1,  weighted  =   TRUE, mode  =  graph.mode)
        Graph.2 <-  graph.adjacency(wt.matrix.2,  weighted  =   TRUE, mode  =  graph.mode)
        
        ind.vec <-  c(rep(TRUE,n),in.sample.ind[n+(1:n)])
        
        Graph.M.1  <-   graph.adjacency(G.comb.w[ind.vec,ind.vec],weighted  =   TRUE ,
                                        mode  =  graph.mode)
        ind.vec <-  c(in.sample.ind[(1:n)],rep(TRUE,n))
        Graph.M.2  <-   graph.adjacency(G.comb.w[ind.vec,ind.vec],weighted  =   TRUE ,
                                        mode  =  graph.mode)
        
      }
    }
    else{  #if wt.matrix.1 is not null
      
      if (sep.graphs){
        #Given weight matrix, generate weighted graph                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
        Graph.1 <-  graph.adjacency(wt.matrix.1 ,weighted  =  TRUE,mode  =   graph.mode)
        Graph.2 <-  graph.adjacency(wt.matrix.2,weighted  =  TRUE,mode  =  graph.mode)
        
        #Don't need to compute Graph.M
        
      }
      else{
        G.comb.w <-  omnibusM(wt.matrix.1,wt.matrix.2,A.M)
        Graph.1 <-  graph.adjacency(wt.matrix.1,weighted  =   TRUE, mode  =  graph.mode)
        Graph.2 <-  graph.adjacency(wt.matrix.2,weighted  =   TRUE,mode  =  graph.mode)
        Graph.M <-  graph.adjacency(G.comb.w,weighted  =   TRUE ,
                                    mode  =  graph.mode)
        ind.vec <-  c(rep(TRUE,n),in.sample.ind[n+(1:n)])
        
        
        Graph.M.1  <-   graph.adjacency(G.comb.w[ind.vec,ind.vec],weighted  =   TRUE ,
                                        mode  =  graph.mode)
        ind.vec <-  c(in.sample.ind[(1:n)],rep(TRUE,n))
        Graph.M.2  <-   graph.adjacency(G.comb.w[ind.vec,ind.vec],weighted  =   TRUE ,
                                        mode  =  graph.mode)
      }
    }
    
  }
  
  #Now that graphs are generated from adjacency or weight matrices,
  # compute dissimilarities using shortest.paths	
  
  if (sep.graphs){
    
    #compute dissimilarities in separate graphs, then impute dissimilarity between different condition
    D.1 <-  shortest.paths(Graph.1)
    D.2 <-  shortest.paths(Graph.2)
    D.1[is.infinite(D.1)] <-  NA
    D.2[is.infinite(D.2)] <-  NA
    D.w <-   (D.1+D.2)/2#mapply(min,D.1,D.2)
    
    D.M <-   omnibusM(D.1,D.2,D.w)
    
  }
  else{
    D.1 <-  shortest.paths(Graph.1)
    D.2 <-  shortest.paths(Graph.2)
    
    
    D.M.1 <-  shortest.paths(Graph.M.1)
    D.M.2 <-  shortest.paths(Graph.M.2)
    
    ind.vec <-  c(rep(TRUE,n),in.sample.ind[n+(1:n)])
    D.M[ind.vec,ind.vec]  <-   D.M.1
    ind.vec <-  c(in.sample.ind[(1:n)],rep(TRUE,n))
    D.M[ind.vec,ind.vec]  <-   D.M.2
    
    
    
    ind.vec <-  c(in.sample.ind[(1:n)],rep(FALSE,n))
    D.M[ind.vec,ind.vec]  <-   D.1[in.sample.ind[(1:n)],in.sample.ind[(1:n)]]
    ind.vec <-  c(rep(FALSE,n),in.sample.ind[n+(1:n)])
    D.M[ind.vec,ind.vec]  <-   D.2[in.sample.ind[(1:n)],in.sample.ind[(1:n)]]
    #compute dissimilarities in joint graph
    
    D.M[is.infinite(D.M)] <-  NA
  }
  #print("max(D.M)")
  #print(max(D.M))
  
  
  Embed.List <-  Embed.Nodes(D.M,  in.sample.ind ,oos ,
                             d  =  d.dim,
                             wt.equalize  =  FALSE,
                             separability.entries.w  =  FALSE,
                             assume.matched.for.oos  =   FALSE,
                             w.vals  =  w.vals.vec)	
  
  
  
  J <-  list()
  for (Y.embed in Embed.List$Y.embeds){
    
    test.samp.size <-  nrow(Y.embed)/2
    Dist  <-  as.matrix(dist(Y.embed))[1:test.samp.size,(1:test.samp.size)+test.samp.size]
    
    J <-  c(J,list(Dist))
    
  }
  
  return(J)
  
}

tMatch.oos<-function(M,oos.log){
  true.pairings <-  0
  pair.names <-  levels(M)
  num.pairs <-  length(pair.names)
  oos.indices <- which(oos.log) 
  num.oospairs <- length(oos.indices)/2
  for (i in oos.indices[1:num.oospairs]){
    paired.ind  <-  which(M == levels(M)[i])
    if (abs(paired.ind[1]-paired.ind[2]) == num.pairs) {
      true.pairings  <-   true.pairings + 1
    }
    # if (!(paired.ind[1]==i)&& (paired.ind[2]==(i+num.pairs)))
    #   stop("An assumption is false")
  }
  print(paste(true.pairings," true matches  out of ", num.inpairs ," pairings"))
  return(true.pairings)


}



Embed.Nodes  <-  function(D.omnibus,
                          in.sample.ind,
                          oos, 
                          d.start,
                          wt.equalize  =  FALSE,
                          separability.entries.w  =  TRUE,
                          assume.matched.for.oos   =   FALSE ,
                          w.vals  =  0.95,
                          oos.embed.n.at.a.time   =   sum(!in.sample.ind)/2,
                          mds.init.method="gower"){
  
  
  oos.use.imputed <-   TRUE
  w.max.index <-  length(w.vals)
  # number of insample pairs
  n <-  sum(in.sample.ind)/2
  # number of oos pairs
  all.m <-  sum(!in.sample.ind)/2
  
  #in.sample.ind <-  which(in.sample.ind)
  
  D.in  <-   D.omnibus[in.sample.ind,in.sample.ind]
  
  # embedding order for oos vertices
  all.oos.indices <-   which(!in.sample.ind)
  oos.sample.indices <-all.oos.indices
  
  
  
  #number  of groups that are  embedded at the same time
  test.m <- all.m
  insample.indices <-   which(in.sample.ind)
  
  # Embed in-sample using different weight matrices (differentw values)
  
  dim.increment<-10
  
  init.conf  <-  NULL
  
  
  # Find the minimum embedding dimension that matches all the seeds correctly
  full.seed.match  <-   FALSE
  
  
 
  
  embed.dim  <-  d.start
  prevTrueMatch = -1
  True.match.last.memory <- rep(-1,3)
  w.vals.l<-w.vals[[1]]
  full.seed.match<-FALSE
 
  Y.embeds.s<-list()
  embed.dim<-embed.dim-dim.increment
  
  while  (!full.seed.match) {
    embed.dim   <-   embed.dim + dim.increment
   
    Y.0 <- Embed.at.dim( D.in,D.omnibus,
                        in.sample.ind,
                        oos, 
                        embed.dim,
                        wt.equalize  =  FALSE,
                        separability.entries.w  =  FALSE,
                        assume.matched.for.oos   =   FALSE ,
                        w.vals  =  w.vals.l,
                        oos.embed.n.at.a.time   =   sum(!in.sample.ind)/2,
                        mds.init.method="gower")
    
    num.pairs <- n+all.m
    Y.0.in.1<-Y.0[1:n,]
    Y.0.in.2<-Y.0[n+(1:n),]
    Y.0.oos.1<-Y.0[(2*n)+(1:all.m),]
    Y.0.oos.2<-Y.0[(2*n+all.m)+(1:all.m),]
    
    pw.dist.insample <- as.matrix(dist(rbind(Y.0.in.1,Y.0.oos.1,Y.0.in.2,Y.0.oos.2)))
    
    cost.mat <- pw.dist.insample[1:num.pairs,num.pairs+(1:num.pairs)]
    rownames(cost.mat) <- 1:num.pairs
    colnames(cost.mat) <- 1:num.pairs
    #insample.match <-   pairmatch()
    
    #numTrueMatch <-  tMatch.insample(insample.match,in.sample.ind)
    #temp.ind<- c(rep(T,n),rep(F,all.m))
    #temp.ind<- c(temp.ind, rep(T,n),rep(F,all.m))
    #numTrueMatch <-  tMatch.insample(insample.match,temp.ind)
    
    matching<- solve_LSAP(cost.mat)
    all.matches <- as.matrix(matching)==1:num.pairs
    numTrueMatch <- sum(all.matches[1:n])
    
    print(paste(numTrueMatch," true matches  out of ", n," pairings"))
    
    
    
    
    if (numTrueMatch/n > 0.95){
      full.seed.match   <-    TRUE
      print(paste("optimal dim is ", embed.dim))
    }
    if (all( True.match.last.memory== numTrueMatch)) {
      full.seed.match   <-    TRUE
      print(paste("optimal dim is ", embed.dim))
    }
    True.match.last.memory[1] <-    True.match.last.memory[2]
    True.match.last.memory[2] <-    True.match.last.memory[3]
    True.match.last.memory[3] <-   numTrueMatch
    
    print("OOS embedding \n")
    
  } #  end while  (!full.seed.match)
  
  full.seed.match<-FALSE
  prevTrueMatch = -1
  True.match.last.memory <- rep(-1,3)
  
  
  Y.embeds<-list()
  embed.dim <- embed.dim-dim.increment
  while  (!full.seed.match) {
    Y.embeds<-list()
    embed.dim   <-   embed.dim + 2
    Y.0 <- Embed.at.dim( D.in,D.omnibus,
                        in.sample.ind,
                        oos, 
                        embed.dim,
                        wt.equalize  =  FALSE,
                        separability.entries.w  =  FALSE,
                        assume.matched.for.oos   =   FALSE ,
                        w.vals  =  w.vals.l,
                        oos.embed.n.at.a.time   =   sum(!in.sample.ind)/2,
                        mds.init.method="gower")
    num.pairs <- n+all.m
    Y.0.in.1<-Y.0[1:n,]
    Y.0.in.2<-Y.0[n+(1:n),]
    Y.0.oos.1<-Y.0[(2*n)+(1:all.m),]
    Y.0.oos.2<-Y.0[(2*n+all.m)+(1:all.m),]
    
    pw.dist.insample <- as.matrix(dist(rbind(Y.0.in.1,Y.0.oos.1,Y.0.in.2,Y.0.oos.2)))
    
    cost.mat <- pw.dist.insample[1:num.pairs,num.pairs+(1:num.pairs)]
    rownames(cost.mat) <- 1:num.pairs
    colnames(cost.mat) <- 1:num.pairs
    #insample.match <-   pairmatch()
    
    #numTrueMatch <-  tMatch.insample(insample.match,in.sample.ind)
    #temp.ind<- c(rep(T,n),rep(F,all.m))
    #temp.ind<- c(temp.ind, rep(T,n),rep(F,all.m))
    #numTrueMatch <-  tMatch.insample(insample.match,temp.ind)
    
    matching<- solve_LSAP(cost.mat)
    all.matches <- as.matrix(matching) == 1:num.pairs
    numTrueMatch <- sum(all.matches[1:n])
    print(paste(numTrueMatch," true matches  out of ", n  ," pairings"))
    
    
    
    if (numTrueMatch/n> 0.95){
      full.seed.match   <-    TRUE
      print(paste("optimal dim is ", embed.dim))
    }
    if (all( True.match.last.memory== numTrueMatch)) {
      full.seed.match   <-    TRUE
      print(paste("optimal dim is ", embed.dim))
    }
    True.match.last.memory[1] <-    True.match.last.memory[2]
    True.match.last.memory[2] <-    True.match.last.memory[3]
    True.match.last.memory[3] <-   numTrueMatch
    
    
    #Y.oos <- Y.0[all.oos.indices,]
    Y.embeds <-  c(Y.embeds,list(rbind(Y.0.oos.1,Y.0.oos.2)))
    print("OOS embedding \n")
   
    
  }
    
  #  }
  #sink()
  
  print("OOS embedding complete")
  #sink()
  return (list(Y.embeds = Y.embeds,matches =matching))
  
}

Embed.at.dim<- function ( D.in,D.omnibus,
                         in.sample.ind,
                         oos, 
                         embed.dim,
                         wt.equalize  =  FALSE,
                         separability.entries.w  =  FALSE,
                         assume.matched.for.oos   =   FALSE ,
                         w.vals  =  0.5,
                         oos.embed.n.at.a.time   =   sum(!in.sample.ind)/2,
                         mds.init.method="gower") {
  oos.use.imputed <-   TRUE
  w.max.index <-  length(w.vals)
  # number of insample pairs
  n <-  sum(in.sample.ind)/2
  # number of oos pairs
  all.m <-  sum(!in.sample.ind)/2
  
  
  insample.indices <-   which(in.sample.ind)
  
  # Embed in-sample using different weight matrices (differentw values)  
  
  init.conf  <- NULL# mds.init.method
  
  # embedding order for oos vertices
  all.oos.indices <-   which(!in.sample.ind)
  oos.sample.indices <-all.oos.indices
  
  
  
  #number  of groups that are  embedded at the same time
  test.m <- all.m
  insample.indices <-   which(in.sample.ind)
  
  

 
  
  l <- 1
  print(paste("OOS embedding for JOFC for w  = ",w.vals[l]))
  
  
  w.val.l  <-   w.vals[l]
  
  Y.0  <-   matrix(0,length(in.sample.ind),embed.dim)
  

  omnibus.oos.D.0  <-   rbind(
    cbind(D.in,D.omnibus[insample.indices,oos.sample.indices]),
    cbind(D.omnibus[oos.sample.indices,insample.indices],
          D.omnibus[oos.sample.indices,oos.sample.indices])
  )
  
  if (sink.number()>0) sink()
  #sink("Embedding.debug.txt")
  
  
  #Compute Weight matrix corresponding in-sample  entries
  # Since we are going to oos-embedding, set the weights  of in-sample embedding of stress
  # We are using previous in-sample embeddings, anyway
  oos.Weight.mat.in <-  matrix(1-w.val.l, 2*n, 2*n)
  
  
  # If assume.matched.for.oos is true, we assume OOS dissimilarities are matched(in reality,
  # they are matched for the matched pairs, but unmatched for the unmatched pairs)
  # If assume.matched.for.oos is true, we ignore the dissimilarities between matched/unmatched 
  # pairs
  
  #oos.Weight.mat.oos <-   matrix(0,2*test.m,2*test.m)
  
  oos.Weight.mat.oos <-   rbind(cbind(matrix(1-w.val.l,test.m,test.m), matrix(0,test.m,test.m) ),
                                cbind(matrix(0,test.m,test.m),matrix(1-w.val.l,test.m,test.m))
                            )
  
  # if (oos.use.imputed is true) we treat the dissimiilarities between  in-sample and out-of-sample measurements
  # from different conditions like fidelity terms
  # otherwise they are ignored
  imp.weight<- 0.5*min(1-w.val.l,w.val.l)
  if (oos.use.imputed){
    oos.Weight.mat.w  <-   matrix(imp.weight,2*n,2*test.m)
  } else{
    oos.Weight.mat.w  <-   rbind(cbind(matrix(imp.weight,n,test.m), matrix(0,n,test.m) ),
                                 cbind(matrix(0,n,test.m),matrix(imp.weight,n,test.m))
                                )
  }
  
  

    
    oos.Weight.mat.w.2  <-   rbind(cbind(matrix(imp.weight,test.m,n), matrix(0,test.m,n) ),
                                   cbind(matrix(0,test.m,n),matrix(imp.weight,test.m,n))
                                  )
    
 # } else {
#      oos.Weight.mat.w.2 <- matrix(0,2*test.m,2*n)
 # }
  
  
  
  
  oos.Weight.mat <-  rbind(cbind(oos.Weight.mat.in,oos.Weight.mat.w),
                           cbind(oos.Weight.mat.w.2,oos.Weight.mat.oos)
                              )
  
  
  
  
  print("NAs supposed to have 0 weights")
  print(sum(oos.Weight.mat[is.na(omnibus.oos.D.0)]))
  oos.Weight.mat[is.na(omnibus.oos.D.0)] <-  0
  omnibus.oos.D.0[is.na(omnibus.oos.D.0)] <-  1
  
  #print("JOFC null omnibus OOS embedding \n")
  
  
  #The argument W foor oosIM function
  #is formed of  weights for in.sample indices in the upper left,
  #weight for oos in the upper right 
  # ALWAYS independent of  isWithin
#   Y.0.embed <-  oosIM(D   =   omnibus.oos.D.0,
#                       X   =   Y.0[in.sample.ind,],
#                       init       =   mds.init.method,
#                       verbose    =   FALSE,
#                       itmax      =   1000,
#                       eps        =   1e-8,
#                       W          =   oos.Weight.mat,
#                       isWithin   =   NULL,
#                       bwOos      =   TRUE)
  
  Y.0 <- smacofM  (omnibus.oos.D.0,
                      ndim    = embed.dim,
                      W       = oos.Weight.mat)
  
  #Y.0.embed
  
  #print("oos.sample.indices")
  #print(oos.sample.indices)
  #Y.0[oos.sample.indices,] <-  Y.0.embed
  #sink()
  
  return(Y.0)
}

Embed.Nodes.one.atat  <-  function(D.omnibus,
                                   in.sample.ind,
                                   oos, 
                                   d.start,
                                   wt.equalize  =  FALSE,
                          separability.entries.w  =  FALSE,
                          assume.matched.for.oos   =   FALSE ,
                          w.vals  =  0.99,
                          oos.embed.n.at.a.time   =   sum(!in.sample.ind)/2,
                          mds.init.method="gower"){
  
  #sink("Embed.debug.txt")
  matching <- 0
  Y.embeds <-  list()
  oos.use.imputed <-   FALSE
  w.max.index <-  length(w.vals)
  # number of insample pairs
  n <-  sum(in.sample.ind)/2
  # number of oos pairs
  all.m <-  sum(!in.sample.ind)/2
  
  #in.sample.ind <-  which(in.sample.ind)
  
  D.in  <-   D.omnibus[in.sample.ind,in.sample.ind]
  
  # embedding order for oos vertices
  all.oos.indices <-   which(!in.sample.ind)
  v.embed.order    <-   sample(all.oos.indices[1:all.m],all.m,replace  =  FALSE)
  v.embed.order.2  <-   sample(all.oos.indices[all.m+(1:all.m)],all.m,replace  =  FALSE)
  embed.order    <-   c(v.embed.order,v.embed.order.2)
  
  #number  of groups that are  embedded at the same time
  num.embed.groups   <-   ceiling(all.m/oos.embed.n.at.a.time)
  insample.indices <-   which(in.sample.ind)
  
  # Embed in-sample using different weight matrices (differentw values)
  
    
    
    
    init.conf  <-  NULL
    
    
    # Find the minimum embedding dimension that matches all the seeds correctly
    full.seed.match  <-   FALSE
    embed.dim  <-  d.start
    prevTrueMatch = -1
  True.match.last.memory <- rep(-1,3)
  
    while  (!full.seed.match) {
    X.embeds.f <- list()
    embed.dim   <-   embed.dim + 5
      
    X.embeds.f <-  JOFC.Insample.Embed(D.in,embed.dim,
                                       w.vals,sep.err.w=TRUE,
                                       init.conf  =  init.conf,
                                       wt.equalize  =  wt.equalize)
    # if (inherits(X.embeds,"try-error")) {
    #   print('Unable to embed via smacof')
    #   X.embeds<-list(cmdscale(D.in, k=d.start))
    #   embed.dim<-d.start
    #   full.seed.match   <-    TRUE
    # }
    
    pw.dist.insample <- as.matrix(dist(X.embeds.f[[1]]))
    
    cost.mat <- pw.dist.insample[1:n,n+(1:n)]
    rownames(cost.mat) <- 1:n
    colnames(cost.mat) <- 1:n
    #insample.match <-   pairmatch()
    
    #numTrueMatch <-  tMatch.insample(insample.match,in.sample.ind)
    #temp.ind<- c(rep(T,n),rep(F,all.m))
    #temp.ind<- c(temp.ind, rep(T,n),rep(F,all.m))
    #numTrueMatch <-  tMatch.insample(insample.match,temp.ind)
    
    matching<- solve_LSAP(cost.mat)
    all.matches <- as.matrix(matching) == 1:n
    numTrueMatch <- sum(all.matches[1:n])
    print(paste(numTrueMatch," true matches  out of ", n  ," pairings"))
    
    if (numTrueMatch /n>0.95){
      full.seed.match   <-    TRUE
      print(paste("optimal dim is ", embed.dim))
    }
    if (all( True.match.last.memory== numTrueMatch)) {
		full.seed.match   <-    TRUE
      print(paste("optimal dim is ", embed.dim))
          }
    True.match.last.memory[1] <-    True.match.last.memory[2]
    True.match.last.memory[2] <-    True.match.last.memory[3]
    True.match.last.memory[3] <-   numTrueMatch
     
  } 
  full.seed.match<-FALSE
  prevTrueMatch = -1
  True.match.last.memory <- rep(-1,3)
  embed.dim <- embed.dim -5
      
      
  while  (!full.seed.match) {
    embed.dim   <-   embed.dim + 5
    
    X.embeds <-  try(JOFC.Insample.Embed(D.in,ndimens=embed.dim,
                                         w.vals,sep.err.w=TRUE,
                                         init.conf  =  init.conf,
                                         wt.equalize  =  wt.equalize))
    if (inherits(X.embeds,"try-error")) {
      print('Unable to embed via smacof')
      embed.dim<-embed.dim-5
      X.embeds<-list(cmdscale(D.in, k=embed.dim))
      
        full.seed.match   <-    TRUE
      }
   # print("Insample embedding complete")
    #    for (l in 1:w.max.index){
      
    l <- 1
      print(paste("OOS embedding for JOFC for w  = ",w.vals[l]))
      
      
      w.val.l  <-   w.vals[l]
      X  <-   X.embeds[[l]]
    dim.X <- dim(X)
    Y.0  <-   matrix(-1,length(in.sample.ind),dim.X[2])
    print(dim(X))
    print(dim(Y.0))
    print(sum(in.sample.ind))
    Y.0[in.sample.ind,] <- X
      
      #Vertices are embedded in groups
      for (embed.group in 1:num.embed.groups){
        #sink()
        #sink("Embedding.debug.txt")
        #embed the next test.m (oos) vertices
        test.m  <-   oos.embed.n.at.a.time 
        if (embed.group == num.embed.groups)
          test.m  <-  all.m-(num.embed.groups-1)*oos.embed.n.at.a.time
        embed.ind <-  embed.order[(oos.embed.n.at.a.time*(embed.group-1))+(1:test.m)]
        #embed.ind <-  sort(embed.ind)
        #embed.order[(oos.embed.n.at.a.time*(embed.group-1))+(1:test.m)] <-  embed.ind
        embed.ind.2 <-  embed.order[all.m+(oos.embed.n.at.a.time*(embed.group-1))+(1:test.m)]
        #embed.ind.2 <-  sort(embed.ind.2)        
        #embed.order[all.m+(oos.embed.n.at.a.time*(embed.group-1))+(1:test.m)] <-  embed.ind.2
        
      #print("test.m and n")
      # print(test.m)
      #print(n)
        
        
        
        oos.sample.indices <-  c(embed.ind,embed.ind.2) 
        
        omnibus.oos.D.0  <-   rbind(
          cbind(D.in,D.omnibus[insample.indices,oos.sample.indices]),
          cbind(D.omnibus[oos.sample.indices,insample.indices],
                D.omnibus[oos.sample.indices,oos.sample.indices])
        )
        if (sink.number()>0) sink()
        #sink("Embedding.debug.txt")
        
        
        #Compute Weight matrix corresponding in-sample  entries
        # Since we are going to oos-embedding, set the weights  of in-sample embedding of stress
        # We are using previous in-sample embeddings, anyway
        oos.Weight.mat.in <-  matrix(0,2*n,2*n)
        
        
        # If assume.matched.for.oos is true, we assume OOS dissimilarities are matched(in reality,
        # they are matched for the matched pairs, but unmatched for the unmatched pairs)
        # If assume.matched.for.oos is true, we ignore the dissimilarities between matched/unmatched 
        # pairs
        
        oos.Weight.mat.oos <-   matrix(0,2*test.m,2*test.m)
        
        
        
        # if (oos.use.imputed is true) we treat the dissimiilarities between  in-sample and out-of-sample measurements
        # from different conditions like fidelity terms
        # otherwise they are ignored
        if (oos.use.imputed){
          oos.Weight.mat.w  <-   matrix(1-w.val.l,2*n,2*test.m)
        } else{
          oos.Weight.mat.w  <-   rbind(cbind(matrix(1-w.val.l,n,test.m), matrix(0,n,test.m) ),
                                       cbind(matrix(0,n,test.m),matrix(1-w.val.l,n,test.m))
          )
        }
        
        
        oos.Weight.mat <-  omnibusM(oos.Weight.mat.in,oos.Weight.mat.oos,oos.Weight.mat.w)
        
        
        
        
        oos.Weight.mat[is.na(omnibus.oos.D.0)] <-  0
        omnibus.oos.D.0[is.na(omnibus.oos.D.0)] <-  1
        
      #print("JOFC null omnibus OOS embedding \n")
        
        
        #The argument W foor oosIM function
        #is formed of  weights for in.sample indices in the upper left,
        #weight for oos in the upper right 
        # ALWAYS independent of  isWithin
        Y.0.embed <-  oosIM(D   =   omnibus.oos.D.0,
                            X   =   X,
                            init       =   mds.init.method,
                            verbose    =   FALSE,
                            itmax      =   1000,
                            eps        =   1e-8,
                            W          =   oos.Weight.mat,
                            isWithin   =   NULL,
                          bwOos      =   TRUE)
      #print("oos.sample.indices")
      #print(oos.sample.indices)
        Y.0[oos.sample.indices,] <-  Y.0.embed
        #sink()
      }
    
    pw.dist.insample <- as.matrix(dist(Y.0))
    num.pairs <- n+all.m
    
    
    cost.mat <- pw.dist.insample[1:num.pairs,num.pairs+(1:num.pairs)]
    rownames(cost.mat) <- 1:num.pairs
    colnames(cost.mat) <- 1:num.pairs
    #insample.match <-   pairmatch()
    
    #numTrueMatch <-  tMatch.insample(insample.match,in.sample.ind)
    #temp.ind<- c(rep(T,n),rep(F,all.m))
    #temp.ind<- c(temp.ind, rep(T,n),rep(F,all.m))
    #numTrueMatch <-  tMatch.insample(insample.match,temp.ind)
    
    matching<- solve_LSAP(cost.mat)
    all.matches <- as.matrix(matching) == 1:num.pairs
    numTrueMatch <- sum(all.matches[in.sample.ind[1:num.pairs]])
    print(paste(numTrueMatch," true matches  out of ", n  ," pairings"))
    
    
    #insample.match <-   pairmatch(pw.dist.insample[1:num.pairs,num.pairs+(1:num.pairs)])
    
    #numTrueMatch <-  tMatch.insample(insample.match,in.sample.ind)
    
    if (numTrueMatch/n > 0.95){
      full.seed.match   <-    TRUE
      print(paste("optimal dim is ", embed.dim))
    }
    if (all( True.match.last.memory== numTrueMatch)) {
      full.seed.match   <-    TRUE
      print(paste("optimal dim is ", embed.dim))
    }
    True.match.last.memory[1] <-    True.match.last.memory[2]
    True.match.last.memory[2] <-    True.match.last.memory[3]
    True.match.last.memory[3] <-   numTrueMatch
    
    
      Y.oos <- Y.0[all.oos.indices,]
      Y.embeds <-  c(Y.embeds,list(Y.oos))
      print("OOS embedding \n")
    
  } #  end while  (!full.seed.match)
  
  #  }
    #sink()
  
  print("OOS embedding complete")
  #sink()
  return (list(Y.embeds=Y.embeds,matches =matching))
  
}

solveMarriage <-   function(Dist){
  matches <-  pairmatch(Dist)    
}

Embed.CMDS  <-   function (D.1,D.2,in.sample.ind,d.dim,oos) {
  n  <-   length(in.sample.ind)/2
  m <-   sum(in.sample.ind)/2
  myD.M   <-   D.M # i use just the object D.M ... none of the original entries!
  myD.M[1:n,1:n]  <-  D.1
  myD.M[(n+1):(2*n),(n+1):(2*n)]  <-  D.2
  myD.M[1:n,(n+1):(2*n)]  <-  (D.1+D.2)/2
  myD.M[(n+1):(2*n),1:n]  <-  (D.1+D.2)/2
  for(i in 1:n) for(j in (n+m+1):(2*n))
    myD.M[i,j]   <-   myD.M[j,i]   <-   c.imp
  for(i in (m+1):n) for(j in (n+1):(2*n)) 
    myD.M[i,j]   <-   myD.M[j,i]   <-   c.imp
  
  
  if (oos){
    
    myD.M.in <-   myD.M[in.sample.ind,in.sample.ind]
    ccc   <-   cmdscale(myD.M.in,k  =  d.dim,eig  =  T)
    #plot(ccc$eig)
    #pairs(ccc$points , col  <-  colvec,pch  <-  c(Ln,Ln))
    #plot(ccc$points[,c(2,3)] , col  <-  colvec,pch  <-  c(Ln,Ln))
    Y.emb <-  oosMDS(myD.M,X = ccc$points,
                     w  =  ifelse(in.sample.ind,1,0),init  =  "gower")
    
    # 			
    U   <-   as.matrix(dist(Y.emb[,c(2,3)]))
  } else {
    
    
    ccc   <-   cmdscale(myD.M,k  = d.dim,eig  =  T)
    #plot(ccc$eig)
    #pairs(ccc$points , col  <-  colvec,pch  <-  c(Ln,Ln))
    #plot(ccc$points[,c(2,3)] , col  <-  colvec,pch  <-  c(Ln,Ln))
    
    U   <-   as.matrix(dist(ccc$points[ c((m+1):(n),(n+m+1):(n+n)) ,c(2,3)]))
  }
  return(U[1:(n-m),(n-m+1):(2*(n-m))])
}


present <-  function(M){
  true.pairings <-  0
  pair.names <-  levels(M)
  num.pairs <-  length(pair.names)
  for (i in 1:num.pairs){
    paired.ind  <-  which(M == levels(M)[i])
    if (abs(paired.ind[1]-paired.ind[2]) == num.pairs){
      true.pairings  <-   true.pairings + 1
    }
  }
  print(paste(true.pairings," true matches  out of ", num.pairs ," pairings"))
  return(true.pairings)
}

tMatch.insample <-  function(M,insample.ind){
  true.pairings <-  0
  pair.names <-  levels(M)
  num.pairs <-  length(pair.names)
  insample.indices <- which(insample.ind) 
  num.inpairs <- length(insample.indices)/2
  for (i in insample.indices[1:num.inpairs]){
    paired.ind  <-  which(M == levels(M)[i])
    if (abs(paired.ind[1]-paired.ind[2]) == num.pairs) {
      true.pairings  <-   true.pairings + 1
    }
    # if (!(paired.ind[1]==i)&& (paired.ind[2]==(i+num.pairs)))
    #   stop("An assumption is false")
  }
  print(paste(true.pairings," true matches  out of ", num.inpairs ," pairings"))
  return(true.pairings)
}

ER   <-   function(n,p)
{
  A   <-   matrix( rbinom(n^2,1,p) , ncol  =  n,byrow  =  T )
  A  <-  A*upper.tri(A)
  A  <-  A+t(A)
  return(A)
}

bitflip   <-   function(G,q10,q01,binary  =  T,symmetric = T,hollow  =  T)
  # takes graph (adjacency matrix) G and flips 1s to 0s with prob q10 and flips 0s to 1s with prob q01
  # assumes binary  =  T,symmetric  =  T,hollow  =  T
{
  n  <-  dim(G)[1]
  for(u in 1:(n-1))
    for(v in (u+1):n)
      G[u,v]  <-  G[v,u]  <-  ifelse(G[u,v],1-rbinom(1,1,q10),rbinom(1,1,q01))
  return(G)
}
adj.Mat.2.P <-  function(A){
  a.avg <-  rowMeans(A)
  for (i in 1:nrow(A)) A[,i] <-  A[,i]/a.avg[i]	
}
diff.dist <-  function(P){
  eig(P)
}


diff.dist.fun <-  function(A,T.diff,dissimilarity=FALSE){
  P <-  transition.matrix(A,dissimilarity  =  dissimilarity)
  D <-  diffusion.distance(P, T.diff, directed   =   FALSE)
  D
}

C_dice_weighted <- function(W){
  n<-nrow(W)
  diag(W)<-0
  D<- matrix(0,n,n)
  for (i in 1:n) {
    for (j in 1:n) {
     
      r_ij= sum((W[i,]* (W[i,]>0 &  W[j,]==0)))+ sum((W[j,]* (W[i,]==0 &  W[j,]>0)))-2*W[i,j]
      a_ij <- sum((W[i,]+W[j,])*((W[i,]>0)&(W[j,]>0)))+2*W[i,j]
      #t <- ifelse(W[i,j]==0,1,0)
      D[i,j] <- (r_ij+2*(W[i,j]==0))/(r_ij+a_ij+2)
    }
    
  }
  diag(D)<-0
  return(5*D)
}


ectime<-function(W){
  n <-nrow(W)
  w<-matrix(1,nrow(W),1)
  t<- W%*%w
  t<- c(t)
  #t<- t+ 1E-3
  T.mat <- diag(t)
  L<- T.mat-W
  
  if ((nrow(L)!=n)| (ncol(L)!=n))
    print(dim(L))
  
  nL<-try((diag(t^(-.5))) %*% L %*% (diag(t^(-.5))))
  if (inherits(nL,"try-error")){
    
    print("nL exception")
    print(str(nL))
    print(str(L))
    print(str(t))
    print(dim(L))
    
    print(dim(t))
    stop()
  }
  
  
  if (any(!is.finite(nL))) {
    print("nL invalid values")
    print(str(nL))
    print(str(L))
    print(str(t))
    print( nL[!is.finite(nL)])
    print( W[!is.finite(W)])
    print( t[!is.finite(t)])
  }
  nG<- ginv(nL)
  G<-(diag(t^(-.5)))%*% nG %*%(diag(t^(.5)))
  v=sum(t)
  ect<- matrix(0,nrow(W),ncol(W))
  for(i in 1:n) {
    for(j in 1:n){
      ect[i,j] <- v* (G[i,i]/t[i]+G[j,j]/t[j]-G[i,j]/t[i]-G[j,i]/t[j])
    }
  }
  ect
}


