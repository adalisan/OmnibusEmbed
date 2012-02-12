
JOFC.Insample.Embed.many.match <-function(D1,D2,D.W,ndimens,w.vals,sep.err.w,init.conf,wt.equalize){
#	if (profile.mode) Rprof("JOFC.FC.out",append=TRUE)
	n.1<- nrow(D1)
	n.2<- nrow(D2)
	D<-omnibusM(D1,D2,D.W)
	
	smacof.embed<-list()
	stress.vec<-c()
	comm.sum.vec<-c()
	fid1.sum.vec<-c()
	fid2.sum.vec<-c()
	
	comm.vec<-c()
	fid1.vec<-c()
	fid2.vec<-c()
	
	
	
	for (w in w.vals){
		W.1<-matrix(1-w,n.1,n.1)
		W.2<-matrix(1-w,n.2,n.2)
		W.W<-matrix(0,n.1,n.2)
		W.W[D.W==0]<-w
		Weight.Mat<-omnibusM(W.1,W.2,W.W)
		Weight.Mat[is.na(D)]<-0
		D[is.na(D)] <-1
		
		new.embed <- smacofM(D,ndimens    ,	W=Weight.Mat        ,
				init    = init.conf,
				verbose = FALSE,
				itmax   = 1000,
				eps     = 1e-6)
		smacof.embed<-c(smacof.embed,list(new.embed ))
		stress.mat <- (as.dist(D) - dist(new.embed))^2
		
			
	}
	
	return(smacof.embed)
}

genG <- function(s=seed,n=1000,np=50,p=0.1,q=0.5) {
	x <- matrix(p, nrow=n, ncol=n)    ## null: (100-np)x98 matrix
	if (np>0) {
		egg <- matrix(q, nrow=np, ncol=np) ## alternative: npx98 matrix
		x[(n-np+1):n,(n-np+1):n] <- egg
	}
	ordering<-sample.int(n)
	diag(x) <- 0
	
	A <- matrix(0,nrow=n,ncol=n)
	for(i in 1:(n-1))
		for(j in i:n)
			A[i,j] <- A[j,i] <- sample(c(0,1),1,prob=c(1-x[i,j],x[i,j]))
	
	return(A)
}


perturbG<-function(G,q){
	n<-nrow(G)
	Flip.mat<-matrix(0,nrow=n,ncol=n)
	for(i in 1:(n-1))
		for(j in i:n){
			Flip.mat[i,j]<-sample(c(0,1),1,prob=c(1-q,q))
			
		}
	for(i in 1:(n-1))
		for(j in i:n){
			if (Flip.mat[i,j]==1){
				G[i,j]<-1-G[i,j]
				G[j,i]<-1-G[j,i]
			}
			
		}
	return(G)
	
}



jofc.many<-function(G,Gp,corr.list,
		in.sample.ind.1,in.sample.ind.2,
		d.dim,
		w.vals.vec,
		graph.is.directed=FALSE,
		oos=TRUE,
        notconnect.wt=10,
		use.weighted.graph=TRUE,
		wt.matrix.1=NULL,
		wt.matrix.2=NULL,
		sep.graphs=TRUE # if TRUE, treat two graphs separately to compute dissimilarities
#and impute W (off-diagonalblock matrix)
# if FALSE, join the graphs and compute dissimilarities from joint graph
	    
){
	n.1<-nrow(G)
	n.2<-nrow(Gp)
	test.m.1<-sum(!in.sample.ind.1)
	test.m.2<-sum(!in.sample.ind.2)
	
	
	graph.mode<- ifelse(graph.is.directed,"directed","undirected")
	
	Graph.1<-graph.empty(directed=graph.is.directed)
	Graph.2<-graph.empty(directed=graph.is.directed)
	Graph.M<-graph.empty(directed=graph.is.directed)
	
	if (!use.weighted.graph) {
		#Given adjacency matrix, generate unweighted graph
		print("Using adjacency for computing dissimilarities")
		Graph.1<-graph.adjacency(G, mode=graph.mode)
		Graph.2<-graph.adjacency(Gp,mode=graph.mode)
    A.M<- matrix(0,n.1,n.2)
	for (corr.i in corr.list){
		for ( i in corr.i[[1]]){
			for ( j in corr.i[[2]]){
				A.M[i,j]<-1
			}
		}
	}
    	G.comb<-omnibusM(G,Gp,A.M)
		Graph.M <- graph.adjacency(G.comb,
				weighted= NULL ,mode=graph.mode)
	} else{
		print("Using weighted graph for computing dissimilarities")
		# If creating  a weighted graph 
		# make the weight matrix from adjacency matrix. Those with same-condition edges have weights of wt.connect/10
		#those  with no edges have weights of wt.connect. Those with "matched edges has weights of matched.cost
		A.M<- matrix(500,n.1,n.2)
		for (corr.i in corr.list){
			for ( i in corr.i[[1]]){
				for ( j in corr.i[[2]]){
					A.M[i,j]<-matched.cost
				}
			}
		}    
		
		
     #A.M[!in.sample.ind[1:n]]<- 0
    if (is.null(wt.matrix.1)){
			#Given adjacency matrix, generate weighted graph
			wt.matrix.1 <- G
			wt.matrix.2 <- Gp
			wt.matrix.1[G==0] <- notconnect.wt
			wt.matrix.2[Gp==0] <- notconnect.wt						
			wt.matrix.1[G==1] <- notconnect.wt/10			
			wt.matrix.2[Gp==1] <- notconnect.wt/10
			G.comb.w<-omnibusM(wt.matrix.1,wt.matrix.2,A.M)
			
			if (sep.graphs){
				Graph.1<-graph.adjacency(wt.matrix.1, weighted= TRUE , mode=graph.mode)
				Graph.2<-graph.adjacency(wt.matrix.2, weighted= TRUE , mode=graph.mode)
			}			else{
        stop("This case is not implemented")
				
			}
		}
		else{  #if wt.matrix.1 is not null
			
			if (sep.graphs){
				#Given weight matrix, generate weighted graph                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
				Graph.1<-graph.adjacency(wt.matrix.1 ,weighted=TRUE,mode= graph.mode)
				Graph.2<-graph.adjacency(wt.matrix.2,weighted=TRUE,mode=graph.mode)
				
				#Don't need to compute Graph.M
				
			}
			else{
				stop("This case not implemented")
			}
		}
		
	}
	
	#Now that graphs are generated from adjacency or weight matrices,
	# compute dissimilarities using shortest.paths	
	
	if (sep.graphs){
		
		#compute dissimilarities in separate graphs, then impute dissimilarity between different condition
		D.1<-shortest.paths(Graph.1)
		D.2<-shortest.paths(Graph.2)
		D.1[is.infinite(D.1)]<-NA
		D.2[is.infinite(D.2)]<-NA
		D.w<- matrix(NA,n.1,n.2)
		for (corr.i in corr.list){
			for ( i in corr.i[[1]]){
				for ( j in corr.i[[2]]){
					D.w[i,j]<-matched.cost
				}
			}
		}
		D.M<- omnibusM(D.1,D.2,D.w)
		
	}
	else{
		stop("This case is not implemented")
	}
	#print("max(D.M)")
	#print(max(D.M))
	
	
Embed.List<-Embed.Nodes.many(D.M[1:n.1,1:n.1],D.M[n.1+(1:n.2),n.1+(1:n.2)],D.M[1:n.1,n.1+(1:n.2)],
		
		in.sample.ind.1,in.sample.ind.2 ,oos ,
		d=d.dim,
		wt.equalize=FALSE,
		separability.entries.w=FALSE,
		assume.matched.for.oos = FALSE,
		w.vals=w.vals.vec)	


J<-list()
for (Y.embed in Embed.List){
	print("dim(Y.embed)")
	print(dim(Y.embed))	
	Dist=as.matrix(dist(Y.embed))[1:test.m.1,(1:test.m.2)+test.m.1]
	J<-c(J,list(Dist))
	
}
return(J)
}



jofc.diffusion.dist.many<-function(G,Gp,corr.list,
		in.sample.ind.1,in.sample.ind.2,
		d.dim,
		w.vals.vec,
		graph.is.directed=FALSE,
    oos=TRUE,
	
		wt.matrix.1=NULL,
		wt.matrix.2=NULL,
		sep.graphs=TRUE # if TRUE, treat two graphs separately to compute dissimilarities
#and impute W (off-diagonalblock matrix)
# if FALSE, join the graphs and compute dissimilarities from joint graph
		
){
	n.1<-nrow(G)
	n.2<-nrow(Gp)
	test.m.1<-sum(!in.sample.ind.1)
	test.m.2<-sum(!in.sample.ind.2)
	
	graph.mode<- ifelse(graph.is.directed,"directed","undirected")
	
	
	if (sep.graphs){
		if (is.null(wt.matrix.1)){
			D.1<-diff.dist.fun(G)
			D.2<-diff.dist.fun(Gp)
		} else{
			D.1<-diff.dist.fun(wt.matrix.1)
			D.2<-diff.dist.fun(wt.matrix.2)
		}
		D.w<- matrix(NA,n.1,n.2)
		for (corr.i in corr.list){
			for ( i in corr.i[[1]]){
				for ( j in corr.i[[2]]){
					D.w[i,j]<-0
				}
			}
		}
		D.M<- omnibusM(D.1,D.2,D.w)
    	}	else{ #combine graphs and compute distances on joint graph
		if (is.null(wt.matrix.1)){
			for (corr.i in corr.list){
				for ( i in corr.i[[1]]){
					for ( j in corr.i[[2]]){
						A.M[i,j]<-1
					}
				}
			}    
			 #A.M[!in.sample.ind[1:n]]<- 0
			G.comb<-omnibusM(G,Gp,A.M)
			#compute dissimilarities in joint graph
			D.M<-diff.dist.fun(G.comb)
			D.M[is.infinite(D.M)]<-NA
		} else{
			wt.W<-matrix(Inf,n.1,n.2)
			for (corr.i in corr.list){
				for ( i in corr.i[[1]]){
					for ( j in corr.i[[2]]){
						wt.W[i,j]<-1
					}
				}
			}    
			
			
			Wt.M<-omnibusM(wt.matrix.1,wt.matrix.2,wt.W)
			D.M<-diff.dist.fun(Wt.M)
            D.M[is.infinite(D.M)]<-NA
		}
		
	}
	
	
	
	Embed.List<-Embed.Nodes.many(D.M[1:n.1,1:n.1],D.M[n.1+(1:n.2),n.1+(1:n.2)],D.M[1:n.1,n.1+(1:n.2)],
			in.sample.ind.1,in.sample.ind.2 ,oos ,
			d=d.dim,
			wt.equalize=FALSE,
			separability.entries.w=FALSE,
			assume.matched.for.oos = FALSE,
			w.vals=w.vals.vec)	
	J<-list()
	for (Y.embed in Embed.List){
		print("dim(Y.embed)")
		print(dim(Y.embed))	
		
		
		Dist=as.matrix(dist(Y.embed))[1:test.m.1,(1:test.m.2)+test.m.1]
		J<-c(J,list(Dist))
		
	}
	return(J)
	
}









Embed.Nodes.many <-function(D.1,D.2,D.W,
                       in.sample.ind.1,
					   in.sample.ind.2,
                       oos, 
                       d,
		wt.equalize=FALSE,
		separability.entries.w=FALSE,
		assume.matched.for.oos = FALSE ,w.vals=0.8){
	n.1<-nrow(D.1)
	n.2<-nrow(D.2)
	n.1.in <- sum(in.sample.ind.1)
	n.2.in <- sum(in.sample.ind.2)
	test.m.1 <- n.1 - n.1.in
	test.m.2 <- n.2 - n.2.in
 	Y.embeds<-list()
	oos.use.imputed<- FALSE
	w.max.index<-length(w.vals)
	  
		
	#in.sample.ind<-which(in.sample.ind)
	
	# Embed in-sample using different weight matrices (differentw values)
	
	
		D.omnibus<-omnibusM(D.1,D.2,D.W)
		D.in.1 <- D.1[in.sample.ind.1,in.sample.ind.1]
		D.in.2 <- D.2[in.sample.ind.2,in.sample.ind.2]
		D.in.W <- D.W[in.sample.ind.1,in.sample.ind.2]
		
		D.in<-omnibusM(D.in.1,D.in.2,D.in.W)
		init.conf=NULL
		if (sum(is.na(D.in))==0) {		
				init.conf<-cmdscale(d=D.in,k=d)
		}	
		
		
		
		
		X.embeds<-JOFC.Insample.Embed.many.match(D.in.1,D.in.2,D.in.W,d,w.vals,separability.entries.w,init.conf=init.conf,
				wt.equalize=wt.equalize)
		for (l in 1:w.max.index){
			if (verbose) print("OOS embedding for JOFC for w= \n")
			
			
			w.val.l <- w.vals[l]
			X <- X.embeds[[l]]
			
			#Compute Weight matrix corresponding in-sample  entries
			# Since we are going to oos-embedding, set the weights  of in-sample embedding of stress
			# We are using previous in-sample embeddings, anyway
			oos.Weight.mat.1<-matrix(0,n.1.in+n.2.in,n.1.in+n.2.in)
			
			
			
			n.oos<- n.1+n.2-(n.1.in+n.2.in)
			#Compute Weight matrix corresponding OOS  entries
			oos.Weight.mat.2<-matrix(1-w.val.l,n.oos,n.oos)
			
			# If assume.matched.for.oos is true, we assume OOS dissimilarities are matched(in reality,
			# they are matched for the matched pairs, but unmatched for the unmatched pairs)
			# If assume.matched.for.oos is true, we ignore the dissimilarities between matched/unmatched 
			# pairs
			
			oos.Weight.mat.2[1:test.m.1,test.m.1+(1:test.m.2)]<-0
			oos.Weight.mat.2[test.m.1+(1:test.m.2),(1:test.m.1)]<-0
			
						
			# if (oos.use.imputed is true) we treat the dissimiilarities between  in-sample and out-of-sample measurements
			# from different conditions like fidelity terms
			# otherwise they are ignored
			if (oos.use.imputed){
				oos.Weight.mat.w <- matrix(1-w.val.l,n.1.in+n.2.in,test.m.1+test.m.2)
			} else{
				oos.Weight.mat.w <- rbind(cbind(matrix(1-w.val.l,n.1.in,test.m.1), matrix(0,n.1.in,test.m.2) ),
						cbind(matrix(0,n.2.in,test.m.1),matrix(1-w.val.l,n.2.in,test.m.2))
				)
			}
			oos.Weight.mat<-omnibusM(oos.Weight.mat.1,oos.Weight.mat.2,oos.Weight.mat.w)
			
			
			in.sample.ind<-c(in.sample.ind.1,in.sample.ind.2)
			omnibus.oos.D.0 <- rbind(
					cbind(D.in,D.omnibus[in.sample.ind,!in.sample.ind]),
					cbind(D.omnibus[!in.sample.ind,in.sample.ind],D.omnibus[!in.sample.ind,!in.sample.ind])
			)
			
			oos.Weight.mat[is.na(omnibus.oos.D.0)]<-0
			omnibus.oos.D.0[is.na(omnibus.oos.D.0)]<-1
			
			if (verbose) print("JOFC null omnibus OOS embedding \n")
			
			Y.0t<-oosIM(D=omnibus.oos.D.0,
					X=X,
					init     = "random",
					verbose  = FALSE,
					itmax    = 1000,
					eps      = 1e-8,
					W        = oos.Weight.mat,
					isWithin = NULL,
					bwOos    = FALSE)
			
			if (verbose) print("JOFC alternative omnibus OOS embedding \n")
			Y.embeds<-c(Y.embeds,list(Y.0t))
		}
		
	
    if (FALSE)
		{
    omnibus.D1.D2<- omnibusM(D1,D2,(D1+D2)/2)
  		init.omnibus<-omnibus.D1.D2  
      init.omnibus[2*n+test.m+(1:test.m),1:(n+test.m)]<-0    
     init.omnibus[1:(n+test.m),2*n+test.m+(1:test.m)]<-0
  
      init.omnibus[n+test.m+(1:n),(n+1):(n+test.m)]<-0    
      init.omnibus[(n+1):(n+test.m),n+test.m+(1:n)]<-0
      
      init.conf=NULL
  	if (sum(is.na(omnibus.D1.D2))==0) {		
				init.conf<-cmdscale(d=init.omnibus,k=d)
		}
    
      W.Mat<-w.val.to.W.mat(w.vals[1],2*(n+test.m),separability.entries.w,wt.equalize)
      W.Mat[2*n+test.m+(1:test.m),1:(n+test.m)]<-0  	
      W.Mat[1:(n+test.m),2*n+test.m+(1:test.m)]<-0
  
      W.Mat[n+test.m+(1:n),(n+1):(n+test.m)]<-0    
      W.Mat[(n+1):(n+test.m),n+test.m+(1:n)]<-0
	
	 
      W.Mat[is.na(omnibus.D1.D2)]<-0
		 omnibus.D1.D2[is.na(omnibus.D1.D2)] <-1
		
    omnibus.D1.D2[W.Mat==0] <- 20
		new.embed <- smacofM (omnibus.D1.D2,
                          ndim=d    ,	W=W.Mat        ,
				init    = init.conf,
				verbose = FALSE,
				itmax   = 1000,
				eps     = 1e-6)
    
    
    Y.embeds<-c(Y.embeds,list(new.embed[!in.sample.ind,]))
	}
  
	
	return(Y.embeds)
	
}

solveMarriage<- function(Dist){
	matches<-fullmatch(Dist)
	return(matches)
	
}

present.many<-function(M,corr.list,in.sample.ind.1,in.sample.ind.2){
	test.m.2.indices <- which(!in.sample.ind.2)
	test.m.1.indices <- which(!in.sample.ind.1)
	
	test.m.1<-length(test.m.1.indices)
	test.m.2<-length(test.m.2.indices)
	precision.list<-rep(0,test.m.2)
	recall.list<-rep(0,test.m.2)
	F.meas.list<-rep(0,test.m.2)
	for (test.i in 1:test.m.2){
		M.i<- test.i+test.m.1
		match.index<- M[M.i]
		matches.i.indic <- M[1:test.m.1]==match.index
		matches.indices.i<- test.m.1.indices[matches.i.indic]
		for (corr.i in corr.list){
			for ( j in corr.i[[2]]){
				if (j==test.m.2.indices[test.i]){
					true.matches<-matches.indices.i%in%corr.i[[1]]
					precision<- sum(true.matches)/length(true.matches)
					recall <-sum(true.matches)/length(corr.i[[1]])
					p.r<-(precision+recall)
					F.meas<-0
					if(p.r!=0)
						F.meas<- 2*precision*recall/p.r
			  		
					precision.list[test.i]<-precision
					recall.list[test.i]<-recall
					F.meas.list[test.i]<-F.meas
					break
				}
			}
		}
	}
	return(list(P=precision.list,R=recall.list,F=F.meas.list))	
}


ER = function(n,p)
{
	A = matrix( rbinom(n^2,1,p) , ncol=n,byrow=T )
	A=A*upper.tri(A)
	A=A+t(A)
	return(A)
}

bitflip = function(G,q10,q01,binary=T,symmetric=T,hollow=T)
# takes graph (adjacency matrix) G and flips 1s to 0s with prob q10 and flips 0s to 1s with prob q01
# assumes binary=T,symmetric=T,hollow=T
{
	n=dim(G)[1]
	for(u in 1:(n-1))
		for(v in (u+1):n)
			G[u,v]=G[v,u]=ifelse(G[u,v],1-rbinom(1,1,q10),rbinom(1,1,q01))
	return(G)
}
adj.Mat.2.P<-function(A){
	a.avg<-rowMeans(A)
	for (i in 1:nrow(A)) A[,i]<-A[,i]/a.avg[i]	
}
diff.dist<-function(P){
	eig(P)
}


diff.dist.fun<-function(A){
	P<-transition.matrix(A,dissimilarity=FALSE)
	D<-diffusion.distance(P, T.diff, directed = FALSE)
	D
}



