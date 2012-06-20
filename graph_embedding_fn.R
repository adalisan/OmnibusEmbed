
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

JOFC.graph.custom.dist<-function(G,Gp,
		in.sample.ind,
		d.dim,
		w.vals.vec,
		graph.is.directed=FALSE
){

	n<-nrow(G)
	graph.mode<- ifelse(graph.is.directed,"directed","undirected")
	
	Graph.1<-graph.empty(directed=graph.is.directed)
	Graph.2<-graph.empty(directed=graph.is.directed)
	Graph.M<-graph.empty(directed=graph.is.directed)
	
		#Given adjacency matrix, generate unweighted graph
		print("Using adjacency for computing dissimilarities")
		Graph.1<-graph.adjacency(G, mode=graph.mode)
		Graph.2<-graph.adjacency(Gp,mode=graph.mode)
		A.M<- diag(n)
		A.M[!in.sample.ind[1:n],!in.sample.ind[1:n]]<- 0 # Make sure vertices that are NOT known to be matched
		#are not connected
		

	#Now that graphs are generated from adjacency or weight matrices,
	# compute dissimilarities using shortest.paths	
	
		
		#compute dissimilarities in separate graphs, then impute dissimilarity between different condition
		D.1<- as.matrix(dist(G,'manhattan'))
		D.2<- as.matrix(dist(Gp,'manhattan'))
		D.1[is.infinite(D.1)]<-NA
		D.2[is.infinite(D.2)]<-NA
		D.w<- (D.1+D.2)/2#mapply(min,D.1,D.2)
		
		D.M<- omnibusM(D.1,D.2,D.w)
		
	
	Embed.List<-Embed.Nodes(D.M,  in.sample.ind ,oos=TRUE ,
			d=d.dim,
			wt.equalize=FALSE,
			separability.entries.w=FALSE,
			assume.matched.for.oos = FALSE,
			w.vals=w.vals.vec)	
	J<-list()
	for (Y.embed in Embed.List){
		
		test.samp.size<-nrow(Y.embed)/2
		Dist=as.matrix(dist(Y.embed))[1:test.samp.size,(1:test.samp.size)+test.samp.size]
		
		J<-c(J,list(Dist))
		
	}
	return(J)


}



JOFC.graph<-function(G,Gp,
		in.sample.ind,
		d.dim,
		w.vals.vec,
		graph.is.directed=FALSE
){

	n<-nrow(G)
	graph.mode<- ifelse(graph.is.directed,"directed","undirected")
	
	Graph.1<-graph.empty(directed=graph.is.directed)
	Graph.2<-graph.empty(directed=graph.is.directed)
	Graph.M<-graph.empty(directed=graph.is.directed)
	
		#Given adjacency matrix, generate unweighted graph
		print("Using adjacency for computing dissimilarities")
		Graph.1<-graph.adjacency(G, mode=graph.mode)
		Graph.2<-graph.adjacency(Gp,mode=graph.mode)
		A.M<- diag(n)
		A.M[!in.sample.ind[1:n],!in.sample.ind[1:n]]<- 0 # Make sure vertices that are NOT known to be matched
		#are not connected
		G.comb<-omnibusM(G,Gp,A.M)

	#Now that graphs are generated from adjacency or weight matrices,
	# compute dissimilarities using shortest.paths	
	
		
		#compute dissimilarities in separate graphs, then impute dissimilarity between different condition
		D.1<-shortest.paths(Graph.1)
		D.2<-shortest.paths(Graph.2)
		D.1[is.infinite(D.1)]<-NA
		D.2[is.infinite(D.2)]<-NA
		D.w<- (D.1+D.2)/2#mapply(min,D.1,D.2)
		
		D.M<- omnibusM(D.1,D.2,D.w)
		
	
	Embed.List<-Embed.Nodes(D.M,  in.sample.ind ,oos=TRUE ,
			d=d.dim,
			wt.equalize=FALSE,
			separability.entries.w=FALSE,
			assume.matched.for.oos = FALSE,
			w.vals=w.vals.vec)	
	J<-list()
	for (Y.embed in Embed.List){
		
		test.samp.size<-nrow(Y.embed)/2
		Dist=as.matrix(dist(Y.embed))[1:test.samp.size,(1:test.samp.size)+test.samp.size]
		
		J<-c(J,list(Dist))
		
	}
	return(J)


}

JOFC.graph.diff<-function(G,Gp,
		in.sample.ind,
		d.dim,
		w.vals.vec,
		graph.is.directed=FALSE,
		T.param=2){

	n<-nrow(G)
	graph.mode<- ifelse(graph.is.directed,"directed","undirected")
	
	Graph.1<-graph.empty(directed=graph.is.directed)
	Graph.2<-graph.empty(directed=graph.is.directed)
	Graph.M<-graph.empty(directed=graph.is.directed)
	
	D.M<-matrix(0,2*n,2*n)	
	D.1<-diff.dist.fun(G,T.param)
	D.2<-diff.dist.fun(Gp,T.param)
	D.w<- (D.1+D.2)/2
		
	D.M<- omnibusM(D.1,D.2,D.w)
		
	
	



			
	
	Embed.List<-Embed.Nodes(D.M,  in.sample.ind ,oos=TRUE ,
			d=d.dim,
			wt.equalize=FALSE,
			separability.entries.w=FALSE,
			assume.matched.for.oos = FALSE,
			w.vals=w.vals.vec)	
	J<-list()
	for (Y.embed in Embed.List){
		
		test.samp.size<-nrow(Y.embed)/2
		Dist=as.matrix(dist(Y.embed))[1:test.samp.size,(1:test.samp.size)+test.samp.size]
		
		J<-c(J,list(Dist))
		
	}
	return(J)


}






jofc<-function(G,Gp,
		in.sample.ind,
		d.dim,
		w.vals.vec,
		graph.is.directed=FALSE,
		oos=TRUE,
		notconnect.wt=10,
		use.weighted.graph=TRUE,
		wt.matrix.1=NULL,
		wt.matrix.2=NULL,
		sep.graphs=TRUE, # if TRUE, treat two graphs separately to compute dissimilarities
#and impute W (off-diagonalblock matrix)
# if FALSE, join the graphs and compute dissimilarities from joint graph
        matched.cost=0.01
){
	
	n<-nrow(G)
	graph.mode<- ifelse(graph.is.directed,"directed","undirected")
	
	Graph.1<-graph.empty(directed=graph.is.directed)
	Graph.2<-graph.empty(directed=graph.is.directed)
	Graph.M<-graph.empty(directed=graph.is.directed)
	
	if (!use.weighted.graph) {
		#Given adjacency matrix, generate unweighted graph
		print("Using adjacency for computing dissimilarities")
		Graph.1<-graph.adjacency(G, mode=graph.mode)
		Graph.2<-graph.adjacency(Gp,mode=graph.mode)
		A.M<- diag(n)
		A.M[!in.sample.ind[1:n],!in.sample.ind[1:n]]<- 0 # Make sure vertices that are NOT known to be matched
		#are not connected
		G.comb<-omnibusM(G,Gp,A.M)
		Graph.M <- graph.adjacency(G.comb,
				weighted= NULL ,mode=graph.mode)
	} else{
		print("Using weighted graph for computing dissimilarities")
		# If creating  a weighted graph 
		# make the weight matrix from adjacency matrix. Those with same-condition edges have weights of wt.connect/10
		#those  with no edges have weights of wt.connect. Those with "matched edges has weights of matched.cost
		A.M<- matrix(1E8,n,n)
		
		diag(A.M) <- matched.cost
		
		A.M[!in.sample.ind[1:n],!in.sample.ind[1:n]]<- 0
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
			}
			else{
				Graph.1<-graph.adjacency(wt.matrix.1,  weighted= TRUE, mode=graph.mode)
				Graph.2<-graph.adjacency(wt.matrix.2,  weighted= TRUE, mode=graph.mode)
				
				ind.vec<-c(rep(TRUE,n),in.sample.ind[n+(1:n)])
				
				Graph.M.1 <- graph.adjacency(G.comb.w[ind.vec,ind.vec],weighted= TRUE ,
						mode=graph.mode)
				ind.vec<-c(in.sample.ind[(1:n)],rep(TRUE,n))
				Graph.M.2 <- graph.adjacency(G.comb.w[ind.vec,ind.vec],weighted= TRUE ,
						mode=graph.mode)
				
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
				G.comb.w<-omnibusM(wt.matrix.1,wt.matrix.2,A.M)
				Graph.1<-graph.adjacency(wt.matrix.1,weighted= TRUE, mode=graph.mode)
				Graph.2<-graph.adjacency(wt.matrix.2,weighted= TRUE,mode=graph.mode)
				Graph.M<-graph.adjacency(G.comb.w,weighted= TRUE ,
						mode=graph.mode)
				ind.vec<-c(rep(TRUE,n),in.sample.ind[n+(1:n)])
				
				
				Graph.M.1 <- graph.adjacency(G.comb.w[ind.vec,ind.vec],weighted= TRUE ,
						mode=graph.mode)
				ind.vec<-c(in.sample.ind[(1:n)],rep(TRUE,n))
				Graph.M.2 <- graph.adjacency(G.comb.w[ind.vec,ind.vec],weighted= TRUE ,
						mode=graph.mode)
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
		D.w<- (D.1+D.2)/2#mapply(min,D.1,D.2)
		
		D.M<- omnibusM(D.1,D.2,D.w)
		
	}
	else{
		D.1<-shortest.paths(Graph.1)
		D.2<-shortest.paths(Graph.2)
		
		
		D.M.1<-shortest.paths(Graph.M.1)
		D.M.2<-shortest.paths(Graph.M.2)
		
		ind.vec<-c(rep(TRUE,n),in.sample.ind[n+(1:n)])
		D.M[ind.vec,ind.vec] <- D.M.1
		ind.vec<-c(in.sample.ind[(1:n)],rep(TRUE,n))
		D.M[ind.vec,ind.vec] <- D.M.2
		
		
		
		ind.vec<-c(in.sample.ind[(1:n)],rep(FALSE,n))
		D.M[ind.vec,ind.vec] <- D.1[in.sample.ind[(1:n)],in.sample.ind[(1:n)]]
		ind.vec<-c(rep(FALSE,n),in.sample.ind[n+(1:n)])
		D.M[ind.vec,ind.vec] <- D.2[in.sample.ind[(1:n)],in.sample.ind[(1:n)]]
		#compute dissimilarities in joint graph
		
		D.M[is.infinite(D.M)]<-NA
	}
	#print("max(D.M)")
	#print(max(D.M))
	
	
	Embed.List<-Embed.Nodes(D.M,  in.sample.ind ,oos ,
			d=d.dim,
			wt.equalize=FALSE,
			separability.entries.w=FALSE,
			assume.matched.for.oos = FALSE,
			w.vals=w.vals.vec)	
	J<-list()
	for (Y.embed in Embed.List){
		
		test.samp.size<-nrow(Y.embed)/2
		Dist=as.matrix(dist(Y.embed))[1:test.samp.size,(1:test.samp.size)+test.samp.size]
		
		J<-c(J,list(Dist))
		
	}
	return(J)
	
}



jofc.diffusion.dist<-function(G,Gp,
		in.sample.ind,
		d.dim,
		w.vals.vec,
		graph.is.directed=FALSE,
		oos=TRUE,
		
		wt.matrix.1=NULL,
		wt.matrix.2=NULL,
		sep.graphs=TRUE ,# if TRUE, treat two graphs separately to compute dissimilarities
#and impute W (off-diagonalblock matrix)
# if FALSE, join the graphs and compute dissimilarities from joint graph
	   T.param=2,
	   matched.cost=0.01
){
	n<-nrow(G)
	graph.mode<- ifelse(graph.is.directed,"directed","undirected")
	D.M<-matrix(0,2*n,2*n)
	
	if (sep.graphs){
		if (is.null(wt.matrix.1)){
			D.1<-diff.dist.fun(G,T.param)
			D.2<-diff.dist.fun(Gp,T.param)
		} else{
			D.1<-diff.dist.fun(wt.matrix.1,T.param)
			D.2<-diff.dist.fun(wt.matrix.2,T.param)
		}
		D.w<- (D.1+D.2)/2
		
		D.M<- omnibusM(D.1,D.2,D.w)
	}	else{
		if (is.null(wt.matrix.1)){
			A.M<- diag(n)
			A.M[!in.sample.ind[1:n],!in.sample.ind[1:n]]<- 0
			G.comb<-omnibusM(G,Gp,A.M)
			#compute dissimilarities in joint graph
			D.M<-diff.dist.fun(G.comb,T.param)
			D.M[is.infinite(D.M)]<-NA
		} else{
			
			Wt.M<-omnibusM(wt.matrix.1,wt.matrix.2,(wt.matrix.1+wt.matrix.2)/2)
			D.M<-diff.dist.fun(Wt.M,T.param)
			D.M[is.infinite(D.M)]<-NA
		}
		
	}
	
	
	
	Embed.List<-Embed.Nodes(D.M,  in.sample.ind ,oos ,
			d=d.dim,
			wt.equalize=FALSE,
			separability.entries.w=FALSE,
			assume.matched.for.oos = FALSE,
			w.vals=w.vals.vec)	
	J<-list()
	for (Y.embed in Embed.List){
		
		test.m<-nrow(Y.embed)/2
		Dist=as.matrix(dist(Y.embed))[1:test.m,(1:test.m)+test.m]
		J<-c(J,list(Dist))
		
	}
	return(J)
	
}









Embed.Nodes <-function(D.omnibus,
		in.sample.ind,
		oos, 
		d,
		wt.equalize=FALSE,
		separability.entries.w=FALSE,
		assume.matched.for.oos = FALSE ,w.vals){
	
	
	Y.embeds<-list()
	oos.use.imputed<- FALSE
	w.max.index<-length(w.vals)
	n<-sum(in.sample.ind)/2
	test.m<-sum(!in.sample.ind)/2
	#in.sample.ind<-which(in.sample.ind)
	
	# Embed in-sample using different weight matrices (differentw values)
	if (oos){	
		sink("Embedding.debug.txt")
    print("test.m and n")
    print(test.m)
    print(n)
		D.in <- D.omnibus[in.sample.ind,in.sample.ind]
		init.conf=NULL
		#if (sum(is.na(D.in))==0) {		
	#		init.conf<-cmdscale(d=D.in,k=d)
	#	}	
		
		
		
		
		X.embeds<-JOFC.Insample.Embed(D.in,d,w.vals,separability.entries.w,init.conf=init.conf,
				wt.equalize=wt.equalize)
    print("Insample embedding complete")
		for (l in 1:w.max.index){
			print("OOS embedding for JOFC for w= \n")
			 print(w.vals[l])
			
			w.val.l <- w.vals[l]
			X <- X.embeds[[l]]
			
			insample.indices<- which(in.sample.ind)
			
			#Compute Weight matrix corresponding in-sample  entries
			# Since we are going to oos-embedding, set the weights  of in-sample embedding of stress
			# We are using previous in-sample embeddings, anyway
			oos.Weight.mat.1<-matrix(0,2*n,2*n)
			
			
			
			
			#Compute Weight matrix corresponding OOS  entries
			oos.Weight.mat.2<-w.val.to.W.mat(w.val.l,(2*test.m),separability.entries.w,wt.equalize)
			
			# If assume.matched.for.oos is true, we assume OOS dissimilarities are matched(in reality,
			# they are matched for the matched pairs, but unmatched for the unmatched pairs)
			# If assume.matched.for.oos is true, we ignore the dissimilarities between matched/unmatched 
			# pairs
			if (!assume.matched.for.oos){
				oos.Weight.mat.2[1:test.m,test.m+(1:test.m)]<-0
				oos.Weight.mat.2[test.m+(1:test.m),(1:test.m)]<-0
			}
			
			
			# if (oos.use.imputed is true) we treat the dissimiilarities between  in-sample and out-of-sample measurements
			# from different conditions like fidelity terms
			# otherwise they are ignored
			if (oos.use.imputed){
				oos.Weight.mat.w <- matrix(1-w.val.l,2*n,2*test.m)
			} else{
				oos.Weight.mat.w <- rbind(cbind(matrix(1-w.val.l,n,test.m), matrix(0,n,test.m) ),
						cbind(matrix(0,n,test.m),matrix(1-w.val.l,n,test.m))
				)
			}
			oos.Weight.mat<-omnibusM(oos.Weight.mat.1,oos.Weight.mat.2,oos.Weight.mat.w)
			
			
			
			omnibus.oos.D.0 <- rbind(
					cbind(D.omnibus[in.sample.ind,in.sample.ind],D.omnibus[in.sample.ind,!in.sample.ind]),
					cbind(D.omnibus[!in.sample.ind,in.sample.ind],D.omnibus[!in.sample.ind,!in.sample.ind])
			)
			
			oos.Weight.mat[is.na(omnibus.oos.D.0)]<-0
			omnibus.oos.D.0[is.na(omnibus.oos.D.0)]<-1
			
			 print("JOFC null omnibus OOS embedding \n")
			
			Y.0t<-oosIM(D=omnibus.oos.D.0,
					X=X,
					init     = "random",
					verbose  = FALSE,
					itmax    = 1000,
					eps      = 1e-8,
					W        = oos.Weight.mat,
					isWithin = NULL,
					bwOos    = FALSE)
			
			 print("JOFC alternative omnibus OOS embedding \n")
			Y.embeds<-c(Y.embeds,list(Y.0t))
     
		}
		 sink()
	}  else{
		
		ind.vec.1<-c(in.sample.ind[1:(n+test.m)],rep(FALSE,n+test.m))  
		ind.vec.2<-c(rep(FALSE,n+test.m),in.sample.ind[n+test.m+(1:(n+test.m))])
		
		ind.vec.3<-c(!in.sample.ind[1:(n+test.m)],rep(FALSE,n+test.m))  
		ind.vec.4<-c(rep(FALSE,n+test.m),!in.sample.ind[n+test.m+(1:(n+test.m))])
		
		
		D1<- rbind(
				cbind(D.omnibus[ind.vec.1,ind.vec.1],D.omnibus[ind.vec.1,ind.vec.3]),
				cbind(D.omnibus[ind.vec.3,ind.vec.1],D.omnibus[ind.vec.3,ind.vec.3])
		)
		
		D2<- rbind(
				cbind(D.omnibus[ind.vec.2,ind.vec.2],D.omnibus[ind.vec.2,ind.vec.4]),
				cbind(D.omnibus[ind.vec.4,ind.vec.2],D.omnibus[ind.vec.4,ind.vec.4])
		)
		
		
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
	
	
	Y.embeds
	
}

solveMarriage<- function(Dist){
	matches<-pairmatch(Dist)
	
	
}

cmds <- function (D.1,D.2,in.sample.ind,d.dim,oos) {
	n <- length(in.sample.ind)/2
	m<- sum(in.sample.ind)/2
	myD.M = D.M # i use just the object D.M ... none of the original entries!
	myD.M[1:n,1:n]=D.1
	myD.M[(n+1):(2*n),(n+1):(2*n)]=D.2
	myD.M[1:n,(n+1):(2*n)]=(D.1+D.2)/2
	myD.M[(n+1):(2*n),1:n]=(D.1+D.2)/2
	for(i in 1:n) for(j in (n+m+1):(2*n)) myD.M[i,j] = myD.M[j,i] = c.imp
	for(i in (m+1):n) for(j in (n+1):(2*n)) myD.M[i,j] = myD.M[j,i] = c.imp
	
	
	if (oos){
		
		myD.M.in<- myD.M[in.sample.ind,in.sample.ind]
		ccc = cmdscale(myD.M.in,k=d.dim,eig=T)
		#plot(ccc$eig)
		#pairs(ccc$points , col=colvec,pch=c(Ln,Ln))
		#plot(ccc$points[,c(2,3)] , col=colvec,pch=c(Ln,Ln))
		Y.emb<-oosMDS(myD.M,X=ccc$points, w=ifelse(in.sample.ind,1,0),init="gower")
		
		# 			
		U = as.matrix(dist(Y.emb[,c(2,3)]))
	} else {
		
		
		ccc = cmdscale(myD.M,k=d.dim,eig=T)
		#plot(ccc$eig)
		#pairs(ccc$points , col=colvec,pch=c(Ln,Ln))
		#plot(ccc$points[,c(2,3)] , col=colvec,pch=c(Ln,Ln))
		
		U = as.matrix(dist(ccc$points[ c((m+1):(n),(n+m+1):(n+n)) ,c(2,3)]))
	}
	return(U[1:(n-m),(n-m+1):(2*(n-m))])
}



present<-function(M){
	true.pairings<-0
	pair.names<-levels(M)
	num.pairs<-length(pair.names)
	for (i in 1:num.pairs){
		paired.ind <-which(M==levels(M)[i])
		if (abs(paired.ind[1]-paired.ind[2])==num.pairs){
			true.pairings <- true.pairings + 1
		}
	}
	print(paste(true.pairings," true matches  out of ", num.pairs ," pairings"))
	return(true.pairings)
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


diff.dist.fun<-function(A,T.diff){
	P<-transition.matrix(A,dissimilarity=FALSE)
	D<-diffusion.distance(P, T.diff, directed = FALSE)
	D
}



