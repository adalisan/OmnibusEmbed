
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



jofc<-function(G,Gp,
		in.sample.ind,
		d.dim,
		graph.is.directed=FALSE,
		notconnect.wt=10,
		use.weighted.graph=TRUE,
		wt.matrix.1=NULL,
		wt.matrix.2=NULL,
		sep.graphs=TRUE # if TRUE, treat two graphs separately to compute dissimilarities
#and impute W (off-diagonalblock matrix)
# if FALSE, join the graphs and compute dissimilarities from joint graph
	    
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
		G.comb<-omnibusM(G,Gp,A.M)
		Graph.M <- graph.adjacency(G.comb,
				weighted= NULL ,mode=graph.mode)
	} else{
		print("Using weighted graph for computing dissimilarities")
		# If creating  a weighted graph 
		# make the weight matrix from adjacency matrix. Those with same-condition edges have weights of wt.connect/10
		#those  with no edges have weights of wt.connect. Those with "matched edges has weights of matched.cost
		A.M<- matrix(500,n,n)
		diag(A.M) <- matched.cost
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
				Graph.2<-graph.adjacency(wt.matrix.2,weighted= TRUE , mode=graph.mode)
			}
			else{
				Graph.1<-graph.adjacency(wt.matrix.1, weighted= TRUE, mode=graph.mode)
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
	
	
	
	oos <- TRUE
	
	Embed.List<-Embed.Nodes(D.M,  in.sample.ind ,oos ,
			d=d.dim,
			wt.equalize=FALSE,
			separability.entries.w=FALSE,
			assume.matched.for.oos = FALSE)	
	J<-matrix()
	for (Y.embed in Embed.List){
		
		test.samp.size<-nrow(Y.embed)/2
		Dist=as.matrix(dist(Y.embed))[1:test.samp.size,(1:test.samp.size)+test.samp.size]
				
		J<-Dist
		
	}
	return(J)
	
}



jofc.diffusion.dist<-function(G,Gp,
		in.sample.ind,
		d.dim,
		graph.is.directed=FALSE,

	
		wt.matrix.1=NULL,
		wt.matrix.2=NULL,
		sep.graphs=TRUE # if TRUE, treat two graphs separately to compute dissimilarities
#and impute W (off-diagonalblock matrix)
# if FALSE, join the graphs and compute dissimilarities from joint graph

){
	n<-nrow(G)
	graph.mode<- ifelse(graph.is.directed,"directed","undirected")
	
	
	if (sep.graphs){
		if (is.null(wt.matrix.1)){
			D.1<-diff.dist.fun(G)
			D.2<-diff.dist.fun(Gp)
		} else{
			D.1<-diff.dist.fun(wt.matrix.1)
			D.2<-diff.dist.fun(wt.matrix.2)
		}
		D.w<- (D.1+D.2)/2
		
		D.M<- omnibusM(D.1,D.2,D.w)
    print(D.M)
	}	else{
		if (is.null(wt.matrix.1)){
			A.M<- diag(n)
			G.comb<-omnibusM(G,Gp,A.M)
			#compute dissimilarities in joint graph
			D.M<-diff.dist.fun(G.comb)
			D.M[is.infinite(D.M)]<-1E10
		} else{
			
			Wt.M<-omnibusM(wt.matrix.1,wt.matrix.2,(wt.matrix.1+wt.matrix.2)/2)
			D.M<-diff.dist.fun(Wt.M)
      D.M[is.infinite(D.M)]<-1E10
		}
		
	}
	
	
	oos <- TRUE
	
	Embed.List<-Embed.Nodes(D.M,  in.sample.ind ,oos ,
			d=d.dim,
			wt.equalize=FALSE,
			separability.entries.w=FALSE,
			assume.matched.for.oos = FALSE)	
	J<-matrix()
	for (Y.embed in Embed.List){

		test.samp.size<-nrow(Y.embed)/2
		Dist=as.matrix(dist(Y.embed))[1:test.samp.size,(1:test.samp.size)+test.samp.size]
		J<-Dist
		
	}
	return(J)
	
}









Embed.Nodes <-function(D.omnibus,in.sample.ind,oos, d,
		wt.equalize=FALSE,
		separability.entries.w=FALSE,
		assume.matched.for.oos = FALSE ){
	
	
	Y.embeds<-list()
	oos.use.imputed<- FALSE
	w.max.index<-length(w.vals)
	
	#in.sample.ind<-which(in.sample.ind)
	
	# Embed in-sample using different weight matrices (differentw values)
	if (oos){	
		n<-sum(in.sample.ind)/2
		test.samp.size<-sum(!in.sample.ind)/2
		D.in <- D.omnibus[in.sample.ind,in.sample.ind]
		init.conf=NULL
		if (sum(is.na(D.in))==0) {		
				init.conf<-cmdscale(d=D.in,k=d)
		}	
		
		
		
		
		X.embeds<-JOFC.Insample.Embed(D.in,d,w.vals,separability.entries.w,init.conf=init.conf,
				wt.equalize=wt.equalize)
		for (l in 1:w.max.index){
			if (verbose) print("OOS embedding for JOFC for w= \n")
			if (verbose) print(w.vals[l])
			
			w.val.l <- w.vals[l]
			X <- X.embeds[[l]]
			
			insample.indices<- which(in.sample.ind)
			
			#Compute Weight matrix corresponding in-sample  entries
			# Since we are going to oos-embedding, set the weights  of in-sample embedding of stress
			# We are using previous in-sample embeddings, anyway
			oos.Weight.mat.1<-matrix(0,2*n,2*n)
			
			
			
			
			#Compute Weight matrix corresponding OOS  entries
			oos.Weight.mat.2<-w.val.to.W.mat(w.val.l,(2*test.samp.size),separability.entries.w,wt.equalize)
			
			# If assume.matched.for.oos is true, we assume OOS dissimilarities are matched(in reality,
			# they are matched for the matched pairs, but unmatched for the unmatched pairs)
			# If assume.matched.for.oos is true, we ignore the dissimilarities between matched/unmatched 
			# pairs
			if (!assume.matched.for.oos){
				oos.Weight.mat.2[1:test.samp.size,test.samp.size+(1:test.samp.size)]<-0
				oos.Weight.mat.2[test.samp.size+(1:test.samp.size),(1:test.samp.size)]<-0
			}
			
			
			# if (oos.use.imputed is true) we treat the dissimiilarities between  in-sample and out-of-sample measurements
			# from different conditions like fidelity terms
			# otherwise they are ignored
			if (oos.use.imputed){
				oos.Weight.mat.w <- matrix(1-w.val.l,2*n,2*test.samp.size)
			} else{
				oos.Weight.mat.w <- rbind(cbind(matrix(1-w.val.l,n,test.samp.size), matrix(0,n,test.samp.size) ),
						cbind(matrix(0,n,test.samp.size),matrix(1-w.val.l,n,test.samp.size))
				)
			}
			oos.Weight.mat<-omnibusM(oos.Weight.mat.1,oos.Weight.mat.2,oos.Weight.mat.w)
			
			
			
			omnibus.oos.D.0 <- rbind(
					cbind(D.omnibus[in.sample.ind,in.sample.ind],D.omnibus[in.sample.ind,!in.sample.ind]),
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
		
	}
	
	Y.embeds
	
}

solveMarriage<- function(Dist){
	matches<-pairmatch(Dist)
	
	
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


diff.dist.fun<-function(A){
	P<-transition.matrix(A,dissimilarity=FALSE)
	D<-diffusion.distance(P, T.diff, directed = FALSE)
	D
}



