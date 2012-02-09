

require(igraph)
require(optmatch)
source("./lib/simulation_math_util_fn.R")
source("./lib/smacofM.R")
source("./lib/oosIM.R")
source("./lib/oosMDS.R")
source("./lib/diffusion_distance.R")
source("./lib/graph_embedding_fn_many.R")

cep=TRUE
verbose= FALSE
oos=TRUE
oos.cep = TRUE

a.cep <-20


n.1 = 40
n.2 = 20
n<- n.2
m = 5 # the first m pairs are known matches ; the last n-m pairs are to-be-matched

npert = 11
nmc = 100
pert=(0:10)/10
nc.jofc.diff.p = matrix(0,npert,nmc)
nc.jofc.weighted.p = matrix(0,npert,nmc)
nc.jofc.unweighted.p = matrix(0,npert,nmc)

nc.jofc.diff.r = matrix(0,npert,nmc)
nc.jofc.weighted.r = matrix(0,npert,nmc)
nc.jofc.unweighted.r = matrix(0,npert,nmc)

nc.jofc.diff.f = matrix(0,npert,nmc)
nc.jofc.weighted.f = matrix(0,npert,nmc)
nc.jofc.unweighted.f = matrix(0,npert,nmc)



nc.cmds = matrix(0,npert,nmc)

matched.cost<-0.01


#w.vals.vec <- c(0.5,0.7,0.9,0.95)
w.vals.vec <- c(0.9)

w.max.index<-length(w.vals.vec)

matched.cost<-0.01 #If matched.cost is equal to 1, consider an unweighted graph, with edges between matched vertices
#If matched.cost is  between 0 and 1, the graph is weighted with edges between matched vertices with weights equal to matched.cost. Edges between 
# vertices of the same condition have either weight 1 or 2 according to whether they're connected according to given adjacency matrix or not.

d.dim<-4
T.diff<-1

dims.for.dist <- 1:d.dim

seed<-123

nmc<-10
for(imc in 1:nmc)
{
	G.orig<-ER(n,0.5)
	diag(G.orig)<-1
	k<-1
	repeat.counts <-1+rgeom(n,0.2)
	new.n <- sum(repeat.counts)
	int.end.indices<-cumsum(repeat.counts)
	int.start.indices<-c(1,int.end.indices+1)
	corr.list<-vector("list",n)
	G<-matrix(0,new.n,new.n)
	for (i in 1:n){
		for (j in 1:repeat.counts[i]){
			G[k,]<-rep(G.orig[i,],times=repeat.counts)
			G[,k]<-rep(G.orig[i,],times=repeat.counts)
			#	G<-perturbG(G,0.1)
			k <- k+1
		}
		corr.list[[i]] <- list(a=(int.start.indices[i]:int.end.indices[i]),b=i)
	}
	
	diag(G.orig)<-0
	diag(G)<-0
	
	for(ipert in 1:npert)
	{
		Y.emb<-NULL
		Gp<-bitflip(G.orig ,pert[ipert],pert[ipert])
		Graph.1<-graph.adjacency(G)
		Graph.2<-graph.adjacency(Gp)
		D.1<-shortest.paths(Graph.1)
		D.2<-shortest.paths(Graph.2)
		oos.sampling<-sample(1:n, size=m, replace=FALSE)
		in.sample.ind.1<-rep(TRUE,new.n)
		for ( s in 1:m){
			a<-int.start.indices[oos.sampling[s]]
			b<-int.end.indices[oos.sampling[s]]
			in.sample.ind.1[a:b]<-FALSE
		}
		
		in.sample.ind.2<-rep(TRUE,n)
		in.sample.ind.2[oos.sampling]<-FALSE
		
		#if (imc==1) print(in.sample.ind)
		
		J.1 =jofc.diffusion.dist.many(G,Gp,corr.list,
				in.sample.ind.1,in.sample.ind.2,
				d.dim=d.dim,
				w.vals.vec,
				graph.is.directed=FALSE,
				oos=oos,
				wt.matrix.1=NULL,
				wt.matrix.2=NULL,
				sep.graphs=TRUE
		)
		
		
		M = solveMarriage(J.1[[1]])
		match.perf.eval <- present.many(M,corr.list,in.sample.ind.1,in.sample.ind.2)
		nc.jofc.diff.p[ipert,imc] = mean(match.perf.eval$P)
		nc.jofc.diff.r[ipert,imc] = mean(match.perf.eval$R)
		nc.jofc.diff.f[ipert,imc] <- mean(match.perf.eval$F)
		
		
		J.2 =jofc.many(G,Gp,corr.list,
				in.sample.ind.1,in.sample.ind.2,
				d.dim=d.dim,
				w.vals.vec,
				graph.is.directed=FALSE,	
				oos=oos,          
				
				wt.matrix.1=NULL,
				wt.matrix.2=NULL,
				sep.graphs=TRUE)
		
		M = solveMarriage(J.2[[1]])
		
		match.perf.eval <- present.many(M,corr.list,in.sample.ind.1,in.sample.ind.2)
		nc.jofc.weighted.p[ipert,imc] = mean(match.perf.eval$P)
		nc.jofc.weighted.r[ipert,imc] = mean(match.perf.eval$R)
		nc.jofc.weighted.f[ipert,imc] = mean(match.perf.eval$F)
		
		
		J.3 =jofc.many(G,Gp,corr.list,
				in.sample.ind.1,in.sample.ind.2,
				d.dim=d.dim,
				w.vals.vec,	
				graph.is.directed=FALSE,
				oos=oos,
				use.weighted.graph=FALSE,
				wt.matrix.1=NULL,
				wt.matrix.2=NULL,
				sep.graphs=TRUE)
		
		M = solveMarriage(J.3[[1]])
		match.perf.eval <- present.many(M,corr.list,in.sample.ind.1,in.sample.ind.2)
		nc.jofc.unweighted.p[ipert,imc] = mean(match.perf.eval$P)
		nc.jofc.unweighted.r[ipert,imc] = mean(match.perf.eval$R)
		nc.jofc.unweighted.f[ipert,imc] =mean(match.perf.eval$F)
		
		
		
		
		
	}
}


### notice that the %correctmatches decreases as perturbation parameter increases! :-)

pdf("plot1.pdf")
colors.vec<-c( "red","blue","orange","green")
colors.vec<-colors.vec[1:4]
plot(pert,apply(nc.jofc.diff.f,1,mean),xlab="perturbation parameter",ylab="Fraction of  correct matches",ylim=c(0,1),col=colors.vec[1])
#points(pert,apply(nc.cmds,1,mean)/(n-m),xlab="perturbation parameter",ylab="Fraction of  correct matches",ylim=c(0,1),pch=2,col=colors.vec[2])

points(pert,apply(nc.jofc.weighted.f,1,mean),xlab="perturbation parameter",ylab="Fraction of  correct matches",ylim=c(0,1),pch=3,col=colors.vec[3])

points(pert,apply(nc.jofc.unweighted.f,1,mean),xlab="perturbation parameter",ylab="Fraction of  correct matches",ylim=c(0,1),pch=4,col=colors.vec[4])

legend.txt<- c("diff dist","weighted.graph","unweighted.graph")
legend(x="topright",legend=legend.txt, col =colors.vec,pch=1:4)
title("marriage o jofc")
abline(h=1/(n-m),lty=2) ### chance?  apparently not!?
abline(v=1/2,lty=2) ### chance?  gotta be!?
dev.off()



pdf("plot2.pdf")
colors.vec<-c( "red","blue","orange","green")
colors.vec<-colors.vec[1:4]
plot(pert,apply(nc.jofc.diff.p,1,mean),xlab="perturbation parameter",ylab="Fraction of  correct matches",ylim=c(0,1),col=colors.vec[1],type="l")
#points(pert,apply(nc.cmds,1,mean)/(n-m),xlab="perturbation parameter",ylab="Fraction of  correct matches",ylim=c(0,1),pch=2,col=colors.vec[2])


lines(pert,apply(nc.jofc.weighted.p,1,mean),xlab="perturbation parameter",ylab="Fraction of  correct matches",ylim=c(0,1),col=colors.vec[3])

lines(pert,apply(nc.jofc.unweighted.p,1,mean),xlab="perturbation parameter",ylab="Fraction of  correct matches",ylim=c(0,1),col=colors.vec[4])

lines(pert,apply(nc.jofc.diff.r,1,mean),xlab="perturbation parameter",ylab="Fraction of  correct matches",ylim=c(0,1),lty=2,col=colors.vec[1])


lines(pert,apply(nc.jofc.weighted.r,1,mean),xlab="perturbation parameter",ylab="Fraction of  correct matches",ylim=c(0,1),lty=2,col=colors.vec[3])

lines(pert,apply(nc.jofc.unweighted.r,1,mean),xlab="perturbation parameter",ylab="Fraction of  correct matches",ylim=c(0,1),lty=2,col=colors.vec[4])


legend.txt<- c("diff dist","weighted.graph","unweighted.graph")
legend(x="topright",legend=legend.txt, col =colors.vec,pch=1:4)
title("marriage o jofc")
abline(v=1/2,lty=2) ### chance?  gotta be!?
dev.off()





