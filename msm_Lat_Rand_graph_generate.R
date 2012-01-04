require(msm)

K<-3 #Num of states

num.V<- 100 #Num. of Vertices

T<-10 #Determines time from start to end for the time series of graphs
Trans.Mat<-cbind(c(0,0.5,0.5),
				c(0.3,0,  0.2))
				c(0.1,0.2, 0 ))

r<-5


gen.Graph.Pair<- function(Traj.Ens,
							inst.ensemble #The index of the trajectory  
					){
P.0<-matrix(0,num.V,num.V)
P.1<-matrix(0,num.V,num.V)
Pi.1<-matrix(0,num.V,K)
Init.P <-
for (v in 1:num.V){
Traj.Ens<-sim.msm(q=Trans.Mat,maxtime=T+r)
start.time<-T-1
window<- (Traj.Ens$times>start.time) & (Traj.Ens$times<=T)

#include starting state also
window.start.ind<-which(window)[1]-1
window[windows.start.ind]<-TRUE   

Traj.in.window<-list(states=Traj.Ens$states[window],times=Traj.Ens$times[window])
num.transition <- length(Traj.in.window$times)-1
state.dur<-Traj.in.window$times[2:(num.transition+1)] -Traj.in.window$times[1:num.transition]
for (k in 1:K){
	times.in.k<-which(Traj.in.window$times[v]==k)
	Pi.1[v,k]<-sum(state.dur[times.in.k])
	
	}
	start.time<-T+r-1
	window<- (Traj.Ens$times>start.time) & (Traj.Ens$times<=T+r)
	window.start.ind<-which(window)[1]-1
	window[windows.start.ind]<-TRUE   
	Traj.in.window<-list(states=Traj.Ens$states[window],times=Traj.Ens$times[window])
	num.transition <- length(Traj.in.window$times)-1
	state.dur<-Traj.in.window$times[2:(num.transition+1)] -Traj.in.window$times[1:num.transition]
	for (k in 1:K){
		times.in.k<-which(Traj.in.window$times[v]==k)
		Pi.1[v,k]<-sum(state.dur[times.in.k])
		
	}
	
	
	
	
	
	
}
P.1[i,j]<-P.1[j,i] <- Pi.0[i,2:K]%.%Pi.0[j,2:K]	
}


