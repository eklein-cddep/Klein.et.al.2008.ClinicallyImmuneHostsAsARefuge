################################################################
#
#			File: run_VectorialCapacityPFPR.R
#			Authors: Eili Klein, Resources for the Future, Dave Smith, National Institutes of Health
#			Purpose: 
#
################################################################


# Declare the files needed
require(odesolve)
source ("functions.R") 
source ("parameters.R") 

t = seq(0,250*365, by=1)

theta = q
c.1=ent.c.1
c.2=ent.c.2/10

ep_inv = c(seq(-6,6,by=.01))
ep = exp(ep_inv)/(ent.a^2*exp(-ent.g*ent.n)/ent.g)
#ep = c(seq(.01,.09,by=.01),seq(.1,9.9,by=.1),seq(10,200,by=1))
#ep = seq(.01,50,by=.005)
fitnessCost_test = seq(1,2.5,by=.25)
Rhold <- matrix(NA,nrow=length(ep),ncol=length(fitnessCost_test)+1,byrow=T,dimnames=list(c(ep),c("V",fitnessCost_test)))
#vecCapacity <- matrix(NA,nrow=length(ep),ncol=20,byrow=TRUE,dimnames=list(c(ep),c("VectorialCapacity","S.1","Iw.1","Ix.1","S.2","Iw.2","Ix.2","Rw","Rx","R0.w1","R0.w2","R0.x1","R0.x2","R0.x","Ix1.New","Ix2.new","Ixnew","Iw1.New","Iw2.New","Iwnew")))

for(j in 1:length(fitnessCost_test)) {
	r.x.1 = r.w.1*fitnessCost_test[j]
	r.x.2 = r.w.2*fitnessCost_test[j]
	for(e in 1:length(ep)) {
		
		params = c(epsilon=ep[e],q=q,g=g,mu=mu,c.1=c.1,c.2=c.2,xi.1=xi.1,xi.2=xi.2,rho.1=rho.1,rho.2=rho.2,r.w.1=r.w.1,r.x.1=r.x.1,r.w.2=r.w.2,r.x.2=r.x.2)
		# First determine the equilibrium of the model with no demand or resistance
		inits = c(S.1=.7, Iw.1=.3,  S.2=0, Iw.2=0)
	  	odeup = data.frame(lsoda(inits,t,runInitODEmodel, c(params,init=TRUE,resDelay=FALSE))	)
	  	
		last = length(odeup[,1])
		V = getVectorialCapacity(ep[e])
		
		R0.w1 = ((ent.b*V*(1-xi.1))/(theta + rho.1 + r.w.1 + mu))*(c.1 + ((theta*c.2)/(rho.2+r.w.2+mu)))
		R0.w2 = ((ent.b*V*c.2*(1-xi.2))/(rho.2 + r.w.2 + mu))
		R0.x1 = ((ent.b*V)/(theta + r.x.1 + mu))*(c.1 + ((theta*c.2)/(r.x.2+mu)))
		R0.x2 = ((ent.b*V*c.2)/(r.x.2 + mu))
		
		S.1 = odeup[last,]$S.1
		Iw.1 = odeup[last,]$Iw.1
		Ix.1 = 0
		S.2 = odeup[last,]$S.2
		Iw.2 = odeup[last,]$Iw.2
		Ix.2 = 0
		fracImmune = Iw.2+S.2
		R0.x = ((1-fracImmune)*R0.x1 + fracImmune*R0.x2)
		R.w = (S.1+S.2)*((1-fracImmune)*R0.w1 + fracImmune*Iw.2*R0.w2)
		R.x = (S.1+S.2)*R0.x
		Rhold[e,1]=V
		Rhold[e,j+1]=R.x
		
	}
	
}


vec = data.frame(Rhold)
Vlog = log(vec$V)
xminmax=c(min(Vlog),max(Vlog))
yminmax=c(0,max(vec$X1)*1.1)
plot(1,1,xlim=xminmax,ylim=yminmax,type="n",xlab="Vectorial Capacity",ylab="Effective Reproductive Rate",main="",cex.main=1,axes=TRUE)
lines(Vlog,vec$X1,lty="solid",lwd=1)
lines(Vlog,vec$X1.25,lty="solid",lwd=1)
lines(Vlog,vec$X1.5,lty="solid",lwd=1)
lines(Vlog,vec$X1.75,lty="solid",lwd=1)
lines(Vlog,vec$X2,lty="solid",lwd=1)
lines(Vlog,vec$X2.5,lty="solid",lwd=1)
legend(-2.85,4.25 , c("Base"), cex=1,ncol=1,bty="n",text.col="black")
legend(-2.85,3.5 , c("20%"), cex=1,ncol=1,bty="n",text.col="black")
legend(-2.85,3 , c("33%"), cex=1,ncol=1,bty="n",text.col="black")
legend(-2.85,2.6 , c("43%"), cex=1,ncol=1,bty="n",text.col="black")
legend(-2.85,2.32 , c("50%"), cex=1,ncol=1,bty="n",text.col="black")
legend(-2.85,1.92 , c("60%"), cex=1,ncol=1,bty="n",text.col="black")


browser()
