################################################################
#
#			File: run_VectorialCapacityPFPR.R
#			Authors: Eili Klein, Resources for the Future, Dave Smith, National Institutes of Health
#			Purpose: 
#
################################################################


# Declare the files needed
#require(odesolve)
require(deSolve) ## edited 2013-05-08
source ("functions.R") 
source ("parameters.R") 


T=40
t = seq(0,T*365, by=1)



theta = q
c.1=ent.c.1
c.2=ent.c.2/10
#ep = c(.01,.02,.03,.04,.05,.06,.07,.08,.075,.08,.09,.1,.2,.5,1,10,20)

ep_inv = c(seq(-4,4,by=.1))
ep = exp(ep_inv)/(ent.a^2*exp(-ent.g*ent.n)/ent.g)

#ep = seq(.01,50,by=.01)
vecCapacity <- matrix(NA,nrow=length(ep),ncol=24,byrow=TRUE,dimnames=list(c(ep),c("VectorialCapacity","EIR","S.1","Iw.1","Ix.1","S.2","Iw.2","Ix.2","PFPR","PFPR.1","PFPR.2","PFPR.W","PFPR.X","Immune","Immune.S","Immune.I","Clinical","Clinical.EIR","Rw","Rx","R0.w1","R0.w2","R0.x1","R0.x2")))

for(e in 1:length(ep)) {
	
	params = c(epsilon=ep[e],q=q,g=g,mu=mu,c.1=c.1,c.2=c.2,xi.1=xi.1,xi.2=xi.2,rho.1=rho.1,rho.2=rho.2,r.w.1=r.w.1,r.x.1=r.x.1,r.w.2=r.w.2,r.x.2=r.x.2)
	odeup = runMainModel(params)
	last = length(odeup[,1])
	EIR = odeup[last,]$EIR	
	V = getVectorialCapacity(ep[e])
	PFPR = odeup[last,]$Iw.1 + odeup[last,]$Iw.2 + odeup[last,]$Ix.1 + odeup[last,]$Ix.2
	PFPR.1 = odeup[last,]$Iw.1 + odeup[last,]$Ix.1
	PFPR.2 = odeup[last,]$Iw.2 + odeup[last,]$Ix.2
	PFPR.W = odeup[last,]$Iw.1 + odeup[last,]$Iw.2
	PFPR.X = odeup[last,]$Ix.1 + odeup[last,]$Ix.2
	Immune = odeup[last,]$S.2 + odeup[last,]$Iw.2 + odeup[last,]$Ix.2
	Immune.S = odeup[last,]$S.2
	Immune.I = odeup[last,]$Iw.2 + odeup[last,]$Ix.2
	
	tm = odeup[,1]
	DISC = exp(-interest*tm)
	
	if(EIR == 0) {
		Clinical = 0
		Clinical.EIR = 0
	} else {
	 	P = c.1*(odeup$Iw.1 + odeup$Ix.1) + c.2*(odeup$Iw.2 + odeup$Ix.2)
	 	h = ent.b*V*P/(1+ent.s*P)
	 	F.w = (c.1*odeup$Iw.1 + c.2*odeup$Iw.2)/P
	 	h.w = h*F.w
	 	#Clinical = sum(DISC*popSize*(odeup$S.1*h.w*xi.1 + odeup$S.2*h.w*xi.2 + odeup$Iw.1*s.1 + odeup$Iw.2*s.2 + odeup$Ix.1*s.1 + odeup$Ix.2*s.2))
	 	Clinical = odeup[last,]$S.1*h.w[last]*xi.1 + odeup[last,]$S.2*h.w[last]*xi.2 + odeup[last,]$Iw.1*s.1 + odeup[last,]$Iw.2*s.2 + odeup[last,]$Ix.1*s.1 + odeup[last,]$Ix.2*s.2
		Clinical.EIR = Clinical/EIR	
		
	}
	
	
	R0.w1 = ((ent.b*V*(1-xi.1))/(theta + rho.1 + r.w.1 + mu))*(c.1 + ((theta*c.2)/(rho.2+r.w.2+mu)))
	R0.w2 = ((ent.b*V*c.2*(1-xi.2))/(rho.2 + r.w.2 + mu))
	R0.x1 = ((ent.b*V)/(theta + r.x.1 + mu))*(c.1 + ((theta*c.2)/(r.x.2+mu)))
	R0.x2 = ((ent.b*V*c.2)/(r.x.2 + mu))
	
	S.1 = odeup[last,]$S.1
	Iw.1 = odeup[last,]$Iw.1
	Ix.1 = odeup[last,]$Ix.1
	S.2 = odeup[last,]$S.2
	Iw.2 = odeup[last,]$Iw.2
	Ix.2 = odeup[last,]$Ix.2
	
	R.w = (S.1+S.2)*(Iw.1*R0.w1 + Iw.2*R0.w2)
	R.x = (S.1+S.2)*(Ix.1*R0.x1 + Ix.2*R0.x2)
	
	vecCapacity[e,] = c(V,EIR,S.1,Iw.1,Ix.1,S.2,Iw.2,Ix.2,PFPR,PFPR.1,PFPR.2,PFPR.W,PFPR.X,Immune,Immune.S,Immune.I,Clinical,Clinical.EIR,R.w,R.x,R0.w1,R0.w2,R0.x1,R0.x2)
	cat(e, " of ", length(ep),"\n") 
}

vec = data.frame(vecCapacity)

par(mar=c(5.1, 4.1, 4.1, 4.1))
Vlog = log(vec$V)
xminmax=c(-3.20,max(Vlog))
yminmax=c(0,1)
plot(1,1,xlim=xminmax,ylim=yminmax,type="n",xlab="Log of Vectorial Capacity",ylab="Plasmodium falciparum Prevalance Rate",main="",cex.main=1,axes=F)
lines(Vlog,vec$PFPR,lty="solid",lwd=2,col="black")
lines(Vlog,vec$PFPR.2,lty="dashed",lwd=1,col="black")
lines(Vlog,vec$PFPR.1,lty="solid",lwd=1,col="black")
legend(.9,.98 , c("Total"), cex=1,ncol=1,bty="n",text.col="black")
legend(.9,.74 , c("Clinically Immune"), cex=1,ncol=1,bty="n",text.col="black")
legend(.9,.2 , c("Nonimmune"), cex=1,ncol=1,bty="n",text.col="black")
axis(2, pretty(yminmax,10))
box()
par(new=T)
yminmax=c(0,max(vec$Clinical)*2)
plot(Vlog, vec$Clinical, axes=F, ylab="", type="b", pch=20, xlab="",xlim=xminmax,ylim=yminmax)
axis(4, pretty(yminmax,10))
axis(1,pretty(xminmax,10)) 
par(mar=c(5.1, 4.1, 4.1, 1.1))
mtext("Clinical Incidence Rate",side=4,cex=1)
par(mar=c(5.1, 4.1, 4.1, 2.1))
legend(.9,.014 , c("Clinical Incidence Rate"), cex=1,ncol=1,bty="n",text.col="black")



