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

ep_inv = c(seq(-6,6,by=.1))
ep = exp(ep_inv)/(ent.a^2*exp(-ent.g*ent.n)/ent.g)
#ep = c(seq(.01,.09,by=.01),seq(.1,9.9,by=.1),seq(10,200,by=1))
#ep = seq(.01,50,by=.005)
vecCapacity <- matrix(NA,nrow=length(ep),ncol=20,byrow=TRUE,dimnames=list(c(ep),c("VectorialCapacity","S.1","Iw.1","Ix.1","S.2","Iw.2","Ix.2","Rw","Rx","R0.w1","R0.w2","R0.x1","R0.x2","R0.x","Ix1.New","Ix2.new","Ixnew","Iw1.New","Iw2.New","Iwnew")))

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
	R0.w = R0.w1
	R.w = (S.1+S.2)*((1-fracImmune)*R0.w1 + fracImmune*Iw.2*R0.w2)
	R.x = (S.1+S.2)*R0.x
	
	inits = c(S.1=S.1,Iw.1=Iw.1,Ix.1=(1/1000000),S.2=S.2,Iw.2=Iw.2,Ix.2=0)
	odeup = data.frame(lsoda(inits,t,runODEmodel, c(params,init=FALSE,resDelay=FALSE))	)
  	last = length(odeup[,1])
		Ix.1new = odeup[last,]$Ix.1
		Ix.2new = odeup[last,]$Ix.2
		Iw.1new = odeup[last,]$Iw.1
		Iw.2new = odeup[last,]$Iw.2
	vecCapacity[e,] = c(V,S.1,Iw.1,Ix.1,S.2,Iw.2,Ix.2,R.w,R.x,R0.w1,R0.w2,R0.x1,R0.x2,R0.x,Ix.1new,Ix.2new,Ix.1new+Ix.2new,Iw.1new,Iw.2new,Iw.1new+Iw.2new)
		
	IXall = Ix.1new + Ix.2new
	IWall = Iw.1new + Iw.2new
	cat(V,": (",IXall+IWall,") Iw:",IWall," | Ix:",IXall," | R0w:",R0.w," | R0x:",R0.x," | Rw:",R.w," | Rx:",R.x," | ","\n")
	#browser()
}

vec = data.frame(vecCapacity)

Vlog = log(vec$V)
line1 = rep(1,length(Vlog))

j = 0
for(i in 1:length(vec$R0.x)) {
	if(vec$R0.x[i] < 1) j = i
}

par(mar=c(5.1, 4.1, 4.1, 4.1))
Vlog_2 = Vlog[-1:-j]
line1 = rep(1,length(Vlog_2))
Rx_2 = vec$Rx[-1:-j]
R0.x_2 = vec$R0.x[-1:-j]
Iw.2_2 = vec$Iw.2[-1:-j]
Ix.2_2 = vec$Ix2.new[-1:-j]
Ix = vec$Ixnew[-1:-j]
Iw = vec$Iwnew[-1:-j]
xminmax=range(Vlog_2)
yminmax=c(0,max(Rx_2)*2)
plot(1,1,xlim=xminmax,ylim=yminmax,type="n",xlab="Log of Vectorial Capacity",ylab="Reproductive Rate",main="",cex.main=1,axes=F)
lines(Vlog_2,line1,lty="dashed",col="black",lwd=1)
lines(Vlog_2,Rx_2,lty="solid",lwd=1,col="black")
lines(Vlog_2,R0.x_2,lty="dotdash",lwd=1,col="black")
legend(-2.5,3 , c(expression(R[0]),"resistant parasite"), cex=1.1,ncol=1,bty="n",text.col="black",xjust=0)
legend(1.7,.88 , c(expression(R[x]),"resistant parasite"), cex=1.1,ncol=1,bty="n",text.col="black")

axis(2, pretty(yminmax,10))
yminmax=c(0,max(Ix.2_2)*2)

par(new=T,lwd=1)
#plot(Vlog_2, Iw.2_2, axes=F, ylab="", type="b", pch=19, xlab="",col="black", ylim=c(0,1))
plot(Vlog_2, Iw.2_2, axes=F, ylab="", type="b", pch=20, xlab="",col="black", ylim=c(0,1))
#par(new=T,lwd=1)
#plot(Vlog_2, Ix, axes=F, ylab="", type="b", pch=21, xlab="",col="black", ylim=c(0,1))
#lines(Vlog_2,Ix.2_2,lty="dashed",col="black",lwd=1)
axis(4, pretty(yminmax,10))
axis(1,pretty(xminmax,10)) 
box()
#legend(2,.7 , c("Rx","R0x","Iw2"),lty=c("solid","dashed","solid"), lwd=c(2,2,1), pch = c(-1, -1, 1), cex=.8,ncol=1)
par(mar=c(5.1, 4.1, 4.1, 1.1),lwd=1)
mtext("Plasmodium falciparum Parasite Rate",side=4,cex=1.1)
par(mar=c(5.1, 4.1, 4.1, 2.1))
legend(2,.9 , c("Proportion clinically immune"), cex=1,ncol=1,bty="n",text.col="black")
#legend(2,.1 , c("all resistant"), cex=1,ncol=1,bty="n",text.col="black")

