################################################################
#
#			File: functions.R
#			Authors: Eili Klein, Resources for the Future, Dave Smith, National Institutes of Health
#			Purpose: functions for use in running model
#
################################################################

############ runInitODEmodel ################### 
# This is the Initial differential equations section that is used to compute the population changes when called by lsoda
# Parameters to be passed include the time, the y values (or different populations) and the parameters
runInitODEmodel = function(t,y,p)
{
	with(as.list(c(p)), {
		
		S.1 = y[1]
		Iw.1 = y[2]
		S.2 = y[3]
		Iw.2 = y[4]
					
		#Determine the vectorial capacity based on the number of mosquitos that emerge each day - the parameter that we use to change the EIR	
		V = getVectorialCapacity(epsilon) 

		H = sum(y)
	  	#browser()
	 	#Determine Transmission and immunity loss rates for each infection
	 	P = c.1*(Iw.1) + c.2*(Iw.2)
	 	h = ent.b*V*P/(1+ent.s*P)
	 	F.w = (c.1*Iw.1 + c.2*Iw.2)/P
	 	h.w = h*F.w

	  	# Births match deaths
  		B = H*mu
  		
		# Stage 1
		dS.1 = mu - S.1*(h.w*(1-xi.1) + mu) + S.2*g + Iw.1*(rho.1 + r.w.1) 
		dIw.1 = S.1*h.w*(1-xi.1) - Iw.1*(rho.1 + r.w.1 + q + mu)
		
		# Stage 2
		dS.2 = -S.2*(h.w*(1-xi.2) + mu) - S.2*g + Iw.2*(rho.2 + r.w.2)
		dIw.2 = S.2*h.w*(1-xi.2) - Iw.2*(rho.2 + r.w.2 + mu) + Iw.1*q
	
		#Return for next calculation
		list(c(dS.1, dIw.1,  dS.2, dIw.2))
	})
}

############ runODEmodel ################### 
# This is the main differential equations section that is used to compute the population changes when called by lsoda
# Parameters to be passed include the time, the y values (or different populations) and the parameters
runODEmodel = function(t,y,p)
{
	with(as.list(c(p)), {
		
		S.1 = y[1]
		Iw.1 = y[2]
		Ix.1 = y[3]
		S.2 = y[4]
		Iw.2 = y[5]
		Ix.2 = y[6]
					
		#Determine the vectorial capacity based on the number of mosquitos that emerge each day - the parameter that we use to change the EIR	
		V = getVectorialCapacity(epsilon) 

		H = sum(y)
	  	#browser()
	 	#Determine Transmission and immunity loss rates for each infection
	 	P = c.1*(Iw.1 + Ix.1) + c.2*(Iw.2 + Ix.2)
	 	h = ent.b*V*P/(1+ent.s*P)
	 	F.w = (c.1*Iw.1 + c.2*Iw.2)/P
	 	h.w = h*F.w
	 	F.x = (c.1*Ix.1 + c.2*Ix.2)/P
	 	h.x = h*F.x

	  	# Births match deaths
  		B = H*mu
  		
		# Stage 1
		dS.1 = mu - S.1*(h.w*(1-xi.1) + h.x + mu) + S.2*g + Iw.1*(rho.1 + r.w.1) + r.x.1*Ix.1
		dIw.1 = S.1*h.w*(1-xi.1) - Iw.1*(rho.1 + r.w.1 + q + mu)
		dIx.1 = h.x*S.1 - Ix.1*(r.x.1 + q + mu)
		
		# Stage 2
		dS.2 = -S.2*(h.w*(1-xi.2) + h.x + mu) - S.2*g + Iw.2*(rho.2 + r.w.2) + r.x.2*Ix.2
		dIw.2 = S.2*h.w*(1-xi.2) - Iw.2*(rho.2 + r.w.2 + mu) + Iw.1*q
		dIx.2 = h.x*S.2 - Ix.2*(r.x.2 + mu) + Ix.1*q
	
		#Return for next calculation
		list(c(dS.1, dIw.1, dIx.1, dS.2, dIw.2, dIx.2))
	})
}



# Determines the transmission rate based on the number of infections present, the transmission between humans and mosquitos and the vectorial capacity
getTransmissionRate = function(V,I.1,I.2,c.1,c.2,H){V*(c.1*I.1 + c.2*I.2)/H}	

# Gets the initial infection levels based on the parameters and runs for the number of years requested
findInitialInfections = function(params,years) {
	T=years
	t = seq(0,T*365, by=1)

	# First determine the equilibrium of the model with no demand or resistance
	inits = c(S.1=.7, Iw.1=.3, S.2=0, Iw.2=0)
  	odeup.init = lsoda(inits,t,runInitODEmodel, c(params))	
  	
  	#This keeps system from blowing up when values are off, but graphs will be all out of whack
  	last = length(odeup.init[,1])
  	inits = c(S.1=odeup.init[last,2],Iw.1=odeup.init[last,3],Ix.1=0,S.2=odeup.init[last,4],Iw.2=odeup.init[last,5],Ix.2=0)
  	inits		
}

# This finds the initial infection level and then runs the model
runMainModel = function(params,delayRes=FALSE,graphFull=FALSE) {
	# Run model to eq and find initial infection levels
	inits = findInitialInfections(params,250)
	paramsList = data.frame(as.list(params))
	# Determine EIR
	S.1 = inits[1]
	names(S.1) <- NULL
	if(isTRUE(all.equal(S.1,1,tolerance=1e-9)) || inits[2] < 0) {
		P = 0
	} else {
		P = paramsList$c.1*inits[2] + paramsList$c.2*inits[5]
	}
	
	
	EIR = (params[1]*ent.a^2*P*exp(-ent.g*ent.n)/(ent.g + ent.a*P))*365
	
	inits[2]=inits[2]-inf.init # create resistance by reducing Iw.1 ....
	inits[3]=inits[3]+inf.init # ... and increasing Ix.1

	# run model, append EIR info and return
	odeup = data.frame(lsoda(inits,t,runODEmodel, c(params,init=FALSE,resDelay=FALSE)))
	odeup.return = cbind(odeup,EIR=rep(EIR,length(t)))
	odeup.return
}

findInitialInfections2 = function(params,years) {
	T=years
	t = seq(0,T*365, by=1)

	# First determine the equilibrium of the model with no demand or resistance
	inits = c(S.1=.7, Iw.1=.3, S.2=0, Iw.2=0)
  	odeup.init = lsoda(inits,t,runInitODEmodel2, c(params))	
  	
  	#This keeps system from blowing up when values are off, but graphs will be all out of whack
  	last = length(odeup.init[,1])
  	inits = c(S.1=odeup.init[last,2],Iw.1=odeup.init[last,3],Ix.1=0,S.2=odeup.init[last,4],Iw.2=odeup.init[last,5],Ix.2=0)
  	inits		
}

# This finds the initial infection level and then runs the model
runMainModel2 = function(params,delayRes=FALSE,graphFull=FALSE) {
	# Run model to eq and find initial infection levels
	inits = findInitialInfections2(params,250)
	paramsList = data.frame(as.list(params))
	# Determine EIR
	S.1 = inits[1]
	names(S.1) <- NULL
	if(isTRUE(all.equal(S.1,1,tolerance=1e-9)) || inits[2] < 0) {
		P = 0
	} else {
		P = paramsList$c.1*inits[2] + paramsList$c.2*inits[5]
	}
	
	
	EIR = (params[1]*ent.a^2*P*exp(-ent.g*ent.n)/(ent.g + ent.a*P))*365
	
	inits[2]=inits[2]-inf.init # create resistance by reducing Iw.1 ....
	inits[3]=inits[3]+inf.init # ... and increasing Ix.1

	# run model, append EIR info and return
	odeup = data.frame(lsoda(inits,t,runODEmodel2, c(params,init=FALSE,resDelay=FALSE)))
	odeup.return = cbind(odeup,EIR=rep(EIR,length(t)))
	odeup.return
}

# graph the results of main model, first graph is all the equations, second is just S,Iw and Ix
graphMainModel = function(odeup,maintitle,S1=TRUE,Iw1=TRUE,Ix1=TRUE,S2=TRUE,Iw2=TRUE,Ix2=TRUE,Sall=TRUE,Iall=TRUE,legX=20,legY=1) {

op <- par(col="black",lty="solid") # the whole list of settable par's.

#par(ask=TRUE)

with(as.list(odeup),{
	time = time/365
	xminmax=c(min(time),max(time))
	yminmax=c(0,1)
	legTex = NA
	legLty = NA
	plot(1,1,xlim=xminmax,ylim=yminmax,type="n",xlab="Time (years)",ylab="",main=maintitle,cex.main=.9,axes=TRUE)
	if (S1) {
		lines(time,S.1,col=colorChoices[1],lty="solid",lwd=.5)
		legTex = "S.1"
		legLty = "solid"
		legLwd = 1
		legCol = colorChoices[1]
	}
	if (Iw1) {
		lines(time,Iw.1,col=colorChoices[2],lty="solid",lwd=.5)
		if (is.na(legTex[1])) {
			legTex = "Iw.1"
			legLty = "solid"
			legLwd = 1
			legCol = colorChoices[2]
		} else {
			legTex = c(legTex,"Iw.1")
			legLty = c(legLty,"solid")
			legLwd = c(legLwd,1)
			legCol = c(legCol,colorChoices[2])
		}
	}
	if (Ix1) {
		lines(time,Ix.1,col=colorChoices[3],lty="solid",lwd=.5)
		if (is.na(legTex[1])) {
			legTex = "Ix.1"
			legLty = "solid"
			legLwd = 1
			legCol = colorChoices[3]	
		} else {
			legTex = c(legTex,"Ix.1")
			legLty = c(legLty,"solid")
			legLwd = c(legLwd,1)
			legCol = c(legCol,colorChoices[3])
		}
	}
	if (S2) {
		lines(time,S.2,col=colorChoices[1],lty="dashed",lwd=.5)
		if (is.na(legTex[1])) {
			legTex = "S.2"
			legLty = "dashed"
			legLwd = 1
			legCol = colorChoices[1]	
		} else {
			legTex = c(legTex,"S.2")
			legLty = c(legLty,"dashed")
			legLwd = c(legLwd,1)
			legCol = c(legCol,colorChoices[1])
		}
	}
	if (Iw2) {
		lines(time,Iw.2,col=colorChoices[2],lty="dashed",lwd=.5)
		if (is.na(legTex[1])) {
			legTex = "Iw.2"
			legLty = "dashed"
			legLwd = 1
			legCol = colorChoices[2]	
		} else {
			legTex = c(legTex,"Iw.2")
			legLty = c(legLty,"dashed")
			legLwd = c(legLwd,1)
			legCol = c(legCol,colorChoices[2])
		}
	}
	if (Ix2) {
		lines(time,Ix.2,col=colorChoices[3],lty="dashed",lwd=.5)
		if (is.na(legTex[1])) {
			legTex = "Ix.2"
			legLty = "dashed"
			legLwd = 1
			legCol = colorChoices[3]	
		} else {
			legTex = c(legTex,"Ix.2")
			legLty = c(legLty,"dashed")
			legLwd = c(legLwd,1)
			legCol = c(legCol,colorChoices[3])
		}
	}
	#plot overall levels as well
	if (Sall) {
		lines(time,S.1+S.2,col="red",lty="solid", lwd=2)
		if (is.na(legTex[1])) {
			legTex = "S - all"
			legLty = "solid"
			legLwd = 2
			legCol = "red"	
		} else {
			legTex = c(legTex,"S - all")
			legLty = c(legLty,"solid")
			legLwd = c(legLwd,2)
			legCol = c(legCol,"red")
		}
	}
	if (Iall) {
		lines(time,Iw.1+Iw.2+Ix.1+Ix.2,col="red",lty="dashed",lwd=2)
		if (is.na(legTex[1])) {
			legTex = "I - all"
			legLty = "dashed"
			legLwd = 2
			legCol = "red"
		} else {
			legTex = c(legTex,"I - all")
			legLty = c(legLty,"dashed")
			legLwd = c(legLwd,2)
			legCol = c(legCol,"red")
		}
	}
		
	par(col="black",lty="solid",lwd=1)
	legend(legX, legY , legTex,lty=legLty, lwd=legLwd, col = legCol, cex=.8,ncol=2)

})

#Reset Par
par(op)

}
