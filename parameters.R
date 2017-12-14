################################################################
#
#			File: parameters.R
#			Author: Eili Klein, Resources for the Future
#			Purpose: Sets main model parameters
#
################################################################


#The economic discounting rate
interest=0.03

#The population size
Pop=10^6

#Time
T=40
T.delay = 10
t = seq(0,T*365, by=1)
t.delay = seq(0,T.delay*365, by=1)

#initial resistant infections
inf.init = 1/1000000 #.000001

################################################################
#
#			Sets entomologic model parameters
#
################################################################

#Entomological Parameters
#Adults that emerge, per human, per day
ent.epsilon=10
#ent.epsilon=0.5
#Human biting rate
ent.a=0.3
#Transmission efficiencies
ent.b=0.8
ent.c.1=0.5
ent.c.2=0.5
#Mosquito Death Rate 
ent.g=1/10
#Days to sporogony
ent.n=10
ent.s = ent.a/ent.g

# Get Vectorial Capacity
# This is a function which generates the vectorial capacity based on model parameters.  
# Note: This actually includes b - the transmission efficiency from mosquitos to humans - which means it is not quite vectorial capacity, but since b would just be multiplied by this afterwords, it was included here
getVectorialCapacity = function(epsilon) {
	epsilon*ent.a^2*exp(-ent.g*ent.n)/ent.g	 
}



# Fraction of infections that are immediatly clinical, treated and don't transmit
xi.1 = .3
xi.2 = .01

# Rate that infections become clinical and are treated
#  rate symptoms arise
s.1 = .025
s.2 = .01
#  fraction treated rate
f.1 = .2
f.2 = .2
# Rho
rho.1 = f.1 * s.1 # = 1/200
rho.2 = f.2 * s.2 / 2# = 1/500

#Background Mortality
mu = 1/(45*365)

# Immunity gain and loss
q=1/(10*365)
g=1/(2*365)

#Disease Induced Mortality
alpha.1 = 180/100000/365
alpha.2 = 180/100000/365

# Natural clearance rate
spon.days = 1/ent.b * 165
fitnessCost = 2.6
r.w.1 = 1/spon.days
r.x.1 = 1/spon.days * fitnessCost
r.w.2 = 1/spon.days
r.x.2 = 1/spon.days * fitnessCost

interest = .03
popSize = 1000000

################################################################
#
#		Graph Options
#
################################################################

colorChoices = c("black","blue","gray","red","pink","yellow","dodgerblue","forestgreen","green","chocolate","gold","deepskyblue","sienna")

