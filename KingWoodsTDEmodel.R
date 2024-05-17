### Kingsolver and Woods Am Nat August 2015:  Basic time-dependent effects model
### In addition, all parameters defined.

# load libraries
library(pracma)				# contains sigmoid function
library(deSolve)			# contains function ode

###### Variables calculated by function

# exported variables
# R		'RNA' for HSPs
# P		'protein' for HSPs
# Mass	body mass of caterpillar
# T		temperature, which is variable within the function, driven by sine function

# internal variables
# Rf		the equilibrium level of R, all else being equal
# Pf		the equilibrium level of P, all else being equal
# Itpc		relative rate of ingestion based on current temperature (taken from the TPC)
# I			actual ingestion rate incorporating Itpc and scaled to body mass
# a, b		parameters of the sigmoidal functions; handled as internal constants

###### Parameters handed to function

# tauR	time constant for decay of R
# tauP	time constant for decay of P
# Rm		maximum level of R possible
# Pm		maximum level of P possible
# Tc		temperature at the inflection point for the  sigmoid function for Rf 
# Rc		level of R at the inflection point for the  sigmoid function for Pf
# C		conversion efficiency for turning ingested food into body mass
# k 		a scaling coefficient that determines how strongly levels of P affect growth

###### External constants for ingestion rate function (Frazier et al 2006 model)

# Topt	    optimal temperature for ingestion curve
# rho,sig   breadth parameters for ingestion curve
# Imax	    maximum rate of ingestion


###### Main function

RsigCG = function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    T = 30 + 10*sin(Time*(pi/12))   		## temperature determined from sine wave function
    Rf = Rm*sigmoid(T,a=0.5,b=Tc)			## calculate equilibrium level of RNA
    Pf = Pm*sigmoid(R,a=1,b=Rc)				## calculate equilibrium level of protein
    Itpc = Imax*exp(-exp((rho*(T-Topt))-6)-sig*(T-Topt)^2)	## ingestion rate from thermal performance curve
    #I = Itpc*aIng*Mass^bIng 
    I = Itpc                         		## size-independent growth version
    dRdt = -(1/tauR)*(R - Rf)				## decay of RNA level toward Rf
    dPdt = -(1/tauP)*(P - Pf)				## decay of protein level toward Pf
    dMassdt = C*I - k*P              		## size-independent growth version
    #dMassdt = C*I - k*P*(Mass^bIng)
    G= C*I - k*P  							## growth 
    return(list(c(dRdt, dPdt, dMassdt),T,I,G))
  })
}

###### 

# set external constants
Topt=34;Imax=5; rho=0.9;sig=0.005

# run
times <- seq(0, 168, by = 0.1)

pars <- c(tauR = 3,  tauP = 20, Rm = 10, Pm=10, Tc = 35, Rc = 5, aIng = 0.02, bIng = 0.8, C = 0.8, k = 0.5)
yini <- c(R = 0, P = 0, Mass = 0.1)
out2 <- ode(func = RsigCG, y = yini, parms = pars, times = times)

#convert to dataframe
out2.df <- data.frame(out2)
colnames(out2.df)[5:7]<- c("T","I","G") 
summary(out2.df)

#----
#plot

#to long format
out.l<- melt(out2.df, id.vars = c("time"), variable.name = "trait")

#take out mass and adjust temperature for plotting
out.plot<- out2.df
out.plot$T<- out.plot$T /10
out.l<- melt(out.plot, id.vars = c("time"), variable.name = "trait")
out.l<-out.l[-which(out.l$trait=="Mass"),]

ggplot(data=out.l, aes(x=time, y =value, color=trait))+geom_line()

#-------------
#Figure 3 constant temp

temps<- c(15, 20, 25, 30, 35, 40)

for(t.ind in 1:6){

RsigCG = function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    T = rep(temps[t.ind], length(Time))   		## temperature determined from sine wave function
    Rf = Rm*sigmoid(T,a=0.5,b=Tc)			## calculate equilibrium level of RNA
    Pf = Pm*sigmoid(R,a=1,b=Rc)				## calculate equilibrium level of protein
    Itpc = Imax*exp(-exp((rho*(T-Topt))-6)-sig*(T-Topt)^2)	## ingestion rate from thermal performance curve
    #I = Itpc*aIng*Mass^bIng 
    I = Itpc                         		## size-independent growth version
    dRdt = -(1/tauR)*(R - Rf)				## decay of RNA level toward Rf
    dPdt = -(1/tauP)*(P - Pf)				## decay of protein level toward Pf
    dMassdt = C*I - k*P              		## size-independent growth version
    #dMassdt = C*I - k*P*(Mass^bIng)
    G= C*I - k*P  							## growth 
    return(list(c(dRdt, dPdt, dMassdt),T,I,G))
  })
}

out <- ode(func = RsigCG, y = yini, parms = pars, times = times)

if(t.ind==1) out3<- out
if(t.ind>1) out3<- rbind(out3, out)

} #end temp loop

#convert to dataframe
out3.df <- data.frame(out3)
colnames(out3.df)[5:7]<- c("T","I","G") 
summary(out3.df)

#plot
#by time
ggplot(data=out3.df, aes(x=time, y =G, color=factor(T)))+geom_line()
#by temp
ggplot(data=out3.df[out3.df$time %in% c(0, 10, 20, 40, 80),], aes(x=T, y =G, color=factor(time)))+geom_line()



