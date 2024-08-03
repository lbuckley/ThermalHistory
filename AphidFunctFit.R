library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)
library(viridisLite)
library(patchwork)
library(ggplot2)
library(deSolve)
library(pracma)				# contains sigmoid function
library(TrenchR)

#Analysis for English grain aphid, Sitobion avenae

setwd("/Users/laurenbuckley/ThermalHistory")
#setwd("/Users/lbuckley/ThermalHistory") #laptop

#Ma et al. 2021. Are extreme high temperatures at low or high latitudes more likely to inhibit the population growth of a globally distributed aphid?
#https://doi.org/10.1016/j.jtherbio.2021.102936
#aphid Rhopalosiphum padi
#Reproductive rate (nymphs/adult/day)
Rr= function(T, a=0.263, b=3.7, T0 =30.1) exp(a*T)-exp(a*T0-(T0-T)/b)

#Kingsolver and Woods model
RsigCG = function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    T = tfun(Time)   		## temperature determined from sine wave function
    Rf = Rm2*pracma::sigmoid(T,a=0.5,b=Tc)			## calculate equilibrium level of RNA
    Pf = Pm*pracma::sigmoid(R,a=1,b=Rc)				## calculate equilibrium level of protein
    Itpc = Rr(T)	## ingestion rate from thermal performance curve
    I = max(Itpc, -1)                         		## size-independent growth version #set minimum but check assumption
    dRdt = -(1/tauR)*(R - Rf)				## decay of RNA level toward Rf
    dPdt = -(1/tauP)*(P - Pf)				## decay of protein level toward Pf
    dMassdt = C*I - k*P              		## size-independent growth version
    G= C*I - k*P  							## growth 
    return(list(c(dRdt, dPdt, dMassdt),T,I,G))
  })
}

pars <- c(tauR = 3,  tauP = 20, Rm2 = 10, Pm=10, Tc = 35, Rc = 5, aIng = 0.02, bIng = 0.8, C = 0.8, k = 0.5)
yini <- c(R = 0, P = 0, Mass = 0.1)

#-------------
#Constant rate TPCs
#development rate
#chinese clones
dr= function(T, Tmax=34.09, a=0.13, b=4.43, c=7.65) {
  dr=exp(a*T)-exp(b-(Tmax-T)/c)
  dr[dr<0]<- 0
  return(dr)
}
  
#European clones
dr.e= function(T, Tmax=32.91, a=0.13, b=4.28, c=7.65) exp(a*T)-exp(b-(Tmax-T)/c)

#survival
sur= function(T, a=345.51, b=0.35) {
  sur=1-T/(a-b*T^2)
  sur[T>30]<- 0
  return(sur)
}

long= function(T, a=32.73, b= -0.91) 
{
  lon= a + b*T
  lon[lon<0]<- 0
  return(lon)
}
  
fec= function(T, a= -69.1, b=12.49, c= -0.34){
  fec=a +b*T +c*T^2
  fec[fec<0]<- 0
  return(fec)
}

#-------------
#estimate performance using constant TPC
#NymphDur, Longevity, Fecundity, BirthDur

perf.fun<- function(ts){
  #nymphal survival
  sur1<- mean(as.numeric(sur(ts)))
  
  #Adult longevity (days)
  long1<- mean(as.numeric(long(ts)))
  
  #Lifetime fecundity (nymphs / adult)
  fecun1<- mean(as.numeric(fec(ts)))
  
  #Developmental rate (1/days)
  dr1<- as.numeric(dr(ts))
  dr1[dr1<0] <- 0
  dr1<- mean(dr1)
  
  return( c(sur1, long1, fecun1, dr1) )
}

#-------------
#Performance detreiment

# pdet_tprev: performace detriment at previous time period
# dt: time interval
# Ea: activation energy
# https://doi.org/10.1080/02656736.2018.1558289: activation energy ΔE=5.78±0.04×10^5 J mole−1, frequency factor A=3.27±11*10^91 sec−1
# Universal gas constant 8.314	J/mol.K; https://doi.org/10.1016/j.rinp.2021.103992
# 10.3109/02656736.2011.580822

# kB: Boltzman 1.380649×10−23 J/K  ##5.673e-08 W_m-2_K-4, 
# 8·62 × 10−5 eV K−1; E = 0·65 eV, Gilooly 
# T: temperature, C
# c1: damage increase over time
# c2: multiplicative change in damage
# c3: linear increase in damage
# Topt: thermal optima

#kB= 1.380649*10^{-23}

pdet<- function(pdet_tprev, dt, Ea=0.65, R=8.62*10^{-5}, T, Topt, c1, c2, c3)  pdet_tprev + dt * exp(-Ea/(R*(T+273.15))) * (c1*pdet_tprev + c2) + c3 *dt

#put in Topt
#exp(-5.78*10^5/(8.314*(5+273.15)))= 2.829373e-109
#exp(-8.62*10^{-5}/(0.65*(5+273.15)))= ~1
  
plot(1:50, pdet(pdet_tprev=0, dt=1, Ea=5.78*10^5, R=8.314, T=1:50, Topt=30, c1=10^100, c2=10.7^100, c3=0.0), type="l")

plot(1:50, pdet(pdet_tprev=0, dt=1, R=8.62*10^{-5}, Ea=0.65, T=1:50, Topt=30, c1=100000000, c2=100000000, c3=0.0), type="l")

#performance with detriment
pf<- function(Ts, Topt, c1=1.02, c2=1.01, c3=0.0){
  
  ts<- 1:length(Ts) #time series corresponding to temperatures
  p<- fec(Ts, a= -69.1, b=12.49, c= -0.34) #time series of performances
  pdets<- c(0,rep(0, length(Ts)))
  
  #change to account for longer exposures
  
  for(time in ts){
  if(Ts[time]>Topt){    
    pdets[time+1]= pdet(pdets[time], dt=1, Ea=0.65, R=8.62*10^{-5}, T=Ts[time], Topt, c1, c2, c3)
    if(pdets[time+1]<0) pdets[time+1]<-0
    if(pdets[time+1]>1) pdets[time+1]<-1
  p[time]<- p[time]* pdets[time+1]
  }
  }
return(cbind(p,pdets[2:length(pdets)]))  
}
#-------------
#TIME SERIES

#construct temperatures
#currently sawtooth, make sine wave?

#day: increased from 27 °C to peak at 35 °C (at 09.00–15.00 h) and then decreased to 27 °C (at 15.00–21.00 h)
#temperatures in the six chambers dropped from 27 °C to different minima (13, 16, 19, 21, 23 and 25 °C) at 21.00–03.00 h and then rose to 27 °C (at 03.00–09.00 h). 

tmin= c(13, 16, 19, 21, 23, 25)

t.hrs= c(22:24,1:8)
tmin.fun= function(tmin) c( 27-1*(27-tmin)/6, 27-2*(27-tmin)/6, 27-3*(27-tmin)/6, 27-4*(27-tmin)/6, 27-5*(27-tmin)/6, tmin, 27-5*(27-tmin)/6, 27-4*(27-tmin)/6, 27-3*(27-tmin)/6, 27-2*(27-tmin)/6, 27-1*(27-tmin)/6)

temps<- matrix(NA, nrow=length(tmin), ncol=24)
for(t.k in 1:length(tmin)){
  ts= tmin.fun(tmin[t.k])
  temps[t.k, 1:8]= ts[match(1:8, t.hrs)]
  temps[t.k, 9:21]= c(27, 27+1*8/6, 27+2*8/6, 27+3*8/6, 27+4*8/6, 27+5*8/6, 35, 27+5*8/6, 27+4*8/6, 27+3*8/6, 27+2*8/6, 27+1*8/6, 27)
  temps[t.k, 22:24]= ts[match(22:24, t.hrs)]
}

temps= as.data.frame(cbind(tmin, temps))
colnames(temps)[2:25]<- 1:24

#to long format
temps.l<- melt(temps, id.vars = c("tmin"), variable.name = "hr")

#plot temperatures, varying night temperatures
ggplot(data=temps.l, aes(x=hr, y =value, color=tmin, group=tmin))+geom_line()

#-------------
temps= rep(temps.l[which(temps.l$tmin==13),"value"],7)

plot(1:50, pdet(pdet_tprev=0, dt=1, R=8.62*10^{-5}, Ea=0.65, T=1:50, Topt=30, c1=1.1, c2=10000000000, c3=-0.3), type="l")

pout<- pf(Ts=temps, Topt=20, c1=1.1, c2=20000000000, c3=0.15)
pTs<- pout[,1]

plot(1:length(temps), pTs, type="l") #, ylim=c(0,10) #performance with damage
points(1:length(temps), fec(temps, a= -69.1, b=12.49, c= -0.34), type="l", col="purple") #just performance
points(1:length(temps), temps, type="l", col="orange") #temperatures
points(1:length(temps), pout[,2], type="l", col="green") #damage


