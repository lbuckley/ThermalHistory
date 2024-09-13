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
#Performance detriment

# pdet_tprev: performance detriment at previous time period
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


pdet<- function(pdet_tprev, dt, Ea=0.65, R=8.62*10^{-5}, T, Topt, c1, c2, c3)  pdet_tprev + dt * exp(-Ea/(R*((T-Topt)+273.15))) * (c1*pdet_tprev + c2) + c3*(Topt-T) *dt

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

#==============================
# Test fits

setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/out/")
#setwd("/Users/lbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/out/") 

temps.all<- read.csv("TempTimeSeries.csv")
PerfDat<- read.csv("PerformanceData.csv")

#fedundity TPC
fec= function(T, a= -69.1, b=12.49, c= -0.34){
  fec=a +b*T +c*T^2
  fec[fec<0]<- 0
  return(fec)
}

#FUNCTIONS
#damage
# c1: damage increase over time
# c2: multiplicative change in damage
# c3: linear increase in damage

damagenew<- function(damage, T, c1, c2, c3, dt, Ea=0.65, R=8.62*10^{-5})  
  {damagenew= damage + dt * exp(-Ea/(R*(T+273.15))) *10^9* (c1*damage + c2) + dt*c3
  if(damagenew<0) damagenew<-0
  if(damagenew>1) damagenew<-1
  return(damagenew)
}

#compute performance

perf<- function(series,c1,c2,c3,scale)  {
  perf=NA
  damage=0
  for(i in 1:length(series)){
    damage=damagenew(damage,T=series[i],c1=c1,c2=c2,c3=c3,dt=1)
    perf= fec(series[i])*(1-damage)
    if(i==1) perf.all=perf
    if(i>1) perf.all=c(perf.all, perf)
  }
  return(perf.all*scale)
}

perf.nodamage<- function(series,c1,c2,c3,scale)  {
  perf=NA
  damage=0
  for(i in 1:length(series)){
    perf= fec(series[i])
    if(i==1) perf.all=perf
    if(i>1) perf.all=c(perf.all, perf)
  }
  return(perf.all*scale)
}

computeperf<- function(series,c1,c2,c3,scale,printdam=FALSE)  {
  perf=0
  damage=0
  for(i in 1:length(series)){
    damage=damagenew(damage,T=series[i],c1=c1,c2=c2,c3=c3,dt=1)
  perf= perf + fec(series[i])*(1-damage)
  }
return(perf*scale)
}

#-----------
#plot parameter values
ts= seq(1, 30, 0.5)
temps= c(ts, rev(ts),ts, rev(ts))

#make parameter combinations
cs<- expand.grid(c1= seq(-0.2, 0.2, 0.1), c2= seq(.6, 1.6, 0.2), c3= seq(-0.004, -0.002, 0.001),
                 scale= seq(0.003, 0.004, 0.001) )

for(k in 1:nrow(cs)){
  p1= perf(temps, c1=cs[k,1], c2=cs[k,2], c3=cs[k,3], scale=cs[k,4])
  ps= cbind(time=1:length(temps), temps, p1, k, cs[k,])
  if(k==1) ps.all<- ps
  if(k>1) ps.all<- rbind(ps.all, ps)
}

ggplot(data=ps.all, aes(x=time, y =p1, color=c1,lty=factor(scale), group=k))+geom_line()+facet_grid(c2~c3) 

plot(1:length(temps), temps, type="l")
p1.nd= perf.nodamage(temps, c1=cs[k,1], c2=cs[k,2], c3=cs[k,3], scale=cs[k,4])
plot(1:length(temps), p1.nd, type="l")

# c1: damage increase over time
# c2: multiplicative change in damage
# c3: linear increase in damage

#-----------
#error
c1fix=0
c2fix=1
c3fix=-0.015

errs<- function(x,temps=temps.all[temps.all$expt==1,], fecundity=fecs[fecs$expt==1,])  {
  totalerror=0
  treats=unique(temps$treatment)
  for(i in 1:length(treats)){
    delta=computeperf(series=temps[temps$treatment==treats[i],"temp"],c1=0,c2=x[1],c3=x[2],scale=x[3],printdam=FALSE)-fecundity[which(fecundity$treatment==treats[i])[1],"value"]
    totalerror=totalerror + delta^2
  }
  return( sqrt(totalerror) )
}

#check fecundity values
fecs<- PerfDat[PerfDat$metric=="fecundity",]

#-----------
#optimize
opt<- optim(p=c(1,-0.015,4e-4), fn=errs, method=c("BFGS") )
opts<- opt$par

#opt estimates
#exp 1
opts<- c(0, 14.0592519, 0.3660127,  0.2241584)
#mbh
opts<- c(0,  1.38711837e+00, -1.79323190e-02,  4.62565842e-04)

#exp 2
opts<- c(0, 0.997508878, -0.026217091,  0.001551206)

#exp 3
opts<- c(0, 1.021437448, -0.025732468,  0.002401354)

#-----------
#plot performance with values
fec1=fec(temps.expt$temp) 

p1= perf(temps.expt$temp, c1=opts[1], c2=opts[2], c3=opts[3], scale=opts[4])
p1.nd= perf.nodamage(temps.expt$temp, c1=opts[1], c2=opts[2], c3=opts[3], scale=opts[4])

d1= rbind(cbind("fec",fec1, temps.expt$time, temps.expt$treatment), 
          cbind("perf.nd",p1.nd, temps.expt$time, temps.expt$treatment), 
          cbind("perf",p1, temps.expt$time, temps.expt$treatment),
          cbind("temp",temps.expt$temp, temps.expt$time, temps.expt$treatment) )
colnames(d1)<- c("metric","value","time","treatment")
d1<- as.data.frame(d1)
d1$value <- as.numeric(d1$value)
d1$time <- as.numeric(d1$time)
d1$treatment <- as.numeric(d1$treatment)

ggplot(data=d1, aes(x=time, y =value, color=factor(treatment)))+geom_line()+facet_wrap(metric~., scale="free_y") #+xlim(0,100)

#plot outcomes
#aggregate
d1.agg= aggregate(d1, list(metric=d1$metric, treatment=d1$treatment), FUN=mean)[,-c(3,6)]
d1.agg$group<- 1
d1.agg$group[d1.agg$metric %in% c("perf","perf.nd")]<- 2

ggplot(data=d1.agg[-which(d1.agg=="temp"),], aes(x=treatment, y =value, color=metric))+geom_line()+facet_wrap(group~., scale="free_y")



