library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)
library(viridisLite)
library(viridis)
library(patchwork)
library(ggplot2)
library(deSolve)
library(pracma)				# contains sigmoid function
library(TrenchR)
library(rvmethod) #gaussian function

#FIT FUNCTION 
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
#find Topt and CTmax
ts=seq(0,40,0.1)
ft<- fec(ts)
topt<- ts[which.max(ft)]
ctmax= ts[which(ft[120:length(ft)]==0)[1]+120]
ctmin= ts[which(ft>0)[1]-1]

#FUNCTIONS
#damage
# tp: threshold for damage between Topt and CTmax; Tdamage= Topt + (CTmax-Topt)*tp
# c1: d_mult: multiplicative change in damage
# c2: d_lin: linear increase in damage
# c3: r_mag: magnitude of repair
# c4: r_breadth: breadth of repair function around Topt

#old damage functions
#damagenew= damage + dt * exp(-Ea/(R*(T+273.15))) *10^9* (c1*damage + c2) + dt*c3

#make repair depend on distance from Topt
#plot(1:40, gaussfunc(1:40, mu = Topt, sigma = 1))
damagenew<- function(damage, T, c1, c2, c3, c4, tp=0, dt, Topt=topt, CTmax=ctmax)  
{ Tdamage= Topt + (CTmax-Topt)*tp
  damagenew= damage + dt*max(T-Tdamage,0)*(c1*damage + c2) 
  damagenew= damagenew*dt*(1-c3*gaussfunc(T, mu = Topt, sigma = c4))
  if(damagenew<0) damagenew<-0
  if(damagenew>1) damagenew<-1
  return(damagenew)
}

#compute performance
perf<- function(series,c1,c2,c3,c4,tp=0,scale)  {
  perf=NA
  damage=0
  for(i in 1:length(series)){
    damage=damagenew(damage,T=series[i],c1=c1,c2=c2,c3=c3,c4=c4,tp=tp,dt=1)
    perf= fec(series[i])*(1-damage)
    if(i==1) perf.all=perf*scale
    if(i>1) perf.all=c(perf.all, perf*scale)
  }
  return(perf.all)
}

perf.nodamage<- function(series,scale)  {
  perf=NA
  for(i in 1:length(series)){
    perf= fec(series[i])
    if(i==1) perf.all=perf*scale
    if(i>1) perf.all=c(perf.all, perf*scale)
  }
  return(perf.all)
}

computeperf<- function(series,c1,c2,c3,c4,tp=0,scale,printdam=FALSE)  {
  p=0
  damage=0
  for(i in 1:length(series)){
    damage=damagenew(damage,T=series[i],c1=c1,c2=c2,c3=c3,c4=c4,tp=tp,dt=1)
  p= p + fec(series[i])*(1-damage)
  }
return(p*scale)
}

#-----------
#plot parameter values
ts= seq(1, 35, 0.5)
temps= c(ts, rev(ts),ts, rev(ts))

#make parameter combinations 
cs<- expand.grid(c1=seq(0, 2, 0.5), c2= seq(0, .01, 0.003), c3= seq(0, 1, 0.25), c4= seq(1, 5, 1),
                 scale= 1 )

#fit values
#cs<- expand.grid(c1=c(1.95,2), c2= c(0.0007, 0.001), c3= c(0.25,0.66), c4= c(1.1, 1.3), scale= 0.01)
cs<- expand.grid(c1=c(1,2), c2= c(0.00001, 0.001), c3= c(0.2,0.9), c4= c(1, 3), scale= 0.01)

for(k in 1:nrow(cs)){
  p1= perf(temps, c1=cs[k,1], c2=cs[k,2], c3=cs[k,3], c4=cs[k,4], scale=cs[k,5])
  ps= cbind(time=1:length(temps), temps, p1, k, cs[k,])
  if(k==1) ps.all<- ps
  if(k>1) ps.all<- rbind(ps.all, ps)
}

#----------------

fec= function(T, a= -69.1, b=12.49, c= -0.34){
  fec=a +b*T +c*T^2
  fec[fec<0]<- 0
  return(fec)
}

{ Tdamage= Topt + (CTmax-Topt)*tp
damagenew= damage + dt*max(T-Tdamage,0)*(c1*damage + c2) 
damagenew= damagenew*dt*(1-c3*gaussfunc(T, mu = Topt, sigma = c4))
if(damagenew<0) damagenew<-0
if(damagenew>1) damagenew<-1
return(damagenew)
}

#----------------
#plot tpc, etc
ts=seq(0,40,0.1)

#fecundity
fs= fec(ts) 
cdat<- as.data.frame(cbind(temp=ts,value=fs,type="fecundity", c1=0, c2=0, c3=0, c4=0, group=1))

#damage
#FIX damage function elsewhere for vector
damage<- function(T, c1, c2, tp=0, damage.p=0, Topt=18.4, CTmax=30.1, dt=1) {
  Tdamage= Topt + (CTmax-Topt)*tp
  Tdif= T-Tdamage
  Tdif[which(Tdif<0)]<- 0
  damagenew= damage.p + dt*Tdif*(c1*damage.p + c2) 
  return(damagenew)
  }

ddat<- as.data.frame(cbind(temp=ts, value=damage(ts, c1=0.1, c2=0.0001, tp=0, damage.p=0), type="damage", c1=0, c2=0, c3=1, c4=1, group=1.1))
ddat2<- as.data.frame(cbind(temp=ts, value=damage(ts, c1=0.9, c2=0.0001, tp=0, damage.p=0), type="damage", c1=0, c2=0, c3=1, c4=1, group=1.2))
ddat3<- as.data.frame(cbind(temp=ts, value=damage(ts, c1=0.1, c2=0.001, tp=0, damage.p=0), type="damage", c1=0, c2=0, c3=1, c4=1, group=1.3))
ddat4<- as.data.frame(cbind(temp=ts, value=damage(ts, c1=0.9, c2=0.001, tp=0, damage.p=0), type="damage", c1=0, c2=0, c3=1, c4=1, group=1.4))

#repair
repair<- function(T, c3, c4, Topt=18.4) c3*gaussfunc(T, mu = Topt, sigma = c4)
rdat<- as.data.frame(cbind(temp=ts, value=repair(ts, c3=1, c4=1), type="repair", c1=0, c2=0, c3=1, c4=1, group=2))
rdat2<- as.data.frame(cbind(temp=ts, value=repair(ts, c3=0.5, c4=1), type="repair", c1=0, c2=0, c3=0.5, c4=1, group=3))
rdat3<- as.data.frame(cbind(temp=ts, value=repair(ts, c3=1, c4=2), type="repair", c1=0, c2=0, c3=1, c4=2, group=4))

#combine
pdat<- rbind(cdat, rdat, rdat2, rdat3, ddat, ddat2, ddat3, ddat4)
pdat$temp<- as.numeric(pdat$temp); pdat$value<- as.numeric(pdat$value);
pdat$c1<- as.numeric(pdat$c1); pdat$c2<- as.numeric(pdat$c2); pdat$c3<- as.numeric(pdat$c3); pdat$c4<- as.numeric(pdat$c4)

#plot
c.fig<- ggplot(data=pdat, aes(x=temp, y =value, color=c3, lty=factor(c4), group=group))+
  geom_line()+facet_wrap(.~type, scales="free_y")+theme_bw()+
  ylab("Performance")+scale_color_viridis()

#========================
#Plots

funct.fig<- ggplot(data=ps.all, aes(x=time, y =p1, color=c3, lty=factor(c4), group=k))+
  geom_line()+facet_grid(c2~c1)+theme_bw()+
  ylab("Performance")+scale_color_viridis()

#plot temps and performance without damage
ts<- as.data.frame(cbind(time=1:length(temps), temp=temps, performance=p1.nd))
#to long format
#ts.l= gather(ts, time, temp:performance, factor_key=TRUE)

temp.fig<- ggplot(data=ts, aes(x=time, y =temp, ))+
  geom_line()+theme_bw()+
  ylab("Temperature (C)")+xlab("time")

p1.nd= perf.nodamage(temps, scale=cs[k,4])
pnd.fig<- ggplot(data=ts, aes(x=time, y =performance, ))+
  geom_line()+theme_bw()+
  ylab("Performance")+xlab("time")
  
#---
layout <- '
A#
B#
CC
'

setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/figures/")
pdf("FunFig.pdf",height = 10, width = 8)
temp.fig/
  pnd.fig +
  funct.fig+
  plot_layout(design = layout)
  #plot_layout(heights = c(1, 1,2))
dev.off()
