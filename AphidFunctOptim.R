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

#make damage depend on distance from Topt
damagenew<- function(damage, T, c1, c2, c3, dt, Topt=18)  
{damagenew= damage + dt*max(T-Topt,0)*(c1*damage + c2) + dt*max(Topt-T,0)*c3
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
    perf= max(0, fec(series[i])*(1-damage))
    if(i==1) perf.all=perf*scale
    if(i>1) perf.all=c(perf.all, perf*scale)
  }
  return(perf.all)
}

perf.nodamage<- function(series,c1,c2,c3,scale)  {
  perf=NA
  damage=0
  for(i in 1:length(series)){
    perf= max(0, fec(series[i]))
    if(i==1) perf.all=perf*scale
    if(i>1) perf.all=c(perf.all, perf*scale)
  }
  return(perf.all)
}

computeperf<- function(series,c1,c2,c3,scale,printdam=FALSE)  {
  perf=0
  damage=0
  for(i in 1:length(series)){
    damage=damagenew(damage,T=series[i],c1=c1,c2=c2,c3=c3,dt=1)
  perf= max(0, (perf + fec(series[i])*(1-damage))*scale)
  }
return(perf)
}

#-----------
#plot parameter values
ts= seq(1, 30, 0.5)
temps= c(ts, rev(ts),ts, rev(ts))

#make parameter combinations
#cs<- expand.grid(c1= seq(-1, 1, 0.5), c2= seq(.6, 1.6, 0.2), c3= seq(-0.005, 0.005, 0.005),
#                 scale= seq(0.6, 1, 0.2) )
cs<- expand.grid(c1= seq(0, 3, 0.5), c2= seq(0.1, 1, 0.3), c3= seq(-0.005, 0, 0.002),
                 scale= seq(0, 1, 0.2) )

for(k in 1:nrow(cs)){
  p1= perf(temps, c1=cs[k,1], c2=cs[k,2], c3=cs[k,3], scale=cs[k,4])
  ps= cbind(time=1:length(temps), temps, p1, k, cs[k,])
  if(k==1) ps.all<- ps
  if(k>1) ps.all<- rbind(ps.all, ps)
}

ggplot(data=ps.all, aes(x=time, y =p1, color=c1, lty=factor(scale), group=k))+geom_line()+facet_grid(c2~c3) 

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

expt<-3

errs<- function(x,temps=temps.all[temps.all$expt==expt,], fecundity=fecs[fecs$expt==expt,])  {
  totalerror=0
  treats=unique(temps$treatment)
  for(i in 1:length(treats)){
    delta=computeperf(series=temps[temps$treatment==treats[i],"temp"],c1=x[1],c2=x[2],c3=x[3],scale=x[4],printdam=FALSE)-fecundity[which(fecundity$treatment==treats[i])[1],"value"]
    totalerror=totalerror + delta^2
  }
  return( sqrt(totalerror) )
}
#scale=x[3]

#check fecundity values
fecs<- PerfDat[PerfDat$metric=="fecundity",]

#-----------
#optimize

#fec.init= mean(fecs[fecs$expt==expt,"value"])/nrow(temps.all[temps.all$expt==expt,])
#opt<- optim(p=c(1,-0.014,0.001), fn=errs, method=c("BFGS") )
#opt<- optim(p=c(1,-0.014,0.001, fec.init), fn=errs, method=c("BFGS") )
#opt<- optim(par=c(0,1,-0.015,0.5), fn=errs, method=c("BFGS") )

opt<- optim(par=c(0,1,-0.01,0.5), fn=errs, NULL, method=c("L-BFGS-B"), 
            lower=c(0,0,-0.015,0), upper=c(3,1,0,1) )

opts<- opt$par
#opts<- c(0, opt$par,0)

# #opt estimates
# #exp 1
# opts<- c(0, 0.9987293, -0.1568119,  0.1257637)
# #mbh
# opts<- c(0,  1.38711837e+00, -1.79323190e-02,  4.62565842e-04)
# 
# #exp 2
# opts<- c(0, 1.74156168 -0.02247459  0.45884653)
# 
# #exp 3
# opts<- c(0, 1.021437448, -0.025732468,  0.002401354)

#-----------
#plot performance with values
temps.expt<- temps.all[temps.all$expt==expt,]

fec1=fec(temps.expt$temp)#/nrow(temps.expt)

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

if(expt==1) d1$treatment <- as.numeric(d1$treatment)
if(expt==2) d1$treatment <- factor(d1$treatment, ordered=TRUE, levels=c(
  "1_3_1", "1_2_1", "1_1_1", "2_3_1", "2_2_1", "2_1_1", "3_3_1", "3_2_1", "3_1_1",
  "1_3_2", "1_2_2", "1_1_2", "2_3_2", "2_2_2", "2_1_2", "3_3_2", "3_2_2", "3_1_2"))
if(expt==3) d1$treatment <- factor(d1$treatment, ordered=TRUE, levels=c(
  "22_0", "22_5", "22_9",  "22_13", "28_0", "30_0", "32_0", "28_5", "30_5", "32_5"))

ggplot(data=d1, aes(x=time, y =value, color=factor(treatment)))+geom_line()+facet_wrap(metric~., scale="free_y") #+xlim(0,100)

#plot outcomes
#aggregate
d1.agg= aggregate(d1, list(metric=d1$metric, treatment=d1$treatment), FUN=mean)[,-c(3,6)]

ggplot(data=d1.agg[-which(d1.agg=="temp"),], aes(x=treatment, y =value, color=metric))+geom_point()+geom_line()



