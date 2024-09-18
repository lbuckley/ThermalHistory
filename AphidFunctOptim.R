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
ts=seq(0,40,0.5)
ft<- fec(ts)
topt<- ts[which.max(ft)]
ctmax= ts[which(ft[20:length(ft)]==0)[1]+20]

#FUNCTIONS
#damage
# c1: damage increase over time
# c2: multiplicative change in damage
# c3: linear increase in damage

# damagenew<- function(damage, T, c1, c2, c3, dt, Ea=0.65, R=8.62*10^{-5})  
#   {damagenew= damage + dt * exp(-Ea/(R*(T+273.15))) *10^9* (c1*damage + c2) + dt*c3
#   if(damagenew<0) damagenew<-0
#   if(damagenew>1) damagenew<-1
#   return(damagenew)
# }

#make damage depend on distance from Topt
damagenew<- function(damage, T, c1, c2, c3, tp, dt, Topt=topt, CTmax= ctmax)  
{ #tdamage= Topt #+ (CTmax-Topt)*tp, fitting suggest threshold close to Topt 
  damagenew= damage + dt*max(T-Topt,0)*(c1*damage + c2) + dt*c3
if(damagenew<0) damagenew<-0
if(damagenew>1) damagenew<-1
return(damagenew)
}

#make repair depend on distance from Topt
#plot(1:40, gaussfunc(1:40, mu = Topt, sigma = 1))
damagenew<- function(damage, T, c1, c2, c3, c4, dt, Topt=topt)  
{ damagenew= damage + dt*max(T-Topt,0)*(c1*damage + c2)
  damagenew= damagenew*dt*(1-c3*gaussfunc(T, mu = Topt, sigma = c4))
  if(damagenew<0) damagenew<-0
  if(damagenew>1) damagenew<-1
  return(damagenew)
}

#compute performance
perf<- function(series,c1,c2,c3,c4,scale)  {
  perf=NA
  damage=0
  for(i in 1:length(series)){
    damage=damagenew(damage,T=series[i],c1=c1,c2=c2,c3=c3,c4=c4,dt=1)
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

computeperf<- function(series,c1,c2,c3,c4,scale,printdam=FALSE)  {
  p=0
  damage=0
  for(i in 1:length(series)){
    damage=damagenew(damage,T=series[i],c1=c1,c2=c2,c3=c3,c4=c4,dt=1)
  p= p + fec(series[i])*(1-damage)
  }
return(p*scale)
}

#-----------
#plot parameter values
ts= seq(1, 35, 0.5)
temps= c(ts, rev(ts),ts, rev(ts))

#make parameter combinations 
#cs<- expand.grid(c1= seq(-1, 1, 0.5), c2= seq(.6, 1.6, 0.2), c3= seq(-0.005, 0.005, 0.005),
#                 scale= seq(0.6, 1, 0.2) )
#cs<- expand.grid(c1=seq(0, 1, 0.5), c2= seq(0.001, .01, 0.005), c3= seq(-0.04, -0.005, 0.003),
#                 scale= 1 )

cs<- expand.grid(c1=seq(0, 2, 0.5), c2= seq(0, .01, 0.003), c3= seq(0, 1, 0.25), c4= seq(1, 5, 1),
                 scale= 1 )

for(k in 1:nrow(cs)){
  p1= perf(temps, c1=cs[k,1], c2=cs[k,2], c3=cs[k,3], c4=cs[k,4], scale=cs[k,5])
  ps= cbind(time=1:length(temps), temps, p1, k, cs[k,])
  if(k==1) ps.all<- ps
  if(k>1) ps.all<- rbind(ps.all, ps)
}

ggplot(data=ps.all, aes(x=time, y =p1, color=c3, lty=factor(c4), group=k))+geom_line()+facet_grid(c2~c1) 

plot(1:length(temps), temps, type="l")
p1.nd= perf.nodamage(temps, scale=cs[k,4])
plot(1:length(temps), p1.nd, type="l")

# c1: damage increase over time
# c2: multiplicative change in damage
# c3: linear increase in damage

#-----------
#error
#extract fecundity values
fecs<- PerfDat[PerfDat$metric=="fecundity",]

expt<-2
#estimate scale
scale.est<- max(fecs[fecs$expt==expt,"value"])/(sum(fec(temps.all[temps.all$expt==expt,"temp"]))/length(unique(fecs[fecs$expt==expt,"treatment"])))

errs<- function(x,temps=temps.all[temps.all$expt==expt,], fecundity=fecs[fecs$expt==expt,]){  
  totalerror=0
  treats=unique(temps$treatment)
  for(i in 1:length(treats)){
    delta=computeperf(series=temps[temps$treatment==treats[i],"temp"],c1=x[1],c2=x[2],c3=x[3],c4=x[4],scale=x[5])-fecundity[which(fecundity$treatment==treats[i])[1],"value"]
    totalerror=totalerror + delta^2
  }
  return( sqrt(totalerror) )
}

#-----------
#optimize

# #Manually optimize
# params<- as.matrix(expand.grid(c1=seq(0, 1, 0.1), c2= seq(0.001, .01, 0.001), c3= seq(-0.04, -0.005, 0.005) ))
# perror=rep(NA, nrow(params))
# 
# temps=temps.all[temps.all$expt==expt,]
# fecundity=fecs[fecs$expt==expt,]
# treats=unique(temps$treatment)
# 
# for(k in 1:nrow(params)){
# totalerror= 0
# for(i in 1:length(treats)){
#   delta=computeperf(series=temps[temps$treatment==treats[i],"temp"],c1=params[k,1],c2=params[k,2],c3=params[k,3],scale=scale.est)-fecundity[which(fecundity$treatment==treats[i])[1],"value"]
#   totalerror=totalerror + delta^2
# }
# perror[k]<- sqrt(totalerror)
# }
# 
# opts= c(params[which.min(perror),], scale.est)

#---------
#opt<- optim(par=c(0,1,-0.01), fn=errs, method=c("BFGS") )
#opt<- optim(par=c(0.02,0.005,-0.002,0.5), fn=errs, NULL, method=c("L-BFGS-B"), 
#            lower=c(0,0,-0.005,0), upper=c(0.05,0.01,-0.0,1) )

c(scale.est, scale.est/5, scale.est*5)

opt<- optim(par=c(2,0.001,0.1,1,scale.est), fn=errs, NULL, method=c("L-BFGS-B"), 
            lower=c(0,0,0,1,scale.est/5), upper=c(4,0.01,1,2,scale.est*5) )

opts<- opt$par
#opts<- c(opt$par, scale.est)

# #opt estimates
# #exp 1
# 0.0000000000  0.0010000000 -0.0400000000  0.0005293701 
# #mbh
# opts<- c(0,  1.38711837e+00, -1.79323190e-02,  4.62565842e-04)
# 80% 0.000000000  0.001000000 -0.007615917
#ctmax  0.0002834419  0.0015512075 -0.0400000000  0.0005293701
# est scale 1.00000000  0.01000000 -0.00500000  0.01006947
# 
# #exp 2
# 0.000000000  0.006677989 -0.040000000  0.003259022
#80% 0.055123653  0.010000000 -0.020471564  0.003259022
#ctmax 1.000000000  0.010000000 -0.017595137  0.003259022
#est scale  1.00000000  0.01000000 -0.00500000  0.00617719
# 
# #exp 3
# 0.002443842  0.004665393 -0.010321043  0.005058577
#80% 0.000000000  0.001000000 -0.005000000  0.005058577
#ctmax 0.000000000  0.001000000 -0.005000000  0.005058577
#est scale 0.4825552  0.0010000 -0.0050000  0.0023397
#-----------
#plot performance with values
temps.expt<- temps.all[temps.all$expt==expt,]

p1= perf(temps.expt$temp, c1=opts[1], c2=opts[2], c3=opts[3], c4=opts[4], scale=opts[5])
p1.nd= perf.nodamage(temps.expt$temp, scale=opts[5])

d1<- data.frame(metric="perf.nd",value=p1.nd, time=temps.expt$time, treatment=temps.expt$treatment) 
d2<- data.frame(metric="perf",value=p1, time=temps.expt$time, treatment=temps.expt$treatment)
d3<- data.frame(metric="temp",value=temps.expt$temp, time=temps.expt$time, treatment=temps.expt$treatment)
d1<- rbind(d1,d2,d3)        
          
d1$value <- as.numeric(d1$value)
d1$time <- as.numeric(d1$time)
d1$metric <- factor(d1$metric)

if(expt==1) d1$treatment <- factor(d1$treatment)
if(expt==2) d1$treatment <- factor(d1$treatment, ordered=TRUE, levels=c(
  "1_3_1", "1_2_1", "1_1_1", "2_3_1", "2_2_1", "2_1_1", "3_3_1", "3_2_1", "3_1_1",
  "1_3_2", "1_2_2", "1_1_2", "2_3_2", "2_2_2", "2_1_2", "3_3_2", "3_2_2", "3_1_2"))
if(expt==3) d1$treatment <- factor(d1$treatment, ordered=TRUE, levels=c(
  "22_0", "22_5", "22_9",  "22_13", "28_0", "30_0", "32_0", "28_5", "30_5", "32_5"))

#ggplot(data=d1, aes(x=time, y =value, color=factor(treatment)))+geom_line()+facet_wrap(metric~., scale="free_y") #+xlim(0,100)

#plot outcomes
#aggregate
d1.agg= aggregate(.~metric+treatment, d1, sum)
d1.agg= d1.agg[-which(d1.agg$metric=="temp"), 1:3]

#add observed
fdat<- fecs[fecs$expt==expt,c("metric","treatment", "value" )]
fdat<- aggregate(.~metric+treatment, fdat, mean)

d1.agg<- rbind(d1.agg, fdat)

ggplot(data=d1.agg, aes(x=treatment, y =value, color=metric, group=metric))+geom_point()+geom_line()

#opts
#expt1  
#expt2 1.9586133275 0.0004568355 0.2706604496 1.1895890090 0.0137586265
#expt3 

