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

#toggle between desktop (y) and laptop (n)
desktop<- "n"

#FIT FUNCTION 
if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/out/")
if(desktop=="n") setwd("/Users/lbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/out/") 

temps.all<- read.csv("TempTimeSeries.csv")
PerfDat<- read.csv("PerformanceData.csv")

#------------
#Constant rate TPCs
#English grain aphid, Sitobion avenae
#Zhao et al 2013
#https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2656.12196

#https://journals.biologists.com/jeb/article/218/14/2289/14375/Daily-temperature-extremes-play-an-important-role

#development rate
#chinese clones
dr.c= function(T, Tmax=34.09, a=0.13, b=4.43, c=7.65) {
  dr=exp(a*T)-exp(b-(Tmax-T)/c)
  dr[dr<0]<- 0
  return(dr)
}

#European clones
dr= function(T, Tmax=32.91, a=0.13, b=4.28, c=7.65){ 
  d=exp(a*T)-exp(b-(Tmax-T)/c)
  d[d<0]<- 0
  return(d)
}

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

#fedundity TPC (nymphs/adult)
fec= function(T, a= -69.1, b=12.49, c= -0.34){
  fec=a +b*T +c*T^2
  fec[fec<0]<- 0
  return(fec)
}

#performance metric
pms<- c("dr", "sur", "long", "fec")
pm.ind<- 4

#find Topt and CTmax
ts=seq(0,40,0.1)

if(pm.ind==1) ft= dr(ts)
if(pm.ind==2) ft= sur(ts)
if(pm.ind==3) ft= long(ts) 
if(pm.ind==4) ft= fec(ts) 

topt<- ts[which.max(ft)]
ctmax= ts[which(ft[120:length(ft)]==0)[1]+120]
ctmin= ts[which(ft>0)[1]-1]

#Functions
#make repair depend on distance from Topt
#plot(1:40, gaussfunc(1:40, mu = Topt, sigma = 1))
damage.rep<- function(damage.p, T, c1, c2, c3, c4, tp=0, dt, Topt=topt, CTmax=ctmax)  
{ Tdamage= Topt + (CTmax-Topt)*tp
Tdif= T-Tdamage
Tdif[which(Tdif<0)]<- 0
damage.p= damage.p + dt*Tdif*(c1*damage.p + c2) 
damage.p[which(damage.p<0)]<-0
damage.p[which(damage.p>1)]<-1
#repair
damage.p= damage.p*(1-c3*gaussfunc(T, mu = Topt, sigma = c4))
return(damage.p)
}

#compute performance
perf<- function(pm, series,c1,c2,c3,c4,tp=0,scale)  {
  p=NA
  damage=0
  for(i in 1:length(series)){
    damage=damage.rep(damage,T=series[i],c1=c1,c2=c2,c3=c3,c4=c4,tp=tp,dt=1)
    
    if(pm==1) p= dr(series[i])*(1-damage)
    if(pm==2) p= sur(series[i])*(1-damage)
    if(pm==3) p= long(series[i])*(1-damage)
    if(pm==4) p= fec(series[i])*(1-damage)
    
    if(i==1) perf.all=p*scale
    if(i>1) perf.all=c(perf.all, p*scale)
  }
  return(perf.all)
}

perf.nodamage<- function(pm, series,scale)  {
  perf=NA
  for(i in 1:length(series)){
    if(pm==1) perf= dr(series[i])
    if(pm==2) perf= sur(series[i])
    if(pm==3) perf= long(series[i])
    if(pm==4) perf= fec(series[i])
    
    if(i==1) perf.all=perf*scale
    if(i>1) perf.all=c(perf.all, perf*scale)
  }
  return(perf.all)
}

#===================
#Fig 1. TPCs, model repair, damage
#Fig 2. model time series

# in Fig 1_Conceptual
#===================
#Fig 3. Temperature schematic
#Expt 1: vary min
#Expt 2: vary max
#Expt 3: vary variance (expt 3, mild means)
#Expt 4: vary means and variance (expt 3, high means)

for(expt in 1:4){
tplot= ggplot(data=temps.all[temps.all$expt %in% expt,], aes(x=time, y =temp, color=factor(treatment)))+geom_line(lwd=1.5)+ xlim(0,40)+
    theme_bw(base_size=16) +theme(legend.position = "bottom")+scale_color_viridis(discrete = TRUE)+labs(color="treatment")

if(expt==1) tplot.e1<- tplot
if(expt==2) tplot.e2<- tplot
if(expt==3) tplot.e3<- tplot
if(expt==4) tplot.e4<- tplot
}

#-------------
#Expt 5: vary duration of heatwave (expt 2)
trt<- expand.grid(hd=1:3, nd=1:3, first= 1) #first: 1 is n, 2 is h
temps<- matrix(1, nrow=nrow(trt), ncol=8)
temps= as.data.frame(cbind(trt, temps))

build.temps<- function(x){
  x[1]<- as.numeric(x[1]); x[2]<- as.numeric(x[2])
  if(x[3]==1) ts<- rep(c(rep(1, x[2]), rep(2, x[1])), ceiling(8/(x[1]+x[2])))
  if(x[3]==2) ts<- rep(c(rep(2, x[1]), rep(1, x[2])), ceiling(8/(x[1]+x[2])))
  return(ts[1:8])
}

#set up temperatures  
tz= t(apply(temps[,1:3], MARGIN=1, FUN=build.temps))
temps[, 4:ncol(temps)]<- tz

#to long format
temps.l<- melt(temps, id.vars = c("hd","nd","first"), variable.name = "day")
temps.l$treat<- paste(temps.l$hd, temps.l$nd, sep="_")
temps.l$treat.nh<- paste(temps.l$hd, temps.l$nd, temps.l$first, sep="_")
temps.l$day<-as.numeric(temps.l$day)

#plot
tplot.e5<- ggplot(temps.l, aes(x = day, y = treat.nh, fill=factor(value))) +
  geom_tile()

#-------------
#Expt 6: vary length of heatwave and timing (expt 5 adult)
#Expt 7: vary length of heatwave and timing (expt 5 nymph)
t6<- as.data.frame(matrix(1, 15, 12))
colnames(t6)<- c("AE1","AE2","AE3","AE4","AE5","AE6","NL1","NL2","NL3","NL4","NL5","NL6")
t6$day<- 1:15
#AE
t6[8:9,"AE1"]<- 2
t6[8:10,"AE2"]<- 2
t6[8:11,"AE3"]<- 2
t6[8:12,"AE4"]<- 2
t6[8:13,"AE5"]<- 2
t6[8:14,"AE6"]<- 2
#NE
t6[1:7,"NL6"]<- 2
t6[2:7,"NL5"]<- 2
t6[3:7,"NL4"]<- 2
t6[4:7,"NL3"]<- 2
t6[5:7,"NL2"]<- 2
t6[6:7,"NL1"]<- 2

#to long format
t6.l<- melt(t6, id.vars = c("day"), variable.name = "treatment")
#order treatments
t6.l$treatment= factor(t6.l$treatment, ordered=T, 
  levels=c("NL6","NL5","NL4","NL3","NL2","NL1","AE1","AE2","AE3","AE4","AE5","AE6"))

tplot.e6<- ggplot(t6.l, aes(x = day, y = treatment, fill=factor(value))) +
  geom_tile()

#write out plot
if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/figures/")
if(desktop=="n") setwd("/Users/lbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/figures/") 

pdf("Fig3_temps.pdf",height = 14, width = 14)
print( tplot.e1+tplot.e2+tplot.e3+tplot.e4+tplot.e5+tplot.e6+
         plot_layout(ncol=3, heights = c(3, 2))+ plot_annotation(tag_levels = 'A') )
dev.off()
#===================
#Fig 4. Performance with and without damage

#===================
#Figs 5. Fecundity comparisons

#===================
#Fig 6. Development comparison


