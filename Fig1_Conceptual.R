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

#Dev rate
dr= function(T, Tmax=32.91, a=0.13, b=4.28, c=7.65){ 
  d=exp(a*T)-exp(b-(Tmax-T)/c)
  d[d<0]<- 0
  return(d)
}

#FUNCTIONS
#damage
# tp: threshold for damage between Topt and CTmax; Tdamage= Topt + (CTmax-Topt)*tp
# c1: d_mult: multiplicative change in damage
# c2: d_linear: linear increase in damage
# c3: r_mag: magnitude of repair
# c4: r_breadth: breadth of repair function around Topt

#damage without repair
#plot(1:40, gaussfunc(1:40, mu = Topt, sigma = 1))
damage<- function(T, c1, c2, tp=0, damage.p=0, Topt=18.4, CTmax=30.1, dt=1) {
  Tdamage= Topt + (CTmax-Topt)*tp
  Tdif= T-Tdamage
  Tdif[which(Tdif<0)]<- 0
  damagenew= damage.p + dt*Tdif*(c1*damage.p + c2) 
  damagenew[which(damagenew<0)]<-0
  damagenew[which(damagenew>1)]<-1
  return(damagenew)
}

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
perf<- function(series,c1,c2,c3,c4,tp=0,scale=1)  {
  p=NA
  damage=0
  for(i in 1:length(series)){
    damage=damage.rep(damage,T=series[i],c1=c1,c2=c2,c3=c3,c4=c4,tp=tp,dt=1)
    p= fec(series[i])*(1-damage)
    if(i==1) perf.all=p*scale
    if(i>1) perf.all=c(perf.all, p*scale)
  }
  return(perf.all)
}

perf.nodamage<- function(series,scale=1)  {
  p=NA
  for(i in 1:length(series)){
    p= fec(series[i])
    if(i==1) perf.all=p*scale
    if(i>1) perf.all=c(perf.all, p*scale)
  }
  return(perf.all)
}

computeperf<- function(series,c1,c2,c3,c4,tp=0,scale=1,printdam=FALSE)  {
  p=0
  damage=0
  for(i in 1:length(series)){
    damage=damage.rep(damage,T=series[i],c1=c1,c2=c2,c3=c3,c4=c4,tp=tp,dt=1)
  p= p + fec(series[i])*(1-damage)
  }
return(p*scale/length(series))
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
cs<- expand.grid(c1=c(1,2), c2= c(0.00001, 0.001), c3= c(0.2,0.9), c4= c(1, 3), scale= 1)

for(k in 1:nrow(cs)){
  p1= perf(temps, c1=cs[k,1], c2=cs[k,2], c3=cs[k,3], c4=cs[k,4], scale=cs[k,5])
  ps= cbind(time=1:length(temps), temps, p1, k, cs[k,])
  if(k==1) ps.all<- ps
  if(k>1) ps.all<- rbind(ps.all, ps)
}

#----------------

#plot tpc, etc
ts=seq(0,40,0.1)

#performance metric
pms<- c("dr", "sur", "long", "fec")
pm.ind<- 1

#fecundity
if(pm.ind==1) fs= dr(ts) 
if(pm.ind==4) fs= fec(ts) 
cdat<- as.data.frame(cbind(temp=ts,value=fs,type="fecundity", c1=0, c2=0, c3=0, c4=0, group=1))

if(pm.ind==1){
perf_wd<- function(T, c1, c2, tp=0, damage.p, Topt=18.4, CTmax=30.1, dt=1) {
  Tdamage= Topt + (CTmax-Topt)*tp
  Tdif= T-Tdamage
  Tdif[which(Tdif<0)]<- 0
  damage.p1= damage.p + dt*Tdif*(c1*damage.p + c2) 
  damage.p1[which(damage.p1<0)]<-0
  damage.p1[which(damage.p1>1)]<-1
  p= dr(T)*(1-damage.p1)
  return(p)
}
}

if(pm.ind==4){
  perf_wd<- function(T, c1, c2, tp=0, damage.p, Topt=18.4, CTmax=30.1, dt=1) {
    Tdamage= Topt + (CTmax-Topt)*tp
    Tdif= T-Tdamage
    Tdif[which(Tdif<0)]<- 0
    damage.p1= damage.p + dt*Tdif*(c1*damage.p + c2) 
    damage.p1[which(damage.p1<0)]<-0
    damage.p1[which(damage.p1>1)]<-1
    p= fec(T)*(1-damage.p1)
    return(p)
  }
}

ddat<- as.data.frame(cbind(temp=ts, value=damage(ts, c1=0.01, c2=0.001, tp=0, damage.p=0.01), type="damage", c1=0.01, c2=0.001, c3=0, c4=0, group=1.1))
ddat2<- as.data.frame(cbind(temp=ts, value=damage(ts, c1=0.1, c2=0.001, tp=0, damage.p=0.01), type="damage", c1=0.1, c2=0.001, c3=0, c4=0, group=1.2))
ddat3<- as.data.frame(cbind(temp=ts, value=damage(ts, c1=0.01, c2=0.03, tp=0, damage.p=0.01), type="damage", c1=0.01, c2=0.03, c3=0, c4=0, group=1.3))
ddat4<- as.data.frame(cbind(temp=ts, value=damage(ts, c1=0.1, c2=0.03, tp=0, damage.p=0.01), type="damage", c1=0.1, c2=0.03,  c3=0, c4=0, group=1.4))

#performance with damage
fdat<- as.data.frame(cbind(temp=ts, value=perf_wd(ts, c1=0.01, c2=0.001, tp=0, damage.p=0.01), type="fecundity", c1=0.01, c2=0.001, c3=0, c4=0, group=3.1))
fdat2<- as.data.frame(cbind(temp=ts, value=perf_wd(ts, c1=0.1, c2=0.001, tp=0, damage.p=0.01), type="fecundity", c1=0.1, c2=0.001, c3=0, c4=0, group=3.2))
fdat3<- as.data.frame(cbind(temp=ts, value=perf_wd(ts, c1=0.01, c2=0.03, tp=0, damage.p=0.01), type="fecundity", c1=0.01, c2=0.03, c3=0, c4=0, group=3.3))
fdat4<- as.data.frame(cbind(temp=ts, value=perf_wd(ts, c1=0.1, c2=0.03, tp=0, damage.p=0.01), type="fecundity", c1=0.1, c2=0.03,  c3=0, c4=0, group=3.4))

#repair
repair<- function(T, c3, c4, Topt=18.4) c3*gaussfunc(T, mu = Topt, sigma = c4)
rdat<- as.data.frame(cbind(temp=ts, value=repair(ts, c3=0.2, c4=1), type="repair", c1=0, c2=0, c3=0.2, c4=1, group=2.1))
rdat2<- as.data.frame(cbind(temp=ts, value=repair(ts, c3=0.5, c4=1), type="repair", c1=0, c2=0, c3=0.05, c4=1, group=2.2))
rdat3<- as.data.frame(cbind(temp=ts, value=repair(ts, c3=0.2, c4=3), type="repair", c1=0, c2=0, c3=0.2, c4=3, group=2.3))
rdat4<- as.data.frame(cbind(temp=ts, value=repair(ts, c3=0.5, c4=3), type="repair", c1=0, c2=0, c3=0.05, c4=3, group=2.4))

#combine
pdat<- rbind(cdat, rdat, rdat2, rdat3, rdat4, ddat, ddat2, ddat3, ddat4,
             fdat, fdat2, fdat3, fdat4)
pdat$temp<- as.numeric(pdat$temp); pdat$value<- as.numeric(pdat$value);
pdat$c1<- as.numeric(pdat$c1); pdat$c2<- as.numeric(pdat$c2); pdat$c3<- as.numeric(pdat$c3); pdat$c4<- as.numeric(pdat$c4)

#-------
#plot

#fecundity
f.fig<- ggplot(data=pdat[which(pdat$type=="fecundity"),], aes(x=temp, y =value, color=factor(c1), lty=factor(c2), group=group))+
  geom_line(size=1.25)+theme_bw()+ theme(text=element_text(size=14))+
  ylab("Performance")+xlab("Temperature (C)") + 
  scale_colour_brewer(palette = "Dark2") +
  theme(legend.position = "bottom",  legend.box = 'vertical')+theme(legend.spacing.y = unit(3, "pt"))+
  labs(colour="d_mult", lty="d_linear")

if(pm.ind==1) {dr.fig<-f.fig; dr.fig= dr.fig+ylab("Development Rate (1/days)")}
if(pm.ind==4) {fec.fig<-f.fig; fec.fig= fec.fig+ylab("Fecundity") }

#damage
d.fig= ggplot(data=pdat[which(pdat$type=="damage"),], aes(x=temp, y =value, color=factor(c1), lty=factor(c2), group=group))+
  geom_line(size=1.25)+theme_bw()+ theme(text=element_text(size=14))+
  ylab("Damage (proportion)") +xlab("Temperature (C)")+ 
  scale_colour_brewer(palette = "Dark2") +
  theme(legend.position = "bottom",  legend.box = 'vertical')+
  labs(colour="d_mult", lty="d_linear")

#repair
r.fig= ggplot(data=pdat[which(pdat$type=="repair"),], aes(x=temp, y =value, color=factor(c3), lty=factor(c4), group=group))+
  geom_line(size=1.25)+theme_bw()+ theme(text=element_text(size=14))+
  ylab("Repair (proportion)") +xlab("Temperature (C)")+ 
  scale_colour_brewer(palette = "Dark2") +
  theme(legend.position = "bottom",  legend.box = 'vertical')+
  labs(colour="r_mag", lty="r_breadth")
  
#========================
#Plot time series

#plot temps and performance without damage
ts<- temps
#ts<- as.data.frame(cbind(time=1:length(temps), temp=temps, performance=p1.nd))
#to long format
#ts.l= gather(ts, time, temp:performance, factor_key=TRUE)

#performance without repair
fdat0<- as.data.frame(cbind(time=1:length(ts), temp=ts, value=perf(ts, c1=0.0, c2=0.0000, c3=0, c4=1, scale=1), type="fecundity", c1=0, c2=0.000, c3=0, c4=1, group=3.0))
fdat<- as.data.frame(cbind(time=1:length(ts), temp=ts, value=perf(ts, c1=0.0001, c2=0.0001, c3=0, c4=1, scale=1), type="fecundity", c1=0.0001, c2=0.0001, c3=0, c4=1, group=3.1))
fdat2<- as.data.frame(cbind(time=1:length(ts), temp=ts, value=perf(ts, c1=0.0001, c2=0.0005, c3=0, c4=1, scale=1), type="fecundity", c1=0.0001, c2=0.0005, c3=0, c4=1, group=3.2))
fdat3<- as.data.frame(cbind(time=1:length(ts), temp=ts, value=perf(ts, c1=0.001, c2=0.0001, c3=0, c4=1, scale=1), type="fecundity", c1=0.001, c2=0.0001,  c3=0, c4=1, group=3.3))
fdat4<- as.data.frame(cbind(time=1:length(ts), temp=ts, value=perf(ts, c1=0.001, c2=0.0005, c3=0, c4=1, scale=1), type="fecundity", c1=0.001, c2=0.0001,  c3=0, c4=1, group=3.4))

fdat5<- as.data.frame(cbind(time=1:length(ts), temp=ts, value=perf(ts, c1=0.005, c2=0.0005, c3=0, c4=1, scale=1), type="fecundity", c1=0.005, c2=0.0005,  c3=0, c4=1, group=3.5))
fdat6<- as.data.frame(cbind(time=1:length(ts), temp=ts, value=perf(ts, c1=0.001, c2=0.01, c3=0, c4=1, scale=1), type="fecundity", c1=0.001, c2=0.01,  c3=0, c4=1, group=3.6))

#add repair
ddat<- as.data.frame(cbind(time=1:length(ts), temp=ts, value=perf(ts, c1=0.001, c2=0.01, c3=0.05, c4=1, scale=1), type="perf repair", c1=0.001, c2=0.01, c3=0.05, c4=1, group=1.1))
ddat2<- as.data.frame(cbind(time=1:length(ts), temp=ts, value=perf(ts, c1=0.001, c2=0.01, c3=0.05, c4=3, scale=1), type="perf repair", c1=0.001, c2=0.01, c3=0.05, c4=3, group=1.2))
ddat3<- as.data.frame(cbind(time=1:length(ts), temp=ts, value=perf(ts, c1=0.001, c2=0.01, c3=0.2, c4=1, scale=1), type="perf repair", c1=0.001, c2=0.01, c3=0.2, c4=1, group=1.3))
ddat4<- as.data.frame(cbind(time=1:length(ts), temp=ts, value=perf(ts, c1=0.001, c2=0.01, c3=0.2, c4=3, scale=1), type="perf repair", c1=0.001, c2=0.01, c3=0.2, c4=3, group=1.4))

#combine
pdat<- rbind(ddat, ddat2, ddat3, ddat4, fdat0, fdat, fdat2, fdat3, fdat4, fdat5, fdat6)
pdat$time<- as.numeric(pdat$time); pdat$temp<- as.numeric(pdat$temp); pdat$value<- as.numeric(pdat$value);
pdat$c1<- as.numeric(pdat$c1); pdat$c2<- as.numeric(pdat$c2); pdat$c3<- as.numeric(pdat$c3); pdat$c4<- as.numeric(pdat$c4)

#-------------
#temp
t.fig.ts<- ggplot(data=pdat[which(pdat$type=="fecundity"),], aes(x=time, y =temp))+
  geom_line(size=1.25)+theme_bw()+ theme(text=element_text(size=14))+
  ylab("Temperature")+
  scale_colour_brewer(palette = "Dark2") +theme(legend.position = "bottom")

#performance without repair
f.fig.ts<- ggplot(data=pdat[which(pdat$type=="fecundity"),], aes(x=time, y =value, color=factor(c1), lty=factor(c2), group=group))+
  geom_line(size=1.25)+theme_bw()+ theme(text=element_text(size=14))+
  ylab("Performance")+
  scale_colour_brewer(palette = "Dark2") +theme(legend.position = "bottom")+
  labs(colour="d_mult", lty="d_linear")

#add repair
pr.fig.ts<- ggplot(data=pdat[which(pdat$type=="perf repair"),], aes(x=time, y =value, color=factor(c3), lty=factor(c4), group=group))+
  geom_line(size=1.25)+theme_bw()+ theme(text=element_text(size=14))+
  ylab("Performance")+
  scale_colour_brewer(palette = "Dark2") +theme(legend.position = "bottom")+
  labs(colour="r_mag", lty="r_breadth")

  
#---
layout <- '
AAA
BBB
BBB
CCC
CCC
'

#plot
setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/figures/")
#setwd("/Users/lbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/figures/") 

pdf("Fig1_Function.pdf",height = 9, width = 9)
(fec.fig + dr.fig) / (d.fig +r.fig) + 
 plot_annotation(tag_levels = 'A')
dev.off()

pdf("Fig2_Function.pdf",height = 9, width = 9)
  t.fig.ts +f.fig.ts +pr.fig.ts +
  plot_layout(design = layout) + plot_annotation(tag_levels = 'A')
dev.off()



