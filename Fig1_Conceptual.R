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

tp1=0.9
pm.ind=4 #run also pm.ind=1

#toggle between desktop (y) and laptop (n)
desktop<- "y"

#FIT FUNCTION 
if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/out/")
if(desktop=="n") setwd("/Users/lbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/out/") 

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
dr.e= function(T, Tmax=32.91, a=0.13, b=4.28, c=7.65){ 
  d=exp(a*T)-exp(b-(Tmax-T)/c)
  d[d<0]<- 0
  return(d)
}

#use Ma et al 2015
dr= function(T, Tmax=32.947, a=0.137, b=4.514, c=7.267){ 
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
damage<- function(T, c1, c2, tp=tp1, damage.p=0, Topt=18.4, CTmax=30.1, dt=1) {
  p=NA
  damage=0
  dur=0 
  
  #plot 3rd hour of given temp
  
  Tdamage= Topt + (CTmax-Topt)*tp
  Tdif= T-Tdamage
  if(length(which(Tdif<0))>0) Tdif[which(Tdif<0)]<- 0
  
  for(i in 1:length(T)){
    #damage
    dur<- 3 #dur + ifelse(Tdif[i]>0, 1, 0)
    #damage.n<- 1- exp(-(c1*dur)-(c2*Tdif[i]))
    damage.n<- c1*dur*ifelse(Tdif[i]>0, 1, 0)+c2*Tdif[i]
    damage= damage + damage.n
    
    if(damage<0) damage<-0
    if(damage>1) damage<-1
    
    if(i==1) damage.all=damage
    if(i>1) damage.all=c(damage.all, damage)
  }
  return(damage.all)
}

#Functions
perf.damage<- function(pm, T,c1,c2,c3,c4,tp=tp1,scale,Topt=topt, CTmax=ctmax)  
{ 
  p=NA
  damage=0
  dur=0 
  
  Tdamage= Topt + (CTmax-Topt)*tp
  Tdif= T-Tdamage
  if(length(which(Tdif<0))>0) Tdif[which(Tdif<0)]<- 0
  
  for(i in 1:length(T)){
    #damage
    dur<- dur + ifelse(Tdif[i]>0, 1, 0)
    #damage.n<- 1- exp(-(c1*dur)-(c2*Tdif[i]))
    damage.n<- c1*dur*ifelse(Tdif[i]>0, 1, 0)+c2*Tdif[i]
    damage= damage + damage.n
    
    if(damage<0) damage<-0
    if(damage>1) damage<-1
    
    #repair
    damage= damage*(1-c3*gaussfunc(T[i], mu = Topt, sigma = c4))
    
    #performance
    if(pm==1) p= dr(T[i])*(1-damage)
    if(pm==2) p= sur(T[i])*(1-damage)
    if(pm==3) p= long(T[i])*(1-damage)
    if(pm==4) p= fec(T[i])*(1-damage)
    
    if(i==1) perf.all=p*scale
    if(i>1) perf.all=c(perf.all, p*scale)
  } #end loop temperature series
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
  p1= perf.damage(pm=4,T=temps, c1=cs[k,1], c2=cs[k,2], c3=cs[k,3], c4=cs[k,4], scale=cs[k,5])
  ps= cbind(time=1:length(temps), temps, p1, k, cs[k,])
  if(k==1) ps.all<- ps
  if(k>1) ps.all<- rbind(ps.all, ps)
}

#----------------

#plot tpc, etc
ts=seq(0,40,0.1)

#performance metric
pms<- c("dr", "sur", "long", "fec")

#fecundity
if(pm.ind==1) fs= dr(ts) 
if(pm.ind==4) fs= fec(ts) 
cdat<- as.data.frame(cbind(temp=ts,value=fs,type="fecundity", c1=0, c2=0, c3=0, c4=0, group=1))

if(pm.ind==1){
  perf_wd<- function(T, c1, c2, tp=tp1, damage.p, Topt=18.4, CTmax=30.1, dt=1) {
    Tdamage= Topt + (CTmax-Topt)*tp
    Tdif= T-Tdamage
    Tdif[which(Tdif<0)]<- 0
    
    dur=3 #use dur=3
    damage.p1<- c1*dur*ifelse(Tdif>0, 1, 0)+c2*Tdif
    damage.p1= damage.p + damage.p1
    damage.p1[which(damage.p1<0)]<-0
    damage.p1[which(damage.p1>1)]<-1
    p= dr(T)*(1-damage.p1)
    return(p)
  }
}

if(pm.ind==4){
  perf_wd<- function(T, c1, c2, tp=tp1, damage.p, Topt=18.4, CTmax=30.1, dt=1) {
    Tdamage= Topt + (CTmax-Topt)*tp
    Tdif= T-Tdamage
    Tdif[which(Tdif<0)]<- 0
    
    dur=3 #use dur=3
    damage.p1<- c1*dur*ifelse(Tdif>0, 1, 0)+c2*Tdif
    damage.p1= damage.p + damage.p1
    damage.p1[which(damage.p1<0)]<-0
    damage.p1[which(damage.p1>1)]<-1
    p= fec(T)*(1-damage.p1)
    return(p)
  }
}

ddat<- as.data.frame(cbind(temp=ts, value=damage(ts, c1=0.00001, c2=0.00001, tp=tp1, damage.p=0.01), type="damage", c1=0.00001, c2=0.00001, c3=0, c4=0, group=1.1))
ddat2<- as.data.frame(cbind(temp=ts, value=damage(ts, c1=0.0001, c2=0.00001, tp=tp1, damage.p=0.01), type="damage", c1=0.0001, c2=0.00001, c3=0, c4=0, group=1.2))
ddat3<- as.data.frame(cbind(temp=ts, value=damage(ts, c1=0.00001, c2=0.0001, tp=tp1, damage.p=0.01), type="damage", c1=0.00001, c2=0.0001, c3=0, c4=0, group=1.3))
ddat4<- as.data.frame(cbind(temp=ts, value=damage(ts, c1=0.0001, c2=0.0001, tp=tp1, damage.p=0.01), type="damage", c1=0.0001, c2=0.0001,  c3=0, c4=0, group=1.4))

#performance with damage
fdat<- as.data.frame(cbind(temp=ts, value=perf_wd(ts, c1=0.00001, c2=0.00001, tp=tp1, damage.p=0), type="fecundity", c1=0.00001, c2=0.00001, c3=0, c4=0, group=3.1))
#fdat2<- as.data.frame(cbind(temp=ts, value=perf_wd(ts, c1=0.1, c2=0.00001, tp=tp1, damage.p=0), type="fecundity", c1=0.1, c2=0.00001, c3=0, c4=0, group=3.2))
#fdat3<- as.data.frame(cbind(temp=ts, value=perf_wd(ts, c1=0.00001, c2=1, tp=tp1, damage.p=0), type="fecundity", c1=0.00001, c2=1, c3=0, c4=0, group=3.3))
#fdat4<- as.data.frame(cbind(temp=ts, value=perf_wd(ts, c1=0.1, c2=1, tp=tp1, damage.p=0), type="fecundity", c1=0.1, c2=1, c3=0, c4=0, group=3.4))

#repair
repair<- function(T, c3, c4, Topt=18.4) c3*gaussfunc(T, mu = Topt, sigma = c4)
rdat<- as.data.frame(cbind(temp=ts, value=repair(ts, c3=0.2, c4=1), type="repair", c1=0, c2=0, c3=0.2, c4=1, group=2.1))
rdat2<- as.data.frame(cbind(temp=ts, value=repair(ts, c3=0.5, c4=1), type="repair", c1=0, c2=0, c3=0.5, c4=1, group=2.2))
rdat3<- as.data.frame(cbind(temp=ts, value=repair(ts, c3=0.2, c4=3), type="repair", c1=0, c2=0, c3=0.2, c4=3, group=2.3))
rdat4<- as.data.frame(cbind(temp=ts, value=repair(ts, c3=0.5, c4=3), type="repair", c1=0, c2=0, c3=0.5, c4=3, group=2.4))

#combine
pdat<- rbind(cdat, rdat, rdat2, rdat3, rdat4, ddat, ddat2, ddat3, ddat4,
             fdat)
pdat$temp<- as.numeric(pdat$temp); pdat$value<- as.numeric(pdat$value);
pdat$c1<- as.numeric(pdat$c1); pdat$c2<- as.numeric(pdat$c2); pdat$c3<- as.numeric(pdat$c3); pdat$c4<- as.numeric(pdat$c4)

#-------
#plot

#fecundity
# f.fig<- ggplot(data=pdat[which(pdat$type=="fecundity"),], aes(x=temp, y =value, color=factor(c1), lty=factor(c2), group=group))+
#   geom_line(size=1.25)+theme_bw()+ theme(text=element_text(size=14))+ 
#   ylab("Performance")+xlab("Temperature (C)") + 
#   scale_colour_brewer(palette = "Dark2") +
#   theme(legend.position = "bottom",  legend.box = 'vertical')+theme(legend.spacing.y = unit(3, "pt"))+
#   labs(colour="d_time", lty="d_temp")

f.fig<- ggplot(data=pdat[which(pdat$type=="fecundity"),], aes(x=temp, y =value, group=group))+
     geom_line(size=1.25)+theme_bw()+ theme(text=element_text(size=14))+ 
     ylab("Performance")+xlab("Temperature (°C)") 

#find Topt and CTmax
ts=seq(0,40,0.1)
ft<- fec(ts)
topt<- ts[which.max(ft)]
topt.p<- max(ft)
ctmax= ts[which(ft[120:length(ft)]==0)[1]+120]

tc= topt + 0.9*(ctmax-topt)

if(pm.ind==4){
  f.fig<- f.fig + geom_point(size=3, aes(x=topt, y =topt.p)) + geom_point(size=3, aes(x=ctmax, y =0))+
    annotate("text", size=6, x = 23, y=topt.p, label = "Topt")+
    annotate("text", size=6, x = 34, y=2, label = "Tmax")+
    geom_segment(aes(x = topt, y = 0, xend = ctmax, yend = 0), linetype="dashed", size=1.25)+
    annotate("text", size=6, x = 28, y=2, label = "c")+
    geom_point(size=3, aes(x=tc, y =fec(tc) ))+
    annotate("text", size=6, x = 32, y=fec(tc), label = "Tc")+
    geom_segment(aes(x = tc, y = -1, xend = tc, yend = 1), linetype="dashed", size=1.25)
}

if(pm.ind==1) {dr.fig<-f.fig; dr.fig= dr.fig+ylab("Development Rate (1/days)")}
if(pm.ind==4) {fec.fig<-f.fig; fec.fig= fec.fig+ylab("Fecundity (nymphs per adult)") }

#damage
d.fig= ggplot(data=pdat[which(pdat$type=="damage"),], aes(x=temp, y =value, color=factor(c1), lty=factor(c2), group=group))+
  geom_line(size=1.25)+theme_bw()+ theme(text=element_text(size=14))+
  ylab("Damage (proportion)") +xlab("Temperature (°C)")+ 
  scale_colour_brewer(palette = "Dark2") +
  theme(legend.position = "bottom",  legend.box = 'vertical')+
  labs(colour= expression(d[time]), lty=expression(d[temp]) )

#repair
r.fig= ggplot(data=pdat[which(pdat$type=="repair"),], aes(x=temp, y =value, color=factor(c3), lty=factor(c4), group=group))+
  geom_line(size=1.25)+theme_bw()+ theme(text=element_text(size=14))+
  ylab("Repair (proportion)") +xlab("Temperature (°C)")+ 
  scale_colour_brewer(palette = "Dark2") +
  theme(legend.position = "bottom",  legend.box = 'vertical')+
  labs(colour=expression(r[mag]), lty=expression(r[breadth]))
  
#========================
#Plot time series

#plot temps and performance without damage
#ts<- temps.all[which(temps.all$expt==6 & temps.all$treatment=="AE6"),"temp"]
ts<- temps.all[which(temps.all$expt==1 & temps.all$treatment=="13"),"temp"]
ts<- ts[1:100]

#ts<- as.data.frame(cbind(time=1:length(temps), temp=temps, performance=p1.nd))
#to long format
#ts.l= gather(ts, time, temp:performance, factor_key=TRUE)

#performance without repair
fdat0<- as.data.frame(cbind(time=1:length(ts), temp=ts, value=perf.damage(pm=4, T=ts, c1=0.0, c2=0.0000, c3=0, c4=1, scale=1), type="fecundity", c1=0, c2=0.000, c3=0, c4=1, group=3.0))
fdat<- as.data.frame(cbind(time=1:length(ts), temp=ts, value=perf.damage(pm=4, T=ts, c1=0.0001, c2=0.0001, c3=0, c4=1, scale=1), type="fecundity", c1=0.0001, c2=0.0001, c3=0, c4=1, group=3.1))
fdat2<- as.data.frame(cbind(time=1:length(ts), temp=ts, value=perf.damage(pm=4, T=ts, c1=0.0001, c2=0.001, c3=0, c4=1, scale=1), type="fecundity", c1=0.0001, c2=0.001, c3=0, c4=1, group=3.2))
fdat3<- as.data.frame(cbind(time=1:length(ts), temp=ts, value=perf.damage(pm=4, T=ts, c1=0.001, c2=0.0001, c3=0, c4=1, scale=1), type="fecundity", c1=0.001, c2=0.0001,  c3=0, c4=1, group=3.3))
fdat4<- as.data.frame(cbind(time=1:length(ts), temp=ts, value=perf.damage(pm=4, T=ts, c1=0.001, c2=0.001, c3=0, c4=1, scale=1), type="fecundity", c1=0.001, c2=0.001,  c3=0, c4=1, group=3.4))

#add repair
ddat<- as.data.frame(cbind(time=1:length(ts), temp=ts, value=perf.damage(pm=4, T=ts, c1=0.001, c2=0.01, c3=0.05, c4=1, scale=1), type="perf repair", c1=0.001, c2=0.01, c3=0.05, c4=1, group=1.1))
ddat2<- as.data.frame(cbind(time=1:length(ts), temp=ts, value=perf.damage(pm=4, T=ts, c1=0.001, c2=0.01, c3=0.05, c4=3, scale=1), type="perf repair", c1=0.001, c2=0.01, c3=0.05, c4=3, group=1.2))
ddat3<- as.data.frame(cbind(time=1:length(ts), temp=ts, value=perf.damage(pm=4, T=ts, c1=0.001, c2=0.01, c3=0.2, c4=1, scale=1), type="perf repair", c1=0.001, c2=0.01, c3=0.2, c4=1, group=1.3))
ddat4<- as.data.frame(cbind(time=1:length(ts), temp=ts, value=perf.damage(pm=4, T=ts, c1=0.001, c2=0.01, c3=0.2, c4=3, scale=1), type="perf repair", c1=0.001, c2=0.01, c3=0.2, c4=3, group=1.4))

#combine
pdat<- rbind(ddat, ddat2, ddat3, ddat4, fdat0, fdat, fdat2, fdat3, fdat4)
pdat$time<- as.numeric(pdat$time); pdat$temp<- as.numeric(pdat$temp); pdat$value<- as.numeric(pdat$value);
pdat$c1<- as.numeric(pdat$c1); pdat$c2<- as.numeric(pdat$c2); pdat$c3<- as.numeric(pdat$c3); pdat$c4<- as.numeric(pdat$c4)

#-------------
#temp
t.fig.ts<- ggplot(data=pdat[which(pdat$type=="fecundity"),], aes(x=time, y =temp))+
  annotate("rect", xmin = 0, xmax = 6, ymin = -Inf, ymax = Inf, alpha = .3)+
  annotate("rect", xmin = 22, xmax = 30, ymin = -Inf, ymax = Inf, alpha = .3)+
  annotate("rect", xmin = 46, xmax = 54, ymin = -Inf, ymax = Inf, alpha = .3)+
  annotate("rect", xmin = 70, xmax = 78, ymin = -Inf, ymax = Inf, alpha = .3)+
  xlim(12,84)+
  geom_line(size=1.25)+theme_bw()+ theme(text=element_text(size=14))+
  ylab("Temperature (°C)")+xlab("Time (hour)")+
  scale_colour_brewer(palette = "Dark2") +theme(legend.position = "bottom")

#performance without repair
f.fig.ts<- ggplot(data=pdat[which(pdat$type=="fecundity"),], aes(x=time, y =value, color=factor(c1), lty=factor(c2), group=group))+
  annotate("rect", xmin = 0, xmax = 6, ymin = -Inf, ymax = Inf, alpha = .3)+
  annotate("rect", xmin = 22, xmax = 30, ymin = -Inf, ymax = Inf, alpha = .3)+
  annotate("rect", xmin = 46, xmax = 54, ymin = -Inf, ymax = Inf, alpha = .3)+
  annotate("rect", xmin = 70, xmax = 78, ymin = -Inf, ymax = Inf, alpha = .3)+
  xlim(12,84)+
  geom_line(size=1.25)+theme_bw()+ theme(text=element_text(size=14))+
  ylab("Performance")+xlab("Time (hour)")+
  scale_colour_brewer(palette = "Dark2") +theme(legend.position = "bottom")+
  labs(colour=expression(d[time]), lty=expression(d[temp]))

#add repair
pr.fig.ts<- ggplot(data=pdat[which(pdat$type=="perf repair"),], aes(x=time, y =value, color=factor(c3), lty=factor(c4), group=group))+
  annotate("rect", xmin = 0, xmax = 6, ymin = -Inf, ymax = Inf, alpha = .3)+
  annotate("rect", xmin = 22, xmax = 30, ymin = -Inf, ymax = Inf, alpha = .3)+
  annotate("rect", xmin = 46, xmax = 54, ymin = -Inf, ymax = Inf, alpha = .3)+
  annotate("rect", xmin = 70, xmax = 78, ymin = -Inf, ymax = Inf, alpha = .3)+
  xlim(12,84)+
  geom_line(size=1.25)+theme_bw()+ theme(text=element_text(size=14))+
  ylab("Performance")+xlab("Time (hour)")+
  scale_colour_brewer(palette = "Dark2") +theme(legend.position = "bottom")+
  labs(colour=expression(r[mag]), lty=expression(r[breadth]) )

  
#---
layout <- '
AAA
BBB
BBB
CCC
CCC
'

#plot
if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/figures/")
if(desktop=="n") setwd("/Users/lbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/figures/") 

pdf("Fig1_Function.pdf",height = 9, width = 9)
(fec.fig + dr.fig) / (d.fig +r.fig) + 
 plot_annotation(tag_levels = 'A')
dev.off()

pdf("Fig3_Function.pdf",height = 9, width = 9)
  t.fig.ts +f.fig.ts +pr.fig.ts +
  plot_layout(design = layout) + plot_annotation(tag_levels = 'A')
dev.off()



