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
library(ggpubr)

#toggle between desktop (y) and laptop (n)
desktop<- "n"

#FIT FUNCTION 
if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/out/")
if(desktop=="n") setwd("/Users/lbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/out/") 

temps.all<- read.csv("TempTimeSeries.csv")
PerfDat<- read.csv("PerformanceData.csv")
out<- read.csv("out_fec.csv")
out.dr<- read.csv("out_dr.csv")

#performance metric
pms<- c("dr", "sur", "long", "fec")
pm.ind<- 4

#scen: #1. baseline fit scale; 2. fix scale; 3. fit tp; 4. drop c1; 5. drop c2 with floor
#scens= c(1,5,5,3,5,5,5)   #tp=1: scens= c(1,3,3,3,3,2,2) 
scens= c(1,1,1,1,1,1,1)

#set up default tp
tp1=0.9

#rename treatment in expt 3
t3<- c("22_0","22_5","22_9","22_13")
t3.lab<- c("00","05","09","13")
temps.all$treatment[temps.all$expt==3]<- t3.lab[match(temps.all$treatment[temps.all$expt==3],t3)]
PerfDat$treatment[PerfDat$expt==3]<- t3.lab[match(PerfDat$treatment[PerfDat$expt==3],t3)]

#drop field treatments
PerfDat <- PerfDat[-which(PerfDat$population=="field"),]
#====================
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
perf.damage<- function(pm, T,c1,c2,c3,c4,tp=tp, tr=tr, scale,Topt=topt, CTmax=ctmax)  
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
    damage.n<- c1*dur*ifelse(Tdif[i]>0, 1, 0)+c2*Tdif[i]
    #damage.n<- 1- exp(-(c1*dur)-(c2*Tdif[i]))
    #damage.n<- c1*exp(dur)*ifelse(Tdif[i]>0, 1, 0)+c2*exp(Tdif[i])
    
    damage= damage + damage.n
    
    if(damage<0) damage<-0
    if(damage>1) damage<-1
    
    #repair
    damage= damage*(1-c3*gaussfunc(T[i], mu = tr, sigma = c4))
    
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
 elab<- paste("expt",expt,sep=" ")
 temps.all$treatment[temps.all$expt==3]<- gsub("22_","", temps.all$treatment[temps.all$expt==3])
 
tplot= ggplot(data=temps.all[temps.all$expt %in% expt,], aes(x=time, y =temp, color=factor(treatment)))+geom_line(lwd=1.5)+ xlim(0,40)+
    theme_bw(base_size=16) +theme(legend.position = "right")+scale_color_viridis(discrete = TRUE)+labs(color="treatment")+
  labs(title=elab, color="treatment (°C)")+xlab("time (hour)")+ylab("temperature")

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
  if(x[3]==1) ts<- rep(c(rep("normal", x[2]), rep("hot", x[1])), ceiling(8/(x[1]+x[2])))
  if(x[3]==2) ts<- rep(c(rep("hot", x[1]), rep("normal", x[2])), ceiling(8/(x[1]+x[2])))
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
temps.l$elab<- paste("expt",temps.l$expt,sep=" ")

#plot
tplot.e5<- ggplot(temps.l, aes(x = day, y = treat.nh, fill=factor(value))) +
  geom_tile()+labs(title="expt 5", fill="day type")+
  theme_bw(base_size=16)+theme(legend.position = "right")+ylab("treatment")

#-------------
#Expt 6: vary length of heatwave and timing (expt 5 adult)
#Expt 7: vary length of heatwave and timing (expt 5 nymph)
t6<- as.data.frame(matrix("normal", 20, 12))
colnames(t6)<- c("AE1","AE2","AE3","AE4","AE5","AE6","NL1","NL2","NL3","NL4","NL5","NL6")
t6$day<- 1:20
#AE
t6[8:9,"AE1"]<- "hot"
t6[8:10,"AE2"]<- "hot"
t6[8:11,"AE3"]<- "hot"
t6[8:12,"AE4"]<- "hot"
t6[8:13,"AE5"]<- "hot"
t6[8:14,"AE6"]<- "hot"
#NE
t6[1:7,"NL6"]<- "hot"
t6[2:7,"NL5"]<- "hot"
t6[3:7,"NL4"]<- "hot"
t6[4:7,"NL3"]<- "hot"
t6[5:7,"NL2"]<- "hot"
t6[6:7,"NL1"]<- "hot"

#to long format
t6.l<- melt(t6, id.vars = c("day"), variable.name = "treatment")
#order treatments
t6.l$treatment= factor(t6.l$treatment, ordered=TRUE, 
  levels=c("NL6","NL5","NL4","NL3","NL2","NL1","AE1","AE2","AE3","AE4","AE5","AE6"))

tplot.e6<- ggplot(t6.l, aes(x = day, y = treatment, fill=factor(value))) +
  geom_tile()+labs(title="expts 6 & 7", fill="day type")+
  theme_bw(base_size=16)+theme(legend.position = "right")+xlim(0.5,16.5)

#write out plot
if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/figures/")
if(desktop=="n") setwd("/Users/lbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/figures/") 

pdf("Fig3_temps.pdf",height = 10, width = 12)
print( tplot.e1+tplot.e2+tplot.e3+tplot.e4+tplot.e5+tplot.e6+
         plot_layout(ncol=2) ) #+ plot_annotation(tag_levels = 'A')
dev.off()
#===================
#Fig 4. Performance with and without damage

for(expt in 1:7){
temps.expt<- temps.all[temps.all$expt==expt,]
#drop first=2 for experiment 5
if(expt==5){
  #put other first in supplement
  treats= matrix(unlist(strsplit(temps.expt$treatment, split = "_")),ncol=3,byrow=T)
  temps.expt= temps.expt[treats[,3]==1,]
}

cs<- as.numeric(out[which(out$expt==expt & out$scenario==scens[expt]),4:10])

p1= perf.damage(pm=pm.ind, temps.expt$temp, c1=cs[1], c2=cs[2], c3=cs[3], c4=cs[4], tp=cs[5], tr=cs[6], scale=cs[7])
p1.nd= perf.nodamage(pm=pm.ind, temps.expt$temp, scale=cs[7])

d1<- data.frame(metric="perf.nd",value=p1.nd, time=temps.expt$time, treatment=temps.expt$treatment) 
d2<- data.frame(metric="perf",value=p1, time=temps.expt$time, treatment=temps.expt$treatment)
d1<- rbind(d1,d2)        

#combine across experiments
d1$expt<- expt
if(expt==1) d1.all<- d1
if(expt>1) d1.all<- rbind(d1.all, d1)

#plot
d1$metric <- revalue(d1$metric, c("perf.nd" = "performance no damage", "perf" = "performance with damage"))
d1$metric <- factor(d1$metric, ordered=TRUE, levels=c("performance no damage", "performance with damage"))
elab= paste("expt",expt,sep=" ")

#plot days 7 to 11
d1$treatment <- factor(d1$treatment)
pplot= ggplot(data=d1, aes(x=time, y =value, color=factor(treatment), lty=metric))+geom_line(lwd=1.5)+xlim(216,312)+
  theme_bw(base_size=16) +theme(legend.position = "right")+scale_color_viridis(discrete = TRUE)+
  labs(color="treatment", title=elab)+guides(lty ="none")+ylab("fecundity")+xlab("time (hour)")

#for adult heatwaves, plot later period
#account for 20 minute data
if(expt==6){
  pplot= ggplot(data=d1, aes(x=time, y =value, color=factor(treatment), lty=metric))+geom_line(lwd=1.5)+xlim(216,312)+
    theme_bw(base_size=16) +theme(legend.position = "right")+scale_color_viridis(discrete = TRUE)+
    labs(color="treatment", title=elab)+guides(lty ="none")+ylab("fecundity")+xlab("time (hour)")
}

if(expt==7){
  pplot= ggplot(data=d1, aes(x=time, y =value, color=factor(treatment), lty=metric))+geom_line(lwd=1.5)+xlim(24,120)+
    theme_bw(base_size=16) +theme(legend.position = "right")+scale_color_viridis(discrete = TRUE)+
    labs(color="treatment", title=elab)+guides(lty ="none")+ylab("fecundity")+xlab("time (hour)")
}

if(expt==1) pplot.e1<- pplot
if(expt==2) pplot.e2<- pplot
if(expt==3) pplot.e3<- pplot
if(expt==4) pplot.e4<- pplot
if(expt==5) pplot.e5<- pplot
if(expt==6) pplot.e6<- pplot
if(expt==7) pplot.e7<- pplot
}

#write out plot
pdf("Fig4_performance.pdf",height = 12, width = 12)
print( pplot.e1+pplot.e2+pplot.e3+pplot.e4+pplot.e5+pplot.e6+pplot.e7+
         plot_layout(ncol=2) ) #+ plot_annotation(tag_levels = 'A')
dev.off()

#===================
#Figs 5. Fecundity comparisons

#extract performance values
if(pm.ind==1) fecs<- PerfDat[PerfDat$metric=="dev_rate",]
if(pm.ind==2) fecs<- PerfDat[PerfDat$metric=="survival",]
if(pm.ind==3) fecs<- PerfDat[PerfDat$metric=="longevity",]
if(pm.ind==4) fecs<- PerfDat[PerfDat$metric=="fecundity",]

#aggregate performance estimates
d1.agg= aggregate(.~metric+treatment+expt, d1.all, mean)

#add observed
fdat<- fecs[,c("metric","treatment", "value","expt")]
fdat<- aggregate(.~metric+treatment+expt, fdat, mean)
d1.agg<- d1.agg[,c("metric","treatment", "value","expt")]
d1.agg<- rbind(d1.agg, fdat)
d1.agg$elab<- paste("expt",d1.agg$expt,sep=" ")

#rename metrics
metrs<- c("fecundity","perf","perf.nd")
metr.lab<- c("observed","performance","performance without damage")
d1.agg$metric= metr.lab[match(d1.agg$metric, metrs)]

for(expt in 1:7){
  d1.agg.e<- d1.agg[which(d1.agg$expt==expt),]
  
#adjust treatments
if(expt==3) {d1.agg.e$treatment<- gsub("22_", "", d1.agg.e$treatment); d1.agg.e$treatment= as.numeric(d1.agg.e$treatment)}
if(expt==4){
  d1.agg.e$tvar<- as.numeric(substr(d1.agg.e$treatment, 4, 5))
  d1.agg.e$treatment<- as.numeric(substr(d1.agg.e$treatment, 1, 2))
  }
if(expt==5){
  #code levels
  treats= matrix(unlist(strsplit(d1.agg.e$treatment, split = "_")),ncol=3,byrow=T)
  colnames(treats)=c("hotdays","normaldays","first") #first: 1 is n, 2 is h
  d1.agg.e= cbind(d1.agg.e, treats)
  #d1.agg.e$hotdays <- revalue(d1.agg.e$hotdays, c("1" = "hotdays: 1", "2" = "hotdays: 2", "3" = "hotdays:3"))
  
  #put other first in supplement
  d1.agg.e= d1.agg.e[d1.agg.e$first==1,]
}
  
#get legend
  if(expt==1){
    fplot= ggplot(data=d1.agg.e, aes(x=treatment, y =value, color=metric, group=metric))+
      geom_point(size=2)+geom_line(lwd=1.5)+
      theme_bw(base_size=16) +theme(legend.position = "right")+scale_color_brewer(palette="Dark2")+guides(colour = guide_legend(nrow = 3))
    
    # Extract the legend. Returns a gtable
    leg <- get_legend(fplot)
    
    # Convert to a ggplot and print
    lplot<- as_ggplot(leg)
  }

#plot
xlabs<-c("Tmin (°C)","Tmax (°C)","Tvar (°C)","Tmean (°C)","# normal days","heatwave length", "heatwave length")  
  
fplot= ggplot(data=d1.agg.e, aes(x=treatment, y =value, color=metric, group=metric))+
  geom_point(size=2)+geom_line(lwd=1.5)+
  theme_bw(base_size=16) +theme(legend.position = "none")+scale_color_brewer(palette="Dark2")+guides(colour = guide_legend(nrow = 3))+
  labs(title=d1.agg.e$elab)+ylab("fecundity")+xlab(xlabs[expt])

if(expt==4){
  d1.agg.e$tmet<- paste(d1.agg.e$metric, d1.agg.e$tvar, sep="_")
  fplot= ggplot(data=d1.agg.e, aes(x=treatment, y =value, color=metric, lty=factor(tvar), group=tmet))+
    geom_point(size=2)+geom_line(lwd=1.5)+
    theme_bw(base_size=16) +theme(legend.position = c(0.9,0.7),legend.background=element_blank())+
    scale_color_brewer(palette="Dark2")+guides(colour ="none")+
    labs(title=d1.agg.e$elab, lty ="Tvar (°C)")+ylab("fecundity")+xlab(xlabs[expt])
}

if(expt==5){
  d1.agg.e$tmet<- paste(d1.agg.e$metric, d1.agg.e$hotdays, sep="_")
  fplot= ggplot(data=d1.agg.e, aes(x=normaldays, y =value, color=metric, group=tmet, lty=hotdays))+
    geom_point(size=2)+geom_line(lwd=1.5)+
    theme_bw(base_size=16) +theme(legend.position = c(0.15,0.7),legend.background=element_blank())+
    scale_color_brewer(palette="Dark2")+guides(colour = "none")+
    labs(title=d1.agg.e$elab, lty ="# hot days")+ylab("fecundity")+xlab(xlabs[expt])
}

if(expt==1) fplot.e1<- fplot
if(expt==2) fplot.e2<- fplot
if(expt==3) fplot.e3<- fplot
if(expt==4) fplot.e4<- fplot
if(expt==5) fplot.e5<- fplot
if(expt==6) fplot.e6<- fplot
if(expt==7) fplot.e7<- fplot
}

#write out plot
pdf("Fig5_fecundity.pdf",height = 10, width = 10)
print( fplot.e1+fplot.e2+fplot.e3+fplot.e4+fplot.e5+fplot.e6+fplot.e7+lplot+
         plot_layout(ncol=2) ) #+ plot_annotation(tag_levels = 'A')
dev.off()

#===================
#Fig 6. Development comparison

#performance metric
pms<- c("dr", "sur", "long", "fec")
pm.ind<- 1

#scen: #1. baseline fit scale; 2. fix scale; 3. fit tp; 4. drop c1; 5. drop c2 with floor
if(pm.ind==1) scens= c(1,1,1,1,1,NA,1)  #c(3,4,5,2,1,NA,5) #dev_rate

#extract performance values
fecs<- PerfDat[PerfDat$metric=="dev_rate",]

#----------
for(expt in c(1:5,7)){
  temps.expt<- temps.all[temps.all$expt==expt,]
  #drop first=2 for experiment 5
  if(expt==5){
    #put other first in supplement
    treats= matrix(unlist(strsplit(temps.expt$treatment, split = "_")),ncol=3,byrow=T)
    temps.expt= temps.expt[treats[,3]==1,]
  }
  
  cs<- as.numeric(out.dr[which(out.dr$expt==expt & out.dr$scenario==scens[expt]),4:10])
  
  p1= perf.damage(pm=pm.ind, temps.expt$temp, c1=cs[1], c2=cs[2], c3=cs[3], c4=cs[4], tp=cs[5], tr=cs[6], scale=cs[7])
  p1.nd= perf.nodamage(pm=pm.ind, temps.expt$temp, scale=cs[7])
  
  d1<- data.frame(metric="perf.nd",value=p1.nd, time=temps.expt$time, treatment=temps.expt$treatment) 
  d2<- data.frame(metric="perf",value=p1, time=temps.expt$time, treatment=temps.expt$treatment)
  d1<- rbind(d1,d2)        
  
  #combine across experiments
  d1$expt<- expt
  if(expt==1) d1.all<- d1
  if(expt>1) d1.all<- rbind(d1.all, d1)
} #end loop experiment

#aggregate performance estimates
d1.agg= aggregate(.~metric+treatment+expt, d1.all, mean)

#---------
#add observed
fdat<- fecs[,c("metric","treatment", "value","expt")]
fdat<- aggregate(.~metric+treatment+expt, fdat, mean)
d1.agg<- d1.agg[,c("metric","treatment", "value","expt")]
d1.agg<- rbind(d1.agg, fdat)
d1.agg$elab<- paste("expt",d1.agg$expt,sep=" ")

#rename metrics
metrs<- c("dev_rate","perf","perf.nd")
metr.lab<- c("observed","performance","performance without damage")
d1.agg$metric= metr.lab[match(d1.agg$metric, metrs)]

for(expt in c(1:5,7)){
  d1.agg.e<- d1.agg[which(d1.agg$expt==expt),]
  
  #adjust treatments
  if(expt==3) {d1.agg.e$treatment<- gsub("22_", "", d1.agg.e$treatment); d1.agg.e$treatment= as.numeric(d1.agg.e$treatment)}
  if(expt==4){
    d1.agg.e$tvar<- as.numeric(substr(d1.agg.e$treatment, 4, 5))
    d1.agg.e$treatment<- as.numeric(substr(d1.agg.e$treatment, 1, 2))
  }
  if(expt==5){
    #code levels
    treats= matrix(unlist(strsplit(d1.agg.e$treatment, split = "_")),ncol=3,byrow=T)
    colnames(treats)=c("hotdays","normaldays","first") #first: 1 is n, 2 is h
    d1.agg.e= cbind(d1.agg.e, treats)  
    
    #put other first in supplement
    d1.agg.e= d1.agg.e[d1.agg.e$first==1,]
  }
  
  #get legend
  if(expt==1){
    fplot= ggplot(data=d1.agg.e, aes(x=treatment, y =value, color=metric, group=metric))+
      geom_point(size=2)+geom_line(lwd=1.5)+
      theme_bw(base_size=16) +theme(legend.position = "right")+scale_color_brewer(palette="Dark2")+guides(colour = guide_legend(nrow = 3))
    
    # Extract the legend. Returns a gtable
    leg <- get_legend(fplot)
    
    # Convert to a ggplot and print
    lplot<- as_ggplot(leg)
  }
  
  #plot
  xlabs<-c("Tmin","Tmax","Tvar","Tmean","# normal days","heatwave length", "heatwave length")  
  
  fplot= ggplot(data=d1.agg.e, aes(x=treatment, y =value, color=metric, group=metric))+
    geom_point(size=2)+geom_line(lwd=1.5)+
    theme_bw(base_size=16) +theme(legend.position = "none")+scale_color_brewer(palette="Dark2")+guides(colour = guide_legend(nrow = 3))+
    labs(title=d1.agg.e$elab)+ylab("development rate (1/day)")+xlab(xlabs[expt])
  
  if(expt==4){
    d1.agg.e$tmet<- paste(d1.agg.e$metric, d1.agg.e$tvar, sep="_")
    fplot= ggplot(data=d1.agg.e, aes(x=treatment, y =value, color=metric, lty=factor(tvar), group=tmet))+
      geom_point(size=2)+geom_line(lwd=1.5)+
      theme_bw(base_size=16) +theme(legend.position = c(0.9,0.70),legend.background=element_blank())+scale_color_brewer(palette="Dark2")+guides(colour ="none")+
      labs(title=d1.agg.e$elab, lty ="Tvar (°C)")+ylab("development rate (1/day)")+xlab(xlabs[expt])
  }
  
  if(expt==5){
    d1.agg.e$tmet<- paste(d1.agg.e$metric, d1.agg.e$hotdays, sep="_")
    fplot= ggplot(data=d1.agg.e, aes(x=normaldays, y =value, color=metric, group=tmet, lty=hotdays))+
      geom_point(size=2)+geom_line(lwd=1.5)+
      theme_bw(base_size=16) +theme(legend.position = c(0.15,0.70),legend.background=element_blank())+scale_color_brewer(palette="Dark2")+guides(colour = "none")+
      labs(title=d1.agg.e$elab, lty ="# hot days")+ylab("development rate (1/day)")+xlab(xlabs[expt])
  }
  
  if(expt==1) fplot.e1<- fplot
  if(expt==2) fplot.e2<- fplot
  if(expt==3) fplot.e3<- fplot
  if(expt==4) fplot.e4<- fplot
  if(expt==5) fplot.e5<- fplot
  if(expt==7) fplot.e7<- fplot
}

#write out plot
pdf("Fig6_dev_rate.pdf",height = 10, width = 10)
print( fplot.e1+fplot.e2+fplot.e3+fplot.e4+fplot.e5+fplot.e7+lplot+
         plot_layout(ncol=2) ) #+ plot_annotation(tag_levels = 'A')
dev.off()
