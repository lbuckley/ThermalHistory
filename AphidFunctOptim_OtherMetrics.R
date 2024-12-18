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
library(dfoptim)

#toggle between desktop (y) and laptop (n)
desktop<- "n"

#performance metric
pms<- c("dr", "sur", "long", "fec")
#pick metric
pm.ind<- 4

#FIT FUNCTION 
if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/out/")
if(desktop=="n") setwd("/Users/lbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/out/") 

temps.all<- read.csv("TempTimeSeries.csv")
PerfDat<- read.csv("PerformanceData.csv")

#-------------
#TPCs

#Ma et al. 2021. Are extreme high temperatures at low or high latitudes more likely to inhibit the population growth of a globally distributed aphid?
#https://doi.org/10.1016/j.jtherbio.2021.102936
#aphid Rhopalosiphum padi
#Reproductive rate (nymphs/adult/day)
Rr= function(T, a=0.263, b=3.7, T0 =30.1) exp(a*T)-exp(a*T0-(T0-T)/b)

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

#====================
#FUNCTIONS
#damage
# tp: threshold for damage between Topt and CTmax; Tdamage= Topt + (CTmax-Topt)*tp
# c1: d_mult: multiplicative change in damage
# c2: d_linear: linear increase in damage
# c3: r_mag: magnitude of repair
# c4: r_breadth: breadth of repair function around Topt

#find Topt and CTmax
ts=seq(0,40,0.1)

if(pm.ind==1) ft= dr(ts)
if(pm.ind==2) ft= sur(ts)
if(pm.ind==3) ft= long(ts) 
if(pm.ind==4) ft= fec(ts) 

topt<- ts[which.max(ft)]
ctmax= ts[which(ft[120:length(ft)]==0)[1]+120]
ctmin= ts[which(ft>0)[1]-1]

# #plot
# plot(ts, dr(ts))
# plot(ts, sur(ts))
# plot(ts, long(ts))
# plot(ts, fec(ts))

#--------------------
perf.damage<- function(pm, T,c1,c2,c3,c4,tp=0,scale,Topt=topt, CTmax=ctmax)  
{ 
  p=NA
  damage=0
  
  Tdamage= Topt + (CTmax-Topt)*tp
  Tdif= T-Tdamage
  if(length(which(Tdif<0))>0) Tdif[which(Tdif<0)]<- 0
  
  for(i in 1:length(T)){
    #damage
    dur<- ifelse(Tdif[i]>0, 1, 0)
    damage.n<- 1- exp(-c1*dur-c2*Tdif[i])
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

#perf.damage(pm=pm.ind, T=temps.all[which(temps.all$expt==1 & temps.all$treatment==13),"temp"], c1=cs[k,1], c2=cs[k,2], c3=cs[k,3], c4=cs[k,4], scale=cs[k,5],Topt=topt, CTmax=ctmax)

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
cs<- expand.grid(c1=seq(0, 1, 0.25), c2= seq(0, .01, 0.003), c3= seq(0, 1, 0.25), c4= seq(1, 5, 1),
                 scale= 1 )

#fit values
#cs<- expand.grid(c1=c(1.95,2), c2= c(0.0007, 0.001), c3= c(0.25,0.66), c4= c(1.1, 1.3), scale= 0.01)
cs<- expand.grid(c1=c(0.001,0.01), c2= c(0.001,0.01), c3= c(0.2,0.9), c4= c(1, 3), scale= 0.01)

for(k in 1:nrow(cs)){
  p1= perf.damage(pm=pm.ind, temps, c1=cs[k,1], c2=cs[k,2], c3=cs[k,3], c4=cs[k,4], scale=cs[k,5])
  ps= cbind(time=1:length(temps), temps, p1, k, cs[k,])
  if(k==1) ps.all<- ps
  if(k>1) ps.all<- rbind(ps.all, ps)
}

funct.fig<- ggplot(data=ps.all, aes(x=time, y =p1, color=c3, lty=factor(c4), group=k))+
  geom_line()+facet_grid(c2~c1)+theme_bw()+
  ylab("Performance")+scale_color_viridis()

#==================
#FIT MODEL, compare AIC of different assumptions
#extract fecundity values
if(pm.ind==1) fecs.all<- PerfDat[PerfDat$metric=="dev_rate",]
if(pm.ind==2) fecs.all<- PerfDat[PerfDat$metric=="survival",]
if(pm.ind==3) fecs.all<- PerfDat[PerfDat$metric=="longevity",]
if(pm.ind==4) fecs.all<- PerfDat[PerfDat$metric=="fecundity",]

#compare AICs of fits
#1. fit 4 params
#2. fit tp
#3. drop c1
#4. drop c2 with floor for damage c2=0.000001

#store output
opts.scale= array(NA, dim=c(7,4,6), dimnames = list(c("e1","e2","e3","e4","e5","e6","e7"), c("s1","s2","s3","s4"), c("c1","c2","c3","c4","tp","scale")))
opts= array(NA, dim=c(7,4,6), dimnames = list(c("e1","e2","e3","e4","e5","e6","e7"), c("s1","s2","s3","s4"), c("c1","c2","c3","c4","tp","scale")))
fit= array(NA, dim=c(7,4,3), dimnames = list(c("e1","e2","e3","e4","e5","e6","e7"), c("s1","s2","s3","s4"), c("sse","convergence","aic")))

#loop through 7 experiments 
for(expt in c(1:7)){ 

  fecs<- fecs.all[fecs.all$expt==expt,]
  tempse<- temps.all[temps.all$expt==expt,]
  
#check that data exist
if(length(unique(fecs[fecs$expt==expt,"treatment"]))>0){
  
  #account for field and lab populations in Figure 4
  #drop treatments with no estimated performance
  if(expt==4){
    fecs<- fecs[which(fecs$population=="lab"),] #field estimates large
  fecs<- fecs[-which(fecs$treatment %in% c("30_0","32_0")),]
  tempse<- tempse[-which(tempse$treatment %in% c("30_0","32_0")),]
  }
  
  #scale
  scale.est= 1
  
  #if performance with out damage is less than observed fecundity, adjust scale
  #performance estimation by treatment
  tempse$p.nd<- perf.nodamage(pm=pm.ind, tempse[,"temp"], scale=1)
  #mean by treatment
  p.nd.t<- aggregate(tempse, list(treat=tempse$treatment), FUN="mean")
  f.t<- aggregate(fecs, list(treat=fecs$treatment), FUN="mean")
  
  match1<- match(p.nd.t$treat, f.t$treat)
  p.nd.t$fec.ratio= f.t$value[match1]/p.nd.t$p.nd
  if(max(p.nd.t$fec.ratio)>1) scale.est<- max(p.nd.t$fec.ratio)
 
  #-----------
#optimize
  #1. fit scale four parameters
  #error function
  errs<- function(x,temps=tempse, fecundity=fecs, scale=scale.est){  
    totalerror=0
    treats=unique(temps$treatment)
    for(i in 1:length(treats)){
      perfs=perf.damage(pm.ind, T=temps[temps$treatment==treats[i],"temp"],c1=x[1],c2=x[2],c3=x[3],c4=x[4],tp=0,scale=scale,Topt=topt, CTmax=ctmax)
      perfs= mean(perfs[96:length(perfs)])
      delta= fecundity[which(fecundity$treatment==treats[i]),"value"]- perfs
      totalerror=totalerror + sum(delta^2)
      }
    return( sqrt(totalerror) )
  }
  
  opt<- nmkb(fn=errs, par=c(0.001,0.001,0.1,1), lower=c(0,0,0,0), upper=c(1,1,1,3)) #c(1.5,1,1,3)
  
  #store output and fits
  opts[expt,1,]<- c(opt$par[1:4], 0, scale.est)
  fit[expt,1,1:2]<- c(opt$value, opt$convergence)
  
  #2. fit tp
  errs<- function(x,temps=tempse, fecundity=fecs, scale=scale.est){
    totalerror=0
    treats=unique(temps$treatment)
    for(i in 1:length(treats)){
      perfs=perf.damage(pm.ind, T=temps[temps$treatment==treats[i],"temp"],c1=x[1],c2=x[2],c3=x[3],c4=x[4],tp=x[5],scale=scale,Topt=topt, CTmax=ctmax)
      perfs= mean(perfs[96:length(perfs)])
      delta= fecundity[which(fecundity$treatment==treats[i]),"value"]- perfs
      totalerror=totalerror + sum(delta^2)
    }
    return( sqrt(totalerror) )
  }

  opt<- nmkb(fn=errs, par=c(1,0.001,0.1,1,0.5), lower=c(0,0.000001,0,0,0), upper=c(1.5,2,1,3,1) )

  #store output and fits
  opts[expt,2,]<- c(opt$par, scale.est)
  fit[expt,2,1:2]<- c(opt$value, opt$convergence)

  #3. drop c1
  errs<- function(x,temps=tempse, fecundity=fecs, scale=scale.est){
    totalerror=0
    treats=unique(temps$treatment)
    for(i in 1:length(treats)){
      perfs=perf.damage(pm.ind, T=temps[temps$treatment==treats[i],"temp"],c1=0,c2=x[1],c3=x[2],c4=x[3],tp=0,scale=scale,Topt=topt, CTmax=ctmax)
      perfs= mean(perfs[96:length(perfs)])
      delta= fecundity[which(fecundity$treatment==treats[i]),"value"]- perfs
      totalerror=totalerror + sum(delta^2)
    }
    return( sqrt(totalerror) )
  }

  opt<- nmkb(fn=errs, par=c(0.001,0.1,1), lower=c(0.000001,0,0), upper=c(2,1,1.5) )

  #store output and fits
  opts[expt,3,]<- c(0, opt$par, 0, scale.est)
  fit[expt,3,1:2]<- c(opt$value, opt$convergence)

  #4. drop c2 with floor for damage c2=0.000001
  errs<- function(x, temps=tempse, fecundity=fecs, scale=scale.est){
    totalerror=0
    treats=unique(temps$treatment)
    for(i in 1:length(treats)){
      perfs=perf.damage(pm.ind, T=temps[temps$treatment==treats[i],"temp"],c1=x[1],c2=0.0005,c3=x[2],c4=x[3],tp=0,scale=scale,Topt=topt, CTmax=ctmax)
      perfs= mean(perfs[96:length(perfs)])
      delta= fecundity[which(fecundity$treatment==treats[i]),"value"]- perfs
      totalerror=totalerror + sum(delta^2)
    }
    return( sqrt(totalerror) )
  }

  opt<- nmkb(fn=errs, par=c(1,0.1,1), lower=c(0,0,0), upper=c(3,1,1.5) )

  #store output and fits
  opts[expt,4,]<- c(opt$par[1], 0.000001, opt$par[2:3], 0, scale.est)
  fit[expt,4,1:2]<- c(opt$value, opt$convergence)

  #---------
  #estimate AICs, https://stackoverflow.com/questions/39999456/aic-on-nls-on-r
  
  scen.params<- c(4,5,3,3)
  
  for(scen in 1:4){
    logL <- 0.5 *(- nrow(fecs) * (log(2*pi)+1-log(nrow(fecs)) + log(fit[expt,scen,1])))
    fit[expt, scen,3]= 2*(scen.params[scen] + 1) - 2 * logL 
    }
   
} #end check data exists
} #end loop experiments
  
  #-----------------
  #Construct table
  expt1<- cbind(expt="1", scenario=1:4, opts[1,,], fit[1,,])
  expt2<- cbind(expt="2", scenario=1:4, opts[2,,], fit[2,,])
  expt3<- cbind(expt="3", scenario=1:4, opts[3,,], fit[3,,])
  expt4<- cbind(expt="4", scenario=1:4, opts[4,,], fit[4,,])
  expt5<- cbind(expt="5", scenario=1:4, opts[5,,], fit[5,,])
  expt6<- cbind(expt="6", scenario=1:4, opts[6,,], fit[6,,])
  expt7<- cbind(expt="7", scenario=1:4, opts[7,,], fit[7,,])
  
  out<- rbind(expt1, expt2, expt3, expt4, expt5, expt6, expt7)
  colnames(out)[3:ncol(out)]<- c("d_mult","d_linear","r_mag","r_breadth","tp","scale","sse","converge?","AIC")
  out<- as.data.frame(out)
  out[,c(2:3,5:8)]<- round(as.numeric(unlist(out[,c(2:3,5:8)])), 4)
  out[,4]<- round(as.numeric(unlist(out[,4])), 6)
  out[,c(9,11)]<- round(as.numeric(unlist(out[,c(9,11)])),0)
  
  out.scale<- rbind(opts.scale[1,,], opts.scale[2,,], opts.scale[3,,])
  
  #save output
  if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/out/")
  if(desktop=="n") setwd("/Users/lbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/out/") 
  
  out_file <- paste("out_", pms[pm.ind], ".csv", sep="")
  write.csv(out, out_file)
  out_file <- paste("opts_scale_", pms[pm.ind], ".csv", sep="")
  write.csv(out.scale, out_file)
  
  #optimization options
  #efficient package: https://cran.r-project.org/web/packages/lbfgs/vignettes/Vignette.pdf
  #https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/
  
  scen1= cbind(opts[,1,], fit[,1,1])
  scen1[,c(3:4,6:7)]= round(scen1[,c(3:4,6:7)],2) 
  