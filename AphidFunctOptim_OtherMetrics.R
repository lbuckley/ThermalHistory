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
desktop<- "y"

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

#performance metric
pms<- c("dr", "sur", "long", "fec")
pm.ind<- 1

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

#old damage functions
#damagenew= damage + dt * exp(-Ea/(R*(T+273.15))) *10^9* (c1*damage + c2) + dt*c3

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

#perf(pm=pm.ind, series=temps[,"temp"], c1=cs[k,1], c2=cs[k,2], c3=cs[k,3], c4=cs[k,4], scale=cs[k,5])

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

computeperf<- function(series,c1,c2,c3,c4,tp=0,scale,printdam=FALSE)  {
  p=0
  damage=0
  for(i in 1:length(series)){
    if(pm==1) perf= dr(series[i])
    if(pm==2) perf= sur(series[i])
    if(pm==3) perf= long(series[i])
    if(pm==4) perf= fec(series[i])
    
    damage=damage.rep(damage,T=series[i],c1=c1,c2=c2,c3=c3,c4=c4,tp=tp,dt=1)
  p= p + perf*(1-damage)
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
  p1= perf(pm=pm.ind, temps, c1=cs[k,1], c2=cs[k,2], c3=cs[k,3], c4=cs[k,4], scale=cs[k,5])
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
if(pm.ind==1) fecs<- PerfDat[PerfDat$metric=="dev_rate",]
if(pm.ind==2) fecs<- PerfDat[PerfDat$metric=="survival",]
if(pm.ind==3) fecs<- PerfDat[PerfDat$metric=="longevity",]
if(pm.ind==4) fecs<- PerfDat[PerfDat$metric=="fecundity",]

#compare AICs of fits
#1. baseline: fit scale
#2. fix scale
#3. fit tp
#4. drop c1
#5. drop c2 with floor for damage c2=0.000001
#compare AIC of scenarios 1 and 2 to determine whether to fix scale

#store output
opts.scale= array(NA, dim=c(3,3,6), dimnames = list(c("expt", "scenario","params")))
opts= array(NA, dim=c(3,6,6), dimnames = list(c("expt", "scenario","params")))
fit= array(NA, dim=c(3,6,2), dimnames = list(c("expt", "scenario","fit"))) #aic and convergence

#loop through 3 experiments
for(expt in 1:3){

#check that data exist
if(length(unique(fecs[fecs$expt==expt,"treatment"]))>0){
  
#estimate scale as max of performance
  if(pm.ind==1) scale.est<- max(fecs[fecs$expt==expt,"value"])/(sum(dr(temps.all[temps.all$expt==expt,"temp"]))/length(unique(temps.all[temps.all$expt==expt,"treatment"])))
  if(pm.ind==2) scale.est<- max(fecs[fecs$expt==expt,"value"])/(sum(sur(temps.all[temps.all$expt==expt,"temp"]))/length(unique(temps.all[temps.all$expt==expt,"treatment"])))
  if(pm.ind==3) scale.est<- max(fecs[fecs$expt==expt,"value"])/(sum(long(temps.all[temps.all$expt==expt,"temp"]))/length(unique(temps.all[temps.all$expt==expt,"treatment"])))
  if(pm.ind==4) scale.est<- max(fecs[fecs$expt==expt,"value"])/(sum(fec(temps.all[temps.all$expt==expt,"temp"]))/length(unique(temps.all[temps.all$expt==expt,"treatment"])))

  #account for field and lab populations in Figure 3
  if(expt==3) fecs<- fecs[which(fecs$population=="field"),]
  
  #-----------
  #Fit scale
  #Fit c1, c2, c3, c4 at scale.est/2, scale.est, scale.est*2
  #Then fix scale for subsequent fits
  
  errs<- function(x,temps=temps.all[temps.all$expt==expt,], fecundity=fecs[fecs$expt==expt,], scale=scale.est){  
    totalerror=0
    treats=unique(temps$treatment)
    for(i in 1:length(treats)){
      delta=computeperf(series=temps[temps$treatment==treats[i],"temp"],c1=x[1],c2=x[2],c3=x[3],c4=x[4],scale=scale)-mean(fecundity[which(fecundity$treatment==treats[i]),"value"])
      #totalerror=totalerror + delta^2
      #return( sqrt(totalerror) )
      
      #try AIC function: https://optimumsportsperformance.com/blog/optimization-algorithms-in-r-returning-model-fit-metrics/
      totalerror= totalerror + length(temps[temps$treatment==treats[i],"temp"])*(log(2*pi)+1+log((sum(delta^2)/length(temps[temps$treatment==treats[i],"temp"])))) + ((length(x)+1)*2)
    }
    return(totalerror)
  }
  
  opt<- optim(par=c(1,0.001,0.1,1), fn=errs, NULL, method=c("L-BFGS-B"), 
              lower=c(0,0.0001,0,0), upper=c(2,1,1,3) )
  
  opts.scale[expt,1,]<- c(opt$par[1:4], 0, scale.est)
  
  #save at fixed scale
  opts[expt,6,]<- c(opt$par[1:4], 0, scale.est)
  fit[expt,6,]<- c(opt$value, opt$convergence)
  
  errs<- function(x,temps=temps.all[temps.all$expt==expt,], fecundity=fecs[fecs$expt==expt,], scale=scale.est*2){  
    totalerror=0
    treats=unique(temps$treatment)
    for(i in 1:length(treats)){
      delta=computeperf(series=temps[temps$treatment==treats[i],"temp"],c1=x[1],c2=x[2],c3=x[3],c4=x[4],scale=scale)-mean(fecundity[which(fecundity$treatment==treats[i]),"value"])
      #totalerror=totalerror + delta^2
      #return( sqrt(totalerror) )
      
      #try AIC function: https://optimumsportsperformance.com/blog/optimization-algorithms-in-r-returning-model-fit-metrics/
      totalerror= totalerror + length(temps[temps$treatment==treats[i],"temp"])*(log(2*pi)+1+log((sum(delta^2)/length(temps[temps$treatment==treats[i],"temp"])))) + ((length(x)+1)*2)
    }
    return(totalerror)
  }
  
  opt<- optim(par=c(1,0.001,0.1,1), fn=errs, NULL, method=c("L-BFGS-B"), 
              lower=c(0,0.0001,0,0), upper=c(2,1,1,3) )
  
  opts.scale[expt,2,]<- c(opt$par[1:4], 0, scale.est*1.5)
  
  errs<- function(x,temps=temps.all[temps.all$expt==expt,], fecundity=fecs[fecs$expt==expt,], scale=scale.est*4){  
    totalerror=0
    treats=unique(temps$treatment)
    for(i in 1:length(treats)){
      delta=computeperf(series=temps[temps$treatment==treats[i],"temp"],c1=x[1],c2=x[2],c3=x[3],c4=x[4],scale=scale)-mean(fecundity[which(fecundity$treatment==treats[i]),"value"])
      #totalerror=totalerror + delta^2
      #return( sqrt(totalerror) )
      
      #try AIC function: https://optimumsportsperformance.com/blog/optimization-algorithms-in-r-returning-model-fit-metrics/
      totalerror= totalerror + length(temps[temps$treatment==treats[i],"temp"])*(log(2*pi)+1+log((sum(delta^2)/length(temps[temps$treatment==treats[i],"temp"])))) + ((length(x)+1)*2)
    }
    return(totalerror)
  }
  
  opt<- optim(par=c(1,0.001,0.1,1), fn=errs, NULL, method=c("L-BFGS-B"), 
              lower=c(0,0.0001,0,0), upper=c(2,1,1,3) )
  
  opts.scale[expt,3,]<- c(opt$par[1:4], 0, scale.est*2)
  
  #Estimate scale using average parameter values
  params<- colMeans(opts.scale[expt,,])
  
  errs<- function(x,temps=temps.all[temps.all$expt==expt,], fecundity=fecs[fecs$expt==expt,], c1=params[1], c2=params[2], c3=params[3], c4=params[4]){  
    totalerror=0
    treats=unique(temps$treatment)
    for(i in 1:length(treats)){
      delta=computeperf(series=temps[temps$treatment==treats[i],"temp"],c1=c1,c2=c2,c3=c3,c4=c4,scale=x[1])-mean(fecundity[which(fecundity$treatment==treats[i]),"value"])
      #try AIC function: https://optimumsportsperformance.com/blog/optimization-algorithms-in-r-returning-model-fit-metrics/
      totalerror= totalerror + length(temps[temps$treatment==treats[i],"temp"])*(log(2*pi)+1+log((sum(delta^2)/length(temps[temps$treatment==treats[i],"temp"])))) + ((length(x)+1)*2)
    }
    return(totalerror)
  }
  
  opt<- optim(par=c(scale.est), fn=errs, NULL, method=c("L-BFGS-B"), 
              lower=c(scale.est/2), upper=c(scale.est*2) )
  
  if(opt$convergence !=0){
    opt<- optim(par=c(scale.est), fn=errs, NULL, method=c("BFGS") )
  }
  
  #update estimate
  scale.est<- opt$par[1]
  
  #-----------
#optimize
  #1. fit scale four parameters
  #error function
  errs<- function(x,temps=temps.all[temps.all$expt==expt,], fecundity=fecs[fecs$expt==expt,], scale=scale.est){  
    totalerror=0
    treats=unique(temps$treatment)
    for(i in 1:length(treats)){
      delta=computeperf(series=temps[temps$treatment==treats[i],"temp"],c1=x[1],c2=x[2],c3=x[3],c4=x[4],scale=scale)-mean(fecundity[which(fecundity$treatment==treats[i]),"value"])
      #totalerror=totalerror + delta^2
      #return( sqrt(totalerror) )
      
      #try AIC function: https://optimumsportsperformance.com/blog/optimization-algorithms-in-r-returning-model-fit-metrics/
      totalerror= totalerror + length(temps[temps$treatment==treats[i],"temp"])*(log(2*pi)+1+log((sum(delta^2)/length(temps[temps$treatment==treats[i],"temp"])))) + ((length(x)+1)*2)
    }
    return(totalerror)
  }
  
  opt<- optim(par=c(1,0.001,0.1,1), fn=errs, NULL, method=c("L-BFGS-B"), 
              lower=c(0,0.000001,0,0), upper=c(5,2,1,5) )
  
  if(opt$convergence !=0){
    opt<- optim(par=c(1,0.001,0.1,1), fn=errs, NULL, method=c("BFGS") )
  }
  
  opts[expt,1,]<- c(opt$par[1:4], 0, scale.est)
  fit[expt,1,]<- c(opt$value, opt$convergence)
  
  #2. fit tp
  errs<- function(x,temps=temps.all[temps.all$expt==expt,], fecundity=fecs[fecs$expt==expt,], scale=scale.est){  
    totalerror=0
    treats=unique(temps$treatment)
    for(i in 1:length(treats)){
      delta=computeperf(series=temps[temps$treatment==treats[i],"temp"],c1=x[1],c2=x[2],c3=x[3],c4=x[4],tp=x[5],scale=scale)-mean(fecundity[which(fecundity$treatment==treats[i]),"value"])
      
      #try AIC function: https://optimumsportsperformance.com/blog/optimization-algorithms-in-r-returning-model-fit-metrics/
      totalerror= totalerror + length(temps[temps$treatment==treats[i],"temp"])*(log(2*pi)+1+log((sum(delta^2)/length(temps[temps$treatment==treats[i],"temp"])))) + ((length(x)+1)*2)
    }
    return(totalerror)
  }
  
  opt<- optim(par=c(1,0.001,0.1,1,0), fn=errs, NULL, method=c("L-BFGS-B"), 
              lower=c(0,0.000001,0,0,0), upper=c(5,2,1,5,1) )
  
  if(opt$convergence !=0){
    opt<- optim(par=c(1,0.001,0.1,1,0), fn=errs, NULL, method=c("BFGS") )
  }
  
  #store output and fits
  opts[expt,2,]<- c(opt$par, scale.est)
  fit[expt,2,]<- c(opt$value, opt$convergence)
  
  #3. drop c1
  errs<- function(x,temps=temps.all[temps.all$expt==expt,], fecundity=fecs[fecs$expt==expt,], scale=scale.est){  
    totalerror=0
    treats=unique(temps$treatment)
    for(i in 1:length(treats)){
      delta=computeperf(series=temps[temps$treatment==treats[i],"temp"],c1=0,c2=x[1],c3=x[2],c4=x[3],scale=scale.est)-mean(fecundity[which(fecundity$treatment==treats[i]),"value"])
      
      #try AIC function: https://optimumsportsperformance.com/blog/optimization-algorithms-in-r-returning-model-fit-metrics/
      totalerror= totalerror + length(temps[temps$treatment==treats[i],"temp"])*(log(2*pi)+1+log((sum(delta^2)/length(temps[temps$treatment==treats[i],"temp"])))) + ((length(x)+1)*2)
    }
    return(totalerror)
  }
  
  opt<- optim(par=c(0.001,0.1,1), fn=errs, NULL, method=c("L-BFGS-B"), 
              lower=c(0.000001,0,1), upper=c(2,1,3) )
  
  if(opt$convergence !=0){
    opt<- optim(par=c(0.001,0.1,1), fn=errs, NULL, method=c("BFGS") )
  }
  
  #store output and fits
  opts[expt,3,]<- c(0, opt$par, 0, scale.est)
  fit[expt,3,]<- c(opt$value, opt$convergence)
  
  #4. drop c2 with floor for damage c2=0.000001
  errs<- function(x,temps=temps.all[temps.all$expt==expt,], fecundity=fecs[fecs$expt==expt,], scale=scale.est){  
    totalerror=0
    treats=unique(temps$treatment)
    for(i in 1:length(treats)){
      delta=computeperf(series=temps[temps$treatment==treats[i],"temp"],c1=x[1],c2=0.0005,c3=x[2],c4=x[3],scale=scale.est)-mean(fecundity[which(fecundity$treatment==treats[i]),"value"])
      
      #try AIC function: https://optimumsportsperformance.com/blog/optimization-algorithms-in-r-returning-model-fit-metrics/
      totalerror= totalerror + length(temps[temps$treatment==treats[i],"temp"])*(log(2*pi)+1+log((sum(delta^2)/length(temps[temps$treatment==treats[i],"temp"])))) + ((length(x)+1)*2)
    }
    return(totalerror)
  }
  
  opt<- optim(par=c(1,0.1,1), fn=errs, NULL, method=c("L-BFGS-B"), 
              lower=c(0,0,0), upper=c(5,1,5) )
  
  if(opt$convergence !=0){
    opt<- optim(par=c(1,0.1,1), fn=errs, NULL, method=c("BFGS"), hessian=TRUE )
  }
  
  #store output and fits
  opts[expt,4,]<- c(opt$par[1], 0.000001, opt$par[2:3], 0, scale.est)
  fit[expt,4,]<- c(opt$value, opt$convergence)
  
  #5 Estimate 4 parameters plus scale
  
  #error function
  errs<- function(x,temps=temps.all[temps.all$expt==expt,], fecundity=fecs[fecs$expt==expt,], scale=scale.est){  
    totalerror=0
    treats=unique(temps$treatment)
    for(i in 1:length(treats)){
      delta=computeperf(series=temps[temps$treatment==treats[i],"temp"],c1=x[1],c2=x[2],c3=x[3],c4=x[4],scale=x[5])-mean(fecundity[which(fecundity$treatment==treats[i]),"value"])
      #totalerror=totalerror + delta^2
      #return( sqrt(totalerror) )
      
      #try AIC function: https://optimumsportsperformance.com/blog/optimization-algorithms-in-r-returning-model-fit-metrics/
      totalerror= totalerror + length(temps[temps$treatment==treats[i],"temp"])*(log(2*pi)+1+log((sum(delta^2)/length(temps[temps$treatment==treats[i],"temp"])))) + ((length(x)+1)*2)
    }
    return(totalerror)
  }
  
  opt<- optim(par=c(1,0.001,0.1,1, scale.est), fn=errs, NULL, method=c("L-BFGS-B"), 
              lower=c(0,0.000001,0,0, scale.est/4), upper=c(5,2,1,5, scale.est*4) )
  
  if(opt$convergence !=0){
    opt<- optim(par=c(1,0.001,0.1,1, scale.est), fn=errs, NULL, method=c("BFGS") )
  }
  
  opts[expt,5,]<- c(opt$par[1:4], 0, opt$par[5])
  fit[expt,5,]<- c(opt$value, opt$convergence)
  
  #95% CI
  #n <- nrow(temps.all[temps.all$expt==expt,])
  #opt$par - 1.96*sqrt(diag(solve(opt$hessian)))/n # lower limit for 95% confint
  #opt$par + 1.96*sqrt(diag(solve(opt$hessian)))/n # upper limit for 95% confint
  
} #end check data exists
} #end loop experiments
  
  #-----------------
  #Construct table
  expt1<- cbind(expt="1", scenario=1:6, opts[1,,], fit[1,,])
  expt2<- cbind(expt="2", scenario=1:6, opts[2,,], fit[2,,])
  expt3<- cbind(expt="3", scenario=1:6, opts[3,,], fit[3,,])
  out<- rbind(expt1, expt2, expt3)
  colnames(out)[3:ncol(out)]<- c("d_mult","d_linear","r_mag","r_breadth","tp","scale","AIC","converge?")
  out<- as.data.frame(out)
  out[,2:8]<- round(as.numeric(unlist(out[,2:8])), 4)
  out[9]<- round(as.numeric(unlist(out[9])),0)
  
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
  
#=====================
#plot performance with values

expt<- 3
#scen: #1. baseline fit scale; 2. fix scale; 3. fit tp; 4. drop c1; 5. drop c2 with floor
scen<- 6

#extract performance values
if(pm.ind==1) fecs<- PerfDat[PerfDat$metric=="dev_rate",]
if(pm.ind==2) fecs<- PerfDat[PerfDat$metric=="survival",]
if(pm.ind==3) fecs<- PerfDat[PerfDat$metric=="longevity",]
if(pm.ind==4) fecs<- PerfDat[PerfDat$metric=="fecundity",]

temps.expt<- temps.all[temps.all$expt==expt,]

p1= perf(pm=pm.ind, temps.expt$temp, c1=opts[expt, scen, 1], c2=opts[expt, scen, 2], c3=opts[expt, scen, 3], c4=opts[expt, scen, 4], tp=opts[expt, scen, 5], scale=opts[expt, scen, 6])
p1.nd= perf.nodamage(pm=pm.ind, temps.expt$temp, scale=opts[expt, scen, 6])

d1<- data.frame(metric="perf.nd",value=p1.nd, time=temps.expt$time, treatment=temps.expt$treatment) 
d2<- data.frame(metric="perf",value=p1, time=temps.expt$time, treatment=temps.expt$treatment)
d3<- data.frame(metric="temp",value=temps.expt$temp, time=temps.expt$time, treatment=temps.expt$treatment)
d1<- rbind(d1,d2,d3)        
          
d1$value <- as.numeric(d1$value)
d1$time <- as.numeric(d1$time)
d1$metric <- revalue(d1$metric, c("temp" = "temperature", "perf.nd" = "performance no damage", "perf" = "performance with damage"))
d1$metric <- factor(d1$metric, ordered=TRUE, levels=c("temperature", "performance no damage", "performance with damage"))

# if(expt==2) d1$treatment <- factor(d1$treatment, ordered=TRUE, levels=c(
#   "1_3_1", "1_2_1", "1_1_1", "2_3_1", "2_2_1", "2_1_1", "3_3_1", "3_2_1", "3_1_1",
#   "1_3_2", "1_2_2", "1_1_2", "2_3_2", "2_2_2", "2_1_2", "3_3_2", "3_2_2", "3_1_2"))
# if(expt==3) d1$treatment <- factor(d1$treatment, ordered=TRUE, levels=c(
#   "22_0", "22_5", "22_9",  "22_13", "28_0", "30_0", "32_0", "28_5", "30_5", "32_5"))

#----------------
#plot time series
#expt 1
if(expt==1){ 
  d1$treatment <- factor(d1$treatment)
  plot1.expt1= ggplot(data=d1, aes(x=time, y =value, color=factor(treatment)))+geom_line(lwd=1.5)+facet_wrap(.~metric, scale="free_y", switch="y", ncol=1)+xlim(0,100)+
  theme_bw(base_size=16) +theme(legend.position = "bottom")+scale_color_viridis(discrete = TRUE)+labs(color="treatment")
}
  
#expt 2
if(expt==2){
  #code levels
  treats= matrix(unlist(strsplit(d1$treatment, split = "_")),ncol=3,byrow=T)
  colnames(treats)=c("hotdays","normaldays","first") #first: 1 is n, 2 is h
  d1= cbind(d1, treats)  
  d1$hotdays <- revalue(d1$hotdays, c("1" = "hotdays: 1", "2" = "hotdays: 2", "3" = "hotdays:3"))
  
  plot1.expt2= ggplot(data=d1[d1$first==1,], aes(x=time, y =value, color=factor(normaldays)))+geom_line(lwd=1.5)+facet_grid(metric~hotdays, scale="free_y", switch="y")+xlim(0,100)+
  theme_bw(base_size=16) +theme(legend.position = "bottom")+scale_color_viridis(discrete = TRUE)+labs(color="normal days")
#put other first in supplement
  }

#expt 3
if(expt==3){ 
  #code levels
  treats= matrix(unlist(strsplit(d1$treatment, split = "_")),ncol=2,byrow=T)
  colnames(treats)=c("mean","variance") 
  d1= cbind(d1, treats)  
  #code experiment
  d1$expt<- "mild means"
  d1$expt[d1$treatment %in% c("28_0", "30_0", "32_0", "28_5", "30_5", "32_5")]<- "high means"
  d1$expt<- factor(d1$expt, ordered=TRUE, levels=c("mild means","high means") )
  
  plot1.expt3= ggplot(data=d1, aes(x=time, y =value, color=factor(treatment), group=treatment))+geom_line(lwd=1.5)+facet_grid(metric~expt, scale="free_y", switch="y")+xlim(0,100)+
    theme_bw(base_size=16) +theme(legend.position = "bottom")+scale_color_viridis(discrete = TRUE)+labs(color="treatment")
}
#----------  
#plot outcomes
if(expt==1){ 
  #aggregate
  d1.agg= aggregate(.~metric+treatment, d1, sum)
  d1.agg= d1.agg[-which(d1.agg$metric=="temperature"), 1:3]
  
  #add observed
  fdat<- fecs[fecs$expt==expt,c("metric","treatment", "value" )]
  fdat<- aggregate(.~metric+treatment, fdat, mean)
  d1.agg<- rbind(d1.agg, fdat)
  
  plot2.expt1= ggplot(data=d1.agg, aes(x=treatment, y =value, color=metric, group=metric))+geom_point(size=2)+geom_line(lwd=1.5)+
  theme_bw(base_size=16) +theme(legend.position = "bottom")+scale_color_brewer(palette="Dark2")+guides(colour = guide_legend(nrow = 3))
}
  
if(expt==2){
  #aggregate
  d1.agg= aggregate(.~metric+treatment, d1[,1:4], sum)
  d1.agg= d1.agg[-which(d1.agg$metric=="temperature"), 1:3]
  
  #add observed
  fdat<- fecs[fecs$expt==expt,c("metric","treatment", "value" )]
  fdat<- aggregate(.~metric+treatment, fdat, mean)
  d1.agg<- rbind(d1.agg, fdat)
  
  #code levels
  treats= matrix(unlist(strsplit(d1.agg$treatment, split = "_")),ncol=3,byrow=T)
  colnames(treats)=c("hotdays","normaldays","first") #first: 1 is n, 2 is h
  d1.agg= cbind(d1.agg, treats)  
  d1.agg$hotdays <- revalue(d1.agg$hotdays, c("1" = "hotdays: 1", "2" = "hotdays: 2", "3" = "hotdays:3"))
  
  plot2.expt2= ggplot(data=d1.agg[d1.agg$first==1,], aes(x=normaldays, y =value, color=metric, group=metric))+geom_point(size=2)+geom_line(lwd=1.5)+
    facet_grid(.~hotdays, scale="free_y", switch="y")+
  theme_bw(base_size=16) +theme(legend.position = "bottom")+scale_color_brewer(palette="Dark2")
  #put other first in supplement
}

if(expt==3){ 
  #aggregate
  d1.agg= aggregate(.~metric+treatment, d1[,1:4], sum)
  d1.agg= d1.agg[-which(d1.agg$metric=="temperature"), 1:3]
  
  #add observed
  fdat<- fecs[fecs$expt==expt,c("metric","treatment", "value" )]
  fdat<- aggregate(.~metric+treatment, fdat, mean)
  d1.agg<- rbind(d1.agg, fdat)
  
  #code experiment
  d1.agg$expt<- "mild means"
  d1.agg$expt[d1.agg$treatment %in% c("28_0", "30_0", "32_0", "28_5", "30_5", "32_5")]<- "high means"
  d1.agg$expt<- factor(d1.agg$expt, ordered=TRUE, levels=c("mild means","high means") )
  
  #code levels
  treats= matrix(unlist(strsplit(d1.agg$treatment, split = "_")),ncol=2,byrow=T)
  colnames(treats)=c("mean","variance") #first: 1 is n, 2 is h
  d1.agg= cbind(d1.agg, treats) 
  
  #code x axis
  d1.agg$treat[d1.agg$expt=="mild means"]<- d1.agg$variance[d1.agg$expt=="mild means"]
  d1.agg$treat[d1.agg$expt=="high means"]<- d1.agg$mean[d1.agg$expt=="high means"]
  d1.agg$treat<- as.numeric(d1.agg$treat)
  
  #code mean
  d1.agg$var<- d1.agg$variance
  d1.agg$var[d1.agg$expt=="mild means"]<- 0
  
  #code group
  d1.agg$group= paste(d1.agg$metric, d1.agg$var, sep="_")
  
  plot2.expt3= ggplot(data=d1.agg, aes(x=treat, y =value, color=metric, group=group, lty=var))+geom_point(size=2)+geom_line(lwd=1.5)+
    facet_grid(.~expt, scale="free", switch="y")+
    theme_bw(base_size=16) +theme(legend.position = "bottom")+scale_color_brewer(palette="Dark2")
}

#write out plot
if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/figures/")
if(desktop=="n") setwd("/Users/lbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/figures/") 

if(expt==1){
  out_file <- paste("AphidsExpt1_", pms[pm.ind], ".pdf", sep="")
  
  pdf(out_file,height = 14, width = 5)
  print(plot1.expt1 +plot2.expt1 +plot_layout(ncol=1, heights = c(3, 1))+ plot_annotation(tag_levels = 'A') )
  dev.off()
  }

if(expt==2){
  out_file <- paste("AphidsExpt2_", pms[pm.ind], ".pdf", sep="")
  
  pdf(out_file,height = 14, width = 14)
  print(plot1.expt2 +plot2.expt2 +plot_layout(ncol=1, heights = c(3, 1)) + plot_annotation(tag_levels = 'A'))
  dev.off()
  }

if(expt==3){
  out_file <- paste("AphidsExpt3_", pms[pm.ind], ".pdf", sep="")
  
  pdf(out_file,height = 14, width = 14)
  print(plot1.expt3 +plot2.expt3 +plot_layout(ncol=1, heights = c(3, 1)) + plot_annotation(tag_levels = 'A'))
  dev.off()
  }

#------------------------------------
#Plot developmental rate comparisons

plot2.expt1= plot2.expt1 + theme(legend.position = "none")
plot2.expt2= plot2.expt2 + theme(legend.position = "none")
plot2.expt3= plot2.expt3 + guides(colour = guide_legend(nrow = 3))

pdf("Fig_DevRate.pdf",height = 10, width = 6)
  print(plot2.expt1 +plot2.expt2 +plot2.expt3 +plot_layout(ncol=1, heights = c(1, 1, 1.2))+ plot_annotation(tag_levels = 'A') )
dev.off()



