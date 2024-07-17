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

#=====
#Zhao et al. 2014. Night warming on hot days produces novel impacts on development, survival and reproduction in a small arthropod
#Dryad data: http://doi.org/10.5061/dryad.q2070 
#English grain aphid, Sitobion avenae

setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/data/aphids/")
#setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/ThermalHistory/data/aphids/")

adat2.dt<- read.csv("Zhaoetal2014/Zhaoetal2014_devtime.csv")
adat2.p<- read.csv("Zhaoetal2014/Zhaoetal2014_AdPerf.csv")
adat2.lt<- read.csv("Zhaoetal2014/Zhaoetal2014_LifeTable.csv")
adat2.sur<- read.csv("Zhaoetal2014/Zhaoetal2014_SurvNymph.csv")

#----
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

#-----
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
  
#----
#dev time
colnames(adat2.dt)[c(1,3,5,7,9)]<-"NTmin"
colnames(adat2.dt)[c(2,4,6,8,10)]<-"dt"
adat2.dt.l<- rbind( cbind(adat2.dt[,1:2], stage="1st"), cbind(adat2.dt[,3:4], stage="2nd"), cbind(adat2.dt[,5:6], stage="3rd"), cbind(adat2.dt[,7:8], stage="4th"),cbind(adat2.dt[,9:10], stage="Nymph") )

#mean metric
adat2.dt.m= adat2.dt.l %>%
  group_by(NTmin, stage) %>%
  summarise(dt= mean(dt))

ggplot(data=adat2.dt.m, aes(x=stage, y =dt, color=factor(NTmin), group=factor(NTmin)))+geom_line()

#as tpc
ggplot(data=adat2.dt.m, aes(x=NTmin, y =1/dt, color=stage, group=stage))+geom_line()

#----
#performance

#mean metric
adat2.p.m= adat2.p %>%
  group_by(NTmin) %>%
  summarise(lon= mean(Longevity, na.rm = TRUE), fec=mean(Fecundtiy, na.rm = TRUE), fec.rate=mean(Fecundity.rate, na.rm = TRUE) )

#to long format
adat2.p.l<- melt(adat2.p.m, id.vars = c("NTmin"), variable.name = "metric")

ggplot(data=adat2.p.l, aes(x=NTmin, y =value))+geom_line()+facet_wrap(.~metric)

#----
#life table
ggplot(data=adat2.lt, aes(x=day, y =L.x., color=factor(NTmin), group=factor(NTmin) ))+geom_line()
ggplot(data=adat2.lt, aes(x=day, y =M.x., color=factor(NTmin), group=factor(NTmin) ))+geom_smooth()

#----
#survival
ggplot(data=adat2.sur, aes(x=NTmin, y =Days, color=factor(NTmin), group=factor(NTmin) ))+geom_point()+geom_smooth()+facet_wrap(.~Status)

#* Night warming compare estimated and observed

#=====
#Ma CS, Wang L, Zhang W, Rudolf V 2018. Resolving biological impacts of multiple heat waves: interaction of hot and recovery days. Oikos 127:622–33
#https://doi.org/10.1111/oik.04699
#data: http://dx.doi.org/10.5061/dryad.5qk4s
#Sitobion avenae

#vary number of successive hot days (1–3 days) and normal interval days (1–3 days)
#either exposed aphids first to hot days or normal days

#Finding: Increasing the duration of hot days in heat waves had a negative effect on various demographic rates and life-time fitness of individuals, but magnitude of this effect was typically contingent on the temporal clustering of hot periods.

#load biological data
adat3.dev<- read.csv("Maetal2017/Development.csv")
adat3.rep<- read.csv("Maetal2017/Reproduction.csv")
adat3.tr<- read.csv("Maetal2017/Traits.csv")

#----
#generate temps
#hot: 35°C and 20°C as the daily maximum and minimum temperatures
#normal: 28°C and 13°C as the daily maximum and minimum temperatures
#sine curve to simulate the diurnal fluctuations with a magnitude of 15°C 

hrs<- 0:23
plot(hrs, diurnal_temp_variation_sine(T_max = 35, T_min = 20, hrs))

#make temperature sequences
hdt= c(diurnal_temp_variation_sine(T_max = 35, T_min = 20, 4:23), diurnal_temp_variation_sine(T_max = 35, T_min = 20, 0:3))
ndt= c(diurnal_temp_variation_sine(T_max = 28, T_min = 13, 4:23), diurnal_temp_variation_sine(T_max = 28, T_min = 13, 0:3))
#fix transition

x<- temps[1,1:3]
tz= rep(c(rep(ndt, x[2]), rep(hdt, x[1])), ceiling(40/(x[1]+x[2])))
plot(1:100, tz[1:100], type="l")

#set up matrix of 30 days of hourly temps
trt<- expand.grid(hd=1:3, nd=1:3, first= 1:2) #first: 1 is n, 2 is h
temps<- matrix(NA, nrow=nrow(trt), ncol=24*30)
temps= as.data.frame(cbind(trt, temps))

build.temps<- function(x){
  x[1]<- as.numeric(x[1]); x[2]<- as.numeric(x[2])
if(x[3]==1) ts<- rep(c(rep(ndt, x[2]), rep(hdt, x[1])), ceiling(40/(x[1]+x[2])))
if(x[3]==2) ts<- rep(c(rep(hdt, x[1]), rep(ndt, x[2])), ceiling(40/(x[1]+x[2])))
  return(ts[1:720])
}

#set up temperatures  
tz= t(apply(temps[,1:3], MARGIN=1, FUN=build.temps))
temps[, 4:ncol(temps)]<- tz

#to long format
temps.l<- melt(temps, id.vars = c("hd","nd","first"), variable.name = "time")
temps.l$treat<- paste(temps.l$hd, temps.l$nd, sep="_")
temps.l$treat.nh<- paste(temps.l$hd, temps.l$nd, temps.l$first, sep="_")
temps.l$time<-as.numeric(temps.l$time)

#plot temperatures
plot(1:100, tz[1,1:100], type="l")
plot(1:200, temps[1,4:203], type="l")
ggplot(data=temps.l, aes(x=time, y =value, color=hd, group=treat.nh))+geom_line()+facet_wrap(.~treat)+xlim(0,100)

#----
#Traits
#to long format
adat3.l<- melt(adat3.tr[,c(2:7,26:29)], id.vars = c("Treatment","ID","H_C","CycleDays","ContinueNormalDay","ContinueHotday"), variable.name = "metric")
adat3.l$ContinueNormalDay <- as.factor(adat3.l$ContinueNormalDay)

#NymphDur, Longevity, Fecundity, BirthDur
fig.trait<- ggplot(data=adat3.l, aes(x=ContinueHotday, y =value, color=ContinueNormalDay, group=ContinueNormalDay))+
  geom_smooth(method='lm') +geom_point()+
  facet_grid(H_C~metric, scales="free_y")

#----
#Other data sets

#development:nymphal period
#find first zero 
dt= apply(adat3.dev[,5:ncol(adat3.dev)], MARGIN=1, FUN=function(x)max(which(x==1)) )
dt[is.infinite(dt)]<-NA
adat3.dev$dev<- dt

#mean metric
adat3.dev.m= adat3.dev %>%
  group_by(TRT, Hot_Cold) %>%
  summarise( dev= mean(dev, na.rm=TRUE) )

ggplot(data=adat3.dev, aes(x=TRT, y =dev, color=Hot_Cold))+geom_point()

ggplot(data=adat3.dev.m, aes(x=TRT, y =dev, color=Hot_Cold))+geom_point()

#----
#Reproduction
adat3.rep$rep<- colSums(adat3.rep[,5:ncol(adat3.rep)])

#mean metric
adat3.rep.m= adat3.rep %>%
  group_by(Treatment, Hot_Cold) %>%
  summarise( rep= mean(rep, na.rm=TRUE) )

ggplot(data=adat3.rep, aes(x=Treatment, y =rep, color=Hot_Cold))+geom_point()

ggplot(data=adat3.rep.m, aes(x=Treatment, y =rep, color=Hot_Cold))+geom_point()

#--------------------------
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

#estimate performance
tz= t(apply(temps[,4:720], MARGIN=1, FUN=perf.fun))
tp= cbind(temps[,1:3], tz)
colnames(tp)[4:7]<- c("sur", "Longevity", "Fecundity", "dr")
tp$NymphDur<- 1/tp$dr

#to long format
tp.l<- melt(tp, id.vars = c("hd","nd","first"), variable.name = "metric")

#--------------------------
#adapt Kingsolver Woods model

tfun<- function(Time) temps1[Time]

#estimate performance
for(treat.k in 1:nrow(temps) ){
  
  temps1<- as.numeric(temps[treat.k,4:723])
  times<- 1:720
  
  out2 <- ode(func = RsigCG, y = yini, parms = pars, times = times)
  
  #convert to dataframe
  out2.df <- data.frame(out2)
  colnames(out2.df)[5:7]<- c("T","I","G") 
  #summary(out2.df)
  
  out2.df$hd<- temps[treat.k,"hd"]
  out2.df$nd<- temps[treat.k,"nd"]
  out2.df$first<- temps[treat.k,"first"]
  
  if(treat.k==1) pout<- out2.df
  if(treat.k>1) pout<- rbind(pout, out2.df)
  
} #end loop treatments

#---------
#Process Kingsolver model estimates
#to long format
pout.l<- melt(pout, id.vars = c("time","T", "hd","nd","first"), variable.name = "metric")

#time series
ggplot(data=pout.l[which(pout.l$first==1),], aes(x=time, y =value, color=factor(nd)))+geom_line()+
  facet_grid(metric~hd, scale="free_y")

#sum performance
perf= pout.l %>%
  group_by(hd, nd, first, metric) %>%
  summarise(value= mean(value))
colnames(perf)[1:3]<- c("ContinueHotday","ContinueNormalDay", "H_C")
perf$H_C<- c("C","H")[perf$H_C]
perf$H_C<- factor(perf$H_C)
perf$ContinueNormalDay <- factor(perf$ContinueNormalDay)

#equate growth to reproduction
perf$metric<- as.character(perf$metric)
perf$metric[which(perf$metric=="G")]<-"Fecundity"
perf$metric<- as.factor(perf$metric)

#plot Kingsolver and Woods estimate
ggplot(perf, aes(x=ContinueHotday, y =value, color=factor(ContinueNormalDay)))+geom_line()+
  facet_grid(metric~H_C, scale="free_y")
#very small differences, scale?

#------------
#plot fitness components together

#combine dataframes
obs<- adat3.l[,c("H_C","ContinueNormalDay","ContinueHotday","metric","value")]
obs$type<- "observed"

tp.l2<- tp.l[tp.l$metric %in% c("Longevity", "Fecundity","NymphDur"),]
colnames(tp.l2)[1:3]<- c("ContinueHotday","ContinueNormalDay", "H_C")
tp.l2$type<- "estimated"
tp.l2$H_C<- c("C","H")[tp.l2$H_C]

est.kw<- perf[which(perf$metric=="Fecundity"),]
est.kw$type<- "king est"
  
fdat<- rbind(obs, tp.l2[,c("H_C","ContinueNormalDay","ContinueHotday","metric","value", "type")],
             est.kw[,c("H_C","ContinueNormalDay","ContinueHotday","metric","value", "type")] )

#Plot NymphDur, Longevity, Fecundity, BirthDur
fplot<- ggplot(data=fdat, aes(x=ContinueHotday, y =value, color=ContinueNormalDay, lty=type))+
  geom_smooth(method='lm') +geom_point()+
  facet_grid(metric~H_C, scales="free_y")

#nymphal duration: more hot days increase dt but predicted to reduce
#longevity: more hot days decrease longevity more than estimated when few normal days; less hot days decreases longevity less than expected when more normal days
#fecundity: more hot days decrease fecundity more than estimated when few normal days; less hot days decreases longevity less than expected when more normal days
#BirthDir

#=====
#Wang, XJ., Ma, CS. Can laboratory-reared aphid populations reflect the thermal performance of field populations in studies on pest science and climate change biology?. J Pest Sci 96, 509–522 (2023). https://doi.org/10.1007/s10340-022-01565-6
#Sitobion avenae

#read data
setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/data/aphids/WangMa2023/")
#setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/ThermalHistory/data/aphids/WangMa2023/")

#mild means
adat4.var<- read.csv("WangMa2023_temp22mean_diffvar.csv")
#high means
adat4.mean<- read.csv("WangMa2023_diffmeans.csv")
#popgrowth
adat4.r<- read.csv("WangMa2023_popgrowth.csv")

#----------
#set up temperatures
#mean: 22 °C; ±0 °C,±5 °C,±9 °C,±13 °C as fluctuating ranges, and designed four #3 high constant (28 °C, 30 °C, 32 °C) and fluctuating regimes (28 ± 5 °C, 30±5 °C, 32±5 °C)

temps= matrix(NA, nrow= 10, ncol=24)
temps[1,]<- 22
temps[2,]<- diurnal_temp_variation_sine(T_max = 27, T_min = 17, 1:24)
temps[3,]<- diurnal_temp_variation_sine(T_max = 31, T_min = 13, 1:24)
temps[4,]<- diurnal_temp_variation_sine(T_max = 35, T_min = 9, 1:24)
temps[5,]<- 28
temps[6,]<- 30
temps[7,]<- 32
#28 ± 5
temps[8,]<- diurnal_temp_variation_sine(T_max = 33, T_min = 23, 1:24)
#30±5
temps[9,]<- diurnal_temp_variation_sine(T_max = 35, T_min = 25, 1:24)
#32±5
temps[10,]<- diurnal_temp_variation_sine(T_max = 37, T_min = 27, 1:24)

temps<- as.data.frame(temps)
colnames(temps)= 1:24

##expand to 50 days
#temps<- as.data.frame(rep(temps, each = 50))

Tmean= c(22, 22, 22, 22, 28, 30, 32, 28, 30, 32)
Tvar= c(0, 5, 9, 13, 0, 0, 0, 5, 5, 5)
treat= c("mild means", "mild means", "mild means", "mild means", "high means", "high means", "high means", "high means", "high means", "high means")
temps= as.data.frame(cbind(Tmean, Tvar, treat, temps))
temps$Tmean_var= paste(temps$Tmean, temps$Tvar, sep="_")

#to long format
temps.l<- melt(temps, id.vars = c("Tmean", "Tvar", "treat", "Tmean_var"), variable.name = "time")
temps.l$time= as.numeric(temps.l$time)
temps.l$value= as.numeric(temps.l$value)

ggplot(data=temps.l, aes(x=time, y =value, color=Tvar, group=Tmean_var))+geom_line()+
  facet_grid(.~treat, scale="free_y")

#----------
#Tmean of 22, difference variance
#sum performance
adat4.var.m= adat4.var %>%
  group_by(Tmean, Tvar, population, metric) %>%
  summarise(value= mean(value))

ggplot(data=adat4.var.m, aes(x=Tvar, y =value, color=population))+geom_point()+ geom_line()+
  facet_grid(metric~Tmean, scale="free_y")

#Different means, redundant with adat4.var
#sum performance
adat4.mean.m= adat4.mean %>%
  group_by(Tmean, Tvar, population, metric) %>%
  summarise(value= mean(value))

ggplot(data=adat4.mean.m, aes(x=Tvar, y =value, color=population))+geom_point()+ geom_line()+
  facet_grid(metric~Tmean, scale="free_y")

#--------------------------
#estimate performance using constant TPC

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

#estimate performance
tz= t(apply(temps[,4:(ncol(temps)-1)], MARGIN=1, FUN=perf.fun))
tp= cbind(temps[,1:3], tz)
colnames(tp)[4:7]<- c("sur", "Longevity", "Fecundity", "dr")
tp$NymphDur<- 1/tp$dr

#to long format
tp.l<- melt(tp, id.vars = c("Tmean","Tvar","treat"), variable.name = "metric")

#plot estimates
ggplot(data=tp.l, aes(x=Tmean, y=value, color=Tvar)) +geom_point()+
  geom_smooth(method='lm') +geom_point()+
  facet_grid(metric~treat, scales="free_y")

#---------------
#plot estimates vs observed

adat4.var.m$type<- "observed"
adat4.mean.m$type<- "observed"
#align metric names
#dat.mv<- rbind(adat4.var.m, adat4.mean.m)
dat.mv<- adat4.var.m
## UPDATE


#add popgrowth
adat4.r[,c("H_C","ContinueNormalDay","ContinueHotday","metric","value")]

Tmean  Tvar population metric     value type     treat 

#----
dat.mv$treat<- "mild means"
dat.mv$treat[which(dat.mv$Tmean>22)]<- "high means"

dat.mv$metric[dat.mv$metric=="dev_rate"]<-"dr"
dat.mv$metric[dat.mv$metric=="fecundity"]<-"Fecundity"
dat.mv$metric[dat.mv$metric=="survival"]<-"sur"

tp.l$population<- NA
tp.l$type<- "estimated"

dat.mv<- rbind(dat.mv, 
               tp.l[,c("Tmean","Tvar","population","metric","value","treat","type")] )
#drop NA
dat.mv<- dat.mv[!is.na(dat.mv$Tmean),]

#plot
ggplot(data=dat.mv, aes(x=Tvar, y =value, color=population, lty=type))+geom_point()+ geom_line()+  facet_grid(metric~Tmean, scale="free_y")

#model with sensitivity to extremes

#=========================================

#Ma et al. 2015. Daily temperature extremes play an important role in predicting thermal effects. The Journal of Experimental Biology 218 (14), 2289-2296
#https://doi.org/10.1242/jeb.122127, no data in paper

#=====
#Zhao et al. The importance of timing of heat events for predicting the dynamics of aphid pest populations. Pest management science, 2019
#https://doi.org/10.1002/ps.5344
#survival and productivity
#no data online
#Sitobion avenae 

