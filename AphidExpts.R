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

#=====
#Expt 1
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

temps= as.data.frame(cbind(tmin, do.call("cbind", rep(list(temps), 30))))
colnames(temps)[2:ncol(temps)]<- 1:(ncol(temps)-1)

#to long format
temps.l<- melt(temps, id.vars = c("tmin"), variable.name = "time")

#plot temperatures, varying night temperatures
ggplot(data=temps.l, aes(x=time, y =value, color=tmin, group=tmin))+geom_line()

#combine temp data
colnames(temps.l)<- c("treatment","time","temp")
temps.l$expt<- 1

temps.all<- temps.l
temps.all$time<- as.numeric(temps.all$time)

#-----
#Assemble performance

#mean metric
adat2.p.m= adat2.p %>%
  group_by(NTmin) %>%
  summarise(lon= mean(Longevity, na.rm = TRUE), fec=mean(Fecundtiy, na.rm = TRUE), fec.rate=mean(Fecundity.rate, na.rm = TRUE) )

#to long format
adat2.p.l<- melt(adat2.p.m, id.vars = c("NTmin"), variable.name = "metric")

ggplot(data=adat2.p.l, aes(x=NTmin, y =value))+geom_line()+facet_wrap(.~metric)

#combine dataframes
obs<- adat2.p.l[,c("NTmin","metric","value")]
colnames(obs)[1]<- c("treatment")
obs$expt<- 1

PerfDat<- obs

#=========================================
#Expt 2
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

#combine temp data
temps.l<- temps.l[,c("treat.nh","time","value")]
colnames(temps.l)<- c("treatment","time","temp")
temps.l$expt<- 2

temps.all<- rbind(temps.all, temps.l)

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
#add performance data
obs<- adat3.l[,c("H_C","ContinueNormalDay","ContinueHotday","metric","value")]

obs$hc<- match(obs$H_C, c("C","H") )
obs$treatment<- paste(obs$ContinueHotday, obs$ContinueNormalDay, obs$hc, sep="_")
obs$expt<- 2

PerfDat<- rbind(PerfDat, obs[,c("treatment","metric","value","expt")])

#=====
#Expt 3
#Wang, XJ., Ma, CS. Can laboratory-reared aphid populations reflect the thermal performance of field populations in studies on pest science and climate change biology?. J Pest Sci 96, 509–522 (2023). https://doi.org/10.1007/s10340-022-01565-6
#Sitobion avenae

#read data
setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/data/aphids/WangMa2023/")
setwd("/Users/lbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/data/aphids/WangMa2023/")
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

##expand to 30 days
temps= do.call("cbind", rep(list(temps), 30))
colnames(temps)[2:ncol(temps)]<- 1:(ncol(temps)-1)

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

#combine temp data
temps.l<- temps.l[,c("Tmean_var","time","value")]
colnames(temps.l)<- c("treatment","time","temp")
temps.l$expt<- 3

temps.all<- rbind(temps.all, temps.l)

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

#---------------
#combine performance metrics

#to long format
r.l<- melt(adat4.r[,c("Tmean","Tvar","population","mean_Ro","mean_rm")], id.vars = c("Tmean","Tvar","population"), variable.name = "metric")

#add popgrowth
dat.mv<- rbind(adat4.var.m, r.l)

dat.mv$treatment= paste(dat.mv$Tmean, dat.mv$Tvar, sep="_")

#add performance data
obs<- dat.mv[,c("treatment","metric","value")]
obs$expt<- 3

#add field to distinguish lab and field
PerfDat$population<- NA

PerfDat<- rbind(PerfDat, obs[,c("treatment","metric","value","expt","population")])

#=========================================

#Ma et al. 2015. Daily temperature extremes play an important role in predicting thermal effects. The Journal of Experimental Biology 218 (14), 2289-2296
#https://doi.org/10.1242/jeb.122127, no data in paper

#=====
#Zhao et al. The importance of timing of heat events for predicting the dynamics of aphid pest populations. Pest management science, 2019
#https://doi.org/10.1002/ps.5344
#survival and productivity
#no data online
#Sitobion avenae 

#====================
#write out data sets

setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/out/")
#setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/ThermalHistory/out/")

#align names
PerfDat$metric[which(PerfDat$metric=="fec")]="fecundity"
PerfDat$metric[which(PerfDat$metric=="Fecundity")]="fecundity"

PerfDat$metric[which(PerfDat$metric=="lon")]="Longevity"

write.csv(temps.all, "TempTimeSeries.csv")
write.csv(PerfDat, "PerformanceData.csv")

#plot temperatures for seven days
ggplot(data=temps.all[which(temps.all$time<169),], aes(x=time, y =temp, color=treatment))+geom_line()+facet_wrap(.~expt)

#plot fecundity
ggplot(data=PerfDat[which(PerfDat$metric=="fecundity"),], aes(x=treatment, y =value))+geom_point()+facet_wrap(.~expt)



