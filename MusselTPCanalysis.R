#load libraries
library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)
library(viridisLite)
library(patchwork)

library(ggplot2)

#--------------------
#load variable TPCs
#https://doi.org/10.1111/1365-2435.13889
#feeding rate, Fig 2d
#shell length growth, Fig 2a
#Mytilus edulis is the dominant species, small fractions of M. galloprovincialis and M. trossulus

#data
#https://doi.pangaea.de/10.1594/PANGAEA.933828
#https://doi.pangaea.de/10.1594/PANGAEA.897938

#load TPC data
st= read.csv("./data/Shortterm_experiment.csv")
lt= read.csv("./data/Longterm_growth.csv")

#load temperature data
temps= read.csv("./data/Longterm_temps.csv")
temps$doy= as.numeric(format(as.POSIXlt(temps$datetime,format="%m/%d/%y %H:%M"),"%j"))
hour = as.numeric(format(as.POSIXlt(temps$datetime,format="%m/%d/%y %H:%M"),"%H"))
minute = as.numeric(format(as.POSIXlt(temps$datetime,format="%m/%d/%y %H:%M"),"%M"))
sec= as.numeric(format(as.POSIXlt(temps$datetime,format="%m/%d/%y %H:%M"),"%S"))
temps$time= hour + minute/60 + sec/3600

#process time
doy= as.numeric(format(as.POSIXlt(st$DateTime,format="%Y-%m-%d %H:%M:%S"),"%j"))
hour = as.numeric(format(as.POSIXlt(st$DateTime,format="%Y-%m-%d %H:%M:%S"),"%H"))
minute = as.numeric(format(as.POSIXlt(st$DateTime,format="%Y-%m-%d %H:%M:%S"),"%M"))
sec= as.numeric(format(as.POSIXlt(st$DateTime,format="%Y-%m-%d %H:%M:%S"),"%S"))
st$time= hour + minute/60 + sec/3600

#distinguish warming and cooling phase
st$temp.phase<- "warming"
st$temp.phase[st$time>=17.5 | st$time<=5] <- "cooling"
st$temp.phase<- factor(st$temp.phase, levels=c("warming","cooling") )
#add 24hr to cooling
st$time[st$time<=5]= st$time[st$time<=5] +24

#short term by temp
st.feed.temp<- ggplot(data=st, aes(x=Temp_C_x, y =WS_feed_J_per_h_S, color=replicate, lty=trial_name))+
  geom_line()+ theme_classic()+
  facet_wrap(~temp.phase, scales="free_x") +theme_classic(base_size = 20)

st.resp.temp<- ggplot(data=st, aes(x=Temp_C_x, y =WS_resp_J_per_h_S, color=replicate, lty=trial_name))+
  geom_point()+ theme_classic()+
  facet_wrap(~temp.phase, scales="free_x") +theme_classic(base_size = 20)

#short term by time
#Feeding rate
# Value used to transform the data
coeff <- 3

st.feed.time<- ggplot(st, aes(x=time, color=replicate, lty=trial_name)) +
  geom_point( aes(y=WS_feed_J_per_h_S)) + 
  geom_point( aes(y=Temp_C_x*coeff)) + # Divide by 10 to get the same range than the temperature
  scale_y_continuous(
    # Features of the first axis
    name = "Feeding rate (J g^-1 h^-1)",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./coeff, name="Temperature")
  )+
  facet_wrap(~temp.phase, scales="free_x") +theme_classic(base_size = 20)

#Respiration rate
st.resp.time<- ggplot(st, aes(x=time, color=replicate, lty=trial_name)) +
  geom_point( aes(y=WS_resp_J_per_h_S)) + 
  geom_point( aes(y=Temp_C_x)) + 
  scale_y_continuous(
    # Features of the first axis
    name = "Respiration rate (J g^-1 h^-1)",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~., name="Temperature")
  )+
  facet_wrap(~temp.phase, scales="free_x") +theme_classic(base_size = 20)

#------------ 
#long term
lt.mean= lt %>%
  group_by(Mean.temperature...C.,Fluctuation.scenario) %>%
  summarise(mean=mean(Length.growth..mm.day.), sd=sd(Length.growth..mm.day.), n=n(), se=sd/sqrt(n) )

lt.fig<- ggplot(data=lt.mean, aes(x=Mean.temperature...C., y =mean, color=Fluctuation.scenario))+
  geom_point()+ geom_line(alpha=0.8, lwd=1) +theme_classic(base_size = 20)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1) +
  ylab("shell length growth (mm/day)")+
  xlab("thermal average (°C)")+
  scale_color_brewer("fluctuation", palette = "Dark2")+
  theme(legend.position = c(0.3, 0.3))
#change to GAM

#--------
#integrate static TPC through time to estimate performance

#process time
doy= as.numeric(format(as.POSIXlt(temps$datetime,format="%m/%d/%y %k:%M"),"%j"))
hour = as.numeric(format(as.POSIXlt(as.character(temps$datetime),format="%m/%d/%y %k:%M"),"%H"))
minute = as.numeric(format(as.POSIXlt(temps$datetime,format="%m/%d/%y %k:%M"),"%M"))
sec= as.numeric(format(as.POSIXlt(temps$datetime,format="%m/%d/%y %k:%M"),"%S"))
temps$time= doy + hour/24 + minute/(24*60) + sec/(24*60*60)

#plot temps
temp.fig<- ggplot(data=temps, aes(x=time, y =Temperature...C., color=Thermal_fluctuation_levels))+
  geom_line()+ theme_classic()+
  facet_grid(Thermal_mean_levels~.) +theme_classic(base_size = 20)

#fit tpc
tpc.fig<- ggplot(data=lt[lt$fluctuation==0,], aes(x=Mean.temperature...C., y =Length.growth..mm.day.))+
  geom_point()+ 
  theme_classic(base_size = 20)

#-----------
#fit tpc
#run MusselTPCfit.R to fit TPCs

tpcs.weibull

#----------
# #estimate performance
# temps$length.growth.static= weibull_1995(temps$Temperature...C., tpcs.weibull[1,1], tpcs.weibull[1,2], tpcs.weibull[1,3], tpcs.weibull[1,4])
# temps$shell.growth.static= weibull_1995(temps$Temperature...C., tpcs.weibull[2,1], tpcs.weibull[2,2], tpcs.weibull[2,3], tpcs.weibull[2,4])
# temps$tissue.growth.static= weibull_1995(temps$Temperature...C., tpcs.weibull[3,1], tpcs.weibull[3,2], tpcs.weibull[3,3], tpcs.weibull[3,4])

##try Bayne data
tpc.beta= c(198.39, 20.00, 32.00, 4.00, 4.00) #2nd parameter is position, 3rd parameter is breadth
temps$length.growth.static= beta_2012(temps$Temperature...C., tpc.beta[1], tpc.beta[2], tpc.beta[3], tpc.beta[4], tpc.beta[5])
temps$shell.growth.static= temps$length.growth.static
temps$tissue.growth.static= temps$length.growth.static

temps$temp <- temps$Temperature...C.
temps$perf= beta_2012(temps$temp, tpc.beta[1], tpc.beta[2], tpc.beta[3], tpc.beta[4], tpc.beta[5])

#configure TPC
wp<- c(0.15, 22, .48*10^5, 2.1*10^4)
plot(10:30, weibull_1995(10:30, wp[1], wp[2], wp[3], wp[4]), type="l")
points(c(18.5,21, 23.5,26), weibull_1995(c(18.5,21, 23.5,26), wp[1], wp[2], wp[3], wp[4]))
points(seq(18.5-4, 18.5+4,1), weibull_1995(seq(18.5-4, 18.5+4,1), wp[1], wp[2], wp[3], wp[4]), type="l",col="green")
points(seq(26-4, 26+4,1), weibull_1995(seq(26-4, 26+4,1), wp[1], wp[2], wp[3], wp[4]), type="l",col="blue")

#a: height; topt:optimal temperature; b:breadth; c:curve shape
temps$perf= weibull_1995(temps$temp, wp[1], wp[2], wp[3], wp[4])

#----------
#plot time series
#temp
t.plot<- ggplot(data=temps[temps$doy=="300",], aes(x=time, y=temp, color=factor(Thermal_fluctuation_levels), lty=factor(Thermal_mean_levels)))+
  geom_line() +theme_classic(base_size = 20)+
  labs(col="fluct", lty="T mean") +theme(legend.position ="none")+
  facet_wrap(.~factor(Thermal_mean_levels), ncol=4) 

#perf
p.plot<- ggplot(data=temps[temps$doy=="300",], aes(x=time, y=perf, color=factor(Thermal_fluctuation_levels), lty=factor(Thermal_mean_levels)))+
  geom_line() +theme_classic(base_size = 20)+
  labs(col="fluct", lty="T mean") +theme(legend.position ="bottom")+
  facet_wrap(.~factor(Thermal_mean_levels), ncol=4)

#setwd for figures
#setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/ThermalHistory/figures/")
setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/figures/")

pdf("FigTPCs_temps.pdf", height = 10, width = 10)
t.plot / p.plot
dev.off()

#long term estimates
lt.est= temps %>%
  group_by(Thermal_mean_levels,Thermal_fluctuation_levels) %>%
  summarise(perf=sum(perf) )

#plot estimates
ggplot(data=lt.est, aes(x=Thermal_mean_levels, y=perf, color=factor(Thermal_fluctuation_levels)))+
  geom_point()+ geom_line(alpha=0.8, lwd=1) +theme_classic(base_size = 20)+theme(legend.position ="bottom")

#==========================
#estimate performance accounting for damage and recovery
#exponential decline at hot temperatures, short term stress happen around 26
#recovery time

#assume stress beyond Topt
#accelerates quadratically
#multiplicative with duration
tcost<- function(Tb, dur, rep) ifelse(Tb < 24, 0, min(.1 * ((Tb - 22)*dur)^(1+0.1*rep),1))
tcost.vect<- function(x) ifelse(x[1] < 24, 0, min(.3 * ((x[1] - 22)*x[2])^(1.2+0.2*x[3]),1))

#cost decays linearly                                
decay <- function(time.int) max(1.0 - time.int*0.2)

#------------------------  
#identify heat stress events
temps$hs<- ifelse(temps$temp <= 24, 0, 1)

#find lengths of heat stress
rles<- rle(temps$hs)
rs<- as.data.frame(cbind(lengths=rles$lengths, values=rles$values, cumsum=cumsum(rles$lengths))) 
rs$inds<- rs$cumsum -rs$lengths +1

#restrict to heat waves
rs<- rs[rs$values==1,]
#heat wave count
rs$rep<- 1:nrow(rs)
#mean heat wave temperature
rs$tmean<-0
for(hs.ind in 1:(nrow(rs))){
  rs$tmean[hs.ind]<- mean(temps$temp[rs$inds[hs.ind]:(rs$inds[hs.ind]+rs$lengths[hs.ind])])
}

#add durations and number
temps$dur<- 0
temps$dur[rs$inds]<- rs$lengths
temps$rep<- 0
temps$rep[rs$inds]<- rs$rep
temps$tmean<- 0
temps$tmean[rs$inds]<- rs$tmean

#cost 
temps$cost<- apply(cbind(temps$tmean, temps$dur, temps$rep), MARGIN=1, FUN=tcost.vect)

#decay
temps$costdec<- 0

inds=which(temps$dur>=1)
for(i in 1:length(inds)){
  dec= temps$cost[inds[i]]*c(1.0, 0.8, 0.6, 0.4, 0.2)
  ncary<- min(5,nrow(temps)-inds[i])
  temps$costdec[inds[i]:(inds[i]+ncary)]<- temps$costdec[inds[i]:(inds[i]+ncary)]+dec
}

#net performance
temps$costdec[temps$costdec>1]<-1
temps$netempserf= temps$perf*(1-temps$costdec)

#=======================
#long term estimates
lt.est= temps %>%
  group_by(Thermal_mean_levels,Thermal_fluctuation_levels) %>%
  summarise(length.growth=sum(length.growth.static), shell.growth=sum(shell.growth.static), tissue.growth=sum(tissue.growth.static) )

#to long format
lt.est.l<- lt.est %>%
  gather("metric", "value", 3:ncol(lt.est))

#plot estimates
lt.est.fig<- ggplot(data=lt.est.l, aes(x=Thermal_mean_levels, y =value, color=factor(Thermal_fluctuation_levels)))+
  geom_point()+ geom_line(alpha=0.8, lwd=1) +theme_classic(base_size = 20)+
  ylab("estimated growth (mm/day)")+
  xlab("thermal average (°C)")+
  scale_color_brewer("fluctuation", palette = "Dark2")+
  theme(legend.position = c(0.1, 0.2))+
  facet_wrap(.~metric, scale="free_y")

#--------
#plot with cost
#long term estimates
lt.est= temps %>%
  group_by(Thermal_mean_levels,Thermal_fluctuation_levels) %>%
  summarise(perf=sum(perf), perf.c=sum(netempserf) )

#to long format
lt.est.l<- lt.est %>%
  gather("metric", "value", 3:ncol(lt.est))

#plot estimates
cost.est.fig<- ggplot(data=lt.est.l, aes(x=Thermal_mean_levels, y =value, color=factor(Thermal_fluctuation_levels), lty=metric))+
  geom_point()+ geom_line(alpha=0.8, lwd=1) +theme_classic(base_size = 20)+
  ylab("estimated growth (mm/day)")+
  xlab("thermal average (°C)")+
  scale_color_brewer("fluctuation", palette = "Dark2")+
  theme(legend.position = c(0.1, 0.2))

#=====================
#observed long term plots

#long term
lt.mean= lt %>%
  group_by(Mean.temperature...C.,fluctuation) %>%
  summarise(length.mean=mean(Length.growth..mm.day.), length.sd=sd(Length.growth..mm.day.),
            shell.mean=mean(Shell.dry.weight.growth..mg.day.), shell.sd=sd(Shell.dry.weight.growth..mg.day.),
            tissue.mean=mean(Tissue.dry.weight.growth..mg.day.), tissue.sd=sd(Tissue.dry.weight.growth..mg.day.),
            n=length(Length.growth..mm.day.), 
            length.se=length.sd/sqrt(n), shell.se=shell.sd/sqrt(n), tissue.se=tissue.sd/sqrt(n) )

#to long format
lt.l<- lt.mean %>%
  gather("metric", "value", c(3,5,7) )

#plot observed
lt.obs.fig<- ggplot(data=lt.l, aes(x=Mean.temperature...C., y =value, color=factor(fluctuation)))+
  geom_point()+ geom_line(alpha=0.8, lwd=1) +theme_classic(base_size = 20)+
  ylab("observed growth (mm/day)")+
  xlab("thermal average (°C)")+
  scale_color_brewer("fluctuation", palette = "Dark2")+
  theme(legend.position = c(0.1, 0.2))+
  facet_wrap(.~metric, scale="free_y")

#-----
#plot observed and estimated

#setwd for figures
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/ThermalHistory/figures/")

#pdf("FigTPCs.pdf", height = 8, width = 10)
#lt.obs.fig / lt.est.fig
#dev.off()

#differences
# estimate higher at low temps than observed (benefit of warm temperatures)
# estimate at 24C higher than observed (stress)














