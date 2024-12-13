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

#===================
#Fig 1. TPCs, model repair, damage

#===================
#Fig 2. model time series

#===================
#Fig 3. Temperature schematic
#Expt 1: vary min

#-------------
#Expt 2: vary max

#-------------
#Expt 3: vary variance (expt 3, mild means)


#-------------
#Expt 4: vary means and variance (expt 3, high means)

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
ggplot(temps.l, aes(x = day, y = treat.nh, fill=factor(value))) +
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

ggplot(t6.l, aes(x = day, y = treatment, fill=factor(value))) +
  geom_tile()

#===================
#Fig 4. Performance with and without damage

#===================
#Figs 5. Fecundity comparisons

#===================
#Fig 6. Development comparison


