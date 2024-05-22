library(ggplot2)

setwd("/Users/laurenbuckley/ThermalHistory")

#Load data
#Ma et al. 2021. Are extreme high temperatures at low or high latitudes more likely to inhibit the population growth of a globally distributed aphid?
#https://doi.org/10.1016/j.jtherbio.2021.102936

#2.1 Temp regimes Fig A2
#During the day, temperature started to increase at 08:00 h, reached and stayed at a high level (35, 37 or 39 ◦C) from 12:00 to 13:00 h, and then decreased to 22 ◦C by 16:00 h. We kept temperature constant at 22 ◦ C for the rest of the day. Importantly, we set both the daily mean (25.2 ◦C) and minimum (22 ◦C) temperatures as fixed

#construct temperatures
temps= read.csv("./data/temps_Maetal2021JTB.csv")

#add other hours
tadd<- cbind(time=c(1:8,16:24),temp=22)
tadd<- rbind( cbind(rep(39, 17),tadd), cbind(rep(37, 17),tadd), cbind(rep(35, 17),tadd) ) 
colnames(tadd)[1]<- "treatment"
temps<- as.data.frame(rbind(temps, tadd))
temps<- temps[order(temps[,1], temps[,2]),]

days<- rep(1:20, nrow(temps))
days<- days[order(days)]
temps.all<- cbind(temps, days)
temps.all$dt<- temps.all$days+temps.all$time/24

#plot temps
ggplot(data=temps.all, aes(x=dt, y=temp, color=factor(treatment)))+geom_line()

#--------------
#2.4 constant temp TPCs
#Params from Table 1, fitting data in appendix

#nymphal survival
Sur= function(T, a=334.1, b= -0.29) 1-T/(a-b*T^2)

#Adult longevity (days)
Long= function(T, a=29.5, b= -0.81) a + b*T

#Lifetime fecundity (nymphs / adult)
Fecun= function(T, a= -0.3, b= 11.1, c= -50.9) a*T^2 + b*T +c

#Reproductive rate (nymphs/adult/day)
Rr= function(T, a=0.263, b=3.7, T0 =30.1) exp(a*T)-exp(a*T0-(T0-T)/b)

#Developmental rate (1/days)
Dr= function(T, a=0.149, b=6.6, T0 =33.1) exp(a*T)-exp(a*T0-(T0-T)/b)

#Intrinsic rate of increase
Rm= function(T, a=0.177, b=5.6, T0 =30.2) exp(a*T)-exp(a*T0-(T0-T)/b)    

#Plot
p1= cbind(t=1:40, p=Sur(1:40), comp="Sur")
p2= cbind(t=1:40, p=Long(1:40), comp="Long")
p3= cbind(t=1:40, p=Fecun(1:40), comp="Fecun")
p4= cbind(t=1:40, p=Rr(1:40), comp="Rr")
p5= cbind(t=1:40, p=Dr(1:40), comp="Dr")
p6= cbind(t=1:40, p=Rm(1:40), comp="Rm")
pall<- as.data.frame(rbind(p1, p2, p3, p4, p5, p6))
pall$t<- as.numeric(pall$t); pall$p<- as.numeric(pall$p)

fig1<- ggplot(data=pall, aes(x=t, y =p))+geom_line()+
  facet_wrap(.~comp, scales="free_y")+
  scale_y_continuous(limits = c(0, NA))

#need to fix fecundity
fec= cbind(temp=c(10,15,20,25,18,22,25,27.5,30),
           fec=c(30,55,60,41,39.5,35.2,42.4,45,1.6))
fec= as.data.frame(fec)
fec$temp2<- fec$temp^2

ggplot(data=fec, aes(x=temp, y =fec))+geom_point()
lm1<- lm(fec~temp + temp2 , data=fec)
#c= -50.9, a=11.1, b= -0.3

#--------------------------
#estimate performance using constant TPC

#nymphal survival
temps.all$sur<- Sur(temps.all$temp)

#Adult longevity (days)
temps.all$long<- Long(temps.all$temp)

#Lifetime fecundity (nymphs / adult)
temps.all$fecun<- Fecun(temps.all$temp)

#Reproductive rate (nymphs/adult/day)
temps.all$rr<- Rr(temps.all$temp)

#Developmental rate (1/days)
temps.all$dr<- Dr(temps.all$temp)

#Intrinsic rate of increase
temps.all$r<- Rm(temps.all$temp)

#set values less than 0 to 0
temps.all$rr[temps.all$rr<0] <- 0
temps.all$dr[temps.all$dr<0] <- 0
temps.all$r[temps.all$r<0] <- 0

#plot constant rates
#sum performance
perf= temps.all %>%
  group_by(treatment) %>%
  summarise(Sur=mean(sur), Long=mean(long), Fecun=mean(fecun), Rr=mean(rr), Dr=mean(dr), Rm=mean(r) )
            
#to long format
perf<- perf %>%
  gather("comp", "p", c(2:ncol(perf)) )
perf$t<- 25 #mean temp

#plot
fig1+  geom_point(data=perf, aes(color=factor(treatment)))

#--------------------------
#adapt Kingsolver Woods model

tfun<- function(Time) temps[Time]

RsigCG = function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    T = tfun(Time)   		## temperature determined from sine wave function
    Rf = Rm2*sigmoid(T,a=0.5,b=Tc)			## calculate equilibrium level of RNA
    Pf = Pm*sigmoid(R,a=1,b=Rc)				## calculate equilibrium level of protein
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

#loop temperatures
treats<- c(35,37,39)

for(treat.k in 1:3){

temps<- c(temps.all[temps.all$treatment==treats[treat.k],"temp"])
times<- 1:length(temps)  #c(temps.all[temps.all$treatment==35,"dt"])

out2 <- ode(func = RsigCG, y = yini, parms = pars, times = times)

#convert to dataframe
out2.df <- data.frame(out2)
colnames(out2.df)[5:7]<- c("T","I","G") 
#summary(out2.df)
out2.df$treatment<- treats[treat.k]

if(treat.k==1) pout<- out2.df
if(treat.k>1) pout<- rbind(pout, out2.df)

} #end loop treatments

#---------
#plot Kingsolver model

#to long format
pout.l<- melt(pout, id.vars = c("time","treatment","T"), variable.name = "trait")

ggplot(data=pout.l, aes(x=time, y =value, color=factor(treatment)))+geom_line()+
  facet_wrap(.~trait, scale="free_y")

#plot Kingsolver model
#sum performance
perf= pout.l %>%
  group_by(treatment) %>%
  summarise(value= mean(value))

#=================================
#Wang and Ma. 2023. Can laboratory‐reared aphid populations reflect the thermal performance of field populations in studies on pest science and climate change biology?  JPS
#https://doi.org/10.1007/s10340-022-01565-6 

#check data formatting
setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/ThermalHistory/data/aphids/")
adat1<- read.csv("WangMa2023_temp22mean.csv")

#mean metric
adat1.mean= adat1 %>%
  group_by(Tmean, Tvar, metric, population) %>%
  summarise(value= mean(value))

ggplot(data=adat1.mean, aes(x=Tvar, y =value, color=population))+geom_point()+
  facet_grid(Tmean~metric, scale="free_y")

#Ma et al. 2015. Daily temperature extremes play an important role in predicting thermal effects. The Journal of Experimental Biology 218 (14), 2289-2296
#https://doi.org/10.1242/jeb.122127, no data in paper

#Zhao et al. 2014. Night warming on hot days produces novel impacts on development, survival and reproduction in a small arthropod
#Dryad data: http://doi.org/10.5061/dryad.q2070 

adat2.dt<- read.csv("Zhaoetal2014_devtime.csv")
adat2.p<- read.csv("Zhaoetal2014_AdPerf.csv")
adat2.lt<- read.csv("Zhaoetal2014_LifeTable.csv")
adat2.sur<- read.csv("Zhaoetal2014_SurvNymph.csv")

names(adat2.dt.l)[1,3,5,7,9]<-"NTmin"
names(adat2.dt.l)[2,4,5,7,10]<-"dt"
adat2.dt.l<- rbind( cbind(adat2.dt[,1:2],"1st"), cbind(adat2.dt[,1:2],"2nd"), cbind(adat2.dt[,1:2],"3rd"),
                    cbind(adat2.dt[,1:2],"4th"),cbind(adat2.dt[,1:2],"Nymph") )
  
  
  

#Zhao et al. The importance of timing of heat events for predicting the dynamics of aphid pest populations. Pest management science, 2019

#Ma CS, Wang L, Zhang W, Rudolf V 2018. Resolving biological impacts of multiple heat waves: interaction of hot and recovery days. Oikos 127:622–33
#https://doi.org/10.1111/oik.04699

#Zhang W, Chang XQ, Hoffmann AA, Zhang S and Ma CS, Impact of hot events at different developmental stages of a moth: the closer to adult stage, the less reproductive output. Sci Rep 5: 1–9 (2015).



