
#Load data
#Ma et al. 2021. Are extreme high temperatures at low or high latitudes more likely to inhibit the population growth of a globally distributed aphid?
#https://doi.org/10.1016/j.jtherbio.2021.102936

#2.1 Temp regimes Fig A2
#During the day, temperature started to increase at 08:00 h, reached and stayed at a high level (35, 37 or 39 ◦C) from 12:00 to 13:00 h, and then decreased to 22 ◦C by 16:00 h. We kept temperature constant at 22 ◦ C for the rest of the day. Importantly, we set both the daily mean (25.2 ◦C) and minimum (22 ◦C) temperatures as fixed

temps= read.csv("./data/temps_Maetal2021JTB.csv")

#construct temperatures

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

ggplot(data=pall, aes(x=t, y =p))+geom_line()+
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


