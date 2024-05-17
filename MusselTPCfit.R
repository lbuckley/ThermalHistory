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
lt= read.csv("./data/Longterm_growth.csv")

#fit tpc
#rTPC, https://github.com/padpadpadpad/rTPC
#https://padpadpadpad.github.io/rTPC/articles/fit_many_models.html

# load packages
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)

#loop through metrics
for(met in 1:3){

if(met==1) d= lt[lt$fluctuation==0,c("Mean.temperature...C.","Length.growth..mm.day.")]
if(met==2) d= lt[lt$fluctuation==0,c("Mean.temperature...C.","Shell.dry.weight.growth..mg.day.")]
if(met==3) d= lt[lt$fluctuation==0,c("Mean.temperature...C.","Tissue.dry.weight.growth..mg.day.")]

colnames(d)=c("temp","rate")

d_fits <- nest(d, data = c(temp, rate)) %>%
  # mutate(beta = map(data, ~nls_multstart(rate~beta_2012(temp = temp, a, b, c, d, e),
  #                                        data = .x,
  #                                        iter = c(6,6,6,6,6),
  #                                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') - 10,
  #                                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') + 10,
  #                                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'beta_2012'),
  #                                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'beta_2012'),
  #                                        supp_errors = 'Y',
  #                                        convergence_count = FALSE)),
  mutate(gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = temp, rmax, topt, a),
                                             data = .x,
                                             iter = c(4,4,4),
                                             start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 10,
                                             start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 10,
                                             lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                             upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)),
         weibull = map(data, ~nls_multstart(rate~weibull_1995(temp = temp, a,topt,b,c),
                                            data = .x,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'weibull_1995') - 10,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'weibull_1995') + 10,
                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'weibull_1995'),
                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'weibull_1995'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)))

# stack models
d_stack <- dplyr::select(d_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', gaussian:weibull) #beta:weibull

# get parameters using tidy
params <- d_stack %>%
  mutate(., est = map(fit, tidy)) %>%
  dplyr::select(-fit) %>%
  unnest(est)

mod= d_fits$weibull[[1]]

# get predictions
preds <- data.frame(temp = seq(min(d$temp), max(d$temp), length.out = 100))
preds <- broom::augment(mod, newdata = preds)

#extract coefficients
tpc.weibull= coef(mod)

#save coefficients
if(met==1) tpcs.weibull= tpc.weibull
if(met>1) tpcs.weibull= rbind(tpcs.weibull, tpc.weibull)

} #end loop metrics

#add names
rownames(tpcs.weibull)= c("length.growth","shell.growth","tissue.growth" )

#------------
# get predictions using augment
newdata <- tibble(temp = seq(min(d$temp), max(d$temp), length.out = 100))
d_preds <- d_stack %>%
  mutate(., preds = map(fit, augment, newdata = newdata)) %>%
  select(-fit) %>%
  unnest(preds)

label_facets_num <- function(string){
  len <- length(string)
  string = paste('(', 1:len, ') ', string, sep = '')
  return(string)
}

# plot
ggplot(d_preds, aes(temp, rate)) +
  geom_point(aes(temp, rate), d) +
  geom_line(aes(temp, .fitted), col = 'blue') +
  facet_wrap(~model_name, labeller = labeller(model_name = label_facets_num), scales = 'free', ncol = 5) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Fits of every model available in rTPC') +
  geom_hline(aes(yintercept = 0), linetype = 2)

#extract model
#mod= d_fits$beta[[1]]
mod= d_fits$weibull[[1]]

# get predictions
preds <- data.frame(temp = seq(min(d$temp), max(d$temp), length.out = 100))
preds <- broom::augment(mod, newdata = preds)

#extract coefficients
tpc.weibull= coef(mod)
plot(1:50, weibull_1995(1:50, tpc.weibull[1], tpc.weibull[2], tpc.weibull[3], tpc.weibull[4]), type="l")

#tpc.beta= coef(mod)
#plot(1:50, beta_2012(1:50, tpc.beta[1], tpc.beta[2], tpc.beta[3], tpc.beta[4], tpc.beta[5]), type="l", ylim=c(0,200))

#------------------------------
#SHORT TERM TPC

#See VariabilityExtremesMussels
#https://journals.plos.org/climate/article?id=10.1371/journal.pclm.0000226
#Source: https://link.springer.com/article/10.1007/s00442-012-2486-6
#use M. edulis clearance rate?

#Bayne et al 1976 data
#Mytilus californianus
#https://www.jstor.org/stable/4215210
temps= c(13, 17.5, 22, 26)
assim.fed= c(134.17, 176.06, 192.39, 113.58) #cals /day)
assim.starv=c(7.58, 9.92, 10.87, 6.42)

tpc.beta= c(198.39, 20.00, 32.00, 4.00, 4.00) 

#-------------------------------
#PLOT FITS
#length.growth
plot(1:40, weibull_1995(1:40, tpcs.weibull[1,1], tpcs.weibull[1,2], tpcs.weibull[1,3], tpcs.weibull[1,4]), type="l", ylim=range(0,.4), col="purple")
#shell.growth
points(1:40, weibull_1995(1:40, tpcs.weibull[2,1], tpcs.weibull[2,2], tpcs.weibull[2,3], tpcs.weibull[2,4]), type="l", col="blue")
#tissue growth
points(1:40, weibull_1995(1:40, tpcs.weibull[3,1], tpcs.weibull[3,2], tpcs.weibull[3,3], tpcs.weibull[3,4]), type="l", col="green")

#Bayne data TPC
tpc.beta= c(198.39, 20.00, 32.00, 4.00, 4.00) 
points(1:40, beta_2012(1:40, tpc.beta[1], tpc.beta[2], tpc.beta[3], tpc.beta[4], tpc.beta[5])/500, type="l", col="red")

#try TPC shape to fit growth
wp<- c(0.15, 22, .48*10^5, 2.1*10^4)
points(10:30, weibull_1995(10:30, wp[1], wp[2], wp[3], wp[4]), type="l")

#performance function
perf<- function(Tb) weibull_1995(Tb, tpcs.weibull[2,1], tpcs.weibull[2,2], tpcs.weibull[2,3], tpcs.weibull[2,4])

#short term feeding rate declines from peak to zero 26 to 28C
#corresponds to shell growth curve?

#STRESS ASSUMPTIONS
#Shape of cost function
#Persistence of cost
#Repair, optimal recovery conditions
#Repeated exposures

#assume stress beyond Topt
#accelerates quadratically
#multiplicative with duration
tcost<- function(Tb, dur, rep) ifelse(Tb < 26, 0, min(.1 * ((Tb - 26)*dur)^(1+0.1*rep),1))
tcost.vect<- function(x) ifelse(x[1] < 26, 0, min(.05 * ((x[1] - 26)*x[2])^(1+0.1*x[3]),1))
#need to use with apply statement

#vary parameters
tcost.vect<- function(x, cost, rep.scale) ifelse(x[1] < 26, 0, min(cost * ((x[1] - 26)*x[2])^(1+rep.scale*x[3]),1))

#cost decays linearly                                
decay <- function(time.int) max(1.0 - time.int*0.2)

#------------------------ 
#cost accelerates 
#vary cost
plot(1:40, apply(cbind(1:40, 1, 1), MARGIN=1, FUN=tcost.vect, cost=0.1, rep.scale=0.1), type="l")
points(1:40, apply(cbind(1:40, 1, 1), MARGIN=1, FUN=tcost.vect, cost=0.05, rep.scale=0.1), type="l", col="green")
points(1:40, apply(cbind(1:40, 1, 1), MARGIN=1, FUN=tcost.vect, cost=0.2, rep.scale=0.1), type="l", col="blue")

#vary rep
plot(1:40, apply(cbind(1:40, 1, 1), MARGIN=1, FUN=tcost.vect, cost=0.1, rep.scale=0.1), type="l")
points(1:40, apply(cbind(1:40, 1, 1), MARGIN=1, FUN=tcost.vect, cost=0.1, rep.scale=0.05), type="l", col="green")
points(1:40, apply(cbind(1:40, 1, 1), MARGIN=1, FUN=tcost.vect, cost=0.1, rep.scale=0.2), type="l", col="blue")

points(1:40, apply(cbind(1:40, 1, 4), MARGIN=1, FUN=tcost.vect, cost=0.1, rep.scale=0.05), type="l", col="green", lty="dashed")
points(1:40, apply(cbind(1:40, 1, 4), MARGIN=1, FUN=tcost.vect, cost=0.1, rep.scale=0.2), type="l", col="blue", lty="dashed")

#vary cost
plot(1:40, apply(cbind(1:40, 1, 1), MARGIN=1, FUN=tcost.vect, cost=0.1, rep.scale=0.1), type="l")
points(1:40, apply(cbind(1:40, 1, 1), MARGIN=1, FUN=tcost.vect, cost=0.05, rep.scale=0.1), type="l", col="green")
points(1:40, apply(cbind(1:40, 1, 1), MARGIN=1, FUN=tcost.vect, cost=0.2, rep.scale=0.1), type="l", col="blue")

points(1:40, apply(cbind(1:40, 2, 1), MARGIN=1, FUN=tcost.vect, cost=0.05, rep.scale=0.1), type="l", col="green", lty="dashed")
points(1:40, apply(cbind(1:40, 2, 1), MARGIN=1, FUN=tcost.vect, cost=0.2, rep.scale=0.1), type="l", col="blue", lty="dashed")

#------------------------  
#temperature time series
# time series Tb, plot
base <- rnorm(50, 24, 5)
time <- seq(1:50)
tp <- as.data.frame(cbind(time=time, temp=base))

#performance
tp$perf<- perf(tp$temp)

#identify heat stress events
tp$hs<- ifelse(tp$temp < 26, 0, 1)

#find lengths of heat stress
rles<- rle(tp$hs)
rs<- as.data.frame(cbind(lengths=rles$lengths, values=rles$values, cumsum=cumsum(rles$lengths))) 
rs$inds<- rs$cumsum -rs$lengths +1
#restrict to heat waves
rs<- rs[rs$values==1,]
#heat wave count
rs$rep<- 1:nrow(rs)
#mean heat wave temperature
rs$tmean<-0
for(hs.ind in 1:(nrow(rs))){
  rs$tmean[hs.ind]<- mean(tp$temp[rs$inds[hs.ind]:(rs$inds[hs.ind]+rs$lengths[hs.ind])])
}

#add durations and number
tp$dur<- 0
tp$dur[rs$inds]<- rs$lengths
tp$rep<- 0
tp$rep[rs$inds]<- rs$rep
tp$tmean<- 0
tp$tmean[rs$inds]<- rs$tmean

#cost 
tp$cost<- apply(cbind(tp$tmean, tp$dur, tp$rep), MARGIN=1, FUN=tcost.vect)
 
#decay
tp$costdec<- 0

inds=which(tp$dur>=1)
for(i in 1:length(inds)){
 dec= tp$cost[inds[i]]*c(1.0, 0.8, 0.6, 0.4, 0.2)
 ncary<- min(5,nrow(tp)-inds[i])
 tp$costdec[inds[i]:(inds[i]+ncary)]<- tp$costdec[inds[i]:(inds[i]+ncary)]+dec
}

#net performance
tp$costdec[tp$costdec>1]<-1
tp$netperf= tp$perf*(1-tp$costdec)

#to long format
tp.l<- tp[,c("time","perf","netperf")] %>%
  gather("metric", "value", 2:3)

#plot
ggplot(tp.l, aes(time, value, color=metric))+geom_line()

#sum performance
tp.sum= tp.l %>%
  group_by(metric) %>%
  summarise(perf=sum(value) )

ggplot(tp.sum, aes(metric, perf))+geom_point()

#--------------------------
#make function for costs

heat.stress <- function(temps, cost, rep.scale){

#performance
tperf<- perf(temps)
tp<- as.data.frame(cbind(temp=temps, perf=tperf))

#identify heat stress events
tp$hs<- ifelse(tp$temp < 26, 0, 1)

#find lengths of heat stress
rles<- rle(tp$hs)
rs<- as.data.frame(cbind(lengths=rles$lengths, values=rles$values, cumsum=cumsum(rles$lengths))) 
rs$inds<- rs$cumsum -rs$lengths +1
#restrict to heat waves
rs<- rs[rs$values==1,]
#heat wave count
rs$rep<- 1:nrow(rs)
#mean heat wave temperature
rs$tmean<-0
for(hs.ind in 1:(nrow(rs))){
  rs$tmean[hs.ind]<- mean(tp$temp[rs$inds[hs.ind]:(rs$inds[hs.ind]+rs$lengths[hs.ind])])
}

#add durations and number
tp$dur<- 0
tp$dur[rs$inds]<- rs$lengths
tp$rep<- 0
tp$rep[rs$inds]<- rs$rep
tp$tmean<- 0
tp$tmean[rs$inds]<- rs$tmean

#cost 
tp$cost<- apply(cbind(tp$tmean, tp$dur, tp$rep), MARGIN=1, FUN=tcost.vect, cost=cost, rep.scale=rep.scale)

#decay
tp$costdec<- 0

inds=which(tp$dur>=1)
for(i in 1:length(inds)){
  dec= tp$cost[inds[i]]*c(1.0, 0.8, 0.6, 0.4, 0.2)
  ncary<- min(5,nrow(tp)-inds[i])
  tp$costdec[inds[i]:(inds[i]+ncary)]<- tp$costdec[inds[i]:(inds[i]+ncary)]+dec
}

#net performance
tp$costdec[tp$costdec>1]<-1
tp$netperf= tp$perf*(1-tp$costdec)

return(tp)
} #end function

hs1<- heat.stress(base, cost=0.1, rep.scale=0.1)
hs1$time= time; hs1$costp= 0.1; hs1$rep.scale= 0.1
hs2<- heat.stress(base, cost=0.05, rep.scale=0.1)
hs2$time= time; hs2$costp= 0.05; hs2$rep.scale= 0.1
hs3<- heat.stress(base, cost=0.2, rep.scale=0.1)
hs3$time= time; hs3$costp= 0.2; hs3$rep.scale= 0.1
hs4<- heat.stress(base, cost=0.1, rep.scale=0.05)
hs4$time= time; hs4$costp= 0.1; hs4$rep.scale= 0.05
hs5<- heat.stress(base, cost=0.1, rep.scale=0.2)
hs5$time= time; hs5$costp= 0.1; hs5$rep.scale= 0.2

hs<- rbind(hs1, hs2, hs3, hs4, hs5)

#time series
ggplot(hs, aes(time, netperf, color=factor(costp), lty=factor(rep.scale)))+geom_line()

#sum performance
p.sum= hs %>%
  group_by(costp, rep.scale) %>%
  summarise(perf=sum(perf), netperf=sum(netperf))

ggplot(p.sum, aes(costp, netperf, color=rep.scale))+geom_point()

#--------------------------

