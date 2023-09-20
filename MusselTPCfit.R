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

