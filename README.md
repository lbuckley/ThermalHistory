---
editor_options: 
  markdown: 
    wrap: sentence
---

# How damage, recovery, and repair alter the fitness impacts of thermal stress

Analysis of how thermal history influences performance

# GENERAL INFORMATION

This README.txt file was updated on February 19, 2025 by Lauren Buckley

## A. Paper associated with this archive

Citation: Buckley LB et al. 2025.
How damage, recovery, and repair alter the fitness impacts of thermal stress.
Integrative and Comparative Biology

Synopsis: Integrating estimates of organisms' thermal sensitivity of performance across time series of body temperatures is a common strategy for estimating the fitness implications of climate variability and change.
The approach is an important first step in accounting for organismal sensitivity when estimating climate responses, but the approach omits the influence of past exposure to thermal stress.
We fit data from experiments exposing aphids to stressful, fluctuating temperatures to a function accounting for thermal history to estimate rates of damage, repair, and carryover effects in response to thermal stress.
Our analysis indicates that thermal stress and damage initiate at temperatures beyond thermal optima and considerably below critical thermal limits; that thermal stress can accumulate multiplicatively with the magnitude and duration of exposure; and that considerable and rapid rates of recovery and repair are possible but that rates decline quickly as temperatures diverge from the thermal optima.
Accounting for damage and repair across thermal history is important to accurately assessing the fitness outcomes of climate change and variability.

## B. Originators

Lauren B. Buckley and Raymond B. Huey, Department of Biology, University of Washington, Seattle, WA 98195-1800, USA

Data orgininator and author: Chun-Sen Ma, School of Life Science, Hebei University, Baoding, 071002, China

## C. Contact information

Lauren Buckley.
Department of Biology, University of Washington, Seattle, WA 98195-1800, USA.
[lbuckley\@uw.edu](mailto:lbuckley@uw.edu){.email}

Chun-Sen Ma.
School of Life Science, Hebei University, Baoding, 071002, China.\
[machunsen\@caas.cn](mailto:machunsen@caas.cn){.email}

## D. Dates of data collection

No new data are collected, but see references for past data utilized.

## E. Geographic Location(s) of data collection

See source publications

## F. Funding Sources

See source publications for funding that enabled data collection.
This work was supported by the National Science Foundation (IOS-2222089 and DEB-1951356 to L.B.B), the Key Program of National Natural Science Foundation of China (32330090), National Key Research and Development Program of China (2022YFD1400401), and Hebei Natural Science Foundation (C20222001042).

# ACCESS INFORMATION

## 1. Licenses/restrictions placed on the data or code

CC0 1.0 Universal (CC0 1.0).
Public Domain Dedication

## 2. Data derived from other sources

Data from 5 publications reporting on aphid experiments are organized into folders within the data folder:

Zhaoetal2014: Expt 1.
vary min in manuscript.\
Zhao F, Zhang W, Hoffmann AA, Ma C.
2014.
Night warming on hot days produces novel impacts on development, survival and reproduction in a small arthropod.
Journal of Animal Ecology 83:769--78.\
Dryad data: <http://doi.org/10.5061/dryad.q2070>

Maetal2015: Expt 2.
vary max in manuscript.\
Ma G, Hoffmann AA, Ma C-S.
2015.
2015.
Daily temperature extremes play an important role in predicting thermal effects.
The Journal of Experimental Biology 218 (14), 2289-2296.
<https://doi.org/10.1242/jeb.122127>.\
Data provided by Chun-Sen Ma

WangMa2023: Expt 3: vary fluctuations, Expt 4: vary means and fluctuations in manuscript.\
Wang X-J, Ma C-S.
2023.
Can laboratory-reared aphid populations reflect the thermal performance of field populations in studies on pest science and climate change biology?
J Pest Sci 96:509--22.\
Data available via paper: <https://doi.org/10.1007/s10340-022-01565-6>

Maetal2017: Expt 5: hot days in manuscript.\
Ma CS, Wang L, Zhang W, Rudolf V 2018.
Resolving biological impacts of multiple heat waves: interaction of hot and recovery days.
Oikos 127:622--33.
<https://doi.org/10.1111/oik.04699>.
Dryad data: <http://dx.doi.org/10.5061/dryad.5qk4s>

Zhaoetal2019: Expt 6: adult heatwaves, Expt 7: nymphal heatwaves in manuscript.\
Zhao F, Xing K, Hoffmann AA, Ma C.
2019.
The importance of timing of heat events for predicting the dynamics of aphid pest populations.
Pest Management Science 75:1866--74.
<https://doi.org/10.1002/ps.5344>.\
Data provided by Chun-Sen Ma

## 3. Recommended citation for this data/code archive

Data: source publications above

Code: Buckley LB et al. 2025.
How damage, recovery, and repair alter the fitness impacts of thermal stress.
<https://github.com/lbuckley/ThermalHistory>.

# DATA & CODE FILE OVERVIEW

This data repository consist of 16 data files, 4 code scripts, and this README document, with the following data and code filenames and variables.
Raw data is in the data folder; data is processed and written to the processed_data folder; parameter estimates are written to the out folder; and figures are written to the figures folder.

## Data files and variables

data/Zhaoetal2014/: [see initial data source for formats].\
1. Zhaoetal2014_devtime.csv: development time data.\
2. Zhaoetal2014_AdPerf.csv: adult performance data.\
3. Zhaoetal2014_LifeTable.csv: life table data.\
4. Zhaoetal2014_SurvNymph.csv: nymphal survival data.

data/Maetal2015/:\
1.Maetal2015_temps.csv: temperature data.\
Hour: hr (0-24).\
Dmax_32: Temperature (degree C) during indicated hour for treatment with 32C daily maximum.\
Dmax_34: Temperature (degree C) during indicated hour for treatment with 34C daily maximum.\
Dmax_36: Temperature (degree C) during indicated hour for treatment with 36C daily maximum.\
Dmax_38: Temperature (degree C) during indicated hour for treatment with 38C daily maximum.\
Dmax_40: Temperature (degree C) during indicated hour for treatment with 40C daily maximum.

2.  Maetal2015_perf.csv: performance data.\
species: Species name.\
Dmax_C: Daily maximum temperature (degree C) for treatment.\
cage No.: Replicate number.\
Dur_I: Developmental time (hours) of first instar.\
Dur_II: Developmental time (hours) of second instar.\
Dur_III: Developmental time (hours) of third instar.\
Dur_IV: Developmental time (hours) of forth instar.\
DurL: Developmental time (hours) until adulthood.\
longevity: longevity (hours).\
fecundity: number of nymphs per adult.

3.  Maetal2015.csv: Inital data. Not used in analysis.

data/WangMa2023/: [see initial data source for formats].\
1. WangMa2023_temp22mean_diffvar.csv: data from mild means experiment.\
2. WangMa2023_diffmeans.csv: data from high means experiment.\
3. WangMa2023_popgrowth.csv: population growth data

data/Maetal2017/: [see initial data source for formats].\
1. Development.csv: Development data.\
2. Reproduction.csv: Reproduction data.\
3. Traits.csv: Trait data

data/Zhaoetal2019/:\
1. Zhaoetal2019_temps.csv: temperature data.\
ind: individual number.\
Time, GMT+08:00: time in "m/d/y h:m" format.\
am_pm: am or pm time.\
Temp, °C: temperature (°C).\
treatment: treatment- normal or high temperatures.

2.  Zhaoetal2019_perf.csv: performance data.\
    Treatments: heatwave treatment.\
    id: id number.\
    Nymph duration: nymphal duration (days).\
    Longevity: longevity (days).\
    Productivity: number of nymphs per adult.\
    Lifeperiod: time until death (days) of those who died immediately.\
    immediately death: No (0) or Yes (1)

3.  Zhaoetal2019.csv: Inital data. Not used in analysis.

## Code scripts and workflow

1.  AphidExpts_TempsFitnesss.R: Assembles aphid data and corresponding environmental data for analysis.
    The data are from studies of the English grain aphid, Sitobion avenae, from the research group led by Chen-Sen Ma.

2.  FunctionOptimization.R: Fits repair and damage function to aphid data

3.  Figs1_3_Conceptual.R: Creates conceptual figures 1 and 3 illustrating the modelling approach

4.  Figs2_4_5_6_Aphids.R: Creates figures 2, which illustrate experimental temperature regimes, and figures 4, 5, and 6, which illustrate analysis results.

# SOFTWARE VERSIONS

R version 4.1.0 (2021-05-18)

Packages:

library(ggplot2) #ggplot2_3.4.4

library(reshape2) #reshape2_1.4.4

library(patchwork) #patchwork_1.2.0

library(plyr) #plyr_1.8.8

library(dplyr) #dplyr_1.1.2

library(tidyr) #tidyr_1.3.0

library(viridis) #viridis_0.6.4

library(viridisLite) #viridisLite_0.4.2

library(deSolve) #deSolve_1.36

library(pracma) #pracma_2.4.4

library(TrenchR) #TrenchR_1.1.1

library(zoo) #zoo_1.8-12

library(rvmethod) #rvmethod_0.1.2

library(dfoptim) #dfoptim_2023.1.0

library(ggpubr) #ggpubr_0.6.0

library(scico) #scico_1.5.0

# REFERENCES

Ma G, Hoffmann AA, Ma C-S.
2015.
2015.
Daily temperature extremes play an important role in predicting thermal effects.
The Journal of Experimental Biology 218 (14), 2289-2296.
<https://doi.org/10.1242/jeb.122127>.

Ma CS, Wang L, Zhang W, Rudolf V 2018.
Resolving biological impacts of multiple heat waves: interaction of hot and recovery days.
Oikos 127:622--33.
<https://doi.org/10.1111/oik.04699> Dryad data: <http://dx.doi.org/10.5061/dryad.5qk4s>

Wang X-J, Ma C-S.
2023.
Can laboratory-reared aphid populations reflect the thermal performance of field populations in studies on pest science and climate change biology?
J Pest Sci 96:509--22.
<https://doi.org/10.1007/s10340-022-01565-6>

Zhao F, Zhang W, Hoffmann AA, Ma C.
2014.
Night warming on hot days produces novel impacts on development, survival and reproduction in a small arthropod.
Journal of Animal Ecology 83:769--78.
Dryad data: <http://doi.org/10.5061/dryad.q2070>

Zhao F, Xing K, Hoffmann AA, Ma C.
2019.
The importance of timing of heat events for predicting the dynamics of aphid pest populations.
Pest Management Science 75:1866--74.https://doi.org/10.1002/ps.5344
