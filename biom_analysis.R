# This analysis measures the effect of reserve protection on Biomass outcomes

rm(list = ls())
library(tidyr)
library(tidyverse)
library(dplyr)
library(lmerTest)
library(viridis)
library(lme4)

# read data
length.dat = read.delim('data/fish_lengths.txt', header = TRUE)
weight.length = read.delim('data/weight_length.txt', header = TRUE)


# Need to get the weighted averages for size and then get biomass in grams of average size fish for each site
w.fish.lengths = length.dat %>% group_by(site.id, res, fish.id) %>% summarise(w.mean.length = weighted.mean(length.cm, prop))
w.fish.lengths.w = w.fish.lengths %>% group_by(site.id, res) %>% spread(fish.id, w.mean.length, fill = NA) 

weight.length.no = weight.length[-19,]  # Remove 'Nema' because we used an estimate of 5cm for all nemacheiline loaches due to their high total counts.
A.cm = weight.length.no$A.cm
B.cm = weight.length.no$B

log.length = log10(w.fish.lengths.w[, -c(1:2)]) # Removes the categorical variables

mass.1 = data.frame(mapply('*', log.length, B.cm))
mass.2 = data.frame(mapply('+', mass.1, A.cm))
mass = 10^mass.2 # This is the mass of fish in grams
Nema.mass = 10^(log10(5)*weight.length[19,3] + weight.length[19,2]) # because we didn't record lengths for loaches, just use 5 cm for all sites; this gets mass
Nema = rep(Nema.mass, length(mass[,1])) # make a vector of weights for 5 cm loach the length of the sites

all.mean.mass = cbind(mass[,1:18], Nema, mass[,19:34]) # stick the loach mass column, Nema, in alphabetical order for the mass matrix
all.mean.mass.data = cbind(as.data.frame(w.fish.lengths.w[,c(1:2)]), all.mean.mass)
all.mean.mass.data.2 = rbind(all.mean.mass.data, all.mean.mass.data) # double the row values for multiplying later by counts.
mean.ind.mass.data = all.mean.mass.data.2[order(all.mean.mass.data.2[,1], all.mean.mass.data.2[,2]),] # reorder the rows by site.id, then res

#write.table(mean.ind.mass.data, file = 'data/Rb/mean_ind_mass_data.txt') # This write a table for average mass per individual to be used in biomass calc.

# read in the counts data
fish_counts = read.delim('data/fish_counts.txt', header = TRUE)
# Now have data in wide format for density analysis
fish_counts$total.fish = rowSums(fish_counts[,6:40])

fish.counts.sum = fish_counts[,-4] %>% group_by(site.id, res, rep) %>% summarise_all(sum) # Combine kmp and aak counts and areas surveyed

# read in the density data
w.fish.dens = read.delim('data/Rd/comb_dens_data.txt', header = TRUE, sep = ' ') # this is fish/m^2

w.fish.biomass.dens.site = w.fish.dens[,-c(1:3,39)] * mean.ind.mass.data[,-c(1:2)] # multiply g/fish by fish/m^2 to get grams per m^2
w.fish.biomass.dens.data = cbind(w.fish.dens[,c(1:3)], w.fish.biomass.dens.site)
w.fish.biomass.dens.data$sum.biom = rowSums(w.fish.biomass.dens.site, na.rm = TRUE) 

total.lm = lmer(log10(sum.biom) ~ res + (1|site.id), data = w.fish.biomass.dens.data)
total.sum = summary(total.lm)
plot(total.lm)

all.biom.CI = confint(total.lm) # Extract 95% confidence intervals for model
coef.Res = total.sum$coefficients[2] # Extracts the Estimate for ReserveY
mean.Res = 10^(coef.Res)
lowCI = 10^(all.biom.CI[4]) # lower 95% CI
upCI = 10^(all.biom.CI[8]) # upper 95% CI

# BIG FISH BIOMASS
fish.biomass.tibl = as_tibble(w.fish.biomass.dens.data)

big.fish.id = c('Ba.de', 'Ch.ba', 'Ch.st', 'Cr.bu', 'Fo.br', 'Ga.na', 
                'Ha.sa', 'He.mi', 'Hy.sa', 'La.ro', 'Ma.ar', 'Ne.st', 
                'No.no', 'Po.sp', 'Pu.gr', 'Pu.re', 'Ra.gu', 'Sc.bu', 
                'Sp.ac', 'To.sp', 'Xe.ca', 'YOY.ns')

big.mass = rowSums(fish.biomass.tibl[, big.fish.id], na.rm = TRUE)
big.mass.data = cbind(fish.biomass.tibl[,1:3], big.mass)

#write.table(big.mass.data, file = 'data/Rb/big_mass_data.txt') # This will write the table for mass.data to be used in other analyses

par(mfrow = c(1,1))
boxplot(log10(big.mass.data$big.mass) ~ big.mass.data$res)
big.biom.lm = lmer(log10(big.mass) ~ res + (1|site.id), data = big.mass.data) # used a log10(N+0.01) transform on the y variables 
big.biom.sum = summary(big.biom.lm)
plot(big.biom.lm)

big.biom.CI = confint(big.biom.lm) # Extract 95% confidence intervals for model
coef.Res = big.biom.sum$coefficients[2] # Extracts the Estimate for ReserveY
mean.Res = 10^(coef.Res) 
lowCI = 10^(big.biom.CI[4]) # lower 95% CI
upCI = 10^(big.biom.CI[8]) # upper 95% CI


# SMALL FISH BIOMASS

fish.biomass.tibl = as_tibble(w.fish.biomass.dens.data)

lil.fish.id = c('Ba.or', 'Ch.bu', 'Da.pu', 'De.sp', 'Gl.sp', 
                'Ho.bu', 'My.ar', 'Nema', 
                'Pa.vo', 'Pe.st', 'Ra.da', 
                'Sy.ru', 'YOY.un')

lil.mass = rowSums(fish.biomass.tibl[, lil.fish.id], na.rm = TRUE)
lil.mass.data = cbind(fish.biomass.tibl[,1:3], lil.mass)
#write.table(lil.mass.data, file = 'data/Rb/lil_mass_data.txt') # This will write the table for mass.data to be used in other analyses

lil.biom.lm = lmer(log10(lil.mass) ~ res + (1|site.id), data = lil.mass.data) # used a log10(N+0.01) transform on the y variables add~ 1 gram to the total
lil.biom.sum = summary(lil.biom.lm) # not significant
plot(lil.biom.lm) 

# ::::: **** No further analysis, as res had no effect on small fish biomass ****

# PRED FISH BIOMASS
fish.biomass.tibl = as_tibble(w.fish.biomass.dens.data)

pred.fish.id = c('Ch.bu', 'Ch.st', 'Ha.sa', 'He.mi', 'No.no', 'Ra.gu', 'Sp.ac', 'Xe.ca') # :::: these fish are TP 3.5+
pred.mass = rowSums(fish.biomass.tibl[, pred.fish.id], na.rm = TRUE)
pred.mass.data = cbind(fish.biomass.tibl[,1:3], pred.mass)

#write.table(pred.mass.data, file = 'data/Rb/pred_mass_data.txt') # This will write the table for mass.data to be used in other analyses

# Pred Model ----

par(mfrow = c(1,1))
pred.biom.lm = lmer(log10(pred.mass+0.01) ~ res + (1|site.id), data = pred.mass.data) # used a log10(N+0.01) transform on the y variables for sites with 0 pred
pred.biom.sum = summary(pred.biom.lm)
plot(pred.biom.lm)

pred.biom.CI = confint(pred.biom.lm) # Extract 95% confidence intervals for model
coef.Res = pred.biom.sum$coefficients[2] # Extracts the Estimate for ReserveY
mean.Res = 10^(coef.Res)
lowCI = 10^(pred.biom.CI[4]) # lower 95% CI
upCI = 10^(pred.biom.CI[8]) # upper 95% CI


#////////////
# Omnivorous Fish
#\\\\\\\\\\\

fish.biomass.tibl = as_tibble(w.fish.biomass.dens.data) 

omni.fish.id = c('Ga.na', 'Nema', 'Ma.ar', 'Hy.sa', 'Ne.st', 'To.sp', 'Ch.ba', 'Pe.st', 'Cr.bu', 'Sy.ru', 'Da.pu', 
                 'De.sp', 'My.ar', 'Ba.or', 'Pa.vo', 'Sy.ru', 'Fo.br', 'Gl.sp', 'Po.sp', 'Pu.gr', 'Pu.re', 'Ra.da', 'YOY.ns') # ::: these fish are TP 2.49-3.5


omni.mass = rowSums(fish.biomass.tibl[, omni.fish.id], na.rm = TRUE)
omni.mass.data = cbind(fish.biomass.tibl[,1:3], omni.mass)

#write.table(omni.mass.data, file = 'data/Rb/omni_mass_data.txt') # This will write the table for mass.data to be used in other analyses

omni.biom.lm = lmer(log10(omni.mass) ~ res + (1|site.id), data = omni.mass.data)
omni.biom.sum = summary(omni.biom.lm)
plot(omni.biom.lm)

omni.biom.CI = confint(omni.biom.lm) # Extract 95% confidence intervals for model
coef.Res = omni.biom.sum$coefficients[2] # Extracts the Estimate for ReserveY
mean.Res = 10^(coef.Res) 
lowCI = 10^(omni.biom.CI[4]) # lower 95% CI
upCI = 10^(omni.biom.CI[8]) # upper 95% CI

#////////////
# Herbivorous Fish
#\\\\\\\\\\\
fish.biomass.tibl = as_tibble(w.fish.biomass.dens.data)

herb.fish.id = c('Ba.de', 'Sc.bu', 'Ho.bu', 'La.ro')

herb.mass = rowSums(fish.biomass.tibl[, herb.fish.id], na.rm = TRUE)
herb.mass.data = cbind(fish.biomass.tibl[,1:3], herb.mass)

#write.table(herb.mass.data, file = 'data/Rb/herb_mass_data.txt') # This will write the table for mass.data to be used in other analyses

herb.biom.lm = lmer(log10(herb.mass + 0.01) ~ res + (1|site.id), data = herb.mass.data) # used a log10(N+0.01) transform on the y variable for site with no herb
herb.biom.sum = summary(herb.biom.lm)
plot(herb.biom.lm)

herb.biom.CI = confint(herb.biom.lm) # Extract 95% confidence intervals for model
coef.Res = herb.biom.sum$coefficients[2] # Extracts the Estimate for ReserveY
mean.Res = 10^(coef.Res)
lowCI = 10^(herb.biom.CI[4]) # lower 95% CI
upCI = 10^(herb.biom.CI[8]) # upper 95% CI

