# This analysis measures the effect of reserve protection on Density outcomes

rm(list = ls())
library(tidyr)
library(tidyverse)
library(dplyr)
library(lme4)
library(lmerTest)
library(MuMIn)
library(viridis)


# read data
fish_counts = read.delim('fish_counts.txt', header = TRUE)
# Now have data in wide format for density analysis
fish_counts$total.fish = rowSums(fish_counts[,6:40])

#write.table(fish_counts, file = 'data/Rd/fish_counts_data.txt') # This will write the table of fish counts for use in other analyses

dens.glmer = glmer(total.fish ~ res + (1|site.id), data = fish_counts, family = 'poisson', offset = log(area))
all.dens.sum = summary(dens.glmer)

dens.glmer.0 = glmer(total.fish ~ 1 + (1|site.id), data = fish_counts, family = 'poisson', offset = log(area))
plot(ranef(dens.glmer.0)) # check random effect qnorm plot

anova(dens.glmer.0, dens.glmer, test = 'Chisq')

plot(dens.glmer) # Check model residuals versus fitted values

dens.CI = confint(dens.glmer) # Extract 95% confidence intervals for model
coef.Res = all.dens.sum$coefficients[2] # Extracts the Estimate for ReserveY
coef.Int = all.dens.sum$coefficients[1]
mean.n.Res = exp(coef.Int)
mean.Res = exp(coef.Res) # this is the # times greater/less per unit area
lowCI = exp(dens.CI[3]) # Extracts the lower CI 
upCI = exp(dens.CI[6]) # Extracts the upper CI

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# Analysis by functional traits ----
#//////////////////////////////

# Divide fish into 2 size categories: 
## 'big.fish.id' == those ids for species having maximum observed size > 20 cm; 
## 'lil.fish.id' == those ids for species having max observed size <= 20cm

big.fish.id = c('Ba.de', 'Ch.ba', 'Ch.st', 'Cr.bu', 'Fo.br', 'Ga.na', 
                'Ha.sa', 'He.mi', 'Hy.sa', 'La.ro', 'Ma.ar', 'Ne.st', 
                'No.no', 'Po.sp', 'Pu.gr', 'Pu.re', 'Ra.gu', 'Sc.bu', 
                'Sp.ac', 'To.sp', 'Xe.ca', 'YOY.ns')

lil.fish.id = c('Ba.or', 'Ch.bu', 'Da.pu', 'De.sp', 'Gl.sp', 
                'Ho.bu', 'My.ar', 'Nema', 'Pa.vo', 'Pe.st', 'Ra.da', 
                'Sy.ru', 'YOY.un')


# Big fish Density Analysis

big.cts = rowSums(fish_counts[ , big.fish.id], na.rm = TRUE) # Pull counts of only the big fish
big.dens = cbind(as.data.frame(fish_counts[,c(1:5)]), big.cts)

#write.table(big.dens, file = 'data/Rd/big_dens_data.txt') # This will write the table for dens.data to be used in other analyses

big.dens.glmer = glmer(big.cts ~ res + (1|site.id), data = big.dens, family = 'poisson', offset = log(area))
big.dens.sum = summary(big.dens.glmer)

# Create the null model, no RESERVE fixed effect
big.dens.glmer.0 = glmer(big.cts ~ 1 + (1|site.id), data = big.dens, family = 'poisson', offset = log(area))
plot(ranef(big.dens.glmer.0)) # check random effect qnorm plot

# Chisq test for fit of Reserve against null model
anova(big.dens.glmer.0, big.dens.glmer, test = 'Chisq')

plot(big.dens.glmer) # Check model residuals versus fitted values

big.dens.CI = confint(big.dens.glmer) # Extract 95% confidence intervals for model
coef.Res = big.dens.sum$coefficients[2] # Extracts the Estimate for ReserveY
coef.Int = big.dens.sum$coefficients[1]
mean.Res = exp(coef.Res)
lowCI = exp(big.dens.CI[3]) # Extracts the lower 95% CI
upCI = exp(big.dens.CI[6]) # Extracts the upper 95% CI

# Little Fish Density Analysis

lil.cts = rowSums(fish_counts[ , lil.fish.id], na.rm = TRUE) # Pull counts of only the lil fish
lil.dens = cbind(as.data.frame(fish_counts[,c(1:5)]), lil.cts)

#write.table(lil.dens, file = 'data/Rd/lil_dens_data.txt') # This will write the table for dens.data to be used in other analyses

lil.dens.glmer = glmer(lil.cts ~ res + (1|site.id), data = lil.dens, family = 'poisson', offset = log(area))
lil.dens.sum = summary(lil.dens.glmer)

# Create the null model, no RESERVE fixed effect
lil.dens.glmer.0 = glmer(lil.cts ~ 1 + (1|site.id), data = lil.dens, family = 'poisson', offset = log(area))
plot(ranef(lil.dens.glmer.0)) # check random effect qnorm plot

# Chisq test for fit of Reserve against null model
anova(lil.dens.glmer.0, lil.dens.glmer, test = 'Chisq')

plot(lil.dens.glmer) # Check model residuals versus fitted values

# ::::: **** No further analysis, as res had no effect on small fish density ****


# Divide fish into 3 trophic categories: 
## 'pred.fish.id' == those ids for species having trophic position > 3.5
## 'omni.fish.id' == those ids for species having trophic position <=3.5 & > 2.5
## 'herb.fish.id' == those ids for species having trophic position <= 2.5
# *** Leave YOY.un out due to unknown trophic position

pred.fish.id = c('Ch.bu', 'Ch.st', 'Ha.sa', 'He.mi', 'No.no', 'Ra.gu', 'Sp.ac', 'Xe.ca')

omni.fish.id = c('Ga.na', 'Nema', 'Ma.ar', 'Hy.sa', 'Ne.st', 'To.sp', 'Ch.ba', 'Pe.st', 'Cr.bu' , 'Da.pu', 
                 'De.sp', 'My.ar', 'Ba.or', 'Pa.vo', 'Sy.ru', 'Fo.br', 'Gl.sp', 'Po.sp', 'Pu.gr', 'Pu.re', 'Ra.da', 'YOY.ns')
herb.fish.id = c('Ba.de', 'Sc.bu', 'Ho.bu', 'La.ro')

# Predatory Fish Density Analysis

pred.cts = rowSums(fish_counts[ , pred.fish.id], na.rm = TRUE) # Pull counts of only the pred fish
pred.dens = cbind(as.data.frame(fish_counts[,c(1:5)]), pred.cts)

#write.table(pred.dens, file = 'data/Rd/pred_dens_data.txt') # This will write the table for dens.data to be used in other analyses

pred.dens.glmer = glmer(pred.cts ~ res + (1|site.id), data = pred.dens, family = 'poisson', offset = log(area))
pred.dens.sum = summary(pred.dens.glmer)

# Create the null model, no RESERVE fixed effect
pred.dens.glmer.0 = glmer(pred.cts ~ 1 + (1|site.id), data = pred.dens, family = 'poisson', offset = log(area))
plot(ranef(pred.dens.glmer.0))

# Chisq test for fit of Reserve against null model
anova(pred.dens.glmer.0, pred.dens.glmer, test = 'Chisq')

plot(pred.dens.glmer) # Check model residuals versus fitted values

pred.dens.CI = confint(pred.dens.glmer) # Extract 95% confidence intervals for model
coef.Res = pred.dens.sum$coefficients[2] # Extracts the Estimate for ReserveY
mean.Res = exp(coef.Res) # 
lowCI = exp(pred.dens.CI[3]) # Extracts the lower 95% CI
upCI = exp(pred.dens.CI[6]) # Extracts the upper 95% CI


# Omni Fish Density Analysis

omni.cts = rowSums(fish_counts[ , omni.fish.id], na.rm = TRUE) # Pull counts of only the omni fish
omni.dens = cbind(as.data.frame(fish_counts[,c(1:5)]), omni.cts)

#write.table(omni.dens, file = 'data/Rd/omni_dens_data.txt') # This will write the table for dens.data to be used in other analyses

omni.dens.glmer = glmer(omni.cts ~ res + (1|site.id), data = omni.dens, family = 'poisson', offset = log(area))
omni.dens.sum = summary(omni.dens.glmer)

# Create the null model, no RESERVE fixed effect
omni.dens.glmer.0 = glmer(omni.cts ~ 1 + (1|site.id), data = omni.dens, family = 'poisson', offset = log(area))
plot(ranef(omni.dens.glmer.0)) # check random effect qnorm plot

# Chisq test for fit of Reserve against null model
anova(omni.dens.glmer.0, omni.dens.glmer, test = 'Chisq')

plot(omni.dens.glmer) # Check model residuals versus fitted values

omni.dens.CI = confint(omni.dens.glmer) # Extract 95% confidence intervals for model
coef.Res = omni.dens.sum$coefficients[2] # Extracts the Estimate for ReserveY
mean.Res = exp(coef.Res) 
lowCI = exp(omni.dens.CI[3]) # Extracts the lower 95% CI
upCI = exp(omni.dens.CI[6]) # Extracts THE upper 95% CI

# Herb Fish Density Analysis

herb.fish.id = c('Ba.de', 'Sc.bu', 'Ho.bu', 'La.ro')

herb.cts = rowSums(fish_counts[ , herb.fish.id], na.rm = TRUE) # Pull counts of only the herb fish
herb.dens = cbind(as.data.frame(fish_counts[,c(1:5)]), herb.cts)

#write.table(herb.dens, file = 'data/Rd/herb_dens_data.txt') # This will write the table for dens.data to be used in other analyses

herb.dens.glmer = glmer(herb.cts ~ res + (1|site.id), data = herb.dens, family = 'poisson', offset = log(area))
herb.dens.sum = summary(herb.dens.glmer)

# Create the null model, no RESERVE fixed effect
herb.dens.glmer.0 = glmer(herb.cts ~ 1 + (1|site.id), data = herb.dens, family = 'poisson', offset = log(area))
plot(ranef(herb.dens.glmer.0)) # check random effect qnorm plot

# Chisq test for fit of Reserve against null model
anova(herb.dens.glmer.0, herb.dens.glmer, test = 'Chisq')

plot(herb.dens.glmer) # Check model residuals versus fitted values

herb.dens.CI = confint(herb.dens.glmer) # Extract 95% confidence intervals for model
coef.Res = herb.dens.sum$coefficients[2] # Extracts the Estimate for ReserveY
mean.Res = exp(coef.Res)
lowCI = exp(herb.dens.CI[3]) # Extracts the lower 95% CI
upCI = exp(herb.dens.CI[6]) # Extracts upper 95% CI

