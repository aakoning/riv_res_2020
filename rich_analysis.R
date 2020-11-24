
# This analysis measures the effect of reserve protection on Richness outcomes

#Load required packages
rm(list = ls())  
library(tidyr)
library(dplyr)
library(vegan)
library(lme4)
library(loggit)

# Read in data
rep.rich.mat = read.delim('richness_matrix.txt', header = TRUE) # These area binary values (0,1) for species presence/absence at a given rep (i.e., observed by either observer)

rep.rich = specnumber(rep.rich.mat[,-c(1:4)]) # calculates the richness per rep (n=2) per site; rep.cts[,-c(1:4)] removes 'site.id', 'res', 'rep', and 'area' columns for analysis

rep.rich.data = cbind(rep.rich.mat[,c(1:4)], rep.rich)
#write.table(rep.rich.data, file = 'data/Rs/rep_rich_data.txt') # This will write the table for rep.rich.data to be used in other analyses

# Poisson regression with a random intercept for 'site.id'
all.rich.glm = glmer(rep.rich ~ res + (1|site.id), offset = log10(area), data = rep.rich.data, family = 'poisson') 
summary(all.rich.glm)
plot(ranef(all.rich.glm)) # This plots Q-Q plot for random effects
plot(all.rich.glm) # Check residuals

all.rich.sum = summary(all.rich.glm) # get summary output
all.rich.confint = confint(all.rich.glm) # generate confidence intervals; default is 95% CI

# Interpretation as incidence rate ratios, need to exponentiate the coefficients and CI.
# rate of X times greater/less for predictor
all.rich.res = exp(all.rich.sum$coefficients[2]) 
all.rich.dn95CI = exp(all.rich.confint[3,1]) # lower 95% CI
all.rich.up95CI = exp(all.rich.confint[3,2]) # upper 95% CI

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# Analysis by functional traits ----
#//////////////////////////////

# Divide fish into 2 size categories: 
## 'big.fish.id' == those ids for species having maximum observed size > 20 cm; 
## 'lil.fish.id' == those ids for species having max observed size <= 20cm

big.fish.id = c('Ba.de', 'Ch.ba', 'Ch.st', 'Cr.bu', 'Fo.br', 'Ga.na', 
                'Ha.sa', 'He.mi', 'Hy.sa', 'La.ro', 'Ma.ar', 'Ne.st', 
                'No.no', 'Po.sp', 'Pu.gr', 'Pu.re', 'Ra.gu', 'Sc.bu', 
                'Sp.ac', 'To.sp', 'Xe.ca')

lil.fish.id = c('Ba.or', 'Ch.bu', 'Da.pu', 'De.sp', 'Gl.sp', 
                'Ho.bu', 'My.ar', 'Nema', 'Pa.vo', 'Pe.st', 'Ra.da', 
                'Sy.ru')  #We don't include YOY.un (just for richness) because we don't know the species


# big.fish Richness Analysis ----

# Create a vector of richness values for only those fish included in the character string set above.
big.rich = specnumber(rep.rich.mat[, big.fish.id])
# Combine richness for 'big.fish' with the site.id, res, rep columns
big.rich.data = cbind(as.data.frame(rep.rich.data[,1:4]), big.rich)
#write.table(big.rich.data, file = 'data/Rs/big_rich_data.txt') # This will write the table for big.rich.data to be used in other analyses

# Poisson Regression
big.rich.glm = glmer(big.rich ~ res + (1|site.id), offset = log10(area), data = big.rich.data, family = 'poisson')
summary(big.rich.glm)
plot(ranef(big.rich.glm))
plot(big.rich.glm)

big.rich.glm.sum = summary(big.rich.glm)
big.rich.confint = confint(big.rich.glm)

big.rich.res = exp(big.rich.glm.sum$coefficients[2])
big.rich.dn95CI = exp(big.rich.confint[3,1])
big.rich.up95CI = exp(big.rich.confint[3,2])

# lil.fish Richness Analysis ----

# Create a vector of richness values for only those fish included in the character string set above.
lil.rich = specnumber(rep.rich.mat[,lil.fish.id])
# Combine richness for lil fish with the site.id, res, rep columns
lil.rich.data = cbind(as.data.frame(rep.rich.data[,1:4]), lil.rich)
#write.table(lil.rich.data, file = 'data/Rs/lil_rich_data.txt') # This will write the table for lil.rich.data to be used in other analyses

# Poisson regression
lil.rich.glm = glmer(lil.rich ~ res + (1|site.id), offset = log10(area), data = lil.rich.data, family = 'poisson')
summary(lil.rich.glm)
plot(ranef(lil.rich.glm))
plot(lil.rich.glm)

# ::::: **** No further analysis, as res had no effect on small fish richness ****


# Divide fish into 3 trophic categories: 
## 'pred.fish.id' == those ids for species having trophic position > 3.5
## 'omni.fish.id' == those ids for species having trophic position <=3.5 & > 2.5
## 'herb.fish.id' == those ids for species having trophic position <= 2.5
# **** we don't include YOY.un in this analysis, as they're trophic position is unidentified
pred.fish.id = c('Ch.bu', 'Ch.st', 'Ha.sa', 'He.mi', 'No.no', 'Ra.gu', 'Sp.ac', 'Xe.ca')

omni.fish.id = c('Ga.na', 'Nema', 'Ma.ar', 'Hy.sa', 'Ne.st', 'To.sp', 'Ch.ba', 'Pe.st', 'Cr.bu' , 'Da.pu', 
                 'De.sp', 'My.ar', 'Ba.or', 'Pa.vo', 'Sy.ru', 'Fo.br', 'Gl.sp', 'Po.sp', 'Pu.gr', 'Pu.re', 'Ra.da')

herb.fish.id = c('Ba.de', 'Sc.bu', 'Ho.bu', 'La.ro')


# pred.fish analysis ----

pred.rich = specnumber(rep.rich.mat[,pred.fish.id])
# Combine richness for pred fish with the site.id, res, rep columns
pred.rich.data = cbind(as.data.frame(rep.rich.data[,1:4]), pred.rich)
#write.table(pred.rich.data, file = 'data/Rs/pred_rich_data.txt')

# Poisson
pred.rich.glm = glmer(pred.rich ~ res + (1|site.id), offset = log10(area), data = pred.rich.data, family = 'poisson')
summary(pred.rich.glm)
plot(ranef(pred.rich.glm))
plot(pred.rich.glm)

pred.rich.glm.sum = summary(pred.rich.glm)
pred.rich.confint = confint(pred.rich.glm)

# ::::: **** No further analysis, as res had no effect on predator fish richness ****

# omni.fish analysis ----

# Create a vector of richness values for only those fish included in the character string set above.
omni.rich = specnumber(rep.rich.mat[,omni.fish.id])
# Combine richness for omni fish with the site.id, res, rep columns
omni.rich.data = cbind(as.data.frame(rep.rich.data[,1:4]), omni.rich)
#write.table(omni.rich.data, file = 'data/Rs/omni_rich_data.txt') # This will write the table for omni.rich.data to be used in other analyses

# Poisson regression
omni.rich.glm = glmer(omni.rich ~ res + (1|site.id), data = omni.rich.data, family = 'poisson')
summary(omni.rich.glm)
plot(ranef(omni.rich.glm))
plot(omni.rich.glm)

omni.rich.glm.sum = summary(omni.rich.glm)
omni.rich.confint = confint(omni.rich.glm)

omni.rich.res = exp(omni.rich.glm.sum$coefficients[2])
omni.rich.dn95CI = exp(omni.rich.confint[3,1]) # lower 95% CI
omni.rich.up95CI = exp(omni.rich.confint[3,2]) # upper 95% CI

# herb.fish analysis ----

herb.rich = specnumber(rep.rich.mat[,herb.fish.id])
# Combine richness for herb fish with the site.id, res, rep columns
herb.rich.data = cbind(as.data.frame(rep.rich.data[,1:4]), herb.rich)
#write.table(herb.rich.data, file = 'data/Rs/herb_rich_data.txt') # This will write the table for herb.rich.data to be used in other analyses

# Poisson
herb.rich.glm = glmer(herb.rich ~ res + (1|site.id), offset = log10(area), data = herb.rich.data, family = 'poisson')
summary(herb.rich.glm)
plot(ranef(herb.rich.glm))
plot(herb.rich.glm)

herb.rich.glm.sum = summary(herb.rich.glm)
herb.rich.confint = confint(herb.rich.glm)

herb.rich.res = exp(herb.rich.glm.sum$coefficients[2]) 
herb.rich.dn95CI = exp(herb.rich.confint[3,1]) # lower 95% CI
herb.rich.up95CI = exp(herb.rich.confint[3,2]) # upper 95% CI

