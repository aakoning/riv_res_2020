# This analysis measures the effect of reserve characteristics on reserve Biomass response (Br)

rm(list = ls())
library(dplyr)
library(tidyr)
library(lme4)
library(MuMIn)
library(arm)

# read data

comb.dens.data = read.delim('data/Rd/comb_dens_data.txt', header = TRUE, sep = ' ') # This reads in the densities of fish (fish/m^2) from reserve_feat_dens_analysis.R, low discharge sites observers are summed, high discharge are weighted means. 
mean.ind.mass.data = read.delim('data/Rb/mean_ind_mass_data.txt', header = TRUE, sep = ' ') # This reads in the average mass per individual for each site, calculated in biomass_analysis_new.R
res.data = read.delim('data/reserve_features.txt', header = TRUE)

# Multiply the density (ind/m^2) by biomass (g/ind) to get (g/m^2)
all.biom.dens = comb.dens.data[,-c(1:3,39)] * mean.ind.mass.data[,-c(1:2)] 
all.biom.dens.data = cbind(comb.dens.data[,c(1:3)], all.biom.dens) # This is g/m^2 for each species
all.biom.dens.data[is.na(all.biom.dens.data)] <- 0 # Replace 'NA' with 0's

all.biom.dens.data$total.biom = rowSums(all.biom.dens.data[,-c(1:3)]) # sum the biomass per rep and create new column 'total.biom'

all.biom.dens.tot = all.biom.dens.data[, c(1:3,39)] # stick the categoricals and total.biom together for all analysis

all.biom.dens.w = all.biom.dens.tot %>% group_by(site.id) %>% pivot_wider(names_from = res, values_from = total.biom) # spreads Reserve column to Y and N columns with density average for observers at each rep (2 per site.id)

all.biom.dens.w.means = all.biom.dens.w %>% 
  group_by(site.id) %>% 
  summarise(N.mean = mean(N), N.sd = sd(N), Y.mean = mean(Y), Y.sd = sd(Y)) # creates dataframe by Village that includes mean and sd of Y and N columns

all.biom.dens.rr = all.biom.dens.w.means %>% 
  group_by(site.id) %>% 
  mutate(rr.biom = log10(Y.mean/N.mean), se.rr.biom = sqrt(((Y.sd)^2/(2*(Y.mean)^2) + (N.sd)^2/(2*(N.mean)^2))), perc.change.dens = ((Y.mean-N.mean)/N.mean)*100) # creates new dataframe for logRR and for difference

all.biom.dens.rr.data = cbind(as.data.frame(all.biom.dens.rr), res.data)

all.biom.dens.rr.lm = lm(rr.biom ~ yrs.prot + houses + pen + mouth.dist + discharge + clo.res.dist 
                    + bet.cent.n + res.area + rd.dist + vil.dist, 
                    data = all.biom.dens.rr.data, na.action = na.pass) # 

summary(all.biom.dens.rr.lm)

# z-score standardizes all predictor variables to account for non-normality
all.biom.dens.rr.lm.std = standardize(all.biom.dens.rr.lm, standardize.y = FALSE) 
summary(all.biom.dens.rr.lm.std)

# dredge function from MuMIn package calculates all possible models and evaluates fit
all.biom.dens.rr.model.set = dredge(all.biom.dens.rr.lm.std)
sum(all.biom.dens.rr.model.set$weight[1:269]) # we use this to calculate the threshold for inclusion in model averaging (ie summed weights ≤ 0.95)
all.biom.dens.rr.95.models = get.models(all.biom.dens.rr.model.set, 1:269)

all.biom.dens.rr.95.models.avg = model.avg(all.biom.dens.rr.95.models, revised.var =TRUE)
summary(all.biom.dens.rr.95.models.avg) # Read the model averaged output

# Create a barplot of top models by importance

all.biom.dens.rr.impt = as.matrix(all.biom.dens.rr.95.models.avg$sw) # extracts the Akaike weights for each term in the top models
all.biom.dens.rr.coef = as.matrix(all.biom.dens.rr.95.models.avg$coefficients[2,-1]) # extracts the coefficients from the averaged model for coloring by positive or negative

all.biom.dens.rr.bar = merge(all.biom.dens.rr.impt, all.biom.dens.rr.coef, by = 'row.names', all.x = TRUE) # combines the Akaike weights data and coefficients 
row.names(all.biom.dens.rr.bar) = all.biom.dens.rr.bar$Row.names
all.biom.dens.rr.bar = all.biom.dens.rr.bar[,-1]

mod.pred = matrix(rep(NA, times = 10), dimnames = list(c('z.yrs.prot', 'z.houses', 'c.pen', 'z.mouth.dist', 'z.discharge', 'z.clo.res.dist', 'z.bet.cent.n', 'z.res.area', 'z.rd.dist', 'z.vil.dist')))
mod.pred.sort = as.matrix(mod.pred[order(row.names(mod.pred)),])

all.biom.dens.rr.all = merge(mod.pred.sort, all.biom.dens.rr.bar, by = 'row.names', all.x = TRUE) # combines the Akaike weights data and coefficients 
all.biom.dens.rr.all = all.biom.dens.rr.all[,-2]
all.biom.dens.rr.all[is.na(all.biom.dens.rr.all)] <- 0
colnames(all.biom.dens.rr.all) = c('p.variable', 'importance', 'coefficient')

target = c('z.res.area', 'z.discharge', 'z.yrs.prot', 'c.pen', 'z.vil.dist', 'z.rd.dist', 'z.houses', 'z.bet.cent.n', 'z.mouth.dist', 'z.clo.res.dist')
all.biom.dens.rr.plot = all.biom.dens.rr.all[match(target, all.biom.dens.rr.all$p.variable),]
par(mfrow = c(1,6))
par(oma = c(5,3,5,3))

all.plot = barplot(all.biom.dens.rr.plot$importance, names.arg = all.biom.dens.rr.plot$p.variable, 
                   col = ifelse(all.biom.dens.rr.plot$coefficient > 0, 'black', 'white'), 
                   las = 2, main = 'All Fishes',
                   ylim = c(0,1)) 


# Get data for partial residual plots
all.biom.full = all.biom.dens.rr.lm.std
all.biom.avg = all.biom.dens.rr.95.models.avg
all.biom.mod = all.biom.full
coefMod = coef(all.biom.avg)

summary(all.biom.avg)

# Changing the coefficents to the averaged coefficient # see https://stackoverflow.com/questions/28633537/partial-residual-plot-based-on-model-average-coefficients-in-r
all.biom.mod$coefficients <- coefMod[names(coef(all.biom.full))]
all.biom.mod$residuals <- all.biom.dens.rr.data$rr.biom - predict(all.biom.mod)
all.biom.pres <- resid(all.biom.mod, type='partial')

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


big.biom.dens = rowSums(all.biom.dens.data[ , big.fish.id], na.rm = TRUE) # Pull counts of only the big fish
big.biom.dens.data = cbind(all.biom.dens.data[,c(1:3)], big.biom.dens) # add back in the categorical variables, and area and size

big.biom.dens.w = big.biom.dens.data %>% group_by(site.id) %>% pivot_wider(names_from = res, values_from = big.biom.dens) # spreads Reserve column to Y and N columns with density average for observers at each rep (2 per site.id)

big.biom.dens.w.means = big.biom.dens.w %>% 
  group_by(site.id) %>% 
  summarise(N.mean = mean(N), N.sd = sd(N), Y.mean = mean(Y), Y.sd = sd(Y)) # creates dataframe by Village that includes mean and sd of Y and N columns

big.biom.dens.rr = big.biom.dens.w.means %>% 
  group_by(site.id) %>% 
  mutate(rr.biom = log10(Y.mean/N.mean), se.rr.biom = sqrt(((Y.sd)^2/(2*(Y.mean)^2) + (N.sd)^2/(2*(N.mean)^2))), perc.change.dens = ((Y.mean-N.mean)/N.mean)*100) # creates new dataframe for logRR and for difference


big.biom.dens.rr.data = cbind(as.data.frame(big.biom.dens.rr), res.data)

big.biom.dens.rr.lm = lm(rr.biom ~ yrs.prot + houses + pen + mouth.dist + discharge + clo.res.dist 
                         + bet.cent.n + res.area + rd.dist + vil.dist, 
                         data = big.biom.dens.rr.data, na.action = na.pass) # 

summary(big.biom.dens.rr.lm)

# z-score standardizes all predictor variables to account for non-normality
big.biom.dens.rr.lm.std = standardize(big.biom.dens.rr.lm, standardize.y = FALSE)

# dredge function from MuMIn package calculates all possible models and evaluates fit
big.biom.dens.rr.model.set = dredge(big.biom.dens.rr.lm.std)
sum(big.biom.dens.rr.model.set$weight[1:299]) # we use this to calculate the threshold for inclusion in model averaging (ie summed weights ≤ 0.95)
big.biom.dens.rr.95.models = get.models(big.biom.dens.rr.model.set, 1:299)

big.biom.dens.rr.95.models.avg = model.avg(big.biom.dens.rr.95.models, revised.var =TRUE) 
summary(big.biom.dens.rr.95.models.avg) # Read the model averaged output

# Create a barplot of top models by importance

big.biom.dens.rr.impt = as.matrix(big.biom.dens.rr.95.models.avg$sw) # extracts the Akaike weights for each term in the top models
big.biom.dens.rr.coef = as.matrix(big.biom.dens.rr.95.models.avg$coefficients[2,-1]) # extracts the coefficients from the averaged model for coloring by positive or negative

big.biom.dens.rr.bar = merge(big.biom.dens.rr.impt, big.biom.dens.rr.coef, by = 'row.names', big.x = TRUE) # combines the Akaike weights data and coefficients 
row.names(big.biom.dens.rr.bar) = big.biom.dens.rr.bar$Row.names
big.biom.dens.rr.bar = big.biom.dens.rr.bar[,-1]

mod.pred = matrix(rep(NA, times = 10), dimnames = list(c('z.yrs.prot', 'z.houses', 'c.pen', 'z.mouth.dist', 'z.discharge', 'z.clo.res.dist', 'z.bet.cent.n', 'z.res.area', 'z.rd.dist', 'z.vil.dist')))
mod.pred.sort = as.matrix(mod.pred[order(row.names(mod.pred)),])

big.biom.dens.rr.all = merge(mod.pred.sort, big.biom.dens.rr.bar, by = 'row.names', big.x = TRUE) # combines the Akaike weights data and coefficients 
big.biom.dens.rr.all = big.biom.dens.rr.all[,-2]
big.biom.dens.rr.all[is.na(big.biom.dens.rr.all)] <- 0
colnames(big.biom.dens.rr.all) = c('p.variable', 'importance', 'coefficient')

target = c('z.res.area', 'z.discharge', 'z.yrs.prot', 'c.pen', 'z.vil.dist', 'z.rd.dist', 'z.houses', 'z.bet.cent.n', 'z.mouth.dist', 'z.clo.res.dist')
big.biom.dens.rr.plot = big.biom.dens.rr.all[match(target, big.biom.dens.rr.all$p.variable),]

big.plot = barplot(big.biom.dens.rr.plot$importance, names.arg = big.biom.dens.rr.plot$p.variable, 
                   col = ifelse(big.biom.dens.rr.plot$coefficient > 0, 'black', 'white'), 
                   las = 2, main = 'Large Fishes',
                   ylim = c(0,1)) 


# Get data for partial residual plots
big.biom.full = big.biom.dens.rr.lm.std
big.biom.avg = big.biom.dens.rr.95.models.avg
big.biom.mod = big.biom.full
coefMod = coef(big.biom.avg)

summary(big.biom.avg)

# Changing the coefficents to the averaged coefficient # see https://stackoverflow.com/questions/28633537/partial-residual-plot-based-on-model-average-coefficients-in-r
big.biom.mod$coefficients <- coefMod[names(coef(big.biom.full))]
# Changing the residuals
big.biom.mod$residuals <- big.biom.dens.rr.data$rr.biom - predict(big.biom.mod)
big.biom.pres <- resid(big.biom.mod, type='partial')


# \\\\\\\\\\
# SMALL FISH
# //////////

lil.biom.dens = rowSums(all.biom.dens.data[ , lil.fish.id], na.rm = TRUE) # Pull counts of only the lil fish
lil.biom.dens.data = cbind(all.biom.dens.data[,c(1:3)], lil.biom.dens) # add back in the categorical variables, and area and size

lil.biom.dens.w = lil.biom.dens.data %>% group_by(site.id) %>% pivot_wider(names_from = res, values_from = lil.biom.dens) # spreads Reserve column to Y and N columns with density average for observers at each rep (2 per site.id)

lil.biom.dens.w.means = lil.biom.dens.w %>% 
  group_by(site.id) %>% 
  summarise(N.mean = mean(N), N.sd = sd(N), Y.mean = mean(Y), Y.sd = sd(Y)) # creates dataframe by Village that includes mean and sd of Y and N columns

lil.biom.dens.rr = lil.biom.dens.w.means %>% 
  group_by(site.id) %>% 
  mutate(rr.biom = log10(Y.mean/N.mean), se.rr.biom = sqrt(((Y.sd)^2/(2*(Y.mean)^2) + (N.sd)^2/(2*(N.mean)^2))), perc.change.dens = ((Y.mean-N.mean)/N.mean)*100) # creates new dataframe for logRR and for difference


lil.biom.dens.rr.data = cbind(as.data.frame(lil.biom.dens.rr), res.data)

lil.biom.dens.rr.lm = lm(rr.biom ~ yrs.prot + houses + pen + mouth.dist + discharge + clo.res.dist 
                         + bet.cent.n + res.area + rd.dist + vil.dist, 
                         data = lil.biom.dens.rr.data, na.action = na.pass) # 

summary(lil.biom.dens.rr.lm) # Overall model not significant

# z-score standardizes all predictor variables to account for non-normality
lil.biom.dens.rr.lm.std = standardize(lil.biom.dens.rr.lm, standardize.y = FALSE) 

# dredge function from MuMIn package calculates all possible models and evaluates fit
lil.biom.dens.rr.model.set = dredge(lil.biom.dens.rr.lm.std)
sum(lil.biom.dens.rr.model.set$weight[1:271]) # we use this to calculate the threshold for inclusion in model averaging (ie summed weights ≤ 0.95)
lil.biom.dens.rr.95.models = get.models(lil.biom.dens.rr.model.set, 1:271)

lil.biom.dens.rr.95.models.avg = model.avg(lil.biom.dens.rr.95.models, revised.var =TRUE) 
summary(lil.biom.dens.rr.95.models.avg) # Read the model averaged output

# Create a barplot of top models by importance

lil.biom.dens.rr.impt = as.matrix(lil.biom.dens.rr.95.models.avg$sw) # extracts the Akaike weights for each term in the top models
lil.biom.dens.rr.coef = as.matrix(lil.biom.dens.rr.95.models.avg$coefficients[2,-1]) # extracts the coefficients from the averaged model for coloring by positive or negative

lil.biom.dens.rr.bar = merge(lil.biom.dens.rr.impt, lil.biom.dens.rr.coef, by = 'row.names', lil.x = TRUE) # combines the Akaike weights data and coefficients 
row.names(lil.biom.dens.rr.bar) = lil.biom.dens.rr.bar$Row.names
lil.biom.dens.rr.bar = lil.biom.dens.rr.bar[,-1]

mod.pred = matrix(rep(NA, times = 10), dimnames = list(c('z.yrs.prot', 'z.houses', 'c.pen', 'z.mouth.dist', 'z.discharge', 'z.clo.res.dist', 'z.bet.cent.n', 'z.res.area', 'z.rd.dist', 'z.vil.dist')))
mod.pred.sort = as.matrix(mod.pred[order(row.names(mod.pred)),])

lil.biom.dens.rr.all = merge(mod.pred.sort, lil.biom.dens.rr.bar, by = 'row.names', lil.x = TRUE) # combines the Akaike weights data and coefficients 
lil.biom.dens.rr.all = lil.biom.dens.rr.all[,-2]
lil.biom.dens.rr.all[is.na(lil.biom.dens.rr.all)] <- 0
colnames(lil.biom.dens.rr.all) = c('p.variable', 'importance', 'coefficient')

target = c('z.res.area', 'z.discharge', 'z.yrs.prot', 'c.pen', 'z.vil.dist', 'z.rd.dist', 'z.houses', 'z.bet.cent.n', 'z.mouth.dist', 'z.clo.res.dist')
lil.biom.dens.rr.plot = lil.biom.dens.rr.all[match(target, lil.biom.dens.rr.all$p.variable),]

lil.plot = barplot(lil.biom.dens.rr.plot$importance, names.arg = lil.biom.dens.rr.plot$p.variable, 
                   col = ifelse(lil.biom.dens.rr.plot$coefficient > 0, 'black', 'white'), 
                   las = 2, main = 'Small Fishes',
                   ylim = c(0,1)) 

# Get data for partial residual plots
lil.biom.full = lil.biom.dens.rr.lm.std
lil.biom.avg = lil.biom.dens.rr.95.models.avg
lil.biom.mod = lil.biom.full
coefMod = coef(lil.biom.avg)

summary(lil.biom.avg)

# Changing the coefficents to the averaged coefficient # see https://stackoverflow.com/questions/28633537/partial-residual-plot-based-on-model-average-coefficients-in-r
lil.biom.mod$coefficients <- coefMod[names(coef(lil.biom.full))]
lil.biom.mod$residuals <- lil.biom.dens.rr.data$rr.biom - predict(lil.biom.mod)
lil.biom.pres <- resid(lil.biom.mod, type='partial')


# Divide fish into 3 trophic categories: 
## 'pred.fish.id' == those ids for species having trophic position > 3.5
## 'omni.fish.id' == those ids for species having trophic position <=3.5 & > 2.5
## 'herb.fish.id' == those ids for species having trophic position <= 2.5
# *** Leave YOY.un out due to unknown trophic position

pred.fish.id = c('Ch.bu', 'Ch.st', 'Ha.sa', 'He.mi', 'No.no', 'Ra.gu', 'Sp.ac', 'Xe.ca')

omni.fish.id = c('Ga.na', 'Nema', 'Ma.ar', 'Hy.sa', 'Ne.st', 'To.sp', 'Ch.ba', 'Pe.st', 'Cr.bu' , 'Da.pu', 
                 'De.sp', 'My.ar', 'Ba.or', 'Pa.vo', 'Sy.ru', 'Fo.br', 'Gl.sp', 'Po.sp', 'Pu.gr', 'Pu.re', 'Ra.da', 'YOY.ns')

herb.fish.id = c('Ba.de', 'Sc.bu', 'Ho.bu', 'La.ro')

# \\\\\\\\\\
# PREDATORY FISH
# //////////

pred.biom.dens = rowSums(all.biom.dens.data[ , pred.fish.id], na.rm = TRUE) # Pull counts of only the pred fish
pred.biom.dens.data = cbind(all.biom.dens.data[,c(1:3)], pred.biom.dens) # add back in the categorical variables, and area and size

pred.biom.dens.w = pred.biom.dens.data %>% group_by(site.id) %>% pivot_wider(names_from = res, values_from = pred.biom.dens) # spreads Reserve column to Y and N columns with density average for observers at each rep (2 per site.id)

pred.biom.dens.w.means = pred.biom.dens.w %>% 
  group_by(site.id) %>% 
  summarise(N.mean = mean(N+0.01), N.sd = sd(N), Y.mean = mean(Y+0.01), Y.sd = sd(Y)) # creates dataframe by Village that includes mean and sd of Y and N columns

pred.biom.dens.rr = pred.biom.dens.w.means %>% 
  group_by(site.id) %>% 
  mutate(rr.biom = log10(Y.mean/N.mean), se.rr.biom = sqrt(((Y.sd)^2/(2*(Y.mean)^2) + (N.sd)^2/(2*(N.mean)^2))), perc.change.dens = ((Y.mean-N.mean)/N.mean)*100) # creates new dataframe for logRR and for difference


pred.biom.dens.rr.data = cbind(as.data.frame(pred.biom.dens.rr), res.data)

pred.biom.dens.rr.lm = lm(rr.biom ~ yrs.prot + houses + pen + mouth.dist + discharge + clo.res.dist 
                         + bet.cent.n + res.area + rd.dist + vil.dist, 
                         data = pred.biom.dens.rr.data, na.action = na.pass) # 

summary(pred.biom.dens.rr.lm) # Overall model not significant

# z-score standardizes all predictor variables to account for non-normality
pred.biom.dens.rr.lm.std = standardize(pred.biom.dens.rr.lm, standardize.y = FALSE) 

# dredge function from MuMIn package calculates all possible models and evaluates fit
pred.biom.dens.rr.model.set = dredge(pred.biom.dens.rr.lm.std)
sum(pred.biom.dens.rr.model.set$weight[1:292]) # we use this to calculate the threshold for inclusion in model averaging (ie summed weights ≤ 0.95)
pred.biom.dens.rr.95.models = get.models(pred.biom.dens.rr.model.set, 1:292)

pred.biom.dens.rr.95.models.avg = model.avg(pred.biom.dens.rr.95.models, revised.var =TRUE) 
summary(pred.biom.dens.rr.95.models.avg) # Read the model averaged output

# Create a barplot of top models by importance

pred.biom.dens.rr.impt = as.matrix(pred.biom.dens.rr.95.models.avg$sw) # extracts the Akaike weights for each term in the top models
pred.biom.dens.rr.coef = as.matrix(pred.biom.dens.rr.95.models.avg$coefficients[2,-1]) # extracts the coefficients from the averaged model for coloring by positive or negative

pred.biom.dens.rr.bar = merge(pred.biom.dens.rr.impt, pred.biom.dens.rr.coef, by = 'row.names', pred.x = TRUE) # combines the Akaike weights data and coefficients 
row.names(pred.biom.dens.rr.bar) = pred.biom.dens.rr.bar$Row.names
pred.biom.dens.rr.bar = pred.biom.dens.rr.bar[,-1]

mod.pred = matrix(rep(NA, times = 10), dimnames = list(c('z.yrs.prot', 'z.houses', 'c.pen', 'z.mouth.dist', 'z.discharge', 'z.clo.res.dist', 'z.bet.cent.n', 'z.res.area', 'z.rd.dist', 'z.vil.dist')))
mod.pred.sort = as.matrix(mod.pred[order(row.names(mod.pred)),])

pred.biom.dens.rr.all = merge(mod.pred.sort, pred.biom.dens.rr.bar, by = 'row.names', pred.x = TRUE) # combines the Akaike weights data and coefficients 
pred.biom.dens.rr.all = pred.biom.dens.rr.all[,-2]
pred.biom.dens.rr.all[is.na(pred.biom.dens.rr.all)] <- 0
colnames(pred.biom.dens.rr.all) = c('p.variable', 'importance', 'coefficient')

target = c('z.res.area', 'z.discharge', 'z.yrs.prot', 'c.pen', 'z.vil.dist', 'z.rd.dist', 'z.houses', 'z.bet.cent.n', 'z.mouth.dist', 'z.clo.res.dist')
pred.biom.dens.rr.plot = pred.biom.dens.rr.all[match(target, pred.biom.dens.rr.all$p.variable),]

pred.plot = barplot(pred.biom.dens.rr.plot$importance, names.arg = pred.biom.dens.rr.plot$p.variable, 
                   col = ifelse(pred.biom.dens.rr.plot$coefficient > 0, 'black', 'white'), 
                   las = 2, main = 'Predators',
                   ylim = c(0,1)) 

# Get data for partial residual plots
pred.biom.full = pred.biom.dens.rr.lm.std
pred.biom.avg = pred.biom.dens.rr.95.models.avg
pred.biom.mod = pred.biom.full
coefMod = coef(pred.biom.avg)

summary(pred.biom.avg)

# Changing the coefficents to the averaged coefficient # see https://stackoverflow.com/questions/28633537/partial-residual-plot-based-on-model-average-coefficients-in-r
pred.biom.mod$coefficients <- coefMod[names(coef(pred.biom.full))]
# Changing the residuals
pred.biom.mod$residuals <- pred.biom.dens.rr.data$rr.biom - predict(pred.biom.mod)
pred.biom.pres <- resid(pred.biom.mod, type='partial')



# \\\\\\\\\\
# OMNIVOROUS FISH
# //////////

omni.biom.dens = rowSums(all.biom.dens.data[ , omni.fish.id], na.rm = TRUE) # Pull counts of only the omni fish
omni.biom.dens.data = cbind(all.biom.dens.data[,c(1:3)], omni.biom.dens) # add back in the categorical variables, and area and size

omni.biom.dens.w = omni.biom.dens.data %>% group_by(site.id) %>% pivot_wider(names_from = res, values_from = omni.biom.dens) # spreads Reserve column to Y and N columns with density average for observers at each rep (2 per site.id)

omni.biom.dens.w.means = omni.biom.dens.w %>% 
  group_by(site.id) %>% 
  summarise(N.mean = mean(N), N.sd = sd(N), Y.mean = mean(Y), Y.sd = sd(Y)) # creates dataframe by Village that includes mean and sd of Y and N columns

omni.biom.dens.rr = omni.biom.dens.w.means %>% 
  group_by(site.id) %>% 
  mutate(rr.biom = log10(Y.mean/N.mean), se.rr.biom = sqrt(((Y.sd)^2/(2*(Y.mean)^2) + (N.sd)^2/(2*(N.mean)^2))), perc.change.dens = ((Y.mean-N.mean)/N.mean)*100) # creates new dataframe for logRR and for difference


omni.biom.dens.rr.data = cbind(as.data.frame(omni.biom.dens.rr), res.data)

omni.biom.dens.rr.lm = lm(rr.biom ~ yrs.prot + houses + pen + mouth.dist + discharge + clo.res.dist 
                          + bet.cent.n + res.area + rd.dist + vil.dist, 
                          data = omni.biom.dens.rr.data, na.action = na.pass) # 

summary(omni.biom.dens.rr.lm)

# z-score standardizes all omniictor variables to account for non-normality
omni.biom.dens.rr.lm.std = standardize(omni.biom.dens.rr.lm, standardize.y = FALSE) 

# dredge function from MuMIn package calculates all possible models and evaluates fit
omni.biom.dens.rr.model.set = dredge(omni.biom.dens.rr.lm.std)
sum(omni.biom.dens.rr.model.set$weight[1:305]) # we use this to calculate the threshold for inclusion in model averaging (ie summed weights ≤ 0.95)
omni.biom.dens.rr.95.models = get.models(omni.biom.dens.rr.model.set, 1:305)

omni.biom.dens.rr.95.models.avg = model.avg(omni.biom.dens.rr.95.models, revised.var =TRUE) 
summary(omni.biom.dens.rr.95.models.avg) # Read the model averaged output

# Create a barplot of top models by importance

omni.biom.dens.rr.impt = as.matrix(omni.biom.dens.rr.95.models.avg$sw) # extracts the Akaike weights for each term in the top models
omni.biom.dens.rr.coef = as.matrix(omni.biom.dens.rr.95.models.avg$coefficients[2,-1]) # extracts the coefficients from the averaged model for coloring by positive or negative

omni.biom.dens.rr.bar = merge(omni.biom.dens.rr.impt, omni.biom.dens.rr.coef, by = 'row.names', omni.x = TRUE) # combines the Akaike weights data and coefficients 
row.names(omni.biom.dens.rr.bar) = omni.biom.dens.rr.bar$Row.names
omni.biom.dens.rr.bar = omni.biom.dens.rr.bar[,-1]

mod.omni = matrix(rep(NA, times = 10), dimnames = list(c('z.yrs.prot', 'z.houses', 'c.pen', 'z.mouth.dist', 'z.discharge', 'z.clo.res.dist', 'z.bet.cent.n', 'z.res.area', 'z.rd.dist', 'z.vil.dist')))
mod.omni.sort = as.matrix(mod.omni[order(row.names(mod.omni)),])

omni.biom.dens.rr.all = merge(mod.omni.sort, omni.biom.dens.rr.bar, by = 'row.names', omni.x = TRUE) # combines the Akaike weights data and coefficients 
omni.biom.dens.rr.all = omni.biom.dens.rr.all[,-2]
omni.biom.dens.rr.all[is.na(omni.biom.dens.rr.all)] <- 0
colnames(omni.biom.dens.rr.all) = c('p.variable', 'importance', 'coefficient')

target = c('z.res.area', 'z.discharge', 'z.yrs.prot', 'c.pen', 'z.vil.dist', 'z.rd.dist', 'z.houses', 'z.bet.cent.n', 'z.mouth.dist', 'z.clo.res.dist')
omni.biom.dens.rr.plot = omni.biom.dens.rr.all[match(target, omni.biom.dens.rr.all$p.variable),]

omni.plot = barplot(omni.biom.dens.rr.plot$importance, names.arg = omni.biom.dens.rr.plot$p.variable, 
                    col = ifelse(omni.biom.dens.rr.plot$coefficient > 0, 'black', 'white'), 
                    las = 2, main = 'Omnivores',
                    ylim = c(0,1)) 

# Get data for partial residual plots
omni.biom.full = omni.biom.dens.rr.lm.std
omni.biom.avg = omni.biom.dens.rr.95.models.avg
omni.biom.mod = omni.biom.full
coefMod = coef(omni.biom.avg)


# Changing the coefficents to the averaged coefficient # see https://stackoverflow.com/questions/28633537/partial-residual-plot-based-on-model-average-coefficients-in-r
omni.biom.mod$coefficients <- coefMod[names(coef(omni.biom.full))]
omni.biom.mod$residuals <- omni.biom.dens.rr.data$rr.biom - predict(omni.biom.mod)
omni.biom.pres <- resid(omni.biom.mod, type='partial')



# \\\\\\\\\\
# HERBIVOROUS FISH
# //////////

herb.biom.dens = rowSums(all.biom.dens.data[ , herb.fish.id], na.rm = TRUE) # Pull counts of only the herb fish
herb.biom.dens.data = cbind(all.biom.dens.data[,c(1:3)], herb.biom.dens) # add back in the categorical variables, and area and size

herb.biom.dens.w = herb.biom.dens.data %>% group_by(site.id) %>% pivot_wider(names_from = res, values_from = herb.biom.dens) # spreads Reserve column to Y and N columns with density average for observers at each rep (2 per site.id)

herb.biom.dens.w.means = herb.biom.dens.w %>% 
  group_by(site.id) %>% 
  summarise(N.mean = mean(N+0.01), N.sd = sd(N), Y.mean = mean(Y+0.01), Y.sd = sd(Y)) # creates dataframe by Village that includes mean and sd of Y and N columns

herb.biom.dens.rr = herb.biom.dens.w.means %>% 
  group_by(site.id) %>% 
  mutate(rr.biom = log10(Y.mean/N.mean), se.rr.biom = sqrt(((Y.sd)^2/(2*(Y.mean)^2) + (N.sd)^2/(2*(N.mean)^2))), perc.change.dens = ((Y.mean-N.mean)/N.mean)*100) # creates new dataframe for logRR and for difference


herb.biom.dens.rr.data = cbind(as.data.frame(herb.biom.dens.rr), res.data)
write.table(herb.biom.dens.rr.data, 'figures/fig_4/data/herb.biom.rr.data.txt')

herb.biom.dens.rr.lm = lm(rr.biom ~ yrs.prot + houses + pen + mouth.dist + discharge + clo.res.dist 
                          + bet.cent.n + res.area + rd.dist + vil.dist, 
                          data = herb.biom.dens.rr.data, na.action = na.pass) # 

summary(herb.biom.dens.rr.lm) # Overall model not significant

# z-score standardizes all herbictor variables to account for non-normality
herb.biom.dens.rr.lm.std = standardize(herb.biom.dens.rr.lm, standardize.y = FALSE) 

# dredge function from MuMIn package calculates all possible models and evaluates fit
herb.biom.dens.rr.model.set = dredge(herb.biom.dens.rr.lm.std)
sum(herb.biom.dens.rr.model.set$weight[1:305]) # we use this to calculate the threshold for inclusion in model averaging (ie summed weights ≤ 0.95)
herb.biom.dens.rr.95.models = get.models(herb.biom.dens.rr.model.set, 1:305)

herb.biom.dens.rr.95.models.avg = model.avg(herb.biom.dens.rr.95.models, revised.var =TRUE) 
summary(herb.biom.dens.rr.95.models.avg) # Read the model averaged output

# Create a barplot of top models by importance

herb.biom.dens.rr.impt = as.matrix(herb.biom.dens.rr.95.models.avg$sw) # extracts the Akaike weights for each term in the top models
herb.biom.dens.rr.coef = as.matrix(herb.biom.dens.rr.95.models.avg$coefficients[2,-1]) # extracts the coefficients from the averaged model for coloring by positive or negative

herb.biom.dens.rr.bar = merge(herb.biom.dens.rr.impt, herb.biom.dens.rr.coef, by = 'row.names', herb.x = TRUE) # combines the Akaike weights data and coefficients 
row.names(herb.biom.dens.rr.bar) = herb.biom.dens.rr.bar$Row.names
herb.biom.dens.rr.bar = herb.biom.dens.rr.bar[,-1]

mod.herb = matrix(rep(NA, times = 10), dimnames = list(c('z.yrs.prot', 'z.houses', 'c.pen', 'z.mouth.dist', 'z.discharge', 'z.clo.res.dist', 'z.bet.cent.n', 'z.res.area', 'z.rd.dist', 'z.vil.dist')))
mod.herb.sort = as.matrix(mod.herb[order(row.names(mod.herb)),])

herb.biom.dens.rr.all = merge(mod.herb.sort, herb.biom.dens.rr.bar, by = 'row.names', herb.x = TRUE) # combines the Akaike weights data and coefficients 
herb.biom.dens.rr.all = herb.biom.dens.rr.all[,-2]
herb.biom.dens.rr.all[is.na(herb.biom.dens.rr.all)] <- 0
colnames(herb.biom.dens.rr.all) = c('p.variable', 'importance', 'coefficient')

target = c('z.res.area', 'z.discharge', 'z.yrs.prot', 'c.pen', 'z.vil.dist', 'z.rd.dist', 'z.houses', 'z.bet.cent.n', 'z.mouth.dist', 'z.clo.res.dist')
herb.biom.dens.rr.plot = herb.biom.dens.rr.all[match(target, herb.biom.dens.rr.all$p.variable),]

herb.plot = barplot(herb.biom.dens.rr.plot$importance, names.arg = herb.biom.dens.rr.plot$p.variable, 
                    col = ifelse(herb.biom.dens.rr.plot$coefficient > 0, 'black', 'white'), 
                    las = 2, main = 'Herbivores',
                    ylim = c(0,1)) 


# Get data for partial residual plots
herb.biom.full = herb.biom.dens.rr.lm.std
herb.biom.avg = herb.biom.dens.rr.95.models.avg
herb.biom.mod = herb.biom.full
coefMod = coef(herb.biom.avg)

summary(herb.biom.avg)

# Changing the coefficents to the averaged coefficient # see https://stackoverflow.com/questions/28633537/partial-residual-plot-based-on-model-average-coefficients-in-r
herb.biom.mod$coefficients <- coefMod[names(coef(herb.biom.full))]
herb.biom.mod$residuals <- herb.biom.dens.rr.data$rr.biom - predict(herb.biom.mod)
herb.biom.pres <- resid(herb.biom.mod, type='partial')


# All averaged model summaries
summary(all.biom.avg) 
summary(big.biom.avg)
summary(lil.biom.avg) 
summary(pred.biom.avg)
summary(omni.biom.avg)
summary(herb.biom.avg)

#\\\\\\\\\\\\\\\\\\\\\\
#Partial residual plots (Fig. 4) ----
#//////////////////////

library(viridis)
col.plot.m = magma(10)
col.plot.v = viridis(10)
big.col = col.plot.m[4]
big.col.a = adjustcolor(big.col, alpha.f = 0.5)
all.col = col.plot.v[6]
all.col.a = adjustcolor(all.col, alpha.f = 0.5)
lil.col = col.plot.m[5]
lil.col.a = adjustcolor(lil.col, alpha.f = 0.5)
pred.col = col.plot.m[6]
pred.col.a = adjustcolor(pred.col, alpha.f = 0.5)
omni.col = col.plot.m[8]
omni.col.a = adjustcolor(omni.col, alpha.f = 0.5)
herb.col = col.plot.v[9]
herb.col.a = adjustcolor(herb.col, alpha.f = 0.5)
all.pch = 21
#axis.col = '#D1D3D4'
axis.col = 'black'

dev.new = TRUE
#pdf('figures/fig_4/single_plots/biom_area_resid_DB.pdf', width = 7, height = 6)
pdf('figures/fig_4/Figure_4_biom_resid_all.pdf', width = 15, height = 4)

#par(mfrow = c(2,1))
#par(mar = c(4,3,2,1), oma = c(0,0,2,0))

par(mfrow = c(1,5), oma = c(0,4,0,1))

# AREA
plot(omni.biom.pres[,"z.res.area"] ~ omni.biom.dens.rr.data$res.area, pch = all.pch, col = omni.col, bg = omni.col.a, cex = 2, bty = 'n',
     ylab = '', xlab = '', axes = FALSE, ylim = c(-2.5,2.5), xlim = c(0, 25000), type = 'n')

abline(h = 0, col = 'grey')

points(herb.biom.pres[,"z.res.area"] ~ herb.biom.dens.rr.data$res.area, pch = all.pch, bg = herb.col.a, cex = 2, col = herb.col)
points(pred.biom.pres[,"z.res.area"] ~ pred.biom.dens.rr.data$res.area, pch = all.pch, bg = pred.col.a, cex = 2, col = pred.col)
points(omni.biom.pres[,"z.res.area"] ~ omni.biom.dens.rr.data$res.area, pch = all.pch, col = omni.col, bg = omni.col.a, cex = 2)
points(all.biom.pres[,"z.res.area"] ~ all.biom.dens.rr.data$res.area, pch = all.pch, col = all.col, cex = 2, bg = all.col.a)

abline(lm(pred.biom.pres[,"z.res.area"]~pred.biom.dens.rr.data$res.area), lty = 2, col = pred.col, lwd = 2) # p < 0.1
abline(lm(herb.biom.pres[,"z.res.area"]~herb.biom.dens.rr.data$res.area), lty = 2, col = herb.col, lwd = 2) # p < 0.1
abline(lm(all.biom.pres[,"z.res.area"]~all.biom.dens.rr.data$res.area), lty = 1, col = all.col, lwd = 2) # p < 0.05
abline(lm(omni.biom.pres[,"z.res.area"]~omni.biom.dens.rr.data$res.area), lty = 1, col = omni.col, lwd = 2) # p < 0.05


axis(side = 1, at = seq(0,25000, 12500), labels = c('0','1.25', '2.5'), cex.axis = 1.5, col = axis.col, col.axis = axis.col)
axis(side = 2, at = c(-2.5, 0, 2.5), labels = c('-2.5', '0', '2.5'), las = 1, cex.axis = 1.5, col = axis.col, col.axis = axis.col)
mtext(side = 1, 'Reserve Area (Ha)', line = 2.5, cex = 1.2, col  = axis.col)
mtext(side = 2, 'Partial Residuals', line = 3, cex = 1.2, col = axis.col)
mtext(expression(paste('B'[italic(r)])), side = 2, line = 5.5, las = 1, cex = 1.7, col = axis.col)

#dev.off()

# Distance to Village

#dev.new = TRUE

#pdf('figures/fig_4/single_plots/biom_vildist_resid_DB.pdf', width = 7, height = 6)


plot(herb.biom.pres[,"z.vil.dist"] ~ all.biom.dens.rr.data$vil.dist, pch = all.pch, col = herb.col, cex = 2, bg = herb.col.a, bty = 'n',
     ylab = '', xlab = '', axes = FALSE, ylim = c(-2,2), xlim = c(0, 2000), type = 'n')

abline(h = 0, col = 'grey')
points(herb.biom.pres[,"z.vil.dist"] ~ herb.biom.dens.rr.data$vil.dist, pch = all.pch, col = herb.col, cex = 2, bg = herb.col.a)
points(big.biom.pres[,"z.vil.dist"] ~ big.biom.dens.rr.data$vil.dist, pch = all.pch, col = big.col, cex = 2, bg = big.col.a)
abline(lm(herb.biom.pres[,"z.vil.dist"]~herb.biom.dens.rr.data$vil.dist), lty = 2, col = herb.col, lwd = 2) # p < 0.05
abline(lm(big.biom.pres[,"z.vil.dist"]~big.biom.dens.rr.data$vil.dist), lty = 2, col = big.col, lwd = 2) # p < 0.05

axis(side = 1, at = c(0, 1000, 2000), labels = c('0','1','2'), cex.axis = 1.5, col = axis.col, col.axis = axis.col)
axis(side = 2, at = c(-2,0,2), las = 1, cex.axis = 1.5, col = axis.col, col.axis = axis.col)
mtext(side = 1, 'Distance to Village (km)', line = 2.5, cex = 1.2, col  = axis.col)


dev.off()

# Create Heatmap for figure 3: for in-text figure, see figure 3 code in figures folder.

all.biom.dens.rr.plot$fish.cat = rep('All', length(all.biom.dens.rr.plot$coefficient))
big.biom.dens.rr.plot$fish.cat = rep('Large', length(big.biom.dens.rr.plot$coefficient))
lil.biom.dens.rr.plot$fish.cat = rep('Small', length(lil.biom.dens.rr.plot$coefficient))
pred.biom.dens.rr.plot$fish.cat = rep('Predators', length(pred.biom.dens.rr.plot$coefficient))
omni.biom.dens.rr.plot$fish.cat = rep('Omnivores', length(omni.biom.dens.rr.plot$coefficient))
herb.biom.dens.rr.plot$fish.cat = rep('Herbivores', length(herb.biom.dens.rr.plot$coefficient))

biom.dens.heat.data = as.data.frame(rbind(all.biom.dens.rr.plot, big.biom.dens.rr.plot, lil.biom.dens.rr.plot, pred.biom.dens.rr.plot, omni.biom.dens.rr.plot, herb.biom.dens.rr.plot))

attach(biom.dens.heat.data)
biom.dens.heat.data$imp.sign <- with(biom.dens.heat.data, ifelse(coefficient <= 0, -importance, importance)) # this adds a new column which makes the sign of the importance match the coefficient for color ramp

biom.dens.heat.wide = as.data.frame(biom.dens.heat.data[,-c(2:3)] %>% group_by(p.variable) %>% pivot_wider(names_from = fish.cat, values_from = imp.sign)) # This removes the importance and coefficient columns, combined into imp.sign, then makes a wide dataframe
rownames(biom.dens.heat.wide) = biom.dens.heat.wide[,1]
biom.dens.heat.wide = biom.dens.heat.wide[,-1]

# need to reorder the variables here
target = rev(c('z.res.area', 'z.discharge', 'z.yrs.prot', 'c.pen', 'z.vil.dist', 'z.rd.dist', 'z.houses', 'z.bet.cent.n', 'z.mouth.dist', 'z.clo.res.dist'))
biom.dens.heat.wide.z = biom.dens.heat.wide[match(target, rownames(biom.dens.heat.wide)),]
biom.dens.heat.wide.z = biom.dens.heat.wide.z[c('All', 'Large', 'Small', 'Predators', 'Omnivores', 'Herbivores')]

# write a table for plotting heatmaps
write.table(biom.dens.heat.wide.z, 'figures/fig_3/data/biom.heatmap.data.txt')

library(RColorBrewer)
col.pal = colorRampPalette(brewer.pal(11,'RdYlBu'))(60)

abc = heatmap(as.matrix(biom.dens.heat.wide.z), Rowv = NA, Colv = NA, col = col.pal, scale = 'none')


 
