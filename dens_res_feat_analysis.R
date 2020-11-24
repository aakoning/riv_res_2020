# This analysis measures the effect of reserve characteristics on reserve Density response (Dr)

#rm(list = ls())

library(tidyr)
library(dplyr)
library(vegan)
library(MuMIn)
library(lme4)
library(arm)
library(glmmTMB)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(pals)


# Set Working Directory
dens.data = read.delim('fish_counts.txt', header = TRUE)
dens.data$total.fish = rowSums(dens.data[6:40])

#dens.data = read.delim('data/Rd/all_dens_data.txt', header = TRUE, sep = ' ')
res.data = read.delim('data/reserve_features.txt', header = TRUE)

# ::::: Get the high.dis DENSITY data for RR calculation ----

#::: Calculate the densities at each survey by survey method

high.dis = dens.data[c(dens.data$site.id == 'ML' | dens.data$site.id == 'ND' | dens.data$site.id == 'SMP' | dens.data$site.id == 'SK'),]
aak.cts = high.dis[high.dis$obs == 'aak', c(6:41)] # selects only the rows for aak observations and only fish counts (no categorical variables)
kmp.cts = high.dis[high.dis$obs == 'kmp', c(6:41)] # selects only the rows for kmp observations and only fish counts (no categorical variables)
high.dis.dens = high.dis[,-c(1:5)] %>% mutate_all(funs(./high.dis$area)) # Calculate the density (fish/m2) by dividing species counts by area 
high.dis.dens.dat = cbind(high.dis[,c(1:5)], high.dis.dens) # add the first 5 columns back for categorical variables
aak.dens = high.dis.dens.dat[high.dis.dens.dat$obs == 'aak', c(6:41)]
kmp.dens = high.dis.dens.dat[high.dis.dens.dat$obs == 'kmp', c(6:41)]

#::: Need to weight the species densities by survey method (obs) to account for different fish communities targeted by each approach
weighted.dens = ((aak.dens*aak.cts)+(kmp.dens*kmp.cts))/(aak.cts + kmp.cts) # This is the weighted combination of these species specific densities

#::: Function to remove NAs and replace with )
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

weighted.dens[is.nan(weighted.dens)] <- 0 

high.dens.dat = cbind(high.dis[high.dis$obs == 'aak', c(1:3)], weighted.dens) # adds the categorical variables back in with the weighted densities for large sites

# ::::: Get the low.dis DENSITY data for RR calculation ----

low.dis = dens.data[c(dens.data$site.id != 'ML' & dens.data$site.id != 'ND' & dens.data$site.id != 'SMP' & dens.data$site.id != 'SK'),]
low.dis.comb = low.dis[,-4] %>% group_by(site.id, res, rep) %>% summarise_all(sum)

low.dis.comb.dens = low.dis.comb[,-c(1:4)] %>% mutate_all(funs(./low.dis.comb$area)) #Calculate the density (fish/m2) by dividing species counts by area 
low.dens.dat = cbind(low.dis[low.dis$obs == 'aak', c(1:3)], low.dis.comb.dens)

#:::::: COMBINE the low.dis and high.dis spec dens data ----
comb.dens.data = rbind(high.dens.dat, low.dens.dat) # Densities are in fish/m^2, high discharge sites are weighted densities
comb.dens.data = comb.dens.data[order(comb.dens.data$site.id),] # Reorder the variables alphabetically by Reserve

#write.table(comb.dens.data, file = 'data/Rd/comb_dens_data.txt') # This will write the table for weighted densities (only 4 large sites) to be used in other analyses

all.fish.dens = comb.dens.data[,c('site.id', 'res', 'rep', 'total.fish')] # get a dataframe of just categorical variables and total.fish density in fish/m^2

all.fish.dens.w = all.fish.dens %>% group_by(site.id) %>% spread(res, total.fish) # spreads Reserve column to Y and N columns with density average for observers at each rep (2 per site.id)

all.fish.dens.w.means = all.fish.dens.w %>% 
  group_by(site.id) %>% 
  summarise(N.mean = mean(N), N.sd = sd(N), Y.mean = mean(Y), Y.sd = sd(Y)) # creates dataframe by Village that includes mean and sd of Y and N columns

all.dens.rr = all.fish.dens.w.means %>% 
  group_by(site.id) %>% 
  mutate(rr.dens = log10(Y.mean/N.mean), se.rr.dens = sqrt(((Y.sd)^2/(2*(Y.mean)^2) + (N.sd)^2/(2*(N.mean)^2))), perc.change.dens = ((Y.mean-N.mean)/N.mean)*100) # creates new dataframe for logRR and for difference

# combine all.dens.rr and res.data; checked alignment of site.ids
###
all.dens.rr.data = cbind(as.data.frame(all.dens.rr), res.data)
###



plot(rr.dens~yrs.prot, data = all.dens.rr.data)
ab = lm(rr.dens~yrs.prot, data = all.dens.rr.data)
summary(ab)

# Create the omnibus model for all Rd
all.dens.rr.lm = lm(rr.dens ~ yrs.prot + houses + pen + mouth.dist + discharge + clo.res.dist 
                    + bet.cent.n + res.area + rd.dist + vil.dist, 
                     data = all.dens.rr.data, na.action = na.pass) # 

summary(all.dens.rr.lm)

# z-score standardizes all predictor variables to account for non-normality
all.dens.rr.lm.std = standardize(all.dens.rr.lm, standardize.y = FALSE) 

# dredge function from MuMIn package calculates all possible models and evaluates fit
all.dens.rr.model.set = dredge(all.dens.rr.lm.std)
sum(all.dens.rr.model.set$weight[1:306]) # we use this to calculate the threshold for inclusion in model averaging (ie summed weights ≤ 0.95)
all.dens.rr.95.models = get.models(all.dens.rr.model.set, 1:306)

all.dens.rr.95.models.avg = model.avg(all.dens.rr.95.models, revised.var =TRUE) 
summary(all.dens.rr.95.models.avg) # Read the model averaged output

# Create a barplot of top models by importance
all.dens.rr.impt = as.matrix(all.dens.rr.95.models.avg$sw) # extracts the Akaike weights for each term in the top models
all.dens.rr.coef = as.matrix(all.dens.rr.95.models.avg$coefficients[2,-1]) # extracts the coefficients from the averaged model for coloring by positive or negative

all.dens.rr.bar = merge(all.dens.rr.impt, all.dens.rr.coef, by = 'row.names', all.x = TRUE) # combines the Akaike weights data and coefficients 
row.names(all.dens.rr.bar) = all.dens.rr.bar$Row.names
all.dens.rr.bar = all.dens.rr.bar[,-1]

mod.pred = matrix(rep(NA, times = 10), dimnames = list(c('z.yrs.prot', 'z.houses', 'c.pen', 'z.mouth.dist', 'z.discharge', 'z.clo.res.dist', 'z.bet.cent.n', 'z.res.area', 'z.rd.dist', 'z.vil.dist')))
mod.pred.sort = as.matrix(mod.pred[order(row.names(mod.pred)),])

all.dens.rr.all = merge(mod.pred.sort, all.dens.rr.bar, by = 'row.names', all.x = TRUE) # combines the Akaike weights data and coefficients 
all.dens.rr.all = all.dens.rr.all[,-2]
all.dens.rr.all[is.na(all.dens.rr.all)] <- 0
colnames(all.dens.rr.all) = c('p.variable', 'importance', 'coefficient')

target = c('z.res.area', 'z.discharge', 'z.yrs.prot', 'c.pen', 'z.vil.dist', 'z.rd.dist', 'z.houses', 'z.bet.cent.n', 'z.mouth.dist', 'z.clo.res.dist')
all.dens.rr.plot = all.dens.rr.all[match(target, all.dens.rr.all$p.variable),]

par(mfrow = c(1,6))
par(oma = c(5,3,5,3))

all.plot = barplot(all.dens.rr.plot$importance, names.arg = all.dens.rr.plot$p.variable, 
                   col = ifelse(all.dens.rr.plot$coefficient > 0, 'black', 'white'), 
                   las = 2, main = 'All Fishes',
                   ylim = c(0,1)) 

mtext('Importance', side = 2, line = 3, cex = 1.2)

# Get data for partial residual plots
all.dens.full = all.dens.rr.lm.std
all.dens.avg = all.dens.rr.95.models.avg
all.dens.mod = all.dens.full
coefMod = coef(all.dens.avg)

summary(all.dens.avg)
summary(all.dens.full)
# Changing the coefficents to the averaged coefficient # see https://stackoverflow.com/questions/28633537/partial-residual-plot-based-on-model-average-coefficients-in-r
all.dens.mod$coefficients <- coefMod[names(coef(all.dens.full))]
all.dens.mod$residuals <- all.dens.rr.data$rr.dens - predict(all.dens.mod)
all.dens.pres <- resid(all.dens.mod, type='partial')

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


# Big fish Rd analysis

big.fish.dens = rowSums(comb.dens.data[ , big.fish.id], na.rm = TRUE) # Pull counts of only the big fish
big.dens.data = cbind(comb.dens.data[,c(1:3)], big.fish.dens) # add back in the categorical variables, and area and size

big.fish.dens.w = big.dens.data %>% group_by(site.id) %>% spread(res, big.fish.dens) # spreads reserve column to Y and N columns with density average for observers at each rep (2 per site.id)

big.fish.dens.w.means = big.fish.dens.w %>% 
  group_by(site.id) %>% 
  summarise(N.mean = mean(N), N.sd = sd(N), Y.mean = mean(Y), Y.sd = sd(Y)) # creates dataframe by site.id that includes mean and sd of Y and N columns


big.dens.rr = big.fish.dens.w.means %>% 
  group_by(site.id) %>% 
  mutate(rr.dens = log10(Y.mean/N.mean), se.rr.dens = sqrt(((Y.sd)^2/(2*(Y.mean)^2) + (N.sd)^2/(2*(N.mean)^2))), perc.change.dens = ((Y.mean-N.mean)/N.mean)*100) # creates new dataframe for logRR and for difference

# combine big.dens.rr and res.data; checked alignment of site.ids
###
big.dens.rr.data = cbind(as.data.frame(big.dens.rr), res.data)
###


# Create the omnibus model for big Rd
big.dens.rr.lm = lm(rr.dens ~ yrs.prot + houses + pen + mouth.dist + discharge + clo.res.dist 
                    + bet.cent.n + res.area + rd.dist + vil.dist, 
                    data = big.dens.rr.data, na.action = na.pass) # 

summary(big.dens.rr.lm)

# z-score standardizes big predictor variables to account for non-normality
big.dens.rr.lm.std = standardize(big.dens.rr.lm, standardize.y = FALSE) 

# dredge function from MuMIn package calculates big possible models and evaluates fit
big.dens.rr.model.set = dredge(big.dens.rr.lm.std)
sum(big.dens.rr.model.set$weight[1:238]) # we use this to calculate the threshold for inclusion in model averaging (ie summed weights ≤ 0.95)
big.dens.rr.95.models = get.models(big.dens.rr.model.set, 1:238)

big.dens.rr.95.models.avg = model.avg(big.dens.rr.95.models, revised.var =TRUE) 
summary(big.dens.rr.95.models.avg) # Read the model averaged output

# Create a barplot of top models by importance

big.dens.rr.impt = as.matrix(big.dens.rr.95.models.avg$sw) # extracts the Akaike weights for each term in the top models
big.dens.rr.coef = as.matrix(big.dens.rr.95.models.avg$coefficients[2,-1]) # extracts the coefficients from the averaged model for coloring by positive or negative

big.dens.rr.bar = merge(big.dens.rr.impt, big.dens.rr.coef, by = 'row.names', big.x = TRUE) # combines the Akaike weights data and coefficients 
row.names(big.dens.rr.bar) = big.dens.rr.bar$Row.names
big.dens.rr.bar = big.dens.rr.bar[,-1]

mod.pred = matrix(rep(NA, times = 10), dimnames = list(c('z.yrs.prot', 'z.houses', 'c.pen', 'z.mouth.dist', 'z.discharge', 'z.clo.res.dist', 'z.bet.cent.n', 'z.res.area', 'z.rd.dist', 'z.vil.dist')))
mod.pred.sort = as.matrix(mod.pred[order(row.names(mod.pred)),])

big.dens.rr.big = merge(mod.pred.sort, big.dens.rr.bar, by = 'row.names', big.x = TRUE) # combines the Akaike weights data and coefficients 
big.dens.rr.big = big.dens.rr.big[,-2]
big.dens.rr.big[is.na(big.dens.rr.big)] <- 0
colnames(big.dens.rr.big) = c('p.variable', 'importance', 'coefficient')

target = c('z.res.area', 'z.discharge', 'z.yrs.prot', 'c.pen', 'z.vil.dist', 'z.rd.dist', 'z.houses', 'z.bet.cent.n', 'z.mouth.dist', 'z.clo.res.dist')
big.dens.rr.plot = big.dens.rr.big[match(target, big.dens.rr.big$p.variable),]

big.plot = barplot(big.dens.rr.plot$importance, names.arg = big.dens.rr.plot$p.variable, 
                   col = ifelse(big.dens.rr.plot$coefficient > 0, 'black', 'white'), 
                   las = 2, main = 'Big Fishes',
                   ylim = c(0,1)) 

big.dens.full = big.dens.rr.lm.std
big.dens.avg = big.dens.rr.95.models.avg
big.dens.mod = big.dens.full
coefMod = coef(big.dens.avg)

summary(big.dens.avg)
summary(big.dens.full)
# Changing the coefficents to the averaged coefficient # see https://stackoverflow.com/questions/28633537/partial-residual-plot-based-on-model-average-coefficients-in-r
big.dens.mod$coefficients <- coefMod[names(coef(big.dens.full))]
big.dens.mod$residuals <- big.dens.rr.data$rr.dens - predict(big.dens.mod)
big.dens.pres <- resid(big.dens.mod, type='partial')


# Small Fish Analysis

lil.fish.dens = rowSums(comb.dens.data[ , lil.fish.id], na.rm = TRUE) # Pull counts of only the lil fish
lil.dens.data = cbind(comb.dens.data[,c(1:3)], lil.fish.dens) # add back in the categorical variables, and area and size

lil.fish.dens.w = lil.dens.data %>% group_by(site.id) %>% spread(res, lil.fish.dens) # spreads reserve column to Y and N columns with density average for observers at each rep (2 per site.id)

lil.fish.dens.w.means = lil.fish.dens.w %>% 
  group_by(site.id) %>% 
  summarise(N.mean = mean(N), N.sd = sd(N), Y.mean = mean(Y), Y.sd = sd(Y)) # creates dataframe by site.id that includes mean and sd of Y and N columns


lil.dens.rr = lil.fish.dens.w.means %>% 
  group_by(site.id) %>% 
  mutate(rr.dens = log10(Y.mean/N.mean), se.rr.dens = sqrt(((Y.sd)^2/(2*(Y.mean)^2) + (N.sd)^2/(2*(N.mean)^2))), perc.change.dens = ((Y.mean-N.mean)/N.mean)*100) # creates new dataframe for logRR and for difference

# combine lil.dens.rr and res.data; checked alignment of site.ids
###
lil.dens.rr.data = cbind(as.data.frame(lil.dens.rr), res.data)
###

# Create the omnibus model for lil Rd
lil.dens.rr.lm = lm(rr.dens ~ yrs.prot + houses + pen + mouth.dist + discharge + clo.res.dist 
                    + bet.cent.n + res.area + rd.dist + vil.dist, 
                    data = lil.dens.rr.data, na.action = na.pass) # 

summary(lil.dens.rr.lm)

# z-score standardizes lil predictor variables to account for non-normality
lil.dens.rr.lm.std = standardize(lil.dens.rr.lm, standardize.y = FALSE) 

# dredge function from MuMIn package calculates lil possible models and evaluates fit
lil.dens.rr.model.set = dredge(lil.dens.rr.lm.std)
sum(lil.dens.rr.model.set$weight[1:212]) # we use this to calculate the threshold for inclusion in model averaging (ie summed weights ≤ 0.95)
lil.dens.rr.95.models = get.models(lil.dens.rr.model.set, 1:212)

lil.dens.rr.95.models.avg = model.avg(lil.dens.rr.95.models, revised.var =TRUE) 
summary(lil.dens.rr.95.models.avg) # Read the model averaged output

# Create a barplot of top models by importance

lil.dens.rr.impt = as.matrix(lil.dens.rr.95.models.avg$sw) # extracts the Akaike weights for each term in the top models
lil.dens.rr.coef = as.matrix(lil.dens.rr.95.models.avg$coefficients[2,-1]) # extracts the coefficients from the averaged model for coloring by positive or negative

lil.dens.rr.bar = merge(lil.dens.rr.impt, lil.dens.rr.coef, by = 'row.names', lil.x = TRUE) # combines the Akaike weights data and coefficients 
row.names(lil.dens.rr.bar) = lil.dens.rr.bar$Row.names
lil.dens.rr.bar = lil.dens.rr.bar[,-1]

mod.pred = matrix(rep(NA, times = 10), dimnames = list(c('z.yrs.prot', 'z.houses', 'c.pen', 'z.mouth.dist', 'z.discharge', 'z.clo.res.dist', 'z.bet.cent.n', 'z.res.area', 'z.rd.dist', 'z.vil.dist')))
mod.pred.sort = as.matrix(mod.pred[order(row.names(mod.pred)),])

lil.dens.rr.lil = merge(mod.pred.sort, lil.dens.rr.bar, by = 'row.names', lil.x = TRUE) # combines the Akaike weights data and coefficients 
lil.dens.rr.lil = lil.dens.rr.lil[,-2]
lil.dens.rr.lil[is.na(lil.dens.rr.lil)] <- 0
colnames(lil.dens.rr.lil) = c('p.variable', 'importance', 'coefficient')

target = c('z.res.area', 'z.discharge', 'z.yrs.prot', 'c.pen', 'z.vil.dist', 'z.rd.dist', 'z.houses', 'z.bet.cent.n', 'z.mouth.dist', 'z.clo.res.dist')
lil.dens.rr.plot = lil.dens.rr.lil[match(target, lil.dens.rr.lil$p.variable),]

lil.plot = barplot(lil.dens.rr.plot$importance, names.arg = lil.dens.rr.plot$p.variable, 
                   col = ifelse(lil.dens.rr.plot$coefficient > 0, 'black', 'white'), 
                   las = 2, main = 'Small Fishes',
                   ylim = c(0,1)) 

lil.dens.full = lil.dens.rr.lm.std
lil.dens.avg = lil.dens.rr.95.models.avg
lil.dens.mod = lil.dens.full
coefMod = coef(lil.dens.avg)

summary(lil.dens.avg)
summary(lil.dens.full)
# Changing the coefficents to the averaged coefficient # see https://stackoverflow.com/questions/28633537/partial-residual-plot-based-on-model-average-coefficients-in-r
lil.dens.mod$coefficients <- coefMod[names(coef(lil.dens.full))]
# Changing the residuals
lil.dens.mod$residuals <- lil.dens.rr.data$rr.dens - predict(lil.dens.mod)
lil.dens.pres <- resid(lil.dens.mod, type='partial')


# Divide fish into 3 trophic categories: 
## 'pred.fish.id' == those ids for species having trophic position > 3.5
## 'omni.fish.id' == those ids for species having trophic position <=3.5 & > 2.5
## 'herb.fish.id' == those ids for species having trophic position <= 2.5
# *** Leave YOY.un out due to unknown trophic position

pred.fish.id = c('Ch.bu', 'Ch.st', 'Ha.sa', 'He.mi', 'No.no', 'Ra.gu', 'Sp.ac', 'Xe.ca')

omni.fish.id = c('Ga.na', 'Nema', 'Ma.ar', 'Hy.sa', 'Ne.st', 'To.sp', 'Ch.ba', 'Pe.st', 'Cr.bu' , 'Da.pu', 
                 'De.sp', 'My.ar', 'Ba.or', 'Pa.vo', 'Sy.ru', 'Fo.br', 'Gl.sp', 'Po.sp', 'Pu.gr', 'Pu.re', 'Ra.da', 'YOY.ns')

herb.fish.id = c('Ba.de', 'Sc.bu', 'Ho.bu', 'La.ro')

# Predatory Fish Analysis

pred.fish.dens = rowSums(comb.dens.data[ , pred.fish.id], na.rm = TRUE) # Pull counts of only the pred fish
pred.dens.data = cbind(comb.dens.data[,c(1:3)], pred.fish.dens) # add back in the categorical variables, and area and size

pred.fish.dens.w = pred.dens.data %>% group_by(site.id) %>% spread(res, pred.fish.dens) # spreads reserve column to Y and N columns with density average for observers at each rep (2 per site.id)

pred.fish.dens.w.means = pred.fish.dens.w %>% 
  group_by(site.id) %>% 
  summarise(N.mean = mean(N + 0.001), N.sd = sd(N), Y.mean = mean(Y + 0.001), Y.sd = sd(Y)) # creates dataframe by site.id that includes mean and sd of Y and N columns


pred.dens.rr = pred.fish.dens.w.means %>% 
  group_by(site.id) %>% 
  mutate(rr.dens = log10(Y.mean/N.mean), se.rr.dens = sqrt(((Y.sd)^2/(2*(Y.mean)^2) + (N.sd)^2/(2*(N.mean)^2))), perc.change.dens = ((Y.mean-N.mean)/N.mean)*100) # creates new dataframe for logRR and for difference

# combine pred.dens.rr and res.data; checked alignment of site.ids
###
pred.dens.rr.data = cbind(as.data.frame(pred.dens.rr), res.data)
###

# Create the predbus model for pred Rd
pred.dens.rr.lm = lm(rr.dens ~ yrs.prot + houses + pen + mouth.dist + discharge + clo.res.dist 
                    + bet.cent.n + res.area + rd.dist + vil.dist, 
                    data = pred.dens.rr.data, na.action = na.pass) # 

summary(pred.dens.rr.lm)

# z-score standardizes pred predictor variables to account for non-normality
pred.dens.rr.lm.std = standardize(pred.dens.rr.lm, standardize.y = FALSE) 

# dredge function from MuMIn package calculates pred possible models and evaluates fit
pred.dens.rr.model.set = dredge(pred.dens.rr.lm.std)
sum(pred.dens.rr.model.set$weight[1:317]) # we use this to calculate the threshold for inclusion in model averaging (ie summed weights ≤ 0.95)
pred.dens.rr.95.models = get.models(pred.dens.rr.model.set, 1:317)

pred.dens.rr.95.models.avg = model.avg(pred.dens.rr.95.models, revised.var =TRUE) 
summary(pred.dens.rr.95.models.avg) # Read the model averaged output

# Create a barplot of top models by importance

pred.dens.rr.impt = as.matrix(pred.dens.rr.95.models.avg$sw) # extracts the Akaike weights for each term in the top models
pred.dens.rr.coef = as.matrix(pred.dens.rr.95.models.avg$coefficients[2,-1]) # extracts the coefficients from the averaged model for coloring by positive or negative

pred.dens.rr.bar = merge(pred.dens.rr.impt, pred.dens.rr.coef, by = 'row.names', pred.x = TRUE) # combines the Akaike weights data and coefficients 
row.names(pred.dens.rr.bar) = pred.dens.rr.bar$Row.names
pred.dens.rr.bar = pred.dens.rr.bar[,-1]

mod.pred = matrix(rep(NA, times = 10), dimnames = list(c('z.yrs.prot', 'z.houses', 'c.pen', 'z.mouth.dist', 'z.discharge', 'z.clo.res.dist', 'z.bet.cent.n', 'z.res.area', 'z.rd.dist', 'z.vil.dist')))
mod.pred.sort = as.matrix(mod.pred[order(row.names(mod.pred)),])

pred.dens.rr.pred = merge(mod.pred.sort, pred.dens.rr.bar, by = 'row.names', pred.x = TRUE) # combines the Akaike weights data and coefficients 
pred.dens.rr.pred = pred.dens.rr.pred[,-2]
pred.dens.rr.pred[is.na(pred.dens.rr.pred)] <- 0
colnames(pred.dens.rr.pred) = c('p.variable', 'importance', 'coefficient')

target = c('z.res.area', 'z.discharge', 'z.yrs.prot', 'c.pen', 'z.vil.dist', 'z.rd.dist', 'z.houses', 'z.bet.cent.n', 'z.mouth.dist', 'z.clo.res.dist')
pred.dens.rr.plot = pred.dens.rr.pred[match(target, pred.dens.rr.pred$p.variable),]

pred.plot = barplot(pred.dens.rr.plot$importance, names.arg = pred.dens.rr.plot$p.variable, 
                   col = ifelse(pred.dens.rr.plot$coefficient > 0, 'black', 'white'), 
                   las = 2, main = 'Predators',
                   ylim = c(0,1)) 

pred.dens.full = pred.dens.rr.lm.std
pred.dens.avg = pred.dens.rr.95.models.avg
pred.dens.mod = pred.dens.full
coefMod = coef(pred.dens.avg)

summary(pred.dens.avg)
summary(pred.dens.full)
# Changing the coefficents to the averaged coefficient # see https://stackoverflow.com/questions/28633537/partial-residual-plot-based-on-model-average-coefficients-in-r
pred.dens.mod$coefficients <- coefMod[names(coef(pred.dens.full))]
pred.dens.mod$residuals <- pred.dens.rr.data$rr.dens - predict(pred.dens.mod)
pred.dens.pres <- resid(pred.dens.mod, type='partial')


# Omnivorous Fish Analysis

omni.fish.dens = rowSums(comb.dens.data[ , omni.fish.id], na.rm = TRUE) # Pull counts of only the omni fish
omni.dens.data = cbind(comb.dens.data[,c(1:3)], omni.fish.dens) # add back in the categorical variables, and area and size

omni.fish.dens.w = omni.dens.data %>% group_by(site.id) %>% spread(res, omni.fish.dens) # spreads reserve column to Y and N columns with density average for observers at each rep (2 per site.id)

omni.fish.dens.w.means = omni.fish.dens.w %>% 
  group_by(site.id) %>% 
  summarise(N.mean = mean(N), N.sd = sd(N), Y.mean = mean(Y), Y.sd = sd(Y)) # creates dataframe by site.id that includes mean and sd of Y and N columns


omni.dens.rr = omni.fish.dens.w.means %>% 
  group_by(site.id) %>% 
  mutate(rr.dens = log10(Y.mean/N.mean), se.rr.dens = sqrt(((Y.sd)^2/(2*(Y.mean)^2) + (N.sd)^2/(2*(N.mean)^2))), perc.change.dens = ((Y.mean-N.mean)/N.mean)*100) # creates new dataframe for logRR and for difference

# combine omni.dens.rr and res.data; checked alignment of site.ids
###
omni.dens.rr.data = cbind(as.data.frame(omni.dens.rr), res.data)
###

# Create the omnibus model for omni Rd
omni.dens.rr.lm = lm(rr.dens ~ yrs.prot + houses + pen + mouth.dist + discharge + clo.res.dist 
                    + bet.cent.n + res.area + rd.dist + vil.dist, 
                    data = omni.dens.rr.data, na.action = na.pass) # 

summary(omni.dens.rr.lm)

# z-score standardizes omni predictor variables to account for non-normality
omni.dens.rr.lm.std = standardize(omni.dens.rr.lm, standardize.y = FALSE) 

# dredge function from MuMIn package calculates omni possible models and evaluates fit
omni.dens.rr.model.set = dredge(omni.dens.rr.lm.std)
sum(omni.dens.rr.model.set$weight[1:340]) # we use this to calculate the threshold for inclusion in model averaging (ie summed weights ≤ 0.95)
omni.dens.rr.95.models = get.models(omni.dens.rr.model.set, 1:340)

omni.dens.rr.95.models.avg = model.avg(omni.dens.rr.95.models, revised.var =TRUE) 
summary(omni.dens.rr.95.models.avg) # Read the model averaged output

# Create a barplot of top models by importance

omni.dens.rr.impt = as.matrix(omni.dens.rr.95.models.avg$sw) # extracts the Akaike weights for each term in the top models
omni.dens.rr.coef = as.matrix(omni.dens.rr.95.models.avg$coefficients[2,-1]) # extracts the coefficients from the averaged model for coloring by positive or negative

omni.dens.rr.bar = merge(omni.dens.rr.impt, omni.dens.rr.coef, by = 'row.names', omni.x = TRUE) # combines the Akaike weights data and coefficients 
row.names(omni.dens.rr.bar) = omni.dens.rr.bar$Row.names
omni.dens.rr.bar = omni.dens.rr.bar[,-1]

mod.pred = matrix(rep(NA, times = 10), dimnames = list(c('z.yrs.prot', 'z.houses', 'c.pen', 'z.mouth.dist', 'z.discharge', 'z.clo.res.dist', 'z.bet.cent.n', 'z.res.area', 'z.rd.dist', 'z.vil.dist')))
mod.pred.sort = as.matrix(mod.pred[order(row.names(mod.pred)),])

omni.dens.rr.omni = merge(mod.pred.sort, omni.dens.rr.bar, by = 'row.names', omni.x = TRUE) # combines the Akaike weights data and coefficients 
omni.dens.rr.omni = omni.dens.rr.omni[,-2]
omni.dens.rr.omni[is.na(omni.dens.rr.omni)] <- 0
colnames(omni.dens.rr.omni) = c('p.variable', 'importance', 'coefficient')

target = c('z.res.area', 'z.discharge', 'z.yrs.prot', 'c.pen', 'z.vil.dist', 'z.rd.dist', 'z.houses', 'z.bet.cent.n', 'z.mouth.dist', 'z.clo.res.dist')
omni.dens.rr.plot = omni.dens.rr.omni[match(target, omni.dens.rr.omni$p.variable),]

omni.plot = barplot(omni.dens.rr.plot$importance, names.arg = omni.dens.rr.plot$p.variable, 
                   col = ifelse(omni.dens.rr.plot$coefficient > 0, 'black', 'white'), 
                   las = 2, main = 'Omnivores',
                   ylim = c(0,1)) 

omni.dens.full = omni.dens.rr.lm.std
omni.dens.avg = omni.dens.rr.95.models.avg
omni.dens.mod = omni.dens.full
coefMod = coef(omni.dens.avg)

summary(omni.dens.avg)
summary(omni.dens.full)
# Changing the coefficents to the averaged coefficient # see https://stackoverflow.com/questions/28633537/partial-residual-plot-based-on-model-average-coefficients-in-r
omni.dens.mod$coefficients <- coefMod[names(coef(omni.dens.full))]
omni.dens.mod$residuals <- omni.dens.rr.data$rr.dens - predict(omni.dens.mod)
omni.dens.pres <- resid(omni.dens.mod, type='partial')

# Herbivorous Fish Analysis

herb.fish.dens = rowSums(comb.dens.data[ , herb.fish.id], na.rm = TRUE) # Pull counts of only the herb fish
herb.dens.data = cbind(comb.dens.data[,c(1:3)], herb.fish.dens) # add back in the categorical variables, and area and size

herb.fish.dens.w = herb.dens.data %>% group_by(site.id) %>% spread(res, herb.fish.dens) # spreads reserve column to Y and N columns with density average for observers at each rep (2 per site.id)

herb.fish.dens.w.means = herb.fish.dens.w %>% 
  group_by(site.id) %>% 
  summarise(N.mean = mean(N + 0.001), N.sd = sd(N), Y.mean = mean(Y + 0.001), Y.sd = sd(Y)) # creates dataframe by site.id that includes mean and sd of Y and N columns


herb.dens.rr = herb.fish.dens.w.means %>% 
  group_by(site.id) %>% 
  mutate(rr.dens = log10(Y.mean/N.mean), se.rr.dens = sqrt(((Y.sd)^2/(2*(Y.mean)^2) + (N.sd)^2/(2*(N.mean)^2))), perc.change.dens = ((Y.mean-N.mean)/N.mean)*100) # creates new dataframe for logRR and for difference

# combine herb.dens.rr and res.data; checked alignment of site.ids
###
herb.dens.rr.data = cbind(as.data.frame(herb.dens.rr), res.data)
###

# Create the herbbus model for herb Rd
herb.dens.rr.lm = lm(rr.dens ~ yrs.prot + houses + pen + mouth.dist + discharge + clo.res.dist 
                     + bet.cent.n + res.area + rd.dist + vil.dist, 
                     data = herb.dens.rr.data, na.action = na.pass) # 

summary(herb.dens.rr.lm)

# z-score standardizes herb herbictor variables to account for non-normality
herb.dens.rr.lm.std = standardize(herb.dens.rr.lm, standardize.y = FALSE) 

# dredge function from MuMIn package calculates herb possible models and evaluates fit
herb.dens.rr.model.set = dredge(herb.dens.rr.lm.std)
sum(herb.dens.rr.model.set$weight[1:349]) # we use this to calculate the threshold for inclusion in model averaging (ie summed weights ≤ 0.95)
herb.dens.rr.95.models = get.models(herb.dens.rr.model.set, 1:349)

herb.dens.rr.95.models.avg = model.avg(herb.dens.rr.95.models, revised.var =TRUE) 
summary(herb.dens.rr.95.models.avg) # Read the model averaged output

# Create a barplot of top models by importance

herb.dens.rr.impt = as.matrix(herb.dens.rr.95.models.avg$sw) # extracts the Akaike weights for each term in the top models
herb.dens.rr.coef = as.matrix(herb.dens.rr.95.models.avg$coefficients[2,-1]) # extracts the coefficients from the averaged model for coloring by positive or negative

herb.dens.rr.bar = merge(herb.dens.rr.impt, herb.dens.rr.coef, by = 'row.names', herb.x = TRUE) # combines the Akaike weights data and coefficients 
row.names(herb.dens.rr.bar) = herb.dens.rr.bar$Row.names
herb.dens.rr.bar = herb.dens.rr.bar[,-1]

mod.herb = matrix(rep(NA, times = 10), dimnames = list(c('z.yrs.prot', 'z.houses', 'c.pen', 'z.mouth.dist', 'z.discharge', 'z.clo.res.dist', 'z.bet.cent.n', 'z.res.area', 'z.rd.dist', 'z.vil.dist')))
mod.herb.sort = as.matrix(mod.herb[order(row.names(mod.herb)),])

herb.dens.rr.herb = merge(mod.herb.sort, herb.dens.rr.bar, by = 'row.names', herb.x = TRUE) # combines the Akaike weights data and coefficients 
herb.dens.rr.herb = herb.dens.rr.herb[,-2]
herb.dens.rr.herb[is.na(herb.dens.rr.herb)] <- 0
colnames(herb.dens.rr.herb) = c('p.variable', 'importance', 'coefficient')

target = c('z.res.area', 'z.discharge', 'z.yrs.prot', 'c.pen', 'z.vil.dist', 'z.rd.dist', 'z.houses', 'z.bet.cent.n', 'z.mouth.dist', 'z.clo.res.dist')
herb.dens.rr.plot = herb.dens.rr.herb[match(target, herb.dens.rr.herb$p.variable),]

herb.plot = barplot(herb.dens.rr.plot$importance, names.arg = herb.dens.rr.plot$p.variable, 
                    col = ifelse(herb.dens.rr.plot$coefficient > 0, 'black', 'white'), 
                    las = 2, main = 'Herbivores',
                    ylim = c(0,1)) 


herb.dens.full = herb.dens.rr.lm.std
herb.dens.avg = herb.dens.rr.95.models.avg
herb.dens.mod = herb.dens.full
coefMod = coef(herb.dens.avg)

summary(herb.dens.avg)
summary(herb.dens.full)
# Changing the coefficents to the averaged coefficient # see https://stackoverflow.com/questions/28633537/partial-residual-plot-based-on-model-average-coefficients-in-r
herb.dens.mod$coefficients <- coefMod[names(coef(herb.dens.full))]
herb.dens.mod$residuals <- herb.dens.rr.data$rr.dens - predict(herb.dens.mod)
herb.dens.pres <- resid(herb.dens.mod, type='partial')

# All averaged model summaries

summary(all.dens.avg) 
summary(big.dens.avg) 
summary(lil.dens.avg) 
summary(pred.dens.avg)
summary(omni.dens.avg)
summary(herb.dens.avg)

#\\\\\\\\\\\\\\\\\\\\\\
#Partial residual plots (Fig. 4)
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
#pdf('figures/fig_4/single_plots/dens_vildist_resid_DB.pdf', width = 7, height = 6)
pdf('figures/fig_4/Figure_4_dens_resid_all.pdf', width = 15, height = 4)

par(mfrow = c(1,5))
#par(oma = c(2,6,2,2))
par(mfrow = c(1,5), oma = c(0,4,0,1))


plot(all.dens.pres[,"z.vil.dist"] ~ all.dens.rr.data$vil.dist, pch = all.pch, col = all.col, cex = 2, bg = all.col.a, bty = 'n',
     ylab = '', xlab = '', axes = FALSE, ylim = c(-2,2), xlim = c(0, 2000), type = 'n')

abline(h = 0, col = 'grey')
points(big.dens.pres[,"z.vil.dist"] ~ big.dens.rr.data$vil.dist, pch = all.pch, col = big.col, cex = 2, bg = big.col.a)
points(all.dens.pres[,"z.vil.dist"] ~ all.dens.rr.data$vil.dist, pch = all.pch, col = all.col, cex = 2, bg = all.col.a)
abline(lm(big.dens.pres[,"z.vil.dist"]~big.dens.rr.data$vil.dist), lty = 1, col = big.col, lwd = 2) # p < 0.05
abline(lm(all.dens.pres[,"z.vil.dist"]~all.dens.rr.data$vil.dist), lty = 1, col = all.col, lwd = 2) # p < 0.05
axis(side = 1, at = c(0, 1000, 2000), labels = c('0','1','2'), cex.axis = 1.5, col.axis = axis.col, col = axis.col)
axis(side = 2, at = c(-2,0,2), las = 1, cex.axis = 1.5, col.axis = axis.col, col = axis.col)
mtext(side = 1, 'Distance to Village (km)', line = 2.5, cex = 1.2, col  = axis.col)
mtext(side = 2, 'Partial Residuals', line = 3, cex = 1.2, col = axis.col)
mtext(expression(paste('D'[italic(r)])), side = 2, line = 5.5, las = 1, cex = 1.7, col = axis.col)


#dev.off()

#Penalty for lil fish
#dev.new = TRUE
#pdf('figures/fig_4/single_plots/dens_pen_resid_DB.pdf', width = 7, height = 6)
box.col = c('gray', lil.col.a)
boxplot(lil.dens.pres[,'c.pen']~lil.dens.rr.data$pen, col = box.col, bty = 'n', outcol = box.col, outpch = 21, outbg = lil.col, outcex = 2, ylab = '', xlab = '', axes = FALSE, ylim = c(-1,1))

axis(side = 2, at = c(-1,0,1), las = 1, cex.axis = 1.5, col.axis = axis.col, col = axis.col)
mtext(side = 1, 'No', at = 1, line = 1.2, col = axis.col, cex = 1.2)
mtext(side = 1, 'Yes', at = 2, line = 1.2, col = axis.col, cex = 1.2)
mtext(side = 1, 'Penalty', at = 1.5, line = 3, cex = 1.2, col = axis.col)

#dev.off()

# Distance to Mouth

plot(herb.dens.pres[,"z.mouth.dist"] ~ herb.dens.rr.data$mouth.dist, pch = all.pch, col = herb.col, cex = 2, bg = herb.col.a, bty = 'n',
     ylab = '', xlab = '', axes = FALSE, ylim = c(-1.5,1.5), xlim = c(10000, 70000), type = 'n')

abline(h = 0, col = 'grey')
points(herb.dens.pres[,"z.mouth.dist"] ~ herb.dens.rr.data$mouth.dist, pch = all.pch, col = herb.col, cex = 2, bg = herb.col.a)

abline(lm(herb.dens.pres[,"z.mouth.dist"]~herb.dens.rr.data$mouth.dist), lty = 2, col = herb.col, lwd = 1.5) # p < 0.05

axis(side = 1, at = seq(10000, 70000, 30000), labels = c('10','40','70'), cex.axis = 1.5)
axis(side = 2, at = c(-1.5,0,1.5), labels = c('-1.5','0', '1.5'), las = 1, cex.axis = 1.5)
mtext(side = 1, 'Distance to Mouth (km)', line = 2.5, cex = 1.2, col  = 'black')

dev.off()


## Get data for heatmap; Fig 3

all.dens.rr.plot$fish.cat = rep('All', length(all.dens.rr.plot$coefficient))
big.dens.rr.plot$fish.cat = rep('Large', length(big.dens.rr.plot$coefficient))
lil.dens.rr.plot$fish.cat = rep('Small', length(lil.dens.rr.plot$coefficient))
pred.dens.rr.plot$fish.cat = rep('Predators', length(pred.dens.rr.plot$coefficient))
omni.dens.rr.plot$fish.cat = rep('Omnivores', length(omni.dens.rr.plot$coefficient))
herb.dens.rr.plot$fish.cat = rep('Herbivores', length(herb.dens.rr.plot$coefficient))

dens.heat.data = as.data.frame(rbind(all.dens.rr.plot, big.dens.rr.plot, lil.dens.rr.plot, pred.dens.rr.plot, omni.dens.rr.plot, herb.dens.rr.plot))

attach(dens.heat.data)
dens.heat.data$imp.sign <- with(dens.heat.data, ifelse(coefficient <= 0, -importance, importance)) # this adds a new column which makes the sign of the importance match the coefficient for color ramp
dens.heat.wide = as.data.frame(dens.heat.data[,-c(2:3)] %>% group_by(p.variable) %>% spread(fish.cat, imp.sign)) # This removes the importance and coefficient columns, combined into imp.sign, then makes a wide dataframe
rownames(dens.heat.wide) = dens.heat.wide[,1]
dens.heat.wide = dens.heat.wide[,-1]

# need to reorder the variables here
target = rev(c('z.res.area', 'z.discharge', 'z.yrs.prot', 'c.pen', 'z.vil.dist', 'z.rd.dist', 'z.houses', 'z.bet.cent.n', 'z.mouth.dist', 'z.clo.res.dist'))
dens.heat.wide.z = dens.heat.wide[match(target, rownames(dens.heat.wide)),]
dens.heat.wide.z = dens.heat.wide.z[c('All', 'Large', 'Small', 'Predators', 'Omnivores', 'Herbivores')]

write.table(dens.heat.wide.z, 'figures/fig_3/data/dens.heatmap.data.txt')
