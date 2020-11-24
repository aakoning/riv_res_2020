
# This analysis measures the effect of reserve characteristics on reserve richness response (Rr)

#rm(list = ls())

library(tidyr)
library(dplyr)
library(arm) # This allows for the standardize function used to standardize coefficients prior to model selection
library(MuMIn) # this package allows for model selection
library(car)

# Set Working Directory

rich.data = read.delim('data/Rs/rep_rich_data.txt', header = TRUE, sep = ' ') # This table is written in richness_analysis_final.R file

# Need to calculate the response ratio (Rs_reserve - Rs_non-reserve)

rep.rich.wide = rich.data[,-4] %>% group_by(site.id) %>% pivot_wider(names_from = res, values_from = rep.rich) # spreads Reserve column to Y and N columns
write.table(rep.rich.wide, 'data/Rs/rep_rich_wide_final.txt')

site.rich.mean.data = rep.rich.wide %>% # This takes the mean of the richness observations at each reserve and non reserve (Average of richness at each rep)
  group_by(site.id) %>% 
  summarise(N.mean = mean(N), N.sd = sd(N), Y.mean = mean(Y), Y.sd = sd(Y)) # creates dataframe by Village that includes mean and sd of Y and N columns

site.rich.mean.rr = site.rich.mean.data %>%
  group_by(site.id) %>%
  mutate(rr.rich = Y.mean - N.mean, perc.change.rich = ((Y.mean-N.mean)/N.mean)*100) # creates new dataframe for RR difference between the mean richness per site


#/////\\\\\\
site.rich.mean.rr
#///////||\\\\\\\\

# read in the variables for reserve features for predicting contribution to Richness Response (Rr).
res.data = read.delim('data/reserve_features.txt', header = TRUE)
all.rich.rr.data = cbind(as.data.frame(site.rich.mean.rr), res.data)

plot(rr.rich~yrs.prot, data = all.rich.rr.data)
ab = lm(rr.rich~yrs.prot, data = all.rich.rr.data)
summary(ab)

all.rich.diff.lm = lm(rr.rich ~ yrs.prot + houses + pen + mouth.dist + discharge + clo.res.dist + bet.cent.n + res.area + rd.dist + vil.dist, data = all.rich.rr.data, na.action = na.fail) # I checked for log10(transform) of discharge, area, yrs.prot, but descreased fit, and raised AIC by 13 pts.
summary(all.rich.diff.lm)

all.rich.diff.lm.std = standardize(all.rich.diff.lm, standardize.y = FALSE) # numeric values are rescaled to have mean = 0 and sd = 0.5 all predictor variables to account for non-normality
summary(all.rich.diff.lm.std)
all.rich.diff.model.set = dredge(all.rich.diff.lm.std, evaluate = TRUE, extra = 'R^2') # dredge function from MuMIn package calculates all possible factorial models and evaluates fit
sum(all.rich.diff.model.set$weight[1:295]) # This shows the summed Aikake weights are < 0.95
all.rich.95.models = get.models(all.rich.diff.model.set, 1:295) # just found the # models < .95 summed Aikake

all.rich.95.models.avg = model.avg(all.rich.95.models, revised.var = TRUE, fit = TRUE, extra = 'R^2') # conducts model averaging based on importance (summed Akaike weight for all models)
summary(all.rich.95.models.avg)

# Create a barplot of top models by importance
all.rich.95.models.avg$sw
all.rich.diff.impt = as.matrix(all.rich.95.models.avg$sw) # extracts the Akaike weights for each term in the top models
all.rich.diff.coef = as.matrix(all.rich.95.models.avg$coefficients[2,-1]) # extracts the coefficients from the averaged model for coloring by positive or negative

all.rich.diff.bar = merge(all.rich.diff.impt, all.rich.diff.coef, by = 'row.names', all.x = TRUE) # combines the Akaike weights data and coefficients 
row.names(all.rich.diff.bar) = all.rich.diff.bar$Row.names
all.rich.diff.bar = all.rich.diff.bar[,-1]

mod.pred = matrix(rep(NA, times = 10), dimnames = list(c('z.yrs.prot', 'z.houses', 'c.pen', 'z.mouth.dist', 'z.discharge', 'z.clo.res.dist', 'z.bet.cent.n', 'z.res.area', 'z.rd.dist', 'z.vil.dist')))
mod.pred.sort = as.matrix(mod.pred[order(row.names(mod.pred)),])

all.rich.diff.all = merge(mod.pred.sort, all.rich.diff.bar, by = 'row.names', all.x = TRUE) # combines the Akaike weights data and coefficients 
all.rich.diff.all = all.rich.diff.all[,-2]
all.rich.diff.all[is.na(all.rich.diff.all)] <- 0
colnames(all.rich.diff.all) = c('p.variable', 'importance', 'coefficient')

target = c('z.res.area', 'z.discharge', 'z.yrs.prot', 'c.pen', 'z.vil.dist', 'z.rd.dist', 'z.houses', 'z.bet.cent.n', 'z.mouth.dist', 'z.clo.res.dist')
all.rich.diff.plot = all.rich.diff.all[match(target, all.rich.diff.all$p.variable),]

par(mfrow = c(1,6))
par(oma = c(5,3,5,3))
all.plot = barplot(all.rich.diff.plot$importance, names.arg = all.rich.diff.plot$p.variable, 
                   col = ifelse(all.rich.diff.plot$coefficient > 0, 'black', 'white'), 
                   las = 2, main = 'All Fishes',
                   ylim = c(0,1)) 

mtext('Importance', side = 2, line = 3, cex = 1.2)

# Get data for partial residual plots
all.rich.full = all.rich.diff.lm.std
all.rich.avg = all.rich.95.models.avg
all.rich.mod = all.rich.full
coefMod = coef(all.rich.avg)


# Changing the coefficents to the averaged coefficient # see https://stackoverflow.com/questions/28633537/partial-residual-plot-based-on-model-average-coefficients-in-r
all.rich.mod$coefficients = coefMod[names(coef(all.rich.full))]
# Changing the residuals
all.rich.mod$residuals = all.rich.rr.data$rr.rich - predict(all.rich.mod)
all.rich.pres = resid(all.rich.mod, type='partial')

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# Analysis by functional traits ----
#//////////////////////////////


# get big fish richness data
big.rich.data = read.delim('data/Rs/big_rich_data.txt', header = TRUE, sep = ' ')

# Get the Richness response (Rs) for Big Fish 
big.rich.wide = big.rich.data[,-4] %>% group_by(site.id) %>% pivot_wider(names_from = res, values_from = big.rich)

big.rich.mean = big.rich.wide %>% # This takes the mean of the richness observations at each res and non res (Average of richness at each rep)
  group_by(site.id) %>% 
  summarise(N.mean = mean(N), N.sd = sd(N), Y.mean = mean(Y), Y.sd = sd(Y)) # creates dataframe by site.id that includes mean and sd of Y and N columns

big.rich.mean.rr = big.rich.mean %>%
  group_by(site.id) %>%
  mutate(rr.rich = Y.mean - N.mean) # creates new dataframe for RR difference between the mean richness per site

#////////////
# big.rich.mean.rr
#\\\\\\\\\\\\

big.rich.rr.data = cbind(res.data, as.data.frame(big.rich.mean.rr))

big.rich.diff.lm = lm(rr.rich ~ yrs.prot + houses + pen + mouth.dist + discharge + clo.res.dist + bet.cent.n + res.area + rd.dist + vil.dist, data = big.rich.rr.data, na.action = na.pass)
summary(big.rich.diff.lm)
big.rich.diff.lm.std = standardize(big.rich.diff.lm, standardize.y = FALSE) # numeric values are rescaled to have mean = 0 and sd = 0.5 big predictor variables to account for non-normality
big.rich.diff.model.set = dredge(big.rich.diff.lm.std) # dredge function from MuMIn package calculates big possible factorial models and evaluates fit
sum(big.rich.diff.model.set$weight[1:228]) # This shows the summed Aikake weights are < 0.95
big.rich.95.models = get.models(big.rich.diff.model.set, 1:228) # just found the # models < .95 summed Aikake weight using sum(big.rich.diff.model.set$weight[1:393])

big.rich.95.models.avg = model.avg(big.rich.95.models, revised.var = TRUE) # conducts model averaging based on importance (summed Akaike weight for big models)
summary(big.rich.95.models.avg)

# Create a barplot of top models by importance
big.rich.95.models.avg$sw
big.rich.diff.impt = as.matrix(big.rich.95.models.avg$sw) # extracts the Akaike weights for each term in the top models
big.rich.diff.coef = as.matrix(big.rich.95.models.avg$coefficients[2,-1]) # extracts the coefficients from the averaged model for coloring by positive or negative

big.rich.diff.bar = merge(big.rich.diff.impt, big.rich.diff.coef, by = 'row.names', big.x = TRUE) # combines the Akaike weights data and coefficients 
row.names(big.rich.diff.bar) = big.rich.diff.bar$Row.names
big.rich.diff.bar = big.rich.diff.bar[,-1]

mod.pred = matrix(rep(NA, times = 10), dimnames = list(c('z.yrs.prot', 'z.houses', 'c.pen', 'z.mouth.dist', 'z.discharge', 'z.clo.res.dist', 'z.bet.cent.n', 'z.res.area', 'z.rd.dist', 'z.vil.dist')))
mod.pred.sort = as.matrix(mod.pred[order(row.names(mod.pred)),])

big.rich.diff.big = merge(mod.pred.sort, big.rich.diff.bar, by = 'row.names', big.x = TRUE) # combines the Akaike weights data and coefficients 
big.rich.diff.big = big.rich.diff.big[,-2]
big.rich.diff.big[is.na(big.rich.diff.big)] <- 0
colnames(big.rich.diff.big) = c('p.variable', 'importance', 'coefficient')

target = c('z.res.area', 'z.discharge', 'z.yrs.prot', 'c.pen', 'z.vil.dist', 'z.rd.dist', 'z.houses', 'z.bet.cent.n', 'z.mouth.dist', 'z.clo.res.dist')
big.rich.diff.plot = big.rich.diff.big[match(target, big.rich.diff.big$p.variable),]

big.plot = barplot(big.rich.diff.plot$importance, names.arg = big.rich.diff.plot$p.variable, 
                   col = ifelse(big.rich.diff.plot$coefficient > 0, 'black', 'white'), 
                   las = 2, main = 'Big Fishes',
                   ylim = c(0,1)) 


# Get the partial residuals for partial residual plots (At end)
big.rich.full = big.rich.diff.lm.std
big.rich.avg = big.rich.95.models.avg
big.rich.mod = big.rich.full

coefMod = coef(big.rich.avg)

# Changing the coefficents to the averaged coefficient # see https://stackoverflow.com/questions/28633537/partial-residual-plot-based-on-model-average-coefficients-in-r
big.rich.mod$coefficients = coefMod[names(coef(big.rich.full))]
big.rich.mod$residuals = big.rich.rr.data$rr.rich - predict(big.rich.mod)
big.rich.pres = resid(big.rich.mod, type='partial')

# Read in lil fish richness data
lil.rich.data = read.delim('data/Rs/lil_rich_data.txt', header = TRUE, sep = ' ')

# Get the Richness response (Rs) for lil Fish 
lil.rich.wide = lil.rich.data[,-4] %>% group_by(site.id) %>% pivot_wider(names_from = res, values_from = lil.rich) 

lil.rich.mean = lil.rich.wide %>% # This takes the mean of the richness observations at each res and non res (Average of richness at each rep)
  group_by(site.id) %>% 
  summarise(N.mean = mean(N), N.sd = sd(N), Y.mean = mean(Y), Y.sd = sd(Y)) # creates dataframe by site.id that includes mean and sd of Y and N columns

lil.rich.mean.rr = lil.rich.mean %>%
  group_by(site.id) %>%
  mutate(rr.rich = Y.mean - N.mean) # creates new dataframe for RR difference between the mean richness per site


#////////////
# lil.rich.rr
#\\\\\\\\\\\\


lil.rich.rr.data = cbind(res.data, as.data.frame(lil.rich.mean.rr))

lil.rich.diff.lm = lm(rr.rich ~ yrs.prot + houses + pen + mouth.dist + discharge + clo.res.dist + bet.cent.n + res.area + rd.dist + vil.dist, data = lil.rich.rr.data, na.action = na.pass)
summary(lil.rich.diff.lm)
lil.rich.diff.lm.std = standardize(lil.rich.diff.lm, standardize.y = FALSE) # numeric values are rescaled to have mean = 0 and sd = 0.5 lil predictor variables to account for non-normality
lil.rich.diff.model.set = dredge(lil.rich.diff.lm.std) # dredge function from MuMIn package calculates lil possible factorial models and evaluates fit
sum(lil.rich.diff.model.set$weight[1:327]) # This shows the summed Aikake weights are < 0.95
lil.rich.95.models = get.models(lil.rich.diff.model.set, 1:327) # just found the # models < .95 summed Aikake weight using sum(lil.rich.diff.model.set$weight[1:393])

lil.rich.95.models.avg = model.avg(lil.rich.95.models, revised.var = TRUE) # conducts model averaging based on importance (summed Akaike weight for lil models)
summary(lil.rich.95.models.avg)

# Create a barplot of top models by importance
lil.rich.95.models.avg$sw
lil.rich.diff.impt = as.matrix(lil.rich.95.models.avg$sw) # extracts the Akaike weights for each term in the top models
lil.rich.diff.coef = as.matrix(lil.rich.95.models.avg$coefficients[2,-1]) # extracts the coefficients from the averaged model for coloring by positive or negative

lil.rich.diff.bar = merge(lil.rich.diff.impt, lil.rich.diff.coef, by = 'row.names', lil.x = TRUE) # combines the Akaike weights data and coefficients 
row.names(lil.rich.diff.bar) = lil.rich.diff.bar$Row.names
lil.rich.diff.bar = lil.rich.diff.bar[,-1]

mod.pred = matrix(rep(NA, times = 10), dimnames = list(c('z.yrs.prot', 'z.houses', 'c.pen', 'z.mouth.dist', 'z.discharge', 'z.clo.res.dist', 'z.bet.cent.n', 'z.res.area', 'z.rd.dist', 'z.vil.dist')))
mod.pred.sort = as.matrix(mod.pred[order(row.names(mod.pred)),])

lil.rich.diff.lil = merge(mod.pred.sort, lil.rich.diff.bar, by = 'row.names', lil.x = TRUE) # combines the Akaike weights data and coefficients 
lil.rich.diff.lil = lil.rich.diff.lil[,-2]
lil.rich.diff.lil[is.na(lil.rich.diff.lil)] <- 0
colnames(lil.rich.diff.lil) = c('p.variable', 'importance', 'coefficient')

target = c('z.res.area', 'z.discharge', 'z.yrs.prot', 'c.pen', 'z.vil.dist', 'z.rd.dist', 'z.houses', 'z.bet.cent.n', 'z.mouth.dist', 'z.clo.res.dist')
lil.rich.diff.plot = lil.rich.diff.lil[match(target, lil.rich.diff.lil$p.variable),]

lil.plot = barplot(lil.rich.diff.plot$importance, names.arg = lil.rich.diff.plot$p.variable, 
                   col = ifelse(lil.rich.diff.plot$coefficient > 0, 'black', 'white'), 
                   las = 2, main = 'Small Fishes',
                   ylim = c(0,1)) 


# Get the partial residuals for partial residual plots (At end)
lil.rich.full = lil.rich.diff.lm.std
lil.rich.avg = lil.rich.95.models.avg
lil.rich.mod = lil.rich.full

# Changing the coefficents to the averaged coefficient # see https://stackoverflow.com/questions/28633537/partial-residual-plot-based-on-model-average-coefficients-in-r
lil.rich.mod$coefficients <- coefMod[names(coef(lil.rich.full))]
lil.rich.mod$residuals <- lil.rich.rr.data$rr.rich - predict(lil.rich.mod)
lil.rich.pres <- resid(lil.rich.mod, type = 'partial')


# Get the Richness response (Rs) for Pred Fish 
# Read in Pred fish richness data
# Read in pred fish richness data
pred.rich.data = read.delim('data/Rs/pred_rich_data.txt', header = TRUE, sep = ' ')

# Get the Richness response (Rs) for pred Fish 
pred.rich.wide = pred.rich.data[,-4] %>% group_by(site.id) %>% pivot_wider(names_from = res, values_from = pred.rich) 

pred.rich.mean = pred.rich.wide %>% # This takes the mean of the richness observations at each res and non res (Average of richness at each rep)
  group_by(site.id) %>% 
  summarise(N.mean = mean(N), N.sd = sd(N), Y.mean = mean(Y), Y.sd = sd(Y)) # creates dataframe by site.id that includes mean and sd of Y and N columns

pred.rich.mean.rr = pred.rich.mean %>%
  group_by(site.id) %>%
  mutate(rr.rich = Y.mean - N.mean) # creates new dataframe for RR difference between the mean richness per site


#////////////
# pred.rich.rr
#\\\\\\\\\\\\

pred.rich.rr.data = cbind(res.data, as.data.frame(pred.rich.mean.rr))

pred.rich.diff.lm = lm(rr.rich ~ yrs.prot + houses + pen + mouth.dist + discharge + clo.res.dist + bet.cent.n + res.area + rd.dist + vil.dist, data = pred.rich.rr.data, na.action = na.pass)
summary(pred.rich.diff.lm)
pred.rich.diff.lm.std = standardize(pred.rich.diff.lm, standardize.y = FALSE) # numeric values are rescaled to have mean = 0 and sd = 0.5 pred predictor variables to account for non-normality
pred.rich.diff.model.set = dredge(pred.rich.diff.lm.std) # dredge function from MuMIn package calculates pred possible factorial models and evaluates fit
sum(pred.rich.diff.model.set$weight[1:212]) # This shows the summed Aikake weights are < 0.95
pred.rich.95.models = get.models(pred.rich.diff.model.set, 1:212) # just found the # models < .95 summed Aikake weight using sum(pred.rich.diff.model.set$weight[1:393])

pred.rich.95.models.avg = model.avg(pred.rich.95.models, revised.var = TRUE) # conducts model averaging based on importance (summed Akaike weight for pred models)
summary(pred.rich.95.models.avg)
# Create a barplot of top models by importance
pred.rich.95.models.avg$sw
pred.rich.diff.impt = as.matrix(pred.rich.95.models.avg$sw) # extracts the Akaike weights for each term in the top models
pred.rich.diff.coef = as.matrix(pred.rich.95.models.avg$coefficients[2,-1]) # extracts the coefficients from the averaged model for coloring by positive or negative

pred.rich.diff.bar = merge(pred.rich.diff.impt, pred.rich.diff.coef, by = 'row.names', pred.x = TRUE) # combines the Akaike weights data and coefficients 
row.names(pred.rich.diff.bar) = pred.rich.diff.bar$Row.names
pred.rich.diff.bar = pred.rich.diff.bar[,-1]

mod.pred = matrix(rep(NA, times = 10), dimnames = list(c('z.yrs.prot', 'z.houses', 'c.pen', 'z.mouth.dist', 'z.discharge', 'z.clo.res.dist', 'z.bet.cent.n', 'z.res.area', 'z.rd.dist', 'z.vil.dist')))
mod.pred.sort = as.matrix(mod.pred[order(row.names(mod.pred)),])

pred.rich.diff.pred = merge(mod.pred.sort, pred.rich.diff.bar, by = 'row.names', pred.x = TRUE) # combines the Akaike weights data and coefficients 
pred.rich.diff.pred = pred.rich.diff.pred[,-2]
pred.rich.diff.pred[is.na(pred.rich.diff.pred)] <- 0
colnames(pred.rich.diff.pred) = c('p.variable', 'importance', 'coefficient')

target = c('z.res.area', 'z.discharge', 'z.yrs.prot', 'c.pen', 'z.vil.dist', 'z.rd.dist', 'z.houses', 'z.bet.cent.n', 'z.mouth.dist', 'z.clo.res.dist')
pred.rich.diff.plot = pred.rich.diff.pred[match(target, pred.rich.diff.pred$p.variable),]

pred.plot = barplot(pred.rich.diff.plot$importance, names.arg = pred.rich.diff.plot$p.variable, 
                   col = ifelse(pred.rich.diff.plot$coefficient > 0, 'black', 'white'), 
                   las = 2, main = 'Predators',
                   ylim = c(0,1)) 


# Get the partial residuals for partial residual plots (At end)
pred.rich.full = pred.rich.diff.lm.std
pred.rich.avg = pred.rich.95.models.avg
pred.rich.mod = pred.rich.full

coefMod = coef(pred.rich.avg)

# Changing the coefficents to the averaged coefficient # see https://stackoverflow.com/questions/28633537/partial-residual-plot-based-on-model-average-coefficients-in-r
pred.rich.mod$coefficients <- coefMod[names(coef(pred.rich.full))]
pred.rich.mod$residuals <- pred.rich.rr.data$rr.rich - predict(pred.rich.mod)
pred.rich.pres <- resid(pred.rich.mod, type='partial')

# Get the Richness response (Rs) for omni Fish 
# Read in omni fish richness data
omni.rich.data = read.delim('data/Rs/omni_rich_data.txt', header = TRUE, sep = ' ')

# Get the Richness response (Rs) for omni Fish 
omni.rich.wide = omni.rich.data[-4] %>% group_by(site.id) %>% pivot_wider(names_from = res, values_from = omni.rich)

omni.rich.mean = omni.rich.wide %>% # This takes the mean of the richness observations at each res and non res (Average of richness at each rep)
  group_by(site.id) %>% 
  summarise(N.mean = mean(N), N.sd = sd(N), Y.mean = mean(Y), Y.sd = sd(Y)) # creates dataframe by site.id that includes mean and sd of Y and N columns

omni.rich.mean.rr = omni.rich.mean %>%
  group_by(site.id) %>%
  mutate(rr.rich = Y.mean - N.mean) # creates new dataframe for RR difference between the mean richness per site

#////////////
# omni.rich.rr
#\\\\\\\\\\\\

omni.rich.rr.data = cbind(res.data, as.data.frame(omni.rich.mean.rr))

omni.rich.diff.lm = lm(rr.rich ~ yrs.prot + houses + pen + mouth.dist + discharge + clo.res.dist + bet.cent.n + res.area + rd.dist + vil.dist, data = omni.rich.rr.data, na.action = na.pass)
summary(omni.rich.diff.lm)
omni.rich.diff.lm.std = standardize(omni.rich.diff.lm, standardize.y = FALSE) # numeric values are rescaled to have mean = 0 and sd = 0.5 omni predictor variables to account for non-normality
omni.rich.diff.model.set = dredge(omni.rich.diff.lm.std) # dredge function from MuMIn package calculates omni possible factorial models and evaluates fit
sum(omni.rich.diff.model.set$weight[1:353]) # This shows the summed Aikake weights are < 0.95
omni.rich.95.models = get.models(omni.rich.diff.model.set, 1:353) # just found the # models < .95 summed Aikake weight using sum(omni.rich.diff.model.set$weight[1:393])

omni.rich.95.models.avg = model.avg(omni.rich.95.models, revised.var = TRUE) # conducts model averaging based on importance (summed Akaike weight for omni models)
summary(omni.rich.95.models.avg)
# Create a barplot of top models by importance
omni.rich.95.models.avg$sw
omni.rich.diff.impt = as.matrix(omni.rich.95.models.avg$sw) # extracts the Akaike weights for each term in the top models
omni.rich.diff.coef = as.matrix(omni.rich.95.models.avg$coefficients[2,-1]) # extracts the coefficients from the averaged model for coloring by positive or negative

omni.rich.diff.bar = merge(omni.rich.diff.impt, omni.rich.diff.coef, by = 'row.names', omni.x = TRUE) # combines the Akaike weights data and coefficients 
row.names(omni.rich.diff.bar) = omni.rich.diff.bar$Row.names
omni.rich.diff.bar = omni.rich.diff.bar[,-1]

mod.pred = matrix(rep(NA, times = 10), dimnames = list(c('z.yrs.prot', 'z.houses', 'c.pen', 'z.mouth.dist', 'z.discharge', 'z.clo.res.dist', 'z.bet.cent.n', 'z.res.area', 'z.rd.dist', 'z.vil.dist')))
mod.pred.sort = as.matrix(mod.pred[order(row.names(mod.pred)),])

omni.rich.diff.omni = merge(mod.pred.sort, omni.rich.diff.bar, by = 'row.names', omni.x = TRUE) # combines the Akaike weights data and coefficients 
omni.rich.diff.omni = omni.rich.diff.omni[,-2]
omni.rich.diff.omni[is.na(omni.rich.diff.omni)] <- 0
colnames(omni.rich.diff.omni) = c('p.variable', 'importance', 'coefficient')

target = c('z.res.area', 'z.discharge', 'z.yrs.prot', 'c.pen', 'z.vil.dist', 'z.rd.dist', 'z.houses', 'z.bet.cent.n', 'z.mouth.dist', 'z.clo.res.dist')
omni.rich.diff.plot = omni.rich.diff.omni[match(target, omni.rich.diff.omni$p.variable),]

omni.plot = barplot(omni.rich.diff.plot$importance, names.arg = omni.rich.diff.plot$p.variable, 
                   col = ifelse(omni.rich.diff.plot$coefficient > 0, 'black', 'white'), 
                   las = 2, main = 'Omnivores',
                   ylim = c(0,1)) 

# Get the partial residuals for partial residual plots (At end)
omni.rich.full = omni.rich.diff.lm.std
omni.rich.avg = omni.rich.95.models.avg
omni.rich.mod = omni.rich.full

coefMod = coef(omni.rich.avg)

# Changing the coefficents to the averaged coefficient # see https://stackoverflow.com/questions/28633537/partial-residual-plot-based-on-model-average-coefficients-in-r
omni.rich.mod$coefficients <- coefMod[names(coef(omni.rich.full))]
omni.rich.mod$residuals <- omni.rich.rr.data$rr.rich - predict(omni.rich.mod)
omni.rich.pres <- resid(omni.rich.mod, type='partial')


# Get the Richness response (Rs) for herb Fish 
# Read in herb fish richness data
herb.rich.data = read.delim('data/Rs/herb_rich_data.txt', header = TRUE, sep = ' ')

# Get the Richness response (Rs) for herb Fish 
herb.rich.wide = herb.rich.data[,-4] %>% group_by(site.id) %>% pivot_wider(names_from = res, values_from = herb.rich)

herb.rich.mean = herb.rich.wide %>% # This takes the mean of the richness observations at each res and non res (Average of richness at each rep)
  group_by(site.id) %>% 
  summarise(N.mean = mean(N), N.sd = sd(N), Y.mean = mean(Y), Y.sd = sd(Y)) # creates dataframe by site.id that includes mean and sd of Y and N columns

herb.rich.mean.rr = herb.rich.mean %>%
  group_by(site.id) %>%
  mutate(rr.rich = Y.mean - N.mean) # creates new dataframe for RR difference between the mean richness per site

#////////////
# herb.rich.rr
#\\\\\\\\\\\\


herb.rich.rr.data = cbind(res.data, as.data.frame(herb.rich.mean.rr))

herb.rich.diff.lm = lm(rr.rich ~ yrs.prot + houses + pen + mouth.dist + discharge + clo.res.dist + bet.cent.n + res.area + rd.dist + vil.dist, data = herb.rich.rr.data, na.action = na.pass)
summary(herb.rich.diff.lm)
herb.rich.diff.lm.std = standardize(herb.rich.diff.lm, standardize.y = FALSE) # numeric values are rescaled to have mean = 0 and sd = 0.5 herb predictor variables to account for non-normality
herb.rich.diff.model.set = dredge(herb.rich.diff.lm.std) # dredge function from MuMIn package calculates herb possible factorial models and evaluates fit
sum(herb.rich.diff.model.set$weight[1:236]) # This shows the summed Aikake weights are < 0.95
herb.rich.95.models = get.models(herb.rich.diff.model.set, 1:236) # just found the # models < .95 summed Aikake weight using sum(herb.rich.diff.model.set$weight[1:393])

herb.rich.95.models.avg = model.avg(herb.rich.95.models, revised.var = TRUE) # conducts model averaging based on importance (summed Akaike weight for herb models)
summary(herb.rich.95.models.avg)

# Create a barplot of top models by importance
herb.rich.95.models.avg$sw
herb.rich.diff.impt = as.matrix(herb.rich.95.models.avg$sw) # extracts the Akaike weights for each term in the top models
herb.rich.diff.coef = as.matrix(herb.rich.95.models.avg$coefficients[2,-1]) # extracts the coefficients from the averaged model for coloring by positive or negative

herb.rich.diff.bar = merge(herb.rich.diff.impt, herb.rich.diff.coef, by = 'row.names', herb.x = TRUE) # combines the Akaike weights data and coefficients 
row.names(herb.rich.diff.bar) = herb.rich.diff.bar$Row.names
herb.rich.diff.bar = herb.rich.diff.bar[,-1]

mod.herb = matrix(rep(NA, times = 10), dimnames = list(c('z.yrs.prot', 'z.houses', 'c.pen', 'z.mouth.dist', 'z.discharge', 'z.clo.res.dist', 'z.bet.cent.n', 'z.res.area', 'z.rd.dist', 'z.vil.dist')))
mod.herb.sort = as.matrix(mod.herb[order(row.names(mod.herb)),])

herb.rich.diff.herb = merge(mod.herb.sort, herb.rich.diff.bar, by = 'row.names', herb.x = TRUE) # combines the Akaike weights data and coefficients 
herb.rich.diff.herb = herb.rich.diff.herb[,-2]
herb.rich.diff.herb[is.na(herb.rich.diff.herb)] <- 0
colnames(herb.rich.diff.herb) = c('p.variable', 'importance', 'coefficient')

target = c('z.res.area', 'z.discharge', 'z.yrs.prot', 'c.pen', 'z.vil.dist', 'z.rd.dist', 'z.houses', 'z.bet.cent.n', 'z.mouth.dist', 'z.clo.res.dist')
herb.rich.diff.plot = herb.rich.diff.herb[match(target, herb.rich.diff.herb$p.variable),]

herb.plot = barplot(herb.rich.diff.plot$importance, names.arg = herb.rich.diff.plot$p.variable, 
                   col = ifelse(herb.rich.diff.plot$coefficient > 0, 'black', 'white'), 
                   las = 2, main = 'Herbivores',
                   ylim = c(0,1)) 

# Get the partial residuals for partial residual plots (At end)
herb.rich.full = herb.rich.diff.lm.std
herb.rich.avg = herb.rich.95.models.avg
herb.rich.mod = herb.rich.full

coefMod = coef(herb.rich.avg)

# Changing the coefficents to the averaged coefficient # see https://stackoverflow.com/questions/28633537/partial-residual-plot-based-on-model-average-coefficients-in-r
herb.rich.mod$coefficients <- coefMod[names(coef(herb.rich.full))]
# Changing the residuals
herb.rich.mod$residuals <- herb.rich.rr.data$rr.rich - predict(herb.rich.mod)
herb.rich.pres <- resid(herb.rich.mod, type='partial')

# Get the averaged model results to check for significance of individual variables for partial residual plots

summary(all.rich.avg) 
summary(big.rich.avg)
summary(lil.rich.avg)
summary(pred.rich.avg)
summary(omni.rich.avg)
summary(herb.rich.avg)

#\\\\\\\\\\\\\\\\\\\\\\
#Partial residual plots (Fig. 4)
#//////////////////////

# set matching color scheme
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



# RESERVE AREA

dev.new = TRUE
pdf('figures/fig_4/Figure_4_rich_resid_all.pdf', width = 15, height = 4)
#pdf('figures/fig_4/single_plots/rich_area_resid_DB.pdf', width = 7, height = 6)


#par(oma = c(2,6,2,2))
par(mfrow = c(1,5), oma = c(0,4,0,1))

plot(omni.rich.pres[,"z.res.area"] ~ omni.rich.rr.data$res.area, pch = all.pch, col = omni.col, bg = omni.col.a, cex = 2, bty = 'n',
     ylab = '', xlab = '', axes = FALSE, ylim = c(-4,6), xlim = c(0, 25000), type = 'n')

abline(h = 0, col = 'grey')
par(xpd = TRUE)
points(omni.rich.pres[,"z.res.area"] ~ omni.rich.rr.data$res.area, pch = all.pch, col = omni.col, bg = omni.col.a, cex = 2)
points(big.rich.pres[,"z.res.area"] ~ big.rich.rr.data$res.area, pch = all.pch, bg = big.col.a, cex = 2, col = big.col, xpd = TRUE)
points(all.rich.pres[,"z.res.area"] ~ all.rich.rr.data$res.area, pch = all.pch, col = all.col, cex = 2, bg = all.col.a)
par(xpd = FALSE)

abline(lm(all.rich.pres[,"z.res.area"]~all.rich.rr.data$res.area), lty = 1, col = all.col, lwd = 2) # p < 0.05
abline(lm(big.rich.pres[,"z.res.area"]~big.rich.rr.data$res.area), lty = 2, col = big.col, lwd = 2) # p < 0.1
abline(lm(omni.rich.pres[,"z.res.area"]~omni.rich.rr.data$res.area), lty = 1, col = omni.col, lwd = 2) # p < 0.05

axis(side = 1, at = seq(0,25000, 12500), labels = c('0','1.25', '2.5'), cex.axis = 1.5, col = axis.col, col.axis = axis.col)
axis(side = 2, at = c(-4,0,6), las = 1, cex.axis = 1.5, col = axis.col, col.axis = axis.col)
mtext(side = 1, 'Reserve Area (Ha)', line = 3, cex = 1.2, col  = axis.col)
mtext(side = 2, 'Partial Residuals', line = 3, cex = 1.2, col = axis.col)
mtext(expression(paste('R'[italic(r)])), side = 2, line = 5.5, las = 1, cex = 1.7, col = axis.col)
#dev.off()

# DISCHARGE

# dev.new = TRUE
# 
# pdf('figures/fig_4/single_plots/rich_discharge_resid_DB.pdf', width = 7, height = 6)
# par(mfrow = c(1,1))
# par(oma = c(0,0,0,0))
# 
# 
plot(pred.rich.pres[,"z.discharge"] ~ pred.rich.rr.data$discharge, pch = all.pch, col = pred.col, cex = 2, bg = pred.col.a, bty = 'n',
     ylab = '', xlab = '', axes = FALSE, ylim = c(-1,2), xlim = c(0,30), type = 'n')

abline(h = 0, col = 'grey')
points(pred.rich.pres[,"z.discharge"] ~ pred.rich.rr.data$discharge, pch = all.pch, col = pred.col, cex = 2, bg = pred.col.a)
points(herb.rich.pres[,"z.discharge"] ~ herb.rich.rr.data$discharge, pch = all.pch, col = herb.col, bg = herb.col.a, cex = 2)


abline(lm(herb.rich.pres[,"z.discharge"]~herb.rich.rr.data$discharge), lty = 1, col = herb.col, lwd = 2) # p < 0.05
abline(lm(pred.rich.pres[,"z.discharge"]~pred.rich.rr.data$discharge), lty = 1, col = pred.col, lwd = 2) # p < 0.05

axis(side = 1, at = seq(0,30, 15), cex.axis = 1.5, col = axis.col, col.axis = axis.col)
axis(side = 2, at = c(-1, 0, 1, 2), las = 1, cex.axis = 1.5, col = axis.col, col.axis = axis.col)
mtext(side = 1, expression(paste('Discharge (m'^'3'*' sec'^'-1'*')')), line = 3, cex = 1.2, col = axis.col)

# dev.off()


# Betweenness Centrality

# dev.new = TRUE

# pdf('figures/fig_4/single_plots/rich_betcent_resid_DB.pdf', width = 7, height = 6)
# par(mfrow = c(1,1))
# par(oma = c(0,0,0,0))
# 
plot(all.rich.pres[,"z.bet.cent.n"] ~ all.rich.rr.data$bet.cent.n, pch = all.pch, col = all.col, bg = all.col.a, cex = 2, bty = 'n',
     ylab = '', xlab = '', axes = FALSE, ylim = c(-3,3), xlim = c(0, 0.5), type = 'n')

abline(h = 0, col = 'grey')
points(all.rich.pres[,"z.bet.cent.n"] ~ all.rich.rr.data$bet.cent.n, pch = all.pch, col = all.col, bg = all.col.a, cex = 2, xpd = NA)
points(big.rich.pres[,"z.bet.cent.n"] ~ big.rich.rr.data$bet.cent.n, pch = all.pch, col = big.col, cex = 2, bg = big.col.a, xpd = NA)

abline(lm(all.rich.pres[,"z.bet.cent.n"]~all.rich.rr.data$bet.cent.n), lty = 1, col = all.col, lwd = 2) # p < 0.05
abline(lm(big.rich.pres[,"z.bet.cent.n"]~big.rich.rr.data$bet.cent.n), lty = 1, col = big.col, lwd = 2) # p < 0.05

axis(side = 1, at = seq(0,0.5, 0.25), cex.axis = 1.5, col = axis.col, col.axis = axis.col)
axis(side = 2, at = c(-3,-1.5,0,1.5,3), labels = c('-3', '-1.5', '0', '1.5', '3'), las = 1, cex.axis = 1.5, col = axis.col, col.axis = axis.col)
mtext(side = 1, 'Betweenness Centrality', line = 3, cex = 1.2, col  = axis.col)

#dev.off()

# Distance to Mouth


plot(herb.rich.pres[,"z.mouth.dist"] ~ herb.rich.rr.data$mouth.dist, pch = all.pch, col = herb.col, cex = 2, bg = herb.col.a, bty = 'n',
     ylab = '', xlab = '', axes = FALSE, ylim = c(-0.5,1), xlim = c(10000, 70000), type = 'n')

abline(h = 0, col = 'grey')
points(herb.rich.pres[,"z.mouth.dist"] ~ herb.rich.rr.data$mouth.dist, pch = all.pch, col = herb.col, cex = 2, bg = herb.col.a)

abline(lm(herb.rich.pres[,"z.mouth.dist"]~herb.rich.rr.data$mouth.dist), lty = 1, col = herb.col, lwd = 1.5) # p < 0.05

axis(side = 1, at = seq(10000, 70000, 30000), labels = c('10','40','70'), cex.axis = 1.5)
axis(side = 2, at = c(-0.5,0,0.5,1), labels = c('-0.5','0', '0.5', '1'), las = 1, cex.axis = 1.5)
mtext(side = 1, 'Distance to Mouth (km)', line = 3, cex = 1.2, col  = 'black')

# Age

# dev.new = TRUE

#pdf('figures/fig_4/single_plots/rich_age_resid_DB.pdf', width = 7, height = 6)
#par(mfrow = c(1,1))
#par(oma = c(0,0,0,0))
plot(all.rich.pres[,"z.yrs.prot"] ~ all.rich.rr.data$yrs.prot, pch = all.pch, col = all.col, bg = all.col.a, cex = 2, bty = 'n',
     ylab = '', xlab = '', axes = FALSE, ylim = c(-10,5), xlim = c(0, 27), type = 'n')

abline(h = 0, col = 'grey')

points(all.rich.pres[,"z.yrs.prot"] ~ all.rich.rr.data$yrs.prot, pch = all.pch, col = all.col, bg = all.col.a, cex = 2)
points(omni.rich.pres[,"z.yrs.prot"] ~ omni.rich.rr.data$yrs.prot, pch = all.pch, col = omni.col, bg = omni.col.a, cex = 2)
points(pred.rich.pres[,"z.yrs.prot"] ~ pred.rich.rr.data$yrs.prot, pch = all.pch, col = pred.col, bg = pred.col.a, cex = 2)
abline(lm(all.rich.pres[,"z.yrs.prot"]~all.rich.rr.data$yrs.prot), lty = 2, col = all.col, lwd = 2) # p < 0.1
abline(lm(pred.rich.pres[,"z.yrs.prot"]~pred.rich.rr.data$yrs.prot), lty = 1, col = pred.col, lwd = 2) # p < 0.05
abline(lm(omni.rich.pres[,"z.yrs.prot"]~omni.rich.rr.data$yrs.prot), lty = 1, col = omni.col, lwd = 2) # p < 0.05
axis(side = 1, at = seq(0,27, 9), cex.axis = 1.5, col = axis.col, col.axis = axis.col)
axis(side = 2, at = c(-10,-5,0,5), las = 1, cex.axis = 1.5, col = axis.col, col.axis = axis.col)

mtext(side = 1, 'Years Protected', line = 3, cex = 1.2, col  = axis.col)

dev.off()


## Get data for heatmap; Fig 3

all.rich.diff.plot$fish.cat = rep('All', length(all.rich.diff.plot$coefficient))
big.rich.diff.plot$fish.cat = rep('Large', length(big.rich.diff.plot$coefficient))
lil.rich.diff.plot$fish.cat = rep('Small', length(lil.rich.diff.plot$coefficient))
pred.rich.diff.plot$fish.cat = rep('Predators', length(pred.rich.diff.plot$coefficient))
omni.rich.diff.plot$fish.cat = rep('Omnivores', length(omni.rich.diff.plot$coefficient))
herb.rich.diff.plot$fish.cat = rep('Herbivores', length(herb.rich.diff.plot$coefficient))

rich.heat.data = as.data.frame(rbind(all.rich.diff.plot, big.rich.diff.plot, lil.rich.diff.plot, pred.rich.diff.plot, omni.rich.diff.plot, herb.rich.diff.plot))

attach(rich.heat.data)
rich.heat.data$imp.sign <- with(rich.heat.data, ifelse(coefficient <= 0, -importance, importance)) # this adds a new column which makes the sign of the importance match the coefficient for color ramp
rich.heat.wide = as.data.frame(rich.heat.data[,-c(2:3)] %>% group_by(p.variable) %>% spread(fish.cat, imp.sign)) # This removes the importance and coefficient columns, combined into imp.sign, then makes a wide dataframe
rownames(rich.heat.wide) = rich.heat.wide[,1]
rich.heat.wide = rich.heat.wide[,-1]

# need to reorder the variables here
target = rev(c('z.res.area', 'z.discharge', 'z.yrs.prot', 'c.pen', 'z.vil.dist', 'z.rd.dist', 'z.houses', 'z.bet.cent.n', 'z.mouth.dist', 'z.clo.res.dist'))
rich.heat.wide.z = rich.heat.wide[match(target, rownames(rich.heat.wide)),]
rich.heat.wide.z = rich.heat.wide.z[c('All', 'Large', 'Small', 'Predators', 'Omnivores', 'Herbivores')]

write.table(rich.heat.wide.z, 'figures/fig_3/data/rich.heatmap.data.txt')


