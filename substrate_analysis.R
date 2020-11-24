# This code analyzes habitat data for all reserve and non-reserve sites
# in order to rule out habitat differences as a potential cause of 
# the patterns observed inside verses outside reserves.

# We also use these data to calculate the mean survey stream width at each location,
# which we used to calculate the survey area for count offsets, included in the input data file 'fish_counts.txt'


rm(list = ls())
library(tidyr)
library(dplyr)
library(vegan)
library(lme4)
library(lmerTest)

substrate = read.delim('habitat_data.txt', header = TRUE)

substrate$dummy = 1

hab.by.plot = substrate %>% group_by(site.id, res, plot) %>% summarise(md.part.size = median(sub.size), m.part.size = mean(sub.size), m.depth = mean(depth.cm), max.depth = mean(max.depth.cm), m.width = mean(width.cm))
hab.by.site = substrate %>% group_by(site.id) %>% summarise(md.part.size = median(sub.size), m.part.size = mean(sub.size), m.depth = mean(depth.cm), max.depth = mean(max.depth.cm), m.width = mean(width.cm))
hab.by.res = substrate %>% group_by(site.id, res) %>% summarise(md.part.size = median(sub.size), m.part.size = mean(sub.size), m.depth = mean(depth.cm, na.rm = TRUE), max.depth = mean(max.depth.cm), m.width = mean(width.cm))


mean.width = as.data.frame(substrate %>% group_by(site.id, res, plot) %>% summarize(mean.width.m = mean(width.cm/100)))
#write.table(mean.width, file = 'data/mean_stream_width.txt', row.names = FALSE)

##########
#::: Statistical Tests for differences by reserve sites ::::::
##########

#::: Test if median substrate particle varies by reserve

# visualize differences, no strongly apparent differences revealed using boxplot or histogram
par(mfrow = c(2,2))

boxplot(log10(hab.by.plot$md.part.size)[hab.by.plot$res == 'Y'], xlab = 'Reserve', ylab = 'log10(Median Particle Size)')
boxplot(log10(hab.by.plot$md.part.size)[hab.by.plot$res == 'N'], xlab = 'Non-Reserve', ylab = 'log10(Median Particle Size)')
hist(log10(hab.by.plot$md.part.size)[hab.by.plot$res == 'Y'], xlab = 'log10(Median Particle Size) \n Reserve', main = '')
hist(log10(hab.by.plot$md.part.size)[hab.by.plot$res == 'N'], xlab = 'log10(Median Particle Size) \n Non-Reserve', main = '')

sub.lm = lmer(log10(md.part.size) ~ res + (1|site.id), data = hab.by.plot) # use site.id as a blocking factor/random intercept
plot(sub.lm)
summary(sub.lm)


#::: Test if mean depth varies by reserve
# visualize differences
boxplot(log10(hab.by.plot$m.depth)[hab.by.plot$res == 'Y'], xlab = 'Reserve', ylab = 'log10(Mean Depth cm)', ylim = c(0.8, 2.25))
boxplot(log10(hab.by.plot$m.depth)[hab.by.plot$res == 'N'], xlab = 'Non-Reserve', ylab = 'log10(Mean Depth cm)', ylim = c(0.8, 2.25))
hist(log10(hab.by.plot$m.depth)[hab.by.plot$res == 'Y'], xlab = 'log10(Mean Depth) \n Reserve', main = '', xlim = c(0.5,2.5), ylim = c(0,15))
hist(log10(hab.by.plot$m.depth)[hab.by.plot$res == 'N'], xlab = 'log10(Mean Depth) \n Non-Reserve', main = '', xlim = c(0.5,2.5), ylim = c(0,15))

depth.lm = lmer(log10(m.depth) ~ res + (1|site.id), data = hab.by.plot)
plot(depth.lm)
summary(depth.lm)

#::: Test if max depth varies by reserve
boxplot(log10(hab.by.plot$max.depth)[hab.by.plot$res == 'Y'], xlab = 'Reserve', ylab = 'log10(Max Depth cm)', ylim = c(1, 3))
boxplot(log10(hab.by.plot$max.depth)[hab.by.plot$res == 'N'], xlab = 'Non-Reserve', ylab = 'log10(Max Depth cm)', ylim = c(1, 3))
hist(log10(hab.by.plot$max.depth)[hab.by.plot$res == 'Y'], xlab = 'log10(Max Depth) \n Reserve', main = '', xlim = c(1,3), ylim = c(0,15))
hist(log10(hab.by.plot$max.depth)[hab.by.plot$res == 'N'], xlab = 'log10(Max Depth) \n Non-Reserve', main = '', xlim = c(1,3), ylim = c(0,15))

max.depth.lm = lmer(log10(max.depth) ~ res + (1|site.id), data = hab.by.plot)
plot(max.depth.lm)
summary(max.depth.lm)

#::: Test if mean width varies by reserve
boxplot(log10(hab.by.plot$m.width)[hab.by.plot$res == 'Y'], xlab = 'Reserve', ylab = 'log10(Mean Width cm)', ylim = c(2, 4))
boxplot(log10(hab.by.plot$m.width)[hab.by.plot$res == 'N'], xlab = 'Non-Reserve', ylab = 'log10(Mean Width cm)', ylim = c(2, 4))
hist(log10(hab.by.plot$m.width)[hab.by.plot$res == 'Y'], xlab = 'log10(Mean Width) \n Reserve', main = '', xlim = c(2,4), ylim = c(0,15))
hist(log10(hab.by.plot$m.width)[hab.by.plot$res == 'N'], xlab = 'log10(Mean Width) \n Non-Reserve', main = '', xlim = c(2,4), ylim = c(0,15))

width.lm = lmer(log10(m.width) ~ res + (1|site.id), data = hab.by.plot)
plot(width.lm)
summary(width.lm)


# :::: Get the Simpson's diversity index for substrate hetergeneity
sub.w = substrate %>% group_by(site.id, res, plot, sub.cat) %>% summarise(count = sum(dummy)) %>% spread(sub.cat, count, fill = 0)
sub.w$site.no = 1:length(sub.w$site.id) # add dummy 'site.id'
sub.w.no = sub.w[,-c(1,2,3)] # take out site, reserve, plot variables
sub.D = diversity(sub.w.no, 'simpson') # calculates the simpson's index 'D' for substrate classes

sub.pca = rda(sub.w.no[,-9]) # Takes out the site.no category and the dummy variable prior to running a PCA
summary(sub.pca)
par(mfrow = c(1,1))
biplot(sub.pca)

hab.scores = scores(sub.pca)
hab.pc1 = as.vector(hab.scores$sites[,1])
hab.pc2 = as.vector(hab.scores$sites[,2])

hab.data = cbind(as.data.frame(hab.by.plot), sub.D, hab.pc1, hab.pc2) # stick the substrate diversity onto the habitat dataframe

#::: Test if substrate diversity differs by reserve
sub.D.lm = lmer(sub.D ~ res + (1|site.id), data = hab.data)
summary(sub.D.lm)

# ::: Test if PC1 differs by reserve
pc1.lm = lmer(hab.pc1 ~ res + (1|site.id), data = hab.data)
summary(pc1.lm)

# ::: Test if PC2 differs by reserve
pc2.lm = lmer(hab.pc2 ~ res + (1|site.id), data = hab.data)
summary(pc2.lm)

