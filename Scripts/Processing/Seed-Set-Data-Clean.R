#### Script to clean up seed set data, includes flowering data from demog data sets ####

rm(list=ls())

#### Load libaries ####

library(plyr)
library(dplyr)
library(ggplot2)
library(Rmisc)
library(lme4)
library(tidyr)
library(reshape)
library(gridExtra)
# plots <- "Plots/"

##### Load datasets ####

treat <- read.csv("Data/Setup/Marina-Treatment-30.csv") # treatment data

# 2016 Data
flo.16 <- read.csv("Data/Census-Data/flo-data-16.csv")  
seed.16 <- read.csv("Data/Seed-set/2016/Seed-set-2016.csv")

# 2017 Data
flo.17 <- read.csv("Data/Census-Data/flo-data-17.csv")
seed.17 <- read.csv("Data/Seed-set/2017/Seed-set-2017.csv") 

# Viability Data
viab <- read.csv("Data/Viability/Viability-Overall.csv") # seed viability data
#viab <- viab[c(1:5,8,12), c(1,8,9)]
viab <- viab[c(1:9,12), c(1,8,9)]

#### 2016: Flowering ####
flo.16$Date <- as.Date(flo.16$Date, "%m/%d/%y")

# zeroes not meaningful here, but in seed set they are, get rid of and calculate various summary stats 
flo.16.long <- melt(flo.16, id.vars = c("Year","Date","Plot","Subplot","Species"), measure.vars = c("flo.1","flo.2","flo.3","flo.4","flo.5"), variable_name = "individual")
names(flo.16.long)[7] <- "n.flo"
flo.16.long <- filter(flo.16.long, n.flo > 0)

# flo.16.long <- filter(flo.16.long, !(Species == "Hemizonia congesta" & Subplot == "No Grass" & Treat.Code == "D" & n.flo > 100))
# 
# flo.16.long <- filter(flo.16.long, !(Species == "Hemizonia congesta" & Subplot == "No Grass" & Treat.Code == "C" & n.flo > 150))

# calculate summary stats
flo.16.summary <- ddply(flo.16.long, .(Year, Date, Plot, Subplot, Species), summarize, avg.flo = mean(n.flo), sd.flo = sd(n.flo), n.ind = length(n.flo))

# ###
# # Option 2: find date with most individuals (only up to 5), then find date with max avg.flo
# ###
# test <- ddply(flo.16.summary, .(Year, Plot, Subplot, Species), summarize, Date.1 = Date[which(n.ind == max(n.ind))][1], Date.2 = Date[which(n.ind == max(n.ind))][2], Date.3 = Date[which(n.ind == max(n.ind))][3], Date.4 = Date[which(n.ind == max(n.ind))][4]) # only 4 dates
# 
# test <- melt(test, id.vars = c("Year", "Plot", "Subplot", "Species"))
# 
# names(test)[6] <- "Date"
# 
# test <- merge(flo.16.summary, test[,c(1:4,6)], by = c("Year", "Date",  "Plot", "Subplot", "Species"), all.x = F)
# 
# test <- ddply(test, .(Year, Plot, Subplot, Species), summarize, Date = max(Date[which(avg.flo == max(avg.flo))]))
# 
# flo.16.summary.2 <- merge(flo.16.summary, test)
# 
# ## SOMEONE TELL ME WHY THIS DOESNT WORK ##
# # 
# # test <- ddply(flo.16.summary, .(Year, Plot, Subplot, Species), summarize, Date = 
# #                 ifelse(length(Date[which(n.ind == max(n.ind))]) > 1, 
# #                        Date[which(n.ind == max(n.ind))][which(avg.flo == max(avg.flo))], 
# #                        Date[which(n.ind == max(n.ind))])) 
# # class(test$Date) <- class(flo.16.summary$Date)
# ###

###
# Option 1: Date of highest number of ind. flowering, break tie with latest date of highest number ind setting seed
###

## known issues: this slightly underestimates flower number per individual but only in 3 cases for AGHE (preferred ones in parentheses): 83N-AGHE (5/20), 93N-AGHE (5/20)

flo.16.summary <- merge(flo.16[,c(1:8)], flo.16.summary, by = c("Year", "Date", "Plot", "Subplot", "Species"))
# flo.dates <- ddply(flo.16.summary, .(Year, Plot, Subplot, Species), summarize, Date =
#                      ifelse(length(Date[which(flo.c == max(flo.c))]) > 1,
#                             max(Date[which(flo.b == max(flo.b))]),
#                             Date[which(flo.c == max(flo.c))]))

flo.dates <- ddply(flo.16.summary, .(Year, Plot, Subplot, Species), summarize, Date =
                     ifelse(length(Date[which(flo.b == max(flo.b))]) > 1,
                            max(Date[which(flo.c == max(flo.c))]),
                            Date[which(flo.b == max(flo.b))]))

class(flo.dates$Date) <- class(flo.16.summary$Date)


# extract flowering data from both long and summary dataset
flo.16.long <- merge(flo.dates, flo.16.long, by = c("Year", "Date", "Plot", "Subplot", "Species"), all.y = F)

flo.16.summary <- merge(flo.dates, flo.16.summary, by = c("Year", "Date", "Plot", "Subplot", "Species"), all.y = F)

# 
# diff <- setdiff(flo.16.summary[,c(1:5,9)], flo.16.summary.2[,c(1:6)]) # 98 diff flowering dates
# 
# # compare options to see if this affects flowering
# flo.16.summary$option <- 1
# flo.16.summary.2$option <- 2
# flo.compare <- rbind(flo.16.summary[,c(1,3,4,5,9,12)], flo.16.summary.2[,c(1,3,4,5,6,9)])
# flo.compare <- merge(flo.compare, unique(treat[,c(1,3)]), by = "Plot")
# 
# flo.compare$int.avg.flo <- as.integer(flo.compare$avg.flo)
# hist(flo.compare$avg.flo)
# flo.compare$option <- as.factor(flo.compare$option)
# m.flo <- glmer.nb(int.avg.flo ~ option + (1|Treat.Code/Subplot) + (1 * Treat.Code|Species), data = flo.compare, glmerControl(calc.derivs=F, optCtrl=list(maxfun=10^4)))
# plot(fitted(m.flo), resid(m.flo))
# qqnorm(resid(m.flo))
# qqline(resid(m.flo), col = 2, lwd = 2, lty = 2)
# summary(m.flo) # option 2 has marginally significantly higher seed set but only driven by AGHE
# 
# # AGHE
# m.flo <- lmer(log(avg.flo) ~ option + (1|Treat.Code/Subplot), data = flo.compare[flo.compare$Species == "Agoseris heterophylla",])
# plot(fitted(m.flo), resid(m.flo))
# qqnorm(resid(m.flo))
# qqline(resid(m.flo), col = 2, lwd = 2, lty = 2)
# summary(m.flo) # option 2 has significantly higher seed set than option two for AGHE
# 
# # HECO
# m.flo <- lmer(log(avg.flo) ~ option + (1|Treat.Code/Subplot), data = flo.compare[flo.compare$Species == "Hemizonia congesta",])
# plot(fitted(m.flo), resid(m.flo))
# qqnorm(resid(m.flo))
# qqline(resid(m.flo), col = 2, lwd = 2, lty = 2)
# summary(m.flo) # nope
# 
# # PLER
# m.flo <- lmer(log(avg.flo) ~ option + (1|Treat.Code/Subplot), data = flo.compare[flo.compare$Species == "Plantago erecta",])
# plot(fitted(m.flo), resid(m.flo))
# qqnorm(resid(m.flo))
# qqline(resid(m.flo), col = 2, lwd = 2, lty = 2)
# summary(m.flo) # nope
# 
# # CLPU
# m.flo <- lmer(log(avg.flo) ~ option + (1|Treat.Code/Subplot), data = flo.compare[flo.compare$Species == "Clarkia purpurea",])
# plot(fitted(m.flo), resid(m.flo))
# qqnorm(resid(m.flo))
# qqline(resid(m.flo), col = 2, lwd = 2, lty = 2)
# summary(m.flo) # OPTION 2 HAS MARGINALLY HIGHER SEED SET IN CLPU
# 
# # LACA
# m.flo <- lmer(log(avg.flo) ~ option + (1|Treat.Code/Subplot), data = flo.compare[flo.compare$Species == "Lasthenia californica",])
# plot(fitted(m.flo), resid(m.flo))
# qqnorm(resid(m.flo))
# qqline(resid(m.flo), col = 2, lwd = 2, lty = 2)
# summary(m.flo) # nope

# final datasets

# all flowering date
flo.16.long <- flo.16.long[,c(1:5,7)]
flo.16.long <- merge(flo.16.long, unique(treat[,c(1,3)]), by = "Plot")

# flo.16.summary = peak flowering
rm(flo.16, flo.dates)

## how does number of flower per individual vary by treatment?

flo.16.long.sum <- summarySE(flo.16.long, measurevar = "n.flo", groupvars = c("Treat.Code", "Subplot", "Species"))

flo.16.long.sum$Treat.Code <- factor(flo.16.long.sum$Treat.Code, levels = c("D", "C", "W"))

ggplot(flo.16.long.sum, aes(x = Treat.Code, y = n.flo, col = Subplot, group = Subplot)) +
  geom_point() +
  geom_errorbar(aes(ymin = n.flo - se, ymax = n.flo + se), width = 0.02) +
  geom_line() +
  facet_wrap(~ Species) +
  theme_classic() # most species dont alter flower number

#### 2016: Seed ####
# Get rid of plots with burrow damage early in the season that ended up decreasing grass cover
## Early burrow damage: 85, 86, 87
## Post mortality damage: 15, not removing
seed.16 <- filter(seed.16, Plot != 87, Plot != 85, Plot != 86)

seed.16$Subplot <- revalue(seed.16$Subplot, c("A" = "Grass", "B" = "No Grass"))
seed.16$Species <- revalue(seed.16$Species, c("AGHE" = "Agoseris heterophylla", "CLGR" = "Clarkia purpurea", "LACA" = "Lasthenia californica", "PLER" = "Plantago erecta", "HECO" = "Hemizonia congesta", "AVFA" = "Avena fatua", "LOMU" = "Lolium multiflorum", "BRHO" = "Bromus hordeaceus", "TACA" = "Taeniatherum caput-medusae", "VUMI" = "Vulpia microstachys"))

seed.16 <- filter(seed.16, used != "no", Species != "LUBI") # remove data that isn't good: this includes seedheads that could not be accurately counted, cases where no heads were found and envelope had seeds but included a note that some had already been dispersed etc
seed.16$Species <- factor(seed.16$Species)
seed.16[is.na(seed.16)] <- 0
seed.16$n.seed.inf <- ifelse(seed.16$Amy.redo > 0, seed.16$Amy.redo, seed.16$n.seed.inf) # replace seed est with those redone by Amy last summer after discovering that the previous counts were way overestimated

seed.16$Year <- 2016
seed.16 <- seed.16[,c(13,2,4,5,3,6)]

#  creating rowID, USE ALL DATA! so get rid of date
seed.16 <- ddply(seed.16, .(Year, Plot, Subplot, Species), mutate, n.ind = 1:NROW(n.seed.inf))
#seed.16.w <- seed.16 %>% spread(n.ind, n.seed.inf)
seed.16 <- merge(seed.16, unique(treat[,c(1,3)]))

# Explore seed set in graphs, check outliers
#boxplot
seed.16$Treat.Code <- factor(seed.16$Treat.Code, levels = c("D", "C", "W"))
ggplot(seed.16, aes(x = Treat.Code, y = n.seed.inf, col = Subplot, by = Subplot)) +
  geom_boxplot() +
  facet_wrap(~ Species) +
  theme_classic()

seed.16.sum <- summarySE(seed.16, measurevar = "n.seed.inf", groupvars = c("Treat.Code", "Subplot", "Species"))

seed.16.sum$Treat.Code <- factor(seed.16.sum$Treat.Code, levels = c("D", "C", "W"))

ggplot(filter(seed.16.sum, Species != "Avena fatua", Species != "Bromus hordeaceus", Species != "Lolium multiflorum", Species != "Taeniatherum caput-medusae", Species != "Vulpia microstachys"), aes(x = Treat.Code, y = n.seed.inf, col = Subplot, group = Subplot)) +
  geom_point() +
  geom_errorbar(aes(ymin = n.seed.inf - se, ymax = n.seed.inf + se), width = 0.02) +
  geom_line() +
  facet_wrap(~ Species) +
  theme_classic() # seeds per inf only changes in SA

#### 2016: Flower + Seed ####
# get avg seed set per individual by combining two average datasets
#seed.avg <- ddply(seed.16, .(Year, Treat.Code, Plot, Subplot, Species), summarize, avg.seed = mean(n.seed.inf), se = sd(n.seed.inf), n.ind.sd = length(n.seed.inf))
seed.avg <- ddply(seed.16, .(Year, Plot, Subplot, Species), summarize, avg.seed = mean(n.seed.inf))


seed.flo.16 <- merge(flo.16.summary, seed.avg, by = c("Year", "Plot", "Subplot", "Species"), all = T)
seed.flo.16 <- filter(seed.flo.16, Species != "Avena fatua", Species != "Bromus hordeaceus", Species != "Lolium multiflorum", Species != "Taeniatherum caput-medusae", Species != "Vulpia microstachys")
seed.flo.16[seed.flo.16$Species == "Plantago erecta",]$avg.seed <- 2 # some missing, only ever has 2 per flower so fill in
seed.flo.16[seed.flo.16$Species == "Hemizonia congesta",]$avg.seed <- ifelse(is.na(seed.flo.16[seed.flo.16$Species == "Hemizonia congesta",]$avg.seed) == T, mean(seed.flo.16[seed.flo.16$Species == "Hemizonia congesta",]$avg.seed, na.rm = T), seed.flo.16[seed.flo.16$Species == "Hemizonia congesta",]$avg.seed) # HECO only ever has 5 seeds per flower, replace with average because sometimes they are eaten
seed.flo.16 <- merge(seed.flo.16, unique(treat[,c(1,3)]))
seed.flo.16$n.seed.ind <- seed.flo.16$avg.seed*seed.flo.16$avg.flo
rm(seed.avg)

# Other weird outliers:
## Plot 19G-LACA: two individuals one with 6 flowers and another with 4 driving the pattern
## Plot 84G-HECO: just weird, no one individual driving the pattern
# seed.wo.out <- filter(seed.flo.16, !(Species == "Lasthenia californica" & Plot == 19 & Subplot == "Grass"))
# seed.wo.out <- filter(seed.wo.out, !(Species == "Hemizonia congesta" & Plot == 84 & Subplot == "Grass"))

seed.16.sum <- summarySE(seed.flo.16, measurevar = "n.seed.ind", groupvars = c("Treat.Code", "Subplot", "Species"), na.rm = T)

seed.16.sum$Treat.Code <- factor(seed.16.sum$Treat.Code, levels = c("D", "C", "W"))

ggplot(seed.16.sum, aes(x = Treat.Code, y = n.seed.ind, col = Subplot, group = Subplot)) +
  geom_point() +
  geom_errorbar(aes(ymin = n.seed.ind - se, ymax = n.seed.ind + se), width = 0.02) +
  geom_line() +
  facet_wrap(~ Species) +
  theme_classic() # still one weird HECO outlier and one weird LACA outlier but overall good, i'd like to remove at least the one LACA but no good reason to at the moment so its staying in for now

#boxplot
seed.flo.16$Treat.Code <- factor(seed.flo.16$Treat.Code, levels = c("D", "C", "W"))
ggplot(seed.flo.16, aes(x = Treat.Code, y = n.seed.ind, col = Subplot, by = Subplot)) +
  geom_boxplot() +
  facet_wrap(~ Species) +
  theme_classic()

set.treat <- summarySE(seed.flo.16, measurevar = "n.seed.ind", groupvars = c("Treat.Code", "Subplot"), na.rm = T)

ggplot(set.treat, aes(x = Treat.Code, y = n.seed.ind, col = Subplot, group = Subplot)) +
  geom_point() +
  geom_errorbar(aes(ymin = n.seed.ind - se, ymax = n.seed.ind + se), width = 0.02) +
  geom_line()

###
# Grass Cover
###
grass <- read.csv("Data/Census-Data/2015-2016/Grass-Cover-2016.csv")
grass$Subplot <- revalue(grass$Subplot, c("A" = "Grass"))
grass <- merge(grass, unique(treat[,c(1,3)]), by = "Plot")
names(grass)[3] <- "cov.g"

final.16 <- merge(seed.flo.16, grass, by = c("Plot", "Subplot", "Treat.Code"), all.x = T)

grass$Treat.Code <- factor(grass$Treat.Code, levels = c("D", "C", "W"))
ggplot(grass, aes(x = Treat.Code, y = cov.g)) +
  geom_boxplot()

grass.sum <- summarySE(grass, measurevar = "cov.g", groupvars = "Treat.Code")
ggplot(grass.sum, aes(x = Treat.Code, y = cov.g)) +
  geom_point() +
  geom_errorbar(aes(ymin = cov.g - se, ymax = cov.g + se), width = 0.02)

final.16$Treat.Code <- factor(final.16$Treat.Code, levels = c("D", "C", "W"))

# # HECO
# heco.g <- ggplot(filter(final.16, Species == "Hemizonia congesta"), aes(x = cov.g, y = n.seed.ind, col = Treat.Code, group = Treat.Code)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = F) +
#   facet_wrap(~ Species) +
#   theme_classic() 
# 
# # AGHE
# aghe.g <- ggplot(filter(final.16, Species == "Agoseris heterophylla"), aes(x = cov.g, y = n.seed.ind, col = Treat.Code, group = Treat.Code)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = F) +
#   facet_wrap(~ Species) +
#   theme_classic() 
# 
# # PLER
# pler.g <- ggplot(filter(final.16, Species == "Plantago erecta"), aes(x = cov.g, y = n.seed.ind, col = Treat.Code, group = Treat.Code)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = F) +
#   facet_wrap(~ Species) +
#   theme_classic() 
# 
# # LACA
# laca.g <- ggplot(filter(final.16, Species == "Lasthenia californica"), aes(x = cov.g, y = n.seed.ind, col = Treat.Code, group = Treat.Code)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = F) +
#   facet_wrap(~ Species) +
#   theme_classic() 
# 
# # CLPU
# clpu.g <- ggplot(filter(final.16, Species == "Clarkia purpurea"), aes(x = cov.g, y = n.seed.ind, col = Treat.Code, group = Treat.Code)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = F) +
#   facet_wrap(~ Species) +
#   theme_classic() 
# 
# grid.arrange(pler.g, heco.g, laca.g, aghe.g, clpu.g, nrow = 2)
# 
# rm(aghe.g, clpu.g, pler.g, heco.g, laca.g, grass, flo.16.summary, set.treat, seed.flo.16, seed.16.sum, seed.16.w, grass.sum)

#### 2017: Flowering ####
flo.17$Date <- as.Date(flo.17$Date, "%m/%d/%y")

# some surveys were conducted really late and species were recorded as "done" instead of "setting seed. replace setting seed wtih # done
flo.17$flo.C <- ifelse(flo.17$flo.C == 0 & rowSums(flo.17[ , 9:13]) != 0, flo.17$N.done, flo.17$flo.C)

# zeroes not meaningful, get rid of and calculate various summary stats 
flo.17.long <- melt(flo.17, id.vars = c("Year","Date","Plot","Subplot","Species"), measure.vars = c("flo.1","flo.2","flo.3","flo.4","flo.5"), variable_name = "individual")
names(flo.17.long)[7] <- "n.flo"
flo.17.long <- filter(flo.17.long, n.flo > 0)

# calculate summary stats
flo.17.summary <- ddply(flo.17.long, .(Year, Date, Plot, Subplot, Species), summarize, avg.flo = mean(n.flo), sd.flo = sd(n.flo), n.ind = length(n.flo))

###
# Option 2: find date with most individuals (only up to 5), then find date with max avg.flo
###
# test <- ddply(flo.17.summary, .(Year, Plot, Subplot, Species), summarize, Date.1 = Date[which(n.ind == max(n.ind))][1], Date.2 = Date[which(n.ind == max(n.ind))][2]) # only 2 dates
# 
# test <- melt(test, id.vars = c("Year", "Plot", "Subplot", "Species"))
# 
# names(test)[6] <- "Date"
# 
# test <- merge(flo.17.summary, test[,c(1:4,6)], by = c("Year", "Date",  "Plot", "Subplot", "Species"), all.x = F)
# 
# test <- ddply(test, .(Year, Plot, Subplot, Species), summarize, Date = max(Date[which(avg.flo == max(avg.flo))]))
# 
# flo.17.summary.2 <- merge(flo.17.summary, test)

###
# Option 1: Date of highest number of ind. flowering, break tie with latest date of highest number ind setting seed
###

## known issues: loads of data not taken because plants were done so there are flowers/ind data but not number of individuals in each stage

flo.17.summary <- merge(flo.17[,c(1:8)], flo.17.summary, by = c("Year", "Date", "Plot", "Subplot", "Species"))

# flo.dates <- ddply(flo.17.summary, .(Year, Plot, Subplot, Species), summarize, Date = ifelse(length(Date[which(flo.C == max(flo.C))]) > 1,
#                             max(Date[which(flo.B == max(flo.B))]),
#                             Date[which(flo.C == max(flo.C))]))

flo.dates <- ddply(flo.17.summary, .(Year, Plot, Subplot, Species), summarize, Date =
                     ifelse(length(Date[which(flo.B == max(flo.B))]) > 1,
                            max(Date[which(flo.C == max(flo.C))]),
                            Date[which(flo.B == max(flo.B))]))

class(flo.dates$Date) <- class(flo.17.summary$Date)


# extract flowering data from both long and summary dataset
flo.17.long <- merge(flo.dates, flo.17.long, by = c("Year", "Date", "Plot", "Subplot", "Species"), all.y = F)

flo.17.summary <- merge(flo.dates, flo.17.summary, by = c("Year", "Date", "Plot", "Subplot", "Species"), all.y = F)
# 
# diff <- setdiff(flo.17.summary[,c(1:5,9,11)], flo.17.summary.2[,c(1:6,8)]) # 45 diffs
# 
# # compare options to see if this affects flowering
# flo.17.summary$option <- 1
# flo.17.summary.2$option <- 2
# flo.compare <- rbind(flo.17.summary[,c(1,3,4,5,9,12)], flo.17.summary.2[,c(1,3,4,5,6,9)])
# flo.compare <- merge(flo.compare, unique(treat[,c(1,3)]), by = "Plot")
# 
# flo.compare$int.avg.flo <- as.integer(flo.compare$avg.flo)
# hist(log(flo.compare$avg.flo))
# flo.compare$option <- as.factor(flo.compare$option)
# m.flo <- glmer.nb(int.avg.flo ~ option + (1|Treat.Code/Subplot) + (1 * Treat.Code|Species), data = flo.compare, glmerControl(calc.derivs=F, optCtrl=list(maxfun=10^4)))
# plot(fitted(m.flo), resid(m.flo))
# qqnorm(resid(m.flo))
# qqline(resid(m.flo), col = 2, lwd = 2, lty = 2)
# summary(m.flo) # no sig difference between options
# 
# # AGHE
# m.flo <- lmer(log(avg.flo) ~ option + (1|Treat.Code/Subplot), data = flo.compare[flo.compare$Species == "Agoseris heterophylla",])
# plot(fitted(m.flo), resid(m.flo))
# qqnorm(resid(m.flo))
# qqline(resid(m.flo), col = 2, lwd = 2, lty = 2)
# summary(m.flo) # nope
# 
# # HECO
# m.flo <- lmer(log(avg.flo) ~ option + (1|Treat.Code/Subplot), data = flo.compare[flo.compare$Species == "Hemizonia congesta",])
# plot(fitted(m.flo), resid(m.flo))
# qqnorm(resid(m.flo))
# qqline(resid(m.flo), col = 2, lwd = 2, lty = 2)
# summary(m.flo) # nope
# 
# # PLER
# m.flo <- lmer(log(avg.flo) ~ option + (1|Treat.Code/Subplot), data = flo.compare[flo.compare$Species == "Plantago erecta",])
# plot(fitted(m.flo), resid(m.flo))
# qqnorm(resid(m.flo))
# qqline(resid(m.flo), col = 2, lwd = 2, lty = 2)
# summary(m.flo) # nope
# 
# # CLPU
# hist(flo.compare[flo.compare$Species == "Clarkia purpurea",]$int.avg.flo)
# m.flo <- lmer(log(avg.flo) ~ option + (1|Treat.Code/Subplot), data = flo.compare[flo.compare$Species == "Clarkia purpurea",])
# plot(fitted(m.flo), resid(m.flo))
# qqnorm(resid(m.flo))
# qqline(resid(m.flo), col = 2, lwd = 2, lty = 2)
# summary(m.flo) # crazy long tail
# 
# # LACA
# m.flo <- lmer(log(avg.flo) ~ option + (1|Treat.Code/Subplot), data = flo.compare[flo.compare$Species == "Lasthenia californica",])
# plot(fitted(m.flo), resid(m.flo))
# qqnorm(resid(m.flo))
# qqline(resid(m.flo), col = 2, lwd = 2, lty = 2)
# summary(m.flo) # nope
# 
# # CAPA
# m.flo <- lmer(log(avg.flo) ~ option + (1|Treat.Code/Subplot), data = flo.compare[flo.compare$Species == "Calycadenia pauciflora",])
# plot(fitted(m.flo), resid(m.flo))
# qqnorm(resid(m.flo))
# qqline(resid(m.flo), col = 2, lwd = 2, lty = 2)
# summary(m.flo) # nope
# 
# ### no difference between options

# final datasets

# all flowering date
flo.17.long <- flo.17.long[,c(1:5,7)]
flo.17.long <- merge(flo.17.long, unique(treat[,c(1,3)]), by = "Plot")

# flo.17.summary = peak flowering
# flo.17.summary2 = peak avg flower per individuals
rm(flo.17, flo.dates)

## how does number of flower per individual vary by treatment?

flo.17.long.sum <- summarySE(flo.17.long, measurevar = "n.flo", groupvars = c("Treat.Code", "Subplot", "Species"))

flo.17.long.sum$Treat.Code <- factor(flo.17.long.sum$Treat.Code, levels = c("D", "C", "W"))

ggplot(flo.17.long.sum, aes(x = Treat.Code, y = n.flo, col = Subplot, group = Subplot)) +
  geom_point() +
  geom_errorbar(aes(ymin = n.flo - se, ymax = n.flo + se), width = 0.02) +
  geom_line() +
  facet_wrap(~ Species) +
  theme_classic() 

ggplot(flo.17.long.sum, aes(x = Treat.Code, y = n.flo, col = Subplot, group = Subplot)) +
  geom_boxplot() +
  facet_wrap(~ Species) +
  theme_classic() # no big outliers

#### 2017: Seed ####
seed.17$Subplot <- revalue(seed.17$Subplot, c("N" = "No Grass", "G" = "Grass", "T" = "Thatch"))

seed.17$Species <- revalue(seed.17$Species, c("AGHE" = "Agoseris heterophylla", "CLPU" = "Clarkia purpurea", "LACA" = "Lasthenia californica", "HECO" = "Hemizonia congesta", "AVFA" = "Avena fatua", "TACA" = "Taeniatherum caput-medusae", "BRHO" = "Bromus hordeaceus", "LOMU" = "Lolium multiflorum", "CAPA" = "Calycadenia pauciflora"))

seed.17 <- filter(seed.17, used != "no")
seed.17 <- seed.17[complete.cases(seed.17),]
colnames(seed.17)[6] <- "n.seed.inf"
seed.17$Year <- 2017
seed.17 <- seed.17[,c(9,1:4,6)] 

#  creating rowID, USE ALL DATA! so get rid of date
seed.17 <- ddply(seed.17, .(Year, Plot, Subplot, Species), mutate, n.ind = 1:NROW(n.seed.inf))
seed.17.w <- seed.17 %>% spread(n.ind, n.seed.inf)
seed.17 <- merge(seed.17, unique(treat[,c(1,3)]))

# Explore seed set in graphs

seed.17.sum <- summarySE(seed.17, measurevar = "n.seed.inf", groupvars = c("Treat.Code", "Subplot", "Species"))

seed.17.sum$Treat.Code <- factor(seed.17.sum$Treat.Code, levels = c("D", "C", "W"))

ggplot(filter(seed.17.sum, Species != "Avena fatua", Species != "Bromus hordeaceus", Species != "Lolium multiflorum", Species != "Taeniatherum caput-medusae"), aes(x = Treat.Code, y = n.seed.inf, col = Subplot, group = Subplot)) +
  geom_point() +
  geom_errorbar(aes(ymin = n.seed.inf - se, ymax = n.seed.inf + se), width = 0.02) +
  geom_line() +
  facet_wrap(~ Species) +
  theme_classic() # no big outliers

ggplot(filter(seed.17.sum, Species != "Avena fatua", Species != "Bromus hordeaceus", Species != "Lolium multiflorum", Species != "Taeniatherum caput-medusae"), aes(x = Treat.Code, y = n.seed.inf, col = Subplot, by = Subplot)) +
  geom_point() +
  facet_wrap(~ Species) +
  theme_classic() # no big outliers

#### 2017: Flower + Seed ####
# get avg seed set per individual by combining two average datasets
seed.avg <- ddply(seed.17, .(Year, Plot, Subplot, Species), summarize, avg.seed = mean(n.seed.inf))

# Didn't collect PLER this year so input 2 for every PLER/Subplot/Treat.Code Combo, then when combined with flowering data those without PLER will be dropped
combos <- data.frame(Year = rep(2017, 90), Plot = rep(unique(treat[,1]), each = 3), Species = rep("Plantago erecta", 90), Subplot = rep(c("No Grass", "Grass", "Thatch"), 30), avg.seed = 2)

combos <- combos[,c(1, 2, 4, 3, 5)]
seed.avg <- rbind(seed.avg, combos)

seed.flo.17 <- merge(flo.17.summary, seed.avg, by = c("Year", "Plot", "Subplot", "Species"), all = T)
seed.flo.17 <- filter(seed.flo.17, Species != "Avena fatua", Species != "Bromus hordeaceus", Species != "Lolium multiflorum", Species != "Taeniatherum caput-medusae", Species != "Vulpia microstachys")
seed.flo.17 <- merge(seed.flo.17, unique(treat[,c(1,3)]))
# HECO only ever has 5 seeds per flower, replace with average because sometimes they are eaten
seed.flo.17[seed.flo.17$Species == "Hemizonia congesta",]$avg.seed <- ifelse(is.na(seed.flo.17[seed.flo.17$Species == "Hemizonia congesta",]$avg.seed) == T, mean(seed.flo.17[seed.flo.17$Species == "Hemizonia congesta",]$avg.seed, na.rm = T), seed.flo.17[seed.flo.17$Species == "Hemizonia congesta",]$avg.seed) 

# CAPA only ever has max 5 seeds per flower, replace with average because sometimes they are eaten
seed.flo.17[seed.flo.17$Species == "Calycadenia pauciflora",]$avg.seed <- ifelse(is.na(seed.flo.17[seed.flo.17$Species == "Calycadenia pauciflora",]$avg.seed) == T, mean(seed.flo.17[seed.flo.17$Species == "Calycadenia pauciflora",]$avg.seed, na.rm = T), seed.flo.17[seed.flo.17$Species == "Calycadenia pauciflora",]$avg.seed)


seed.flo.17$n.seed.ind <- seed.flo.17$avg.seed*seed.flo.17$avg.flo

rm(seed.avg, combos)


seed.17.sum <- summarySE(seed.flo.17, measurevar = "n.seed.ind", groupvars = c("Treat.Code", "Subplot", "Species"), na.rm = T)

seed.17.sum$Treat.Code <- factor(seed.17.sum$Treat.Code, levels = c("D", "C", "W"))

ggplot(seed.17.sum, aes(x = Treat.Code, y = n.seed.ind, col = Subplot, group = Subplot)) +
  geom_point() +
  geom_errorbar(aes(ymin = n.seed.ind - se, ymax = n.seed.ind + se), width = 0.02) +
  geom_line() +
  facet_wrap(~ Species) +
  theme_classic() # maybe a weird outlier for HECO 

#boxplot - check outliers
seed.flo.17$Treat.Code <- factor(seed.flo.17$Treat.Code, levels = c("D", "C", "W"))
ggplot(seed.flo.17, aes(x = Treat.Code, y = n.seed.ind, col = Subplot, by = Subplot)) +
  geom_boxplot() +
  facet_wrap(~ Species) +
  theme_classic() # nothing really suspicious, a few AGHE and HECO, esp drought plot 85, possibly burrow damage for 85 but hard to tell when it happened

###
# Grass Cover
###
grass <- read.csv("Data/Census-Data/2016-2017/Grass-Cover-2017.csv")[,1:4]
grass$Subplot <- revalue(grass$Subplot, c("G" = "Grass", "T" = "Thatch"))
grass <- merge(grass, unique(treat[,c(1,3)]), by = "Plot")
grass <- ddply(grass, .(Plot, Subplot, Treat.Code), summarize, cov.g = sum(Cover))

final.17 <- merge(seed.flo.17, grass, by = c("Plot", "Subplot", "Treat.Code"), all.x = T)

# grass$Treat.Code <- factor(grass$Treat.Code, levels = c("D", "C", "W"))
# ggplot(grass, aes(x = Treat.Code, y = cov.g, by = Subplot, fill = Subplot)) +
#   geom_boxplot() # in general lower grass cover in thatch treatments
# 
# grass.sum <- summarySE(grass, measurevar = "cov.g", groupvars = c("Treat.Code", "Subplot"))
# ggplot(grass.sum, aes(x = Treat.Code, y = cov.g, group = Subplot, col = Subplot)) +
#   geom_point() +
#   geom_errorbar(aes(ymin = cov.g - se, ymax = cov.g + se), width = 0.02) +
#   geom_line() # weird interaction, no differences in grass cover between watering treatments, but in watered plots larger dif between thatch and grass? 
# 
# # HECO
# (heco.g <- ggplot(filter(final.17, Species == "Hemizonia congesta", Subplot != "No Grass"), aes(x = cov.g, y = n.seed.ind, col = Treat.Code, group = Treat.Code)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = F) +
#   facet_wrap(~ Subplot) +
#   theme_classic() ) # one weird drought outlier
# 
# # AGHE
# (aghe.g <- ggplot(filter(final.17, Species == "Agoseris heterophylla", Subplot != "No Grass"), aes(x = cov.g, y = n.seed.ind, col = Treat.Code, group = Treat.Code)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = F) +
#   facet_wrap(~ Subplot) +  
#   theme_classic() ) # same trend with agoseris: strong negative effect of grass that disappears in watered plots
# 
# # PLER - not as affect by grass
# (pler.g <- ggplot(filter(final.17, Species == "Plantago erecta", Subplot != "No Grass"), aes(x = cov.g, y = n.seed.ind, col = Treat.Code, group = Treat.Code)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = F) +
#   facet_wrap(~ Subplot) +
#   theme_classic() )
# 
# # LACA
# (laca.g <- ggplot(filter(final.17, Species == "Lasthenia californica", Subplot != "No Grass"), aes(x = cov.g, y = n.seed.ind, col = Treat.Code, group = Treat.Code)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = F) +
#   facet_wrap(~ Subplot) +
#   theme_classic() )
# 
# # CLPU
# (clpu.g <- ggplot(filter(final.17, Species == "Clarkia purpurea", Subplot != "No Grass"), aes(x = cov.g, y = n.seed.ind, col = Treat.Code, group = Treat.Code)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = F) +
#   facet_wrap(~ Subplot) +
#   theme_classic() ) # not enough data, one huge outlier in thatch
# 
# # CAPA
# (capa.g <- ggplot(filter(final.17, Species == "Calycadenia pauciflora", Subplot != "No Grass"), aes(x = cov.g, y = n.seed.ind, col = Treat.Code, group = Treat.Code)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = F) +
#   facet_wrap(~ Subplot) +
#   theme_classic() ) # drought and control not enough data for thatch, watering has positive effect on seeds per ind that increases with grass, huge neg effect of grass but only in control plots? 

#grid.arrange(pler.g, heco.g, laca.g, aghe.g, clpu.g, capa.g, nrow = 2)

#rm(aghe.g, clpu.g, pler.g, heco.g, laca.g, capa.g, grass, grass.sum, flo.17.long.sum, flo.17.summary, seed.17.sum, seed.17.w, seed.flo.17, treat)

# final cleaned flower/seedset/grass dataset: final.17, flo.17.long, seed.17

#### Combine years together and add viability ####
names(final.16)[13] <- "avg.seed.inf"
names(final.17)[13] <- "avg.seed.inf"

final.16 <- merge(final.16, viab[1:5, c(1,3)])
final.17 <- merge(final.17, viab[c(1, 6:10), c(1,3)])

names(final.17)[7:9] <- c("flo.a", "flo.b", "flo.c")
final <- rbind(final.16, final.17)

flo.long <- rbind(flo.16.long, flo.17.long)

seed.16 <- merge(seed.16, viab[1:5, c(1,3)])
seed.17 <- merge(seed.17, viab[c(1, 6:12), c(1,3)])
seed <- rbind(seed.16, seed.17)

rm(final.16, final.17, flo.16.long, flo.17.long, seed.16, seed.17)

write.table(final, "Data/Cleaned Data for Analysis/final-flo-seed.csv", sep = ",", row.names = F)
write.table(flo.long, "Data/Cleaned Data for Analysis/final-flo-long.csv", sep = ",", row.names = F)
write.table(seed, "Data/Cleaned Data for Analysis/final-seed-long.csv", sep = ",", row.names = F)

#### Treatment effects on grasses ####
ggplot(filter(seed.16.sum, Species == "Avena fatua" | Species == "Bromus hordeaceus"| Species == "Lolium multiflorum"| Species == "Taeniatherum caput-medusae"), aes(x = Treat.Code, y = n.seed.inf)) +
  geom_point() +
  geom_errorbar(aes(ymin = n.seed.inf - se, ymax = n.seed.inf + se), width = 0.02) +
  geom_line() +
  facet_wrap(~ Species) +
  theme_classic() 

grass.seed <- ddply(filter(seed.16, Species == "Avena fatua" | Species == "Bromus hordeaceus"| Species == "Lolium multiflorum"| Species == "Taeniatherum caput-medusae"), .(Plot, Year, Species, Treat.Code), summarize, n.seed.inf = mean(n.seed.inf))
grass.seed$Treat.Code <- factor(grass.seed$Treat.Code, levels = c("C", "D", "W"))
hist(grass.seed$n.seed.inf)
hist(log(grass.seed$n.seed.inf))
mg.16 <- lmer(log(n.seed.inf) ~ Treat.Code + (1|Plot), data = grass.seed)
plot(fitted(mg.16), resid(mg.16))
qqnorm(resid(mg.16)) 
qqline(resid(mg.16), col = 2, lwd = 2, lty = 2) 
summary(mg.16) # no difs

ggplot(filter(seed.17.sum, Species == "Avena fatua" | Species == "Bromus hordeaceus"| Species == "Lolium multiflorum"| Species == "Taeniatherum caput-medusae", Subplot == "Grass"), aes(x = Treat.Code, y = n.seed.inf)) +
  geom_point() +
  geom_errorbar(aes(ymin = n.seed.inf - se, ymax = n.seed.inf + se), width = 0.02) +
  geom_line() +
  facet_wrap(~ Species) +
  theme_classic() 

grass.seed.17 <- ddply(filter(seed.17, Species == "Avena fatua" | Species == "Bromus hordeaceus"| Species == "Lolium multiflorum"| Species == "Taeniatherum caput-medusae", Subplot == "Grass"), .(Plot, Year, Species, Treat.Code), summarize, n.seed.inf = mean(n.seed.inf))
grass.seed.17$Treat.Code <- factor(grass.seed.17$Treat.Code, levels = c("C", "D", "W"))
grass.seed <- rbind(grass.seed, grass.seed.17)
hist(grass.seed$n.seed.inf)
hist(log(grass.seed$n.seed.inf))

mg <- lmer(log(n.seed.inf) ~ Treat.Code + (1|Species) + (1|Year), data = grass.seed)
plot(fitted(mg), resid(mg))
qqnorm(resid(mg)) 
qqline(resid(mg), col = 2, lwd = 2, lty = 2) 
summary(mg) # no difs