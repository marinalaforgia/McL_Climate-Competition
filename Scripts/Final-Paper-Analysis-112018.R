### FINAL PAPER SCRIPT ###

# Chapter 2 analysis 11-20-2018 #

rm(list=ls()) 

#### Load Libraries ####
library(Rmisc)
library(ggplot2)
library(plyr)
library(dplyr)
library(lme4)
library(lmerTest)
library(arm)
library(coefplot)
library(reshape2)
library(tidyr)
library(cowplot)
library(car)

#### Load: Bootstraps ####
## load previous boots
# CI boots
load("Data/Analysis-output/CI/CI-1.Rdata") # mortality
load("Data/Analysis-output/CI/CI-2.Rdata") # seed set
load("Data/Analysis-output/CI/CI-3.Rdata") # seed bank carryover
load("Data/Analysis-output/CI/CI-4.Rdata") # germination

# Parameter boots
load("Data/Analysis-output/Param-sims/BS-m.Rdata") # mortality
load("Data/Analysis-output/Param-sims/BS-F.Rdata") # seed set
load("Data/Analysis-output/Param-sims/BS-g.Rdata") # germination
load("Data/Analysis-output/20190206-lambda-sim.Rdata")

#### Prep: Vital Rates ####

###
# Grass
###
grass <- read.csv("Data/Post-Processing/grass-cover.csv")
names(grass)[3] <- "cover.g"

treat <- read.csv("Data/Setup/Marina-Treatment-30.csv")

###
# Germ and mortality data. plots 85, 86, 87 have been removed due to burrow damage
###
dem.16 <- read.csv("Data/Post-Processing/dem-data-16.csv")
dem.16 <- filter(dem.16, Plot != 87, Plot != 85, Plot != 86) # removed plots 85, 86, and 87 due to burrow damage early in the season

dem.17 <- read.csv("Data/Post-Processing/dem-data-17.csv")

dem <- rbind(dem.16, dem.17)
dem <- merge(dem, grass, by = c("Year", "Plot", "Subplot", "Treat.Code"), all.x = T)
dem <- dem[, -c(8,9,12)] # get rid of thatch for now and extra columns
dem <- dem[dem$Subplot != "Thatch", ] 
dem <- filter(dem, !(Year == 2017 & Treat.Code == "D")) # shelters in 2017 did not apply a drought
dem$Subplot <- factor(dem$Subplot, levels = c("No Grass", "Grass"))
dem$p.mort <- dem$tot.mort/dem$germ.proj
dem$p.germ <- dem$germ.tot/dem$viable
rm(dem.16, dem.17)

###
# Flowering/seed set data 
###

# 10/28/2018 NOTE: changed code for final-flo-seed.csv to take number of flowers per individual from date of highest number flowering 
#flo.long <- read.csv("Data/Post-Processing/final-flo-long.csv")
#seed.long <- read.csv("Data/Post-Processing/final-seed-long.csv")
flo.seed <- read.csv("Data/Post-Processing/final-flo-seed.csv")
flo.seed <- flo.seed[flo.seed$Subplot != "Thatch",]
flo.seed <- filter(flo.seed, !(Year == 2017 & Treat.Code == "D")) 
flo.seed$Subplot <- factor(flo.seed$Subplot, levels = c("No Grass", "Grass"))
#flo.seed <- filter(flo.seed, n.seed.ind != 0) # only removes 3 scenarios
flo.seed$n.seed.ind <- flo.seed$n.seed.ind * flo.seed$p.viable # adjust for viability
flo.seed <- flo.seed[,-16]


###
# Seed survival
###
sb <- read.csv("Data/Post-Processing/seed-carryover-plot.csv")[,c(1:9)]
sb <- filter(sb, !(Year == 2017 & Treat.Code == "D"))
names(sb)[1] <- "Species"
sb <- sb[,-c(5,6)]

sb.spp <- ddply(sb, .(Species), summarize, p.sd.surv = mean(p.surv))

#### Prep: Traits ####
trait.w <- read.csv("Data/Post-Processing/final-traits-w.csv")
trait.w <- trait.w[,-c(5,7,9)]

sla.13 <- read.csv("Data/Post-Processing/final-sla-13c.csv")[,c(1,3,5)]

PCA.G <- prcomp(trait.w[,c(18, 20)], scale = T) # Greenhouse trait
PCA.F <- prcomp(trait.w[, c(7, 8)], scale = T) # Field trait
PCA.s13 <- prcomp(sla.13[, c(2, 3)], scale = T) # Full SLA + D13C field

# biplot(PCA.G)
# biplot(PCA.F)
# biplot(PCA.s13)

# summary(PCA.G)
# summary(PCA.F)
# summary(PCA.s13)

trait.w$PC.G <- PCA.G$x[,1] 
trait.w$PC.F <- PCA.F$x[,1] 
trait.w$PC.s13 <- PCA.s13$x[c(2,3,5,6,7,9),1] 

# merge datasets with traits
dem <- merge(dem, trait.w[,c(1,23:25)], by = "Species")
sb <- merge(sb, trait.w[,c(1,23:25)], by = "Species")
flo.seed <- merge(flo.seed, trait.w[,c(1,23:25)], by = "Species")

rm(PCA.F, PCA.G, PCA.s13, sla.13, grass)

# merge datasets with grass cover
#dem <- merge(dem, grass, all.x = T)

#### Prep: Plot-level Lambda ####

full <- merge(dem, flo.seed[,-c(6:13)], by = c("Year","Plot","Treat.Code", "Subplot","Species", "PC.F", "PC.G", "PC.s13"), all = T)

full <- merge(full, sb, by = c("Year","Plot","Treat.Code","Species","PC.F", "PC.G", "PC.s13"), all.x = T)

full$p.germ <- full$germ.tot/full$viable
full$L.sb <- full$p.surv*(1-full$p.germ)
full$L.sa <- ifelse(full$germ.proj > 0 , full$p.germ * (1 - (full$tot.mort/full$germ.proj)), 0) 
full$L.seeds <- ifelse(full$n.seed.ind >= 0, full$L.sa * full$n.seed.ind, 0)

full$L <- full$L.sb + full$L.seeds 
full <- filter(full, !(Species == "Lasthenia californica" & Year == 2016), !(Species == "Plantago erecta" & Year == 2016)) # get rid of lasthenia and plantago in 2016 beacuse of germination mix-up

#### M1: Mortality ####
dem$Year <- as.factor(dem$Year)

m1.trait <- glmer(cbind(tot.mort, germ.proj - tot.mort) ~ Treat.Code * Subplot * PC.F + (1|Year:Plot:Species), family = binomial, dem, glmerControl(calc.derivs = F))
plot(fitted(m1.trait), resid(m1.trait))
qqnorm(resid(m1.trait)) 
qqline(resid(m1.trait), col = 2, lwd = 2, lty = 2) 
summary(m1.trait) 

#### M2: Seed Set ####

# Model without 3 way interaction previously determined to be best model; both models are fine
hist(flo.seed$n.seed.ind)
hist(log(flo.seed$n.seed.ind + 1))

m2.trait <- lmer(log(n.seed.ind + 1) ~ Treat.Code * Subplot * PC.F - Treat.Code : Subplot : PC.F + (1|Year:Plot:Species), flo.seed) 
plot(fitted(m2.trait), resid(m2.trait))
qqnorm(resid(m2.trait)) 
qqline(resid(m2.trait), col = 2, lwd = 2, lty = 2) 
summary(m2.trait)

#### M3: Seed Carryover ####

m3.trait <- glmer(cbind(Count, avg.num.per.sub - Count) ~  PC.F * Treat.Code + (1|Year:Plot:Species), family = binomial, data = sb)
plot(fitted(m3.trait), resid(m3.trait))
qqnorm(resid(m3.trait))
qqline(resid(m3.trait), col = 2, lwd = 2, lty = 2)
summary(m3.trait)

#### M4: Germination ####
# removing PLER 2016 and LACA 2016 because of mis-ID
m4.trait <- glmer(cbind(germ.tot, viable-germ.tot) ~ PC.F + (1|Year:Plot:Subplot:Species), family = binomial, data = filter(dem, !(Species == "Lasthenia californica" & Year == 2016), !(Species == "Plantago erecta" & Year == 2016)))
plot(fitted(m4.trait), resid(m4.trait))
qqnorm(resid(m4.trait))
qqline(resid(m4.trait), col = 2, lwd = 2, lty = 2)
summary(m4.trait)

#### M5: Plot Lambda ####
hist(log(full$L + .5))
m5.trait <- lmer(log(L + 0.5) ~ Treat.Code * Subplot * PC.F + (1|Year:Plot:Species), data = full)
plot(fitted(m5.trait), resid(m5.trait))
qqnorm(resid(m5.trait))
qqline(resid(m5.trait), col = 2, lwd = 2, lty = 2)
summary(m5.trait)

#### Prep: Bootstrap CIs ####

df <- expand.grid(Treat.Code = c("D", "C", "W"), Subplot = c("Grass", "No Grass"), PC.F = seq(from = min(trait.w$PC.F), to = max(trait.w$PC.F), length.out = 1000)) 

# ###
# # Mortality
# ###
# 
# boot1 <- bootMer(m1.trait, function(x) predict(x, newdata = df, type = 'response', allow.new.levels = T, re.form = NA), nsim = 1000, ncpus = 2, verbose = T) # 34 did not converge (3.4%); according to Ben Bolker, this number isn't that bad (https://stackoverflow.com/questions/30235758/bootstrap-failed-using-mixed-model-in-lme4-package)
# 
# PI.boot <- data.frame(
#            p.mort = apply(boot1$t, 2, function(x) as.numeric(quantile(x, probs = .5, na.rm=TRUE))),
#            p.mort.lo = apply(boot1$t, 2, function(x) as.numeric(quantile(x, probs = .025, na.rm=TRUE))),
#            p.mort.hi = apply(boot1$t, 2, function(x) as.numeric(quantile(x, probs = .975, na.rm=TRUE))))
# 
# CI.1 <- merge(PI.boot, df, by = "row.names")
# save(CI.1, file = "Data/Analysis-output/CI/CI-1.Rdata")
# 
# ###
# # Seed Set
# ###
# 
# boot2 <- bootMer(m2.trait, function(x) predict(x, newdata = df, type = 'response', allow.new.levels = T, re.form = NA), nsim = 1000, verbose = T)
# 
# PI.boot <- data.frame(
#            n.seed.ind = apply(boot2$t, 2, function(x) as.numeric(quantile(x, probs = .5, na.rm=TRUE))),
#            seed.lo = apply(boot2$t, 2, function(x) as.numeric(quantile(x, probs = .025, na.rm=TRUE))),
#            seed.hi = apply(boot2$t, 2, function(x) as.numeric(quantile(x, probs = .975, na.rm=TRUE))))
# 
# CI.2 <- merge(PI.boot, df, by = "row.names")
# save(CI.2, file = "Data/Analysis-output/CI/CI-2.Rdata")
# 
# ###
# # Seed Carryover
# ###
# 
# boot3 <- bootMer(m3.trait, function(x) predict(x, newdata = df, type = 'response', allow.new.levels = T, re.form = NA), nsim = 1000, verbose = T) # 16/1000 did not converge
# 
# PI.boot <- data.frame(p.surv = apply(boot3$t, 2, function(x) mean(x)),
#            sb.surv.lo = apply(boot3$t, 2, function(x) as.numeric(quantile(x, probs=.05, na.rm=TRUE))),
#            sb.surv.hi = apply(boot3$t, 2, function(x) as.numeric(quantile(x, probs=.95, na.rm=TRUE))))
# 
# CI.3 <- merge(PI.boot, df, by = "row.names")
# save(CI.3, file = "Data/Analysis-output/CI/CI-3.Rdata")
# 
###
# Germination
###
#
# boot4 <- bootMer(m4.trait, function(x) predict(x, newdata = df, type = 'response', allow.new.levels = T, re.form = NA), nsim = 1000, verbose = T) #??
# 
# PI.boot <- data.frame(germ = apply(boot4$t, 2, function(x) mean(x)),
#            germ.lo = apply(boot4$t, 2, function(x) as.numeric(quantile(x, probs=.05, na.rm=TRUE))),
#            germ.hi = apply(boot4$t, 2, function(x) as.numeric(quantile(x, probs=.95, na.rm=TRUE))))
# 
# CI.4 <- merge(PI.boot, df, by = "row.names")
# save(CI.4, file = "Data/Analysis-output/CI/CI-4.Rdata")
# # 
# ###
# # Plot Lambda
# ###
# 
# boot6 <- bootMer(m6.trait, function(x) predict(x, newdata = df, type = 'response', allow.new.levels = T, re.form = NA), nsim = 1000, ncpus = 2, verbose = T)
# 
# PI.boot <- data.frame(
#            L = apply(boot6$t, 2, function(x) as.numeric(quantile(x, probs = .5, na.rm=TRUE))),
#            lower = apply(boot6$t, 2, function(x) as.numeric(quantile(x, probs = .025, na.rm=TRUE))),
#            upper = apply(boot6$t, 2, function(x) as.numeric(quantile(x, probs = .975, na.rm=TRUE))))
# 
# CI.6 <- merge(PI.boot, df, by = "row.names")
# save(CI.6, file = "Data/Analysis-output/CI/CI-6.Rdata")
# 
# rm(boot1, boot2, boot3, boot4, boot6, PI.boot)

#### Prep: Bootstrap Vital Rates ####

# Dataframe for parameter boots
df.bt <- expand.grid(Treat.Code = c("D", "C", "W"), Subplot = c("Grass", "No Grass"), PC.F = unique(trait.w$PC.F)) 


# ###
# # Mortality
# ###
# BS.m <- bootMer(m1.trait, function(x) predict(x, newdata = df.bt, type = 'response', allow.new.levels = T, re.form = NA), nsim = 1000, ncpus = 2, verbose = T) #32/1000 = 3.2% convergence warnings
# 
# save(BS.m, file = "Data/Analysis-output/Param-sims/BS-m.Rdata")
# 
# ###
# # Seed set
# ###
# BS.F <- bootMer(m2.trait, function(x) predict(x, newdata = df.bt, type = 'response', allow.new.levels = T, re.form = NA), nsim = 1000, ncpus = 2, verbose = T)
# 
# save(BS.F, file = "Data/Analysis-output/Param-sims/BS-F.Rdata")
#
# ###
# # Seed carryover
# ###
# 
# BS.sb <- bootMer(m3.trait, function(x) predict(x, newdata = df.bt, type = 'response', allow.new.levels = T, re.form = NA), nsim = 1000, ncpus = 2, verbose = T) #11/1000 did not converge
# 
# save(BS.sb, file = "Data/Analysis-output/Param-sims/BS-sb.Rdata")
# 
###
# Germination
###
# 
# BS.g <- bootMer(m4.trait, function(x) predict(x, newdata = df.bt, type = 'response', allow.new.levels = T, re.form = NA), nsim = 1000, ncpus = 2, verbose = T) #1/1000 failed to converge
# 
# save(BS.g, file = "Data/Analysis-output/Param-sims/BS-g.Rdata")

#### Prep: Boots seed avg ####

BS <- list(BS.m, BS.F, BS.g)
# Create dataframe of parameter means and sds
params <- c("p.mort.b", "seed.set.b", "p.germ.b", "p.mort.sd",  "seed.set.sd", "p.germ.sd")
pm <- as.data.frame(matrix(NA, 36, length(params)))
names(pm) = params
pm <- cbind(df.bt, pm)

for(i in 1:3){
    pm[, i+3] = apply(BS[[i]]$t, 2, function(x) mean(x))
    pm[, i+6] = apply(BS[[i]]$t, 2, function(x) sd(x))
}

# merge with seed survival averages
pm <- merge(pm, trait.w[,c(1,24)], by = "PC.F")
pm <- merge(pm, sb.spp, by = "Species")

# Revalue levels for creating scenarios and matching later
pm$Species <- revalue(pm$Species, c("Agoseris heterophylla" =  "AGHE", "Calycadenia pauciflora" = "CAPA", "Clarkia purpurea" = "CLPU", "Hemizonia congesta" = "HECO", "Lasthenia californica" = "LACA", "Plantago erecta" = "PLER"))
pm$Subplot <- revalue(pm$Subplot, c("No Grass" = "N", "Grass" = "G"))
pm$Scenario <- paste(pm$Treat.Code, pm$Subplot, pm$Species, sep = ".")
Scenario <- pm$Scenario
later <- pm[,c(3,4,1,2,12)]

# Create empty Lambda dataframe to fill
L.sim <- as.data.frame(Scenario)
sims = 10001
L.sim[,c(2:sims)] <- NA

# Other dfs for storage and later merger
g.sim <- L.sim
m.sim <- L.sim
s.sim <- L.sim
F.sim <- L.sim

# Calculate L from random draws of parameter distributions
for(j in 2:sims){
  for(i in 1:length(Scenario)){
    g.sim[i,j] <- rbinom(n = 1, size = 100, prob = pm[i, "p.germ.b"])/100 #germ
    s.sim[i,j] <- rbinom(n = 1, size = 100, prob = pm[i, "p.sd.surv"])/100 #seed surv
    m.sim[i,j] <- rnorm(n = 1, mean = pm[i, "p.mort.b"], sd = pm[i, "p.mort.sd"])
    F.sim[i,j] <- exp(rnorm(n = 1, mean = pm[i, "seed.set.b"], sd = pm[i, "seed.set.sd"])) - 1 #fecund
    L.sim[i,j] = s.sim[i,j]*(1 - g.sim[i,j]) + g.sim[i,j]*(1 - m.sim[i,j])*F.sim[i,j]
  }
}

# melt all dataframes
sim.list <- list(L.sim, g.sim, s.sim, m.sim, F.sim)
for(i in 1:5){ 
  sim.list[[i]] <- melt(sim.list[[i]])
}

# merge all dataframes
names(sim.list[[1]])[3] <- "L"
names(sim.list[[2]])[3] <- "g"
names(sim.list[[3]])[3] <- "s"
names(sim.list[[4]])[3] <- "m"
names(sim.list[[5]])[3] <- "Fe"

L.sim <- merge(sim.list[[1]], sim.list[[2]], by = c("Scenario", "variable"))
L.sim <- merge(L.sim, sim.list[[3]], by = c("Scenario", "variable"))
L.sim <- merge(L.sim, sim.list[[4]], by = c("Scenario", "variable"))
L.sim <- merge(L.sim, sim.list[[5]], by = c("Scenario", "variable"))
L.sim <- merge(L.sim, later, by = "Scenario")

B.L <- L.sim[,c(8:11,3:7)]

B.L.sum <- ddply(B.L, .(Species, Treat.Code, Subplot, PC.F), summarize, L = mean(L), g = mean(g), sb = mean(s), m = mean(m), Fe = mean(Fe))

save(B.L, file = "Data/Analysis-output/20190206-lambda-sim.Rdata")
save(B.L.sum, file = "Data/Analysis-output/20190206-lambda-sim-sum.Rdata")

rm(g.sim, s.sim, m.sim, F.sim, params, pm, i, j, later, Scenario, sim.list, sims, BS, L.sim)


#### M7: Boot Lambda ####

###
# Lambda with seed survival averages
###
B.L.sum$Treat.Code <- factor(B.L.sum$Treat.Code, levels = c("C", "D", "W"))
B.L.sum$Subplot <- factor(B.L.sum$Subplot, levels = c("N", "G"))
B.L.sum$Subplot <- revalue(B.L.sum$Subplot, c("N" = "No Grass", "G" = "Grass"))
hist(B.L.sum$L)
hist(log(B.L.sum$L))
hist(log(B.L.sum$L + 1))
B.L.sum$log.L <- log(B.L.sum$L + 1)
m7.trait <- lm(log.L ~ Treat.Code * Subplot * PC.F, B.L.sum)
plot(fitted(m7.trait), resid(m7.trait))
qqnorm(resid(m7.trait))
qqline(resid(m7.trait), col = 2, lwd = 2, lty = 2) 
summary(m7.trait)
hist(resid(m7.trait))

# Get confidence intervals for lambda

B.L.CI <- predict(m7.trait, newdata = df, interval = "confidence")
B.L.CI <- merge(B.L.CI, df, by = "row.names")

#### Fig 1: RGR v WUE ####
trait.w$Species.short <- revalue(trait.w$Species, c("Agoseris heterophylla" =  "AGHE", "Calycadenia pauciflora" = "CAPA", "Clarkia purpurea" = "CLPU", "Hemizonia congesta" = "HECO", "Lasthenia californica" = "LACA", "Plantago erecta" = "PLER"))


plot.trait <- ggplot(trait.w, aes(x = D13C.F, y = RGR.la.F)) +
  theme_classic() +
  geom_smooth(method = "lm", se = F, col = "black", size = 1) +
  geom_text(aes(label = Species.short), hjust = .5, vjust = .5) +
  labs(y = "Relative Growth Rate", x = "Water Use Efficiency") +
  theme(axis.text = element_text(size = 10), 
        plot.title = element_text(size = 30, face="bold", vjust = 2),
        axis.title = element_text(size = 13), 
        strip.text = element_text(size = 15),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_x_reverse(lim = c(25, 20)) +
  labs(x = "Carbon isotope discrimination", y = "Relative Growth Rate")
  #labs(x = "Carbon isotope discrimination (∆, \u2030)", y = expression(paste("Relative Growth Rate (", cm^{2}, ")"%.%"(", cm^{2}, ")"^{-1}%.%day^{-1})))

ggsave(plot.trait, filename = "Figures/trait.tiff", width = 82, height = 75, units = "mm", dpi = 600)

#### Fig 2: Mortality v PC ####
dem.sum <- summarySE(dem, measurevar = "p.mort", groupvars = c("Treat.Code", "Subplot", "Species", "PC.F"), na.rm = T)

dem.sum$Treat.Code <- factor(dem.sum$Treat.Code, levels = c("D", "C", "W"))

# control and shelter
plot.m1.trait.d <- ggplot(dem.sum[dem.sum$Treat.Code != "W",], aes(y = p.mort, x = PC.F, col = Treat.Code, group = Treat.Code)) +
  geom_ribbon(data = CI.1[CI.1$Treat.Code != "W",], aes(ymin = p.mort.lo, ymax = p.mort.hi), alpha = 0.1, linetype = 0) +
  geom_point() +
  geom_errorbar(aes(ymin = p.mort - se, ymax = p.mort + se)) +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 10), 
        plot.title = element_text(size=30, face="bold", vjust = 2),
        axis.title = element_text(size = 15), 
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = c(0.12, 0.82),
        legend.key.size = unit(1.5, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_color_manual(values = c("red1", "orange1"), labels = c("Drought", "Control")) +
  scale_y_continuous(labels = scales::percent, limits = c(0,.85)) +
  labs(y = "Mortality", x = "Drought Tolerance") +
  facet_wrap(~ Subplot) +
  geom_line(data = CI.1[CI.1$Treat.Code != "W",], aes(y = p.mort))

ggsave(plot.m1.trait.d, filename = "Figures/mortality-trait-d.tiff", width = 173, height = 90, units = "mm", dpi = 600)

# control and water  
plot.m1.trait.w <- ggplot(dem.sum[dem.sum$Treat.Code != "D",], aes(y = p.mort, x = PC.F, col = Treat.Code, group = Treat.Code)) +
  geom_ribbon(data = CI.1[CI.1$Treat.Code != "D",], aes(ymin = p.mort.lo, ymax = p.mort.hi), alpha = 0.1, linetype = 0) +
  geom_point() +
  geom_errorbar(aes(ymin = p.mort - se, ymax = p.mort + se)) +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 10), 
        plot.title = element_text(size=30, face="bold", vjust = 2),
        axis.title = element_text(size = 15), 
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = c(0.12, 0.82),
        legend.key.size = unit(1.5, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_color_manual(values = c("orange1", "mediumblue"), labels = c("Control", "Watered")) +
  scale_y_continuous(labels = scales::percent, limits = c(0,.8)) +
  labs(y = "Mortality", x = "Drought Tolerance") +
  facet_wrap(~ Subplot) +
  geom_line(data = CI.1[CI.1$Treat.Code != "D",], aes(y = p.mort))

ggsave(plot.m1.trait.w, filename = "Figures/mortality-trait-w.tiff", width = 173, height = 90, units = "mm", dpi = 600)

#### Fig 3: Fecund v PC ####

flo.seed.sum <- summarySE(flo.seed, measurevar = "n.seed.ind", groupvars = c("Treat.Code", "Subplot", "Species", "PC.F"), na.rm = T)
flo.seed.sum$Treat.Code <- factor(flo.seed.sum$Treat.Code, levels = c("D", "C", "W"))

## drought
plot.m2.trait.d <- ggplot(flo.seed.sum[flo.seed.sum$Treat.Code != "W",], aes(y = n.seed.ind, x = PC.F, col = Treat.Code, group = Treat.Code)) +
  geom_ribbon(data = CI.2[CI.2$Treat.Code != "W",], aes(ymin = exp(seed.lo) - 1, ymax = exp(seed.hi) - 1), alpha = 0.1, linetype = 0) +
  geom_point() +
  geom_errorbar(aes(ymin = n.seed.ind - se, ymax = n.seed.ind + se)) +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 15), 
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 11),
        legend.position = c(.85, .78),
        legend.key.size = unit(1.4, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_color_manual(values = c("red1", "orange1"), labels = c("Drought", "Control")) +
  labs(y = "Seeds per Individual", x = "Drought Tolerance") +
  facet_wrap(~ Subplot) +
  geom_line(data = CI.2[CI.2$Treat.Code != "W",], aes(y = exp(n.seed.ind) - 1))

ggsave(plot.m2.trait.d, filename = "Figures/seedset-trait-d.tiff", width = 173, height = 90, units = "mm", dpi = 600)

## Watering
plot.m2.trait.w <- ggplot(flo.seed.sum[flo.seed.sum$Treat.Code != "D",], aes(y = n.seed.ind, x = PC.F, col = Treat.Code, group = Treat.Code)) +
  geom_ribbon(data = CI.2[CI.2$Treat.Code != "D",], aes(ymin = exp(seed.lo) - 1, ymax = exp(seed.hi) - 1), alpha = 0.1, linetype = 0) +
  geom_point() +
  geom_errorbar(aes(ymin = n.seed.ind - se, ymax = n.seed.ind + se)) +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 15), 
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 11),
        legend.position = c(.85, .78),
        legend.key.size = unit(1.4, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_color_manual(values = c("orange1", "mediumblue"), labels = c("Control", "Watered")) +
  labs(y = "Seeds per Individual", x = "Drought Tolerance") +
  facet_wrap(~ Subplot) +
  geom_line(data = CI.2[CI.2$Treat.Code != "D",], aes(y = exp(n.seed.ind) - 1))

ggsave(plot.m2.trait.w, filename = "Figures/seedset-trait-w.tiff", width = 173, height = 90, units = "mm", dpi = 600)

#### Fig: Germ v PC ####
germ.sum <- summarySE(filter(dem, !(Species == "Lasthenia californica" & Year == 2016), !(Species == "Plantago erecta" & Year == 2016)), measurevar = "p.germ", groupvars = c("Species", "PC.F"), na.rm = T)

ggplot(germ.sum, aes(y = p.germ, x = PC.F)) +
  geom_point() +
  geom_errorbar(aes(ymin = p.germ - se, ymax = p.germ + se)) +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 20), 
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.position = "right",
        legend.key.size = unit(2, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  labs(y = "Germination", x = "Drought Tolerance") +
  geom_line(data = CI.4, aes(y = germ)) +
  geom_ribbon(data = CI.4, aes(y = NULL, ymin = germ.lo, ymax = germ.hi), alpha = 0.2, linetype = 0) +
  ylim(0,1) +
  facet_wrap(~Year)


#### Fig: Seed Surv v Species ####
sb.sum.trt <- summarySE(sb, groupvars = c("Species", "PC.F", "Treat.Code"), measurevar = "p.surv")
sb.sum.trt$Treat.Code <- factor(sb.sum.trt$Treat.Code, levels = c("D", "C", "W"))

ggplot(sb.sum.trt, aes(y = p.surv, x = Treat.Code)) +
  geom_bar(stat = "identity", position = position_dodge(width=0.9), col = "black") +
  geom_errorbar(aes(ymin = p.surv - se, ymax = p.surv + se), width = 0.2, position = position_dodge(width=0.9)) +
   theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 17),
        legend.text = element_text(size = 12),
        legend.position = "right",
        legend.key.size = unit(1.4, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  labs(y = "Forb Proportion Seed Survival") +
  facet_wrap(~Species)


#### Fig 4: L v PC ####
plot.lambda <- ggplot(B.L.sum, aes(y = L, x = PC.F, col = Treat.Code, group = Treat.Code)) +
  geom_point() +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 15), 
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 11),
        legend.position = c(.85, .78),
        legend.key.size = unit(1.4, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_color_manual(values = c("orange1", "red1", "mediumblue"), labels = c("Control", "Drought", "Watered")) +
  labs(y = "Population Growth Rate", x = "Drought Tolerance") +
  facet_wrap(~ Subplot) +
  geom_line(data = B.L.CI, aes(y = exp(fit)-1))

ggsave(plot.lambda , filename = "Figures/lambda.tiff", width = 173, height = 90, units = "mm", dpi = 600)
