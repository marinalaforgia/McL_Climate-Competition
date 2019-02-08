# Script for Seed Carryover Data Analysis
rm(list=ls())


#### Load Libraries ####
library(plyr)
library(lme4)
library(Rmisc)
library(ggplot2)
library(pscl)

### Read in Data ####
sb16 <- read.csv("Data/Seed bag/Seed_bags_2016.csv")
sb17 <- read.csv("Data/Seed bag/Seed_bags_2017.csv")
treat <- read.csv("Data/Setup/Marina-Treatment-30.csv")
added <- read.csv("Data/Setup/Exp-Seed-Species.csv")
viab <- read.csv("Data/Viability/Viability-Overall.csv") # seed viability data
trait.w <- read.csv("Data/Post-Processing/final-traits-w.csv")
trait.w <- trait.w[,-7]

#PCA.t <- prcomp(trait.w[,c(15,17)], scale = T) # greenhouse trait
PCA.t <- prcomp(trait.w[,c(7,8)], scale = T)
biplot(PCA.t)
summary(PCA.t)
trait.w$PC1 <- PCA.t$x[,1] 

colnames(sb17)[3] <- "n.viable"

sb16$Year <- 2016
sb17$Year <- 2017

# 2016 viability
added16 <- filter(added, Species %in% unique(sb17$Species), exp.year == 2016)[,c(2,4)] # avg number of seed added per subplot in 2016
viab16 <- viab[c(1:5), c(1,9)] # viability data
added16$Species <- revalue(added16$Species, c("PLER" = "Plantago erecta", "AGHE" = "Agoseris heterophylla", "LACA" = "Lasthenia californica", "CLPU" = "Clarkia purpurea", "HECO" = "Hemizonia congesta", "CAPA" = "Calycadenia pauciflora", "AVFA" = "Avena fatua", "BRHO" = "Bromus hordeaceus", "TACA" = "Taeniatherum caput-medusae", "LOMU" = "Lolium multiflorum"))
added16 <- merge(added16, viab16, by = "Species")
added16$viable <- added16$avg.num.per.sub*added16$p.viable
added16$viable <- as.integer(added16$viable)
added16$Year <- 2016


sb16$Species <- revalue(sb16$Species, c("AGHE" = "Agoseris heterophylla", "CLPU" = "Clarkia purpurea", "LACA" = "Lasthenia californica", "PLER" = "Plantago erecta", "HECO" = "Hemizonia congesta"))
sb16 <- merge(sb16, added16[,c(1,4,5)], by = c("Species", "Year"))

# 2017 viability
added17 <- filter(added, Species %in% unique(sb17$Species), exp.year == 2017)[,c(2,4)] # avg number of seed added per subplot in 2017
added17$Species <- revalue(added17$Species, c("AGHE" = "Agoseris heterophylla", "CLPU" = "Clarkia purpurea", "LACA" = "Lasthenia californica", "PLER" = "Plantago erecta", "HECO" = "Hemizonia congesta", "CAPA" = "Calycadenia pauciflora"))
viab17 <- viab[c(1:3,7,8,12), c(1,9)] 
added17 <- merge(added17, viab17, by = "Species")
added17$viable <- added17$avg.num.per.sub*added17$p.viable
added17$viable <- as.integer(added17$viable)
added17$Year <- 2017

sb17$Species <- revalue(sb17$Species, c("AGHE" = "Agoseris heterophylla", "CLPU" = "Clarkia purpurea", "LACA" = "Lasthenia californica", "PLER" = "Plantago erecta", "HECO" = "Hemizonia congesta", "CAPA" = "Calycadenia pauciflora"))
sb17 <- merge(sb17, added17[,c(1,4,5)], by = c("Species", "Year"))

sb <- rbind(sb16, sb17)

sb <- merge(sb, unique(treat[,c(1,3)]), by = "Plot", all.y = F)
colnames(sb)[2] <- "Species_Name"

colnames(sb)[4] <- "Count"

sb$avg.num.per.sub <- ifelse(sb$Count > sb$viable, sb$Count, sb$viable) # cap at high
sb$p.surv <- sb$Count/sb$avg.num.per.sub

sb <- merge(sb, trait.w, by.x = "Species_Name", by.y = "Species")
write.table(sb, "Data/Cleaned Data for Analysis/seed-carryover-plot.csv", sep = ",", row.names = F)
#### Model: traits on seed carryover ####
sb$Year <- as.factor(sb$Year)

sb.sum <- ddply(sb, .(Species_Name, PC1), summarize, Count = sum(Count), avg.num.per.sub = sum(avg.num.per.sub))
sb.sum$p.surv <- sb.sum$Count/sb.sum$avg.num.per.sub

m.sb.t <- glm(cbind(Count, avg.num.per.sub - Count) ~ PC1, family = binomial, data = sb)
plot(fitted(m.sb.t), resid(m.sb.t))
qqnorm(resid(m.sb.t)) 
qqline(resid(m.sb.t), col = 2,lwd=2,lty=2) 
summary(m.sb.t) 

m.sb.t <- lm(p.surv ~ PC1, data = sb.sum)
plot(fitted(m.sb.t), resid(m.sb.t))
qqnorm(resid(m.sb.t)) 
qqline(resid(m.sb.t), col = 2,lwd=2,lty=2) 
summary(m.sb.t) # using percentage gives more what i would expect

m.sb.t <- glmer(cbind(Count, avg.num.per.sub-Count) ~ PC1 + (1|Year) + (1|Plot), family = binomial, data = sb)
plot(fitted(m.sb.t), resid(m.sb.t))
qqnorm(resid(m.sb.t)) 
qqline(resid(m.sb.t), col = 2,lwd=2,lty=2) 
summary(m.sb.t)

ggplot(sb.sum, aes(x = PC1, y = Count/avg.num.per.sub)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  theme_classic() + 
  geom_text(aes(label = Species_Name), hjust = .5, vjust = .5) +
  labs(x = "Drought tolerance (high to low)", y = "Seed Carryover")# noooo

sb$p.surv <- sb$Count/sb$avg.num.per.sub
sb.summary <- summarySE(sb, measurevar = "p.surv", groupvars = c("Species_Name", "PC1"))

ggplot(sb.summary, aes(x = PC1, y = p.surv)) +
  geom_point() +
  geom_errorbar(aes(ymin = p.surv - se, ymax = p.surv + se), width = 0.02) +
  geom_smooth(method = "lm", se = F) +
  theme_classic()

#### Model: Year effect on Forbs ####
sb$Year <- as.factor(sb$Year)
m0.nb <- glmer.nb(Count ~ Year + (1|Species_Name), data = sb[sb$Group == "Forb",])
m0.p <- glmer(Count ~ Year + (1|Species_Name), family = poisson, data = sb[sb$Group == "Forb",])
lrtest(m0.nb, m0.p) # neg binom better
plot(fitted(m0.nb), resid(m0.nb))
qqnorm(resid(m0.nb)) 
qqline(resid(m0.nb), col = 2,lwd=2,lty=2) 
summary(m0.nb) # carryover is lower in general in 2017, don't lump together 

# Aggregated binomial model
m0.ab <- glmer(cbind(Count, avg.num.per.sub - Count) ~ Year + (1|Species_Name), family = binomial, data = sb[sb$Group == "Forb",])
plot(fitted(m0.ab), resid(m0.ab))
qqnorm(resid(m0.ab)) 
qqline(resid(m0.ab), col = 2, lwd = 2, lty = 2) 
summary(m0.ab) 

#### Model: TACA ####
m.TACA.d <- glm(Count ~ Treat.Code, family = poisson, data = sb[sb$Species_Name == "Taeniatherum caput-medusae" & sb$Treat.Code != "W",])
dispersiontest(m.TACA.d, trafo = 1) # not sig greater than 0, stick to poisson
plot(fitted(m.TACA.d), resid(m.TACA.d))
qqnorm(resid(m.TACA.d)) 
qqline(resid(m.TACA.d), col = 2,lwd=2,lty=2) 
summary(m.TACA.d) # no effect of drought on carryover

m.TACA.w <- glm(Count ~ Treat.Code, family = poisson, data = sb[sb$Species_Name == "Taeniatherum caput-medusae" & sb$Treat.Code != "D",])
dispersiontest(m.TACA.w, trafo = 1) # not sig greater than 0, stick to poisson
plot(fitted(m.TACA.w), resid(m.TACA.w))
qqnorm(resid(m.TACA.w)) 
qqline(resid(m.TACA.w), col = 2,lwd=2,lty=2) 
summary(m.TACA.w) # no effect of watering on carryover

# sum TACA carryover across treatments

#### Model: AVFA ####
m.AVFA.d <- glm(Count ~ Treat.Code, family = poisson, data = sb[sb$Species_Name == "Avena fatua" & sb$Treat.Code != "W",])
dispersiontest(m.AVFA.d, trafo = 1) #poisson good
plot(fitted(m.AVFA.d), resid(m.AVFA.d))
qqnorm(resid(m.AVFA.d)) 
qqline(resid(m.AVFA.d), col = 2,lwd=2,lty=2) 
summary(m.AVFA.d) # no effect of drought on carryover

m.AVFA.w <- glm(Count ~ Treat.Code, family = poisson, data = sb[sb$Species_Name == "Avena fatua" & sb$Treat.Code != "D",])
dispersiontest(m.AVFA.w, trafo = 1) #poisson good
plot(fitted(m.AVFA.w), resid(m.AVFA.w))
qqnorm(resid(m.AVFA.w)) 
qqline(resid(m.AVFA.w), col = 2,lwd=2,lty=2) 
summary(m.AVFA.w) # no effect of watering on carryover

# sum AVFA carryover across treatments

#### Model: BRHO ####
m.BRHO.d <- glm(Count ~ Treat.Code, family = poisson, data = sb[sb$Species_Name == "Bromus hordeaceus" & sb$Treat.Code != "W",])
plot(fitted(m.BRHO.d), resid(m.BRHO.d))
qqnorm(resid(m.BRHO.d)) 
qqline(resid(m.BRHO.d), col = 2,lwd=2,lty=2) 
summary(m.BRHO.d) # no effect of drought on carryover

m.BRHO.w <- glm(Count ~ Treat.Code, family = poisson, data = sb[sb$Species_Name == "Bromus hordeaceus" & sb$Treat.Code != "D",])
plot(fitted(m.BRHO.w), resid(m.BRHO.w))
qqnorm(resid(m.BRHO.w)) 
qqline(resid(m.BRHO.w), col = 2,lwd=2,lty=2) 
summary(m.BRHO.w) # no effect of watering on carryover

# sum BRHO carryover across treatments

#### Model: LOMU ####
m.LOMU.d <- glm(Count ~ Treat.Code, family = poisson, data = sb[sb$Species_Name == "Lolium multiflorum" & sb$Treat.Code != "W",])
dispersiontest(m.LOMU.d, trafo = 1) # overdispersion present, use neg. bin model
m.LOMU.d <- glm.nb(Count ~ Treat.Code, data = sb[sb$Species_Name == "Lolium multiflorum" & sb$Treat.Code != "W",])
plot(fitted(m.LOMU.d), resid(m.LOMU.d))
qqnorm(resid(m.LOMU.d)) 
qqline(resid(m.LOMU.d), col = 2,lwd=2,lty=2) 
summary(m.LOMU.d) # no effect of drought on carryover

m.LOMU.w.nb <- glm.nb(Count ~ Treat.Code, data = sb[sb$Species_Name == "Lolium multiflorum" & sb$Treat.Code != "D",])
m.LOMU.w <- glm(Count ~ Treat.Code, family = poisson, data = sb[sb$Species_Name == "Lolium multiflorum" & sb$Treat.Code != "D",])
dispersiontest(m.LOMU.w, trafo = 1) # overdispersion present, use neg. bin model

plot(fitted(m.LOMU.w.nb), resid(m.LOMU.w.nb))
qqnorm(resid(m.LOMU.w.nb)) 
qqline(resid(m.LOMU.w.nb), col = 2,lwd=2,lty=2) 
summary(m.LOMU.w.nb) # no effect of watering on carryover

# sum LOMU carryover across treatments

#### Model: Differences in grass carryover ####
# If we let intercept vary by species, that parameter should absorb the differences in seed added per species, no?
m.grass.w.nb <- glmer.nb(Count ~ Treat.Code + (1|Species_Name), data = sb[sb$Group == "Grass" & sb$Treat.Code != "D",])
m.grass.w <- glmer(Count ~ Treat.Code + (1|Species_Name), family = poisson, data = sb[sb$Group == "Grass" & sb$Treat.Code != "D",])
lrtest(m.grass.w.nb, m.grass.w) # nb model better

plot(fitted(m.grass.w.nb), resid(m.grass.w.nb))
qqnorm(resid(m.grass.w.nb)) 
qqline(resid(m.grass.w.nb), col = 2,lwd=2,lty=2) 
summary(m.grass.w.nb) # no effect of watering on carryover

# aggregated binomial model better?
# Watering
m.grass.w.ab <- glmer(cbind(Count, avg.num.per.sub - Count) ~ Treat.Code + (1|Species_Name), family = binomial, data = sb[sb$Group == "Grass" & sb$Treat.Code != "D",])
plot(fitted(m.grass.w.ab), resid(m.grass.w.ab))
qqnorm(resid(m.grass.w.ab)) 
qqline(resid(m.grass.w.ab), col = 2,lwd=2,lty=2) 
summary(m.grass.w.ab) 

# Drought
m.grass.d.ab <- glmer(cbind(Count, avg.num.per.sub - Count) ~ Treat.Code + (1|Species_Name), family = binomial, data = sb[sb$Group == "Grass" & sb$Treat.Code != "W",])
plot(fitted(m.grass.d.ab), resid(m.grass.d.ab))
qqnorm(resid(m.grass.d.ab)) 
qqline(resid(m.grass.d.ab), col = 2,lwd=2,lty=2) 
summary(m.grass.d.ab) # no species level effects but drought marginally increased carryover and watering marginally decreased carryover

# marginal effects overall, lump together
sb.sum.grass <- summarySE(sb[sb$Group == "Grass",], measurevar = "p.surv", groupvars = "Species_Name")

#### Model 1: Drought effect on seed carryover ####
m1 <- glmer(cbind(Count, avg.num.per.sub - Count) ~ Treat.Code + (1|Year) + (1|Species_Name), family = binomial, data = sb)
plot(fitted(m1), resid(m1))
qqnorm(resid(m1)) 
qqline(resid(m1), col = 2, lwd = 2, lty = 2) 
summary(m1) # no effect of drought or of watering

#### Model 2: Watering effect on seed carryover ####
m2 <- glmer(cbind(Count, avg.num.per.sub - Count) ~ Treat.Code + (1|Year) + (1 + Treat.Code|Species_Name), family = binomial, data = sb[sb$Treat.Code != "D" & sb$Group == "Forb",])
plot(fitted(m2), resid(m2))
qqnorm(resid(m2)) 
qqline(resid(m2), col = 2,lwd=2,lty=2) 
summary(m2) # no effect of water


ranef.m1 <- data.frame(Species_Name = row.names(ranef(m1)[[1]]), Diff.avg = ranef(m1)[[1]][,1])

m1.trait <- merge(sb.summary, ranef.m1, by = "Species_Name")

ggplot(m1.trait, aes(x = PC1, y = Diff.avg)) +
  geom_smooth(method = "lm", se = F) +
  theme_classic() +
  geom_text(aes(label = Species_Name), hjust = .5, vjust = .5)
  
#### Model: HECO ####
hist(sb[sb$Species_Name == "Hemizonia congesta",]$Count)

m4.h.d <- glmer.nb(Count ~ Treat.Code + (1|Year), data = sb[sb$Species == "Hemizonia congesta",])
#m4.h.d2 <- glmer(Count ~ Treat.Code + (1|Year), family = poisson, data = sb[sb$Species_Name == "Hemizonia congesta" & sb$Treat.Code != "W",])
#anova(m4.h.d2, m4.h.d) # neg binomial better
plot(fitted(m4.h.d), resid(m4.h.d))
qqnorm(resid(m4.h.d)) 
qqline(resid(m4.h.d), col = 2, lwd = 2, lty = 2) 
summary(m4.h.d) 

# Watering
#m4.h.w.2 <- glmer(Count ~ Treat.Code + (1|Year), family = "poisson", data = sb[sb$Species_Name == "Hemizonia congesta" & sb$Treat.Code != "D",])
m4.h.w <- glmer.nb(Count ~ Treat.Code + (1|Year), data = sb[sb$Species_Name == "Hemizonia congesta" & sb$Treat.Code != "D",])
#anova(m4.h.w.2, m4.h.w) # negbinom better
plot(fitted(m4.h.w), resid(m4.h.w))
qqnorm(resid(m4.h.w)) 
qqline(resid(m4.h.w), col = 2,lwd=2,lty=2) 
summary(m4.h.w)  

#### Model: AGHE ####

hist(sb[sb$Species_Name == "Agoseris heterophylla",]$Count) 

# Drought
m4.a.d <- glmer.nb(Count ~ Treat.Code + (1|Year), data = sb[sb$Species_Name == "Agoseris heterophylla" & sb$Treat.Code != "W",])
m4.a.d.2 <- glmer(Count ~ Treat.Code + (1|Year), family = "poisson", data = sb[sb$Species_Name == "Agoseris heterophylla" & sb$Treat.Code != "W",])
anova(m4.a.d.2, m4.a.d) #neg binom better
plot(fitted(m4.a.d), resid(m4.a.d))
qqnorm(resid(m4.a.d)) 
qqline(resid(m4.a.d), col = 2, lwd = 2, lty = 2) 
summary(m4.a.d) 

# Watering
m4.a.w.2 <- glmer(Count ~ Treat.Code + (1|Year), family = "poisson", data = sb[sb$Species_Name == "Agoseris heterophylla" & sb$Treat.Code != "D",])
m4.a.w <- glmer.nb(Count ~ Treat.Code + (1|Year), data = sb[sb$Species_Name == "Agoseris heterophylla" & sb$Treat.Code != "D",])
#anova(m4.a.w.2, m4.a.w) #neg binom better
plot(fitted(m4.a.w), resid(m4.a.w))
qqnorm(resid(m4.a.w)) 
qqline(resid(m4.a.w), col = 2,lwd=2,lty=2) 
summary(m4.a.w)  


#### Model: PLER ####

hist(sb[sb$Species_Name == "Plantago erecta",]$Count)


# Drought
m4.p.d <- glmer.nb(Count ~ Treat.Code + (1|Year), data = sb[sb$Species_Name == "Plantago erecta" & sb$Treat.Code != "W",])
#m4.p.d.2 <- glmer(Count ~ Treat.Code + (1|Year), family = "poisson", data = sb[sb$Species_Name == "Plantago erecta" & sb$Treat.Code != "W",])
#anova(m4.p.d.2, m4.p.d) #neg binom better
plot(fitted(m4.p.d), resid(m4.p.d))
qqnorm(resid(m4.p.d)) 
qqline(resid(m4.p.d), col = 2, lwd = 2, lty = 2) 
summary(m4.p.d) 

# Watering
#m4.p.w.2 <- glmer(Count ~ Treat.Code + (1|Year), family = "poisson", data = sb[sb$Species_Name == "Plantago erecta" & sb$Treat.Code != "D",])
m4.p.w <- glmer.nb(Count ~ Treat.Code + (1|Year), data = sb[sb$Species_Name == "Plantago erecta" & sb$Treat.Code != "D",])
#anova(m4.p.w.2, m4.p.w) #neg binom better
plot(fitted(m4.p.w), resid(m4.p.w))
qqnorm(resid(m4.p.w)) 
qqline(resid(m4.p.w), col = 2,lwd=2,lty=2) 
summary(m4.p.w)  

#### Model: CLPU ####
hist(sb[sb$Species_Name == "Clarkia purpurea",]$Count)

# Drought
m4.cl.d <- glmer.nb(Count ~ Treat.Code + (1|Year), data = sb[sb$Species_Name == "Clarkia purpurea" & sb$Treat.Code != "W",])
#m4.cl.d.2 <- glmer(Count ~ Treat.Code + (1|Year), family = "poisson", data = sb[sb$Species_Name == "Clarkia purpurea" & sb$Treat.Code != "W",])
#anova(m4.cl.d.2, m4.cl.d) #neg binom better
plot(fitted(m4.cl.d), resid(m4.cl.d))
qqnorm(resid(m4.cl.d)) 
qqline(resid(m4.cl.d), col = 2, lwd = 2, lty = 2) 
summary(m4.cl.d) 

# Watering
#m4.cl.w.2 <- glmer(Count ~ Treat.Code + (1|Year), family = "poisson", data = sb[sb$Species_Name == "Clarkia purpurea" & sb$Treat.Code != "D",])
m4.cl.w <- glmer.nb(Count ~ Treat.Code + (1|Year), data = sb[sb$Species_Name == "Clarkia purpurea" & sb$Treat.Code != "D",])
#anova(m4.cl.w.2, m4.cl.w) #neg binom better
plot(fitted(m4.cl.w), resid(m4.cl.w))
qqnorm(resid(m4.cl.w)) 
qqline(resid(m4.cl.w), col = 2,lwd=2,lty=2) 
summary(m4.cl.w)  


#### Model: LACA ####
hist(sb[sb$Species_Name == "Lasthenia californica",]$Count)

# Drought
m4.l.d <- glmer.nb(Count ~ Treat.Code + (1|Year), data = sb[sb$Species_Name == "Lasthenia californica" & sb$Treat.Code != "W",])
# m4.l.d.2 <- glmer(Count ~ Treat.Code + (1|Year), family = "poisson", data = sb[sb$Species_Name == "Lasthenia californica" & sb$Treat.Code != "W",])
# anova(m4.l.d.2, m4.l.d) #neg binom better
plot(fitted(m4.l.d), resid(m4.l.d))
qqnorm(resid(m4.l.d)) 
qqline(resid(m4.l.d), col = 2, lwd = 2, lty = 2) 
summary(m4.l.d) 

# Watering
#m4.l.w.2 <- glmer(Count ~ Treat.Code + (1|Year), family = "poisson", data = sb[sb$Species_Name == "Lasthenia californica" & sb$Treat.Code != "D",])
m4.l.w <- glmer.nb(Count ~ Treat.Code + (1|Year), data = sb[sb$Species_Name == "Lasthenia californica" & sb$Treat.Code != "D",])
#anova(m4.l.w.2, m4.l.w) #neg binom better
plot(fitted(m4.l.w), resid(m4.l.w))
qqnorm(resid(m4.l.w)) 
qqline(resid(m4.l.w), col = 2,lwd=2,lty=2) 
summary(m4.l.w) 

#### Model: CAPA ####
# Drought
m4.ca.d <- glm.nb(Count ~ Treat.Code, data = sb[sb$Species_Name == "Calycadenia pauciflora" & sb$Treat.Code != "W",])
#m4.ca.d.2 <- glm(Count ~ Treat.Code, family = "poisson", data = sb[sb$Species_Name == "Calycadenia pauciflora" & sb$Treat.Code != "W",])
#dispersiontest(m4.ca.d.2, trafo = 1) #neg binom better
plot(fitted(m4.ca.d), resid(m4.ca.d))
qqnorm(resid(m4.ca.d)) 
qqline(resid(m4.ca.d), col = 2, lwd = 2, lty = 2) 
summary(m4.ca.d) 

# Watering
#m4.ca.w.2 <- glm(Count ~ Treat.Code, family = "poisson", data = sb[sb$Species_Name == "Calycadenia pauciflora" & sb$Treat.Code != "D",])
m4.ca.w <- glm.nb(Count ~ Treat.Code, data = sb[sb$Species_Name == "Calycadenia pauciflora" & sb$Treat.Code != "D",])
#dispersiontest(m4.ca.w.2, trafo = 1) #negbinom better
plot(fitted(m4.ca.w), resid(m4.ca.w))
qqnorm(resid(m4.ca.w)) 
qqline(resid(m4.ca.w), col = 2,lwd=2,lty=2) 
summary(m4.ca.w) 

#### Graphs and summaries ####

sb.sum <- summarySE(sb, measurevar = "p.surv", groupvars = c("Species_Name"))
write.table(sb.sum, "Seed-survival-soil.csv", row.names = F, sep = ",")

sb.sum.yr <- summarySE(sb, measurevar = "p.surv", groupvars = c("Species_Name", "Year"))

#### extra code ####
# sb.sum$Treat.Code <- revalue(sb.sum$Treat.Code, c("C" = "Control", "D" = "Drought", "W" = "Watered"))
# 
# sb.sum$Species_Name <- factor(sb.sum$Species_Name, levels = c("Lasthenia californica", "Clarkia purpurea", "Agoseris heterophylla","Hemizonia congesta","Plantago erecta"))
# 
# # Bar plot, species by treatment
# graph1a <- ggplot(sb.sum, aes(x = Treat.Code, y = Count, fill = Treat.Code)) +
#   geom_bar(position = "dodge", stat = "identity", color = "black") +
#   geom_errorbar(aes(ymin = Count - se, ymax = Count + se), width = 0.5, size = 0.9, position = position_dodge(width = 0.9)) +
#   facet_wrap(~Species_Name) +
#   theme_bw() +
#   theme(legend.title = element_blank(), 
#         axis.text=element_text(size = 20), 
#         plot.title = element_text(size = 30, face = "bold", margin = margin(0, 0, 20, 0)), 
#         axis.title = element_text(size = 30), 
#         strip.text = element_text(size = 25), 
#         axis.title.y = element_text(margin = margin(0, 20, 0, 0)),
#         axis.title.x = element_blank(),
#         legend.position="none") +
#   ylim(0,100) + 
#   labs(y = "Number of Viable Seeds")
#   
# ggsave(graph1a, filename = "barplot-sppxtrt.png", path = plots, width = 12, height = 8, dpi = 300)
# 
# graph1b <- ggplot(sb.sum, aes(x = Treat.Code, y = Count, fill = Treat.Code)) +
#   geom_bar(position = "dodge", stat = "identity", color = "black") +
#   geom_errorbar(aes(ymin = Count - se, ymax = Count + se), width = 0.5, size = 0.9, position = position_dodge(width = 0.9)) +
#   facet_wrap(~Species_Name, nrow = 1) +
#   theme_bw() +
#   theme(legend.title = element_blank(), 
#         axis.text=element_text(size = 15), 
#         plot.title = element_text(size = 30, face = "bold", margin = margin(0, 0, 20, 0)), 
#         axis.title = element_text(size = 18), 
#         strip.text = element_text(size = 15), 
#         axis.title.y = element_text(margin = margin(0, 20, 0, 0)),
#         axis.title.x = element_blank(),
#         legend.position="none") +
#   ylim(0,100) + 
#   labs(y = "Number of Viable Seeds")
# 
# ggsave(graph1b, filename = "barplot-sppxtrt.png", path = plots, width = 15, height = 3, dpi = 300)
# 
# # Hypothetical Graph
# st.c <- data.frame(Count = 30, Treat = "Control", Type = "Stress Tolerator")
# sa.c <- data.frame(Count = 60, Treat = "Control", Type = "Stress Avoider")
# st.d <- data.frame(Count = 30, Treat = "Drought", Type = "Stress Tolerator")
# sa.d <- data.frame(Count = 70, Treat = "Drought", Type = "Stress Avoider")
# st.w <- data.frame(Count = 20, Treat = "Watered", Type = "Stress Tolerator")
# sa.w <- data.frame(Count = 50, Treat = "Watered", Type = "Stress Avoider")
# hyp <- rbind(st.c, sa.c, st.d, st.w, sa.d, sa.w)
# 
# graph2 <- ggplot(hyp, aes(x = Treat, y = Count, fill = Type)) +
#   geom_bar(position = "dodge", stat = "identity", color = "black") +
#   theme_bw() +
#   theme(legend.title = element_blank(), 
#         axis.text=element_text(size = 20), 
#         plot.title = element_text(size = 30, face = "bold", margin = margin(0, 0, 20, 0)), 
#         axis.title = element_text(size = 30), 
#         strip.text = element_text(size = 25), 
#         axis.title.y = element_text(margin = margin(0, 20, 0, 0)),
#         axis.title.x = element_blank(),
#         legend.text = element_text(size = 20),
#         legend.key.size = unit(3, 'lines')) +
#   ylim(0,100) + 
#   labs(y = "Number of Viable Seeds", title = "Hypothetical Results")
# 
# ggsave(graph2, filename = "hypbarplot-typxtrt.png", path = plots, width = 12, height = 8, dpi = 300)
# 
# # treatment vs. average count
# graph3 <- ggplot(sb.sum, aes(x = Treat.Code, y = Count)) +
#   geom_point(aes(color = Species_Name, shape = Species_Name), size = 6) +
#   geom_errorbar(aes(ymin = Count - se, ymax = Count + se), width = 0.1, size = 0.6) +
#   theme_bw() +
#   theme(legend.title = element_blank(),  
#         legend.text = element_text(size = 20),
#         legend.key.size = unit(3, 'lines'),
#         axis.text=element_text(size = 20), 
#         plot.title = element_text(size = 30, face = "bold", margin = margin(0, 0, 20, 0)), 
#         axis.title = element_text(size = 30), 
#         strip.text = element_text(size = 25), 
#         axis.title.y = element_text(margin = margin(0, 20, 0, 0)),
#         axis.title.x = element_blank()) +
#   ylim(0, 100) + 
#   labs(y = "Number of Viable Seeds")
# 
# ggsave(graph3, filename = "dotplot-sppxtrt.png", path = plots, width = 12, height = 8, dpi = 300)
# 
