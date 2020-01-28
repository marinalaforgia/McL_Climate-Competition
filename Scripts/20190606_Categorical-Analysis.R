# Paper Analysis 2019-06-06: Categorical Analysis of Functional Groups

rm(list=ls()) 

#### Load Libraries ####
library(Rmisc)
library(ggplot2)
library(plyr)
library(dplyr)
library(lme4)
library(lmerTest)
library(reshape2)
library(cowplot)
library(car)
library(multcomp)
library(tidyverse)
library(emmeans)
library(sjstats)

#### Prep: Vital Rates ####

###
# Grass
###
grass <- read.csv("Data/Post-Processing/grass-cover.csv")
names(grass)[3] <- "cover.g"

treat <- read.csv("Data/Setup/Marina-Treatment-30.csv")

###
# Germ and mortality data. 
###

dem.16 <- read.csv("Data/Post-Processing/dem-data-16.csv")

# removed plots 85, 86, and 87 due to burrow damage early in the season
# removed grass plots in 82 and 89 because less than 10% grass cover

dem.16 <- filter(dem.16, Plot != 87, Plot != 85, Plot != 86, !(Plot == 82 & Subplot == "Grass"), !(Plot == 89 & Subplot == "Grass"))

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

# flower number = number of flowers per individual from date of highest number flowering 
#flo.long <- read.csv("Data/Post-Processing/final-flo-long.csv")
#seed.long <- read.csv("Data/Post-Processing/final-seed-long.csv")
flo.seed <- read.csv("Data/Post-Processing/final-flo-seed.csv")
flo.seed <- flo.seed[flo.seed$Subplot != "Thatch",]
flo.seed <- filter(flo.seed, !(Year == 2017 & Treat.Code == "D")) 
flo.seed$Subplot <- factor(flo.seed$Subplot, levels = c("No Grass", "Grass"))
flo.seed$n.seed.ind <- flo.seed$n.seed.ind * flo.seed$p.viable # adjust for viability

###
# Seed survival
###
sb <- read.csv("Data/Post-Processing/seed-carryover-plot.csv")[,c(1:9)]
sb <- filter(sb, !(Year == 2017 & Treat.Code == "D"), !(Year == 2016 & Plot == 85), !(Year == 2016 & Plot == 86), !(Year == 2016 & Plot == 87))
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
trait.w$strat <- ifelse(trait.w$PC.F > 0, "ST", "SA")

# merge datasets with traits
dem <- merge(dem, trait.w[,c(1,26)], by = "Species")
sb <- merge(sb, trait.w[,c(1,26)], by = "Species")
flo.seed <- merge(flo.seed, trait.w[,c(1,26)], by = "Species")

rm(PCA.F, PCA.G, PCA.s13, sla.13)

# merge datasets with grass cover
# dem <- merge(dem, grass, all.x = T)
# dem[is.na(dem$cover.g),]$cover.g <- 0

#### Prep: Plot-level Lambda ####

### Merge with flowering data ###

full <- merge(dem, flo.seed[,-c(6:9,11:12)], by = c("Year","Plot","Treat.Code", "Subplot","Species", "strat"), all = T)

# replace NAs in n.seed.inf with species averages
flo.seed.sum <- summarySE(flo.seed, measurevar = "avg.seed.inf", groupvars = c("Species"), na.rm = T)

for(j in unique(full[is.na(full$avg.seed.inf),]$Species)) {
    full[is.na(full$avg.seed.inf) & full$Species == j,]$avg.seed.inf <- flo.seed.sum[flo.seed.sum$Species == j,]$avg.seed.inf
  }

# replace NAs in avg.flo with species averages
flo.seed.sum <- summarySE(flo.seed, measurevar = "avg.flo", groupvars = c("Species"), na.rm = T)

for(j in unique(full[is.na(full$avg.flo),]$Species)) {
    full[is.na(full$avg.flo) & full$Species == j,]$avg.flo <- flo.seed.sum[flo.seed.sum$Species == j,]$avg.flo
  }

# Put viability estimates back in 
flo.seed.sum <- summarySE(flo.seed, measurevar = "p.viable", groupvars = c("Species", "Year"), na.rm = T)

for(j in unique(full[is.na(full$p.viable) & full$Year == 2017,]$Species)) {
    full[is.na(full$p.viable) & full$Species == j & full$Year == 2017,]$p.viable <- flo.seed.sum[flo.seed.sum$Species == j & flo.seed.sum$Year == 2017,]$p.viable
  }

# Multiply together to get n.seed.ind
full$n.seed.ind <- full$avg.seed.inf * full$avg.flo * full$p.viable

# merge with seed bank data
full <- merge(full, sb, by = c("Year","Plot","Treat.Code","Species", "strat"), all = T)

# Replace seed surv NA with species averages
full[full$Species == "Agoseris heterophylla" & is.na(full$p.surv),]$p.surv <- sb.spp[sb.spp$Species == "Agoseris heterophylla",]$p.sd.surv

full[full$Species == "Clarkia purpurea" & is.na(full$p.surv),]$p.surv <- sb.spp[sb.spp$Species == "Clarkia purpurea",]$p.sd.surv

# avg germ for PLER and LACA from 2017 due to identification mix-up
germ.sum <- summarySE(filter(dem, Year != 2016), measurevar = "p.germ", groupvars = c("Species"), na.rm = T)

full[full$Species == "Lasthenia californica" & full$Year == 2016,]$p.germ <- germ.sum[germ.sum$Species == "Lasthenia californica",]$p.germ

full[full$Species == "Plantago erecta" & full$Year == 2016,]$p.germ <- germ.sum[germ.sum$Species == "Plantago erecta",]$p.germ

# calculate L according: L = s*(1-g) + g*(1-m)*F
full$L.sb <- full$p.surv*(1-full$p.germ)
full$L.sa <- full$p.germ * (1 -  full$p.mort)
full$L.seeds <- full$L.sa * full$n.seed.ind
full$L <- full$L.sb + full$L.seeds
full$L <- ifelse(full$L.sa == 0, full$L.sb, full$L)
full$L <- ifelse(is.na(full$L.sa) == T, full$L.sb, full$L)

#### M0: Grass Cover ####
grass <- filter(grass, !(Year == 2017 & Treat.Code == "D"), Subplot == "Grass") 

grass <- filter(grass, !(Plot == 87 & Year == 2016), !(Plot == 85 & Year == 2016), !(Plot == 86 & Year == 2016), !(Plot == 82 & Year == 2016), !(Plot == 89 & Year == 2016))

m0 <- lmer(log(cover.g) ~ Treat.Code + (1|Year), grass)
plot(fitted(m0), resid(m0))
qqnorm(resid(m0)) 
qqline(resid(m0), col = 2, lwd = 2, lty = 2) 
summary(m0) 

grass.sum <- summarySE(grass, groupvars = c("Treat.Code"), measurevar = "cover.g")

## compare to grass cover in core plots
# grass.core <- read.csv("~/Documents/UC-Davis/Projects/McL_SH_Watering-Exp/Annual-Processing/Data/Core_Plot_2019.csv")
# treat2 <- read.csv("~/Documents/UC-Davis/Projects/McL_SH_Watering-Exp/Annual-Processing/Data/Treatment.csv")
# treat2 <- filter(treat2, Serpentine == "S")
# #grass.core <- merge(grass.core, unique(treat[,c(1,3)]))
# grass.core <- merge(grass.core, treat2[,c(1:2)])
# grass.core <- filter(grass.core[,c(1,29,2,9)], Year == 2016 | Year == 2017)
# 
# grass.core.sum <- summarySE(grass.core, groupvars = c("Treatment", "Year"), measurevar = "CoverTotal.Exotic.Annual.Grass") # same as in my subplots

#### M1: Mortality ####
dem$Year <- as.factor(dem$Year)

m1.trait <- glmer(cbind(tot.mort, germ.proj - tot.mort) ~ Treat.Code * Subplot * strat + (1|Year:Plot:Species), family = binomial, dem, glmerControl(calc.derivs = F))
plot(fitted(m1.trait), resid(m1.trait))
qqnorm(resid(m1.trait)) 
qqline(resid(m1.trait), col = 2, lwd = 2, lty = 2) 
summary(m1.trait) 

# contrasts
D.N.ST <- c(1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0)
D.G.ST <- c(1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0)
W.N.ST <- c(1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0)
W.G.ST <- c(1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1)
C.N.ST <- c(1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0)
C.G.ST <- c(1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0)

D.N.SA <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
D.G.SA <- c(1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0)
W.N.SA <- c(1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
W.G.SA <- c(1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0)
C.N.SA <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
C.G.SA <- c(1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0)

K <- rbind("D.N.ST - C.N.ST" = D.N.ST - C.N.ST,
           "W.N.ST - C.N.ST" = W.N.ST - C.N.ST,
           "D.G.ST - C.G.ST" = D.G.ST - C.G.ST,
           "W.G.ST - C.G.ST" = W.G.ST - C.G.ST,
           
           "D.G.ST - D.N.ST" = D.G.ST - D.N.ST,
           "W.G.ST - W.N.ST" = W.G.ST - W.N.ST,
           "C.G.ST - C.N.ST" = C.G.ST - C.N.ST,
           
           "D.N.SA - C.N.SA" = D.N.SA - C.N.SA,
           "W.N.SA - C.N.SA" = W.N.SA - C.N.SA,
           "D.G.SA - C.G.SA" = D.G.SA - C.G.SA,
           "W.G.SA - C.G.SA" = W.G.SA - C.G.SA,
           
           "D.G.SA - D.N.SA" = D.G.SA - D.N.SA,
           "W.G.SA - W.N.SA" = W.G.SA - W.N.SA,
           "C.G.SA - C.N.SA" = C.G.SA - C.N.SA,
          
           "D.G.ST - D.G.SA" = D.G.ST - D.G.SA,
           "W.G.ST - W.G.SA" = W.G.ST - W.G.SA,
           "C.G.ST - C.G.SA" = C.G.ST - C.G.SA,

           "D.N.ST - D.N.SA" = D.N.ST - D.N.SA,
           "W.N.ST - W.N.SA" = W.N.ST - W.N.SA,
           "C.N.ST - C.N.SA" = C.N.ST - C.N.SA)
           
summary(glht(m1.trait, linfct = K), test = adjusted("BH"))           


#### M2: Seed Set ####

hist(flo.seed$n.seed.ind)
hist(log(flo.seed$n.seed.ind + 1))

m2.trait0 <- lmer(log(n.seed.ind + 1) ~ Treat.Code * Subplot * strat + (1|Year:Plot:Species), flo.seed)

m2.trait <- lmer(log(n.seed.ind + 1) ~ Treat.Code * Subplot * strat - Treat.Code : Subplot : strat + (1|Year:Plot:Species), flo.seed)
anova(m2.trait0, m2.trait)

plot(fitted(m2.trait), resid(m2.trait))
qqnorm(resid(m2.trait)) 
qqline(resid(m2.trait), col = 2, lwd = 2, lty = 2) 
summary(m2.trait)

# contrasts
D.N.ST <- c(1, 1, 0, 0, 1, 0, 0, 1, 0, 0)
D.G.ST <- c(1, 1, 0, 1, 1, 1, 0, 1, 0, 1)
W.N.ST <- c(1, 0, 1, 0, 1, 0, 0, 0, 1, 0)
W.G.ST <- c(1, 0, 1, 1, 1, 0, 1, 0, 1, 1)
C.N.ST <- c(1, 0, 0, 0, 1, 0, 0, 0, 0, 0)
C.G.ST <- c(1, 0, 0, 1, 1, 0, 0, 0, 0, 1)

D.N.SA <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0)
D.G.SA <- c(1, 1, 0, 1, 0, 1, 0, 0, 0, 0)
W.N.SA <- c(1, 0, 1, 0, 0, 0, 0, 0, 0, 0)
W.G.SA <- c(1, 0, 1, 1, 0, 0, 1, 0, 0, 0)
C.N.SA <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
C.G.SA <- c(1, 0, 0, 1, 0, 0, 0, 0, 0, 0)

K3 <- rbind("D.N.ST - C.N.ST" = D.N.ST - C.N.ST,
           "W.N.ST - C.N.ST" = W.N.ST - C.N.ST,
           "D.G.ST - C.G.ST" = D.G.ST - C.G.ST,
           "W.G.ST - C.G.ST" = W.G.ST - C.G.ST,

           "D.N.SA - C.N.SA" = D.N.SA - C.N.SA,
           "W.N.SA - C.N.SA" = W.N.SA - C.N.SA,
           "D.G.SA - C.G.SA" = D.G.SA - C.G.SA,
           "W.G.SA - C.G.SA" = W.G.SA - C.G.SA,

           "D.G.SA - D.N.SA" = D.G.SA - D.N.SA,
           "W.G.SA - W.N.SA" = W.G.SA - W.N.SA,
           "C.G.SA - C.N.SA" = C.G.SA - C.N.SA,

           "D.G.ST - D.N.ST" = D.G.ST - D.N.ST,
           "W.G.ST - W.N.ST" = W.G.ST - W.N.ST,
           "C.G.ST - C.N.ST" = C.G.ST - C.N.ST,
           
           "D.G.ST - D.G.SA" = D.G.ST - D.G.SA,
           "W.G.ST - W.G.SA" = W.G.ST - W.G.SA,
           "C.G.ST - C.G.SA" = C.G.ST - C.G.SA,

           "D.N.ST - D.N.SA" = D.N.ST - D.N.SA,
           "W.N.ST - W.N.SA" = W.N.ST - W.N.SA,
           "C.N.ST - C.N.SA" = C.N.ST - C.N.SA)

            
summary(glht(m2.trait, linfct = K3), test = adjusted("BH"))    


#### M3: Seed Carryover ####
m3.trait <- glmer(cbind(Count, avg.num.per.sub - Count) ~  Treat.Code * strat + (1|Year:Plot:Species), family = binomial, data = sb)
plot(fitted(m3.trait), resid(m3.trait))
qqnorm(resid(m3.trait))
qqline(resid(m3.trait), col = 2, lwd = 2, lty = 2)
summary(m3.trait) 

#### M4: Germination ####
# removing PLER 2016 and LACA 2016 because of mis-ID
m4.trait <- glmer(cbind(germ.tot, viable-germ.tot) ~ strat + (1|Year:Plot:Species), family = binomial, data = filter(dem, !(Species == "Lasthenia californica" & Year == 2016), !(Species == "Plantago erecta" & Year == 2016)))
plot(fitted(m4.trait), resid(m4.trait))
qqnorm(resid(m4.trait))
qqline(resid(m4.trait), col = 2, lwd = 2, lty = 2)
summary(m4.trait)

#### M6: Lambda - NAF ####

hist(log(full$L + 0.5)) # add small constant to normalize 
m5.trait2 <- lmer(log(L + .5) ~ Treat.Code * Subplot + (1|Year:Plot:Species), data = full)
plot(fitted(m5.trait2), resid(m5.trait2))
qqnorm(resid(m5.trait2))
qqline(resid(m5.trait2), col = 2, lwd = 2, lty = 2)
summary(m5.trait2)

#### M5: Lambda - Grass ####
hist(log(full$L + .5)) # add small constant to normalize 
m5.trait3 <- lmer(log(L + .5) ~ strat * Subplot + (1|Year:Plot:Species), data = full)
plot(fitted(m5.trait3), resid(m5.trait3))
qqnorm(resid(m5.trait3))
qqline(resid(m5.trait3), col = 2, lwd = 2, lty = 2)
summary(m5.trait3) 

#### M6: Lambda - strat ####
m5.trait <- lmer(log(L + .5) ~ Treat.Code * Subplot * strat + (1|Year:Plot:Species), data = full)
plot(fitted(m5.trait), resid(m5.trait))
qqnorm(resid(m5.trait))
qqline(resid(m5.trait), col = 2, lwd = 2, lty = 2)
summary(m5.trait)
summary(glht(m5.trait, linfct = K), test = adjusted("BH"))  
  
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
  labs(x = "Carbon isotope discrimination", y = "Relative growth rate")
  #labs(x = "Carbon isotope discrimination (∆, \u2030)", y = expression(paste("Relative Growth Rate (", cm^{2}, ")"%.%"(", cm^{2}, ")"^{-1}%.%day^{-1})))

ggsave(plot.trait, filename = "Figures/trait.tiff", width = 3, height = 2.7, units = "in", dpi = 600)

#### Fig 2: Lambda NAF ####
(full.all <- summarySE(full, groupvars = c("Treat.Code", "Subplot"), measurevar = "L", na.rm = T))
full.all$Subplot <- revalue(full.all$Subplot, c("No Grass" = "No grass"))

fig.lam.all <- ggplot(full.all, aes(y = L, x = Subplot, col = Treat.Code, group = Treat.Code)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = L - se, ymax = L + se, width = 0.06)) +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(size = 10), 
        plot.title = element_text(size=30, face="bold", vjust = 2),
        axis.title.y = element_text(size = 15), 
        axis.text.x = element_text(size = 13, colour = "black"),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 10),
        legend.position = c(0.77, 0.74),
        legend.key.size = unit(1.5, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_color_manual(values = c("grey40", "red1", "blue"), labels = c("Control", "Drought", "Watered")) +
  labs(y = "Lambda")

ggsave(fig.lam.all, filename = "Figures/Fig-lam-all.tiff", width = 3, height = 2.7, units = "in", dpi = 600)

#### Fig 3: Lambda Grass ####
(full.grass <- summarySE(full, groupvars = c("strat", "Subplot"), measurevar = "L", na.rm = T))
full.grass$Subplot <- revalue(full.grass$Subplot, c("No Grass" = "No grass"))

fig.lam.g <- ggplot(full.grass, aes(y = L, x = Subplot, col = strat, group = strat)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = L - se, ymax = L + se, width = 0.06)) +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(size = 10), 
        plot.title = element_text(size=30, face="bold", vjust = 2),
        axis.title.y = element_text(size = 15), 
        axis.text.x = element_text(size = 13, colour = "black"),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 10),
        legend.position = c(0.66, 0.81),
        legend.key.size = unit(1.5, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_color_manual(values = c("darkcyan", "#440154FF"), labels = c("Acquisitive", "Conservative")) +
  #scale_color_viridis_d() +
  ylim(0,25) +
  labs(y = "Lambda")

ggsave(fig.lam.g, filename = "Figures/Fig-lam-grass.tiff", width = 3, height = 2.7, units = "in", dpi = 600)

#### Fig 4: Lambda ####
full.sum <- summarySE(full, measurevar = "L", groupvars = c("Treat.Code", "Subplot", "strat"), na.rm = T)
full.sum$strat <- revalue(full.sum$strat, c("SA" = "Acquisitive", "ST" = "Conservative"))
full.sum$Subplot <- revalue(full.sum$Subplot, c("No Grass" = "No grass"))

Fig.lam <- ggplot(full.sum, aes(y = L, x = Subplot, col = Treat.Code, group = Treat.Code)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = L - se, ymax = L + se, width = 0.06)) +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(size = 10), 
        plot.title = element_text(size=30, face="bold", vjust = 2),
        axis.title.y = element_text(size = 15), 
        axis.text.x = element_text(size = 13, colour = "black"),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = c(0.87, 0.75),
        legend.key.size = unit(1.5, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_color_manual(values = c("grey40", "red1", "blue"), labels = c("Control", "Drought", "Watered")) +
  labs(y = "Lambda") +
  facet_wrap(~ strat) 

ggsave(Fig.lam, filename = "Figures/Fig-lam.tiff", width = 6, height = 3, units = "in", dpi = 600)

#### Fig 5: Mortality ####
dem.sum2 <- summarySE(dem, measurevar = "p.mort", groupvars = c("Treat.Code", "Subplot", "strat"), na.rm = T)
dem.sum2$strat <- revalue(dem.sum2$strat, c("SA" = "Acquisitive", "ST" = "Conservative"))
dem.sum2$Subplot <- revalue(dem.sum2$Subplot, c("No Grass" = "No grass"))

Fig.Mort <- ggplot(dem.sum2, aes(y = p.mort, x = Subplot, col = Treat.Code, group = Treat.Code)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = p.mort - se, ymax = p.mort + se, width = 0.06)) +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(size = 10), 
        plot.title = element_text(size=30, face="bold", vjust = 2),
        axis.title.y = element_text(size = 15), 
        axis.text.x = element_text(size = 13, colour = "black"),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = c(0.12, 0.75),
        legend.key.size = unit(1.5, 'lines'),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_color_manual(values = c("grey40", "red1", "blue"), labels = c("Control", "Drought", "Watered")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,.85)) +
  labs(y = "Mortality") +
  facet_wrap(~ strat) 

ggsave(Fig.Mort, filename = "Figures/Fig-mort.tiff", width = 6, height = 3, units = "in", dpi = 600)

#### Fig 6: Seed Set ####
flo.seed.sum2 <- summarySE(flo.seed, measurevar = "n.seed.ind", groupvars = c("Treat.Code", "Subplot", "strat"), na.rm = T)
flo.seed.sum2$strat <- revalue(flo.seed.sum2$strat, c("SA" = "Acquisitive", "ST" = "Conservative"))
flo.seed.sum2$Subplot <- revalue(flo.seed.sum2$Subplot, c("No Grass" = "No grass"))

Fig.seed <- ggplot(flo.seed.sum2, aes(y = n.seed.ind, x = Subplot, col = Treat.Code, group = Treat.Code)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = n.seed.ind - se, ymax = n.seed.ind + se, width = 0.06)) +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(size = 10), 
        plot.title = element_text(size = 30, face="bold", vjust = 2),
        axis.title.y = element_text(size = 15), 
        axis.text.x = element_text(size = 13, colour = "black"),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = c(0.87, 0.75),
        legend.key.size = unit(1.5, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_color_manual(values = c("grey40", "red1", "blue"), labels = c("Control", "Drought", "Watered")) +
  labs(y = "Seed Set") +
  facet_wrap(~ strat) 

ggsave(Fig.seed, filename = "Figures/Fig-seed.tiff", width = 6, height = 3, units = "in", dpi = 600)


#### Fig. S1: Germination ####
germ.sum2 <- summarySE(filter(dem, !(Species == "Lasthenia californica" & Year == 2016), !(Species == "Plantago erecta" & Year == 2016)), measurevar = "p.germ", groupvars = c("strat"), na.rm = T)
germ.sum2$strat <- revalue(germ.sum2$strat, c("SA" = "Acquisitive", "ST" = "Conservative"))

Fig.germ <- ggplot(germ.sum2, aes(y = p.germ, x = strat)) +
  geom_bar(stat="identity", fill = "darkcyan") +
  geom_errorbar(aes(ymin = p.germ - se, ymax = p.germ + se, width = 0.06)) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size=30, face="bold", vjust = 2),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 9),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = "none",
        legend.key.size = unit(1.5, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, .7)) +
  labs(y = "Germination")

ggsave(Fig.germ, filename = "Figures/Fig-germ.tiff", width = 3, height = 2.5, units = "in", dpi = 600)

#### Fig. S2: Seed survival ####
surv.sum <- summarySE(full, measurevar = "p.surv", groupvars = c("strat", "Treat.Code"), na.rm = T)
surv.sum$Treat.Code <- factor(surv.sum$Treat.Code, c("D", "C", "W"))
surv.sum$strat <- revalue(surv.sum$strat, c("SA" = "Acquisitive", "ST" = "Conservative"))

Fig.surv <- ggplot(surv.sum, aes(y = p.surv, x = strat, fill = Treat.Code, group = Treat.Code)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  geom_errorbar(position = position_dodge(0.9), aes(ymin = p.surv - se, ymax = p.surv + se, width = 0.1)) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size=30, face="bold", vjust = 2),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 9),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 9),
        legend.position = c(.23, .77),
        legend.key.size = unit(1, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
  labs(y = "Seed survival") +
  scale_fill_manual(values = c("grey40", "red1", "blue"), labels = c("Control", "Drought", "Watered"))

ggsave(Fig.surv, filename = "Figures/Fig-surv.tiff", width = 3, height = 2.5, units = "in", dpi = 600)

#### Fig. S3: Spp level lambda ####
full.sum2 <- summarySE(full, measurevar = "L", groupvars = c("Treat.Code", "Subplot", "Species", "strat", "Year"), na.rm = T)
full.sum2$strat <- revalue(full.sum2$strat, c("SA" = "Drought avoider", "ST" = "Drought tolerator"))
full.sum2$Subplot <- revalue(full.sum2$Subplot, c("No Grass" = "No grass"))

Fig.SI6.a <- ggplot(full.sum2[full.sum2$Year == 2017,], aes(y = L, x = Subplot, col = Treat.Code, group = Treat.Code)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = L - se, ymax = L + se, width = 0.06)) +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(size = 10), 
        #plot.title = element_text(size=30, face="bold", vjust = 2),
        axis.title.y = element_text(size = 12), 
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        #legend.position = c(0.9, 0.68),
        legend.key.size = unit(1.5, 'lines'),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  #scale_color_manual(values = c("grey40", "red1", "blue"), labels = c("Control", "Drought", "Watered")) +
  labs(y = "Lambda") +
  #scale_y_continuous(limits = c(0,45)) +
  facet_wrap(~ Species) 

ggsave(Fig.SI6.a, filename = "Figures/Fig-SI6-a.tiff", width = 6, height = 2.5, units = "in", dpi = 600)

Fig.SI6.b <- ggplot(full.sum2[full.sum2$strat == "Drought tolerator",], aes(y = L, x = Subplot, col = Treat.Code, group = Treat.Code)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = L - se, ymax = L + se, width = 0.06)) +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(size = 10), 
        #plot.title = element_text(size=30, face="bold", vjust = 2),
        axis.title.y = element_text(size = 12), 
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = "none",
        legend.key.size = unit(1.5, 'lines'),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_color_manual(values = c("grey40", "red1", "blue"), labels = c("Control", "Drought", "Watered")) +
  labs(y = "Lambda") +
  scale_y_continuous(limits = c(0,45)) +
  facet_wrap(~ Species) 

ggsave(Fig.SI6.b, filename = "Figures/Fig-SI6-b.tiff", width = 6, height = 2.5, units = "in", dpi = 600)

#### Fig. S4: Spp mortality ####
dem.sum3 <- summarySE(dem, measurevar = "p.mort", groupvars = c("Treat.Code", "Subplot", "Species", "strat"), na.rm = T)
dem.sum3$strat <- revalue(dem.sum3$strat, c("SA" = "Drought avoider", "ST" = "Drought tolerator"))
dem.sum3$Subplot <- revalue(dem.sum3$Subplot, c("No Grass" = "No grass"))

Fig.MortA.SI <- ggplot(dem.sum3[dem.sum3$strat == "Drought avoider",], aes(y = p.mort, x = Subplot, col = Treat.Code, group = Treat.Code)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = p.mort - se, ymax = p.mort + se, width = 0.06)) +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(size = 10), 
        #plot.title = element_text(size=30, face="bold", vjust = 2),
        axis.title.y = element_text(size = 12), 
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = c(0.9, 0.7),
        legend.key.size = unit(1.5, 'lines'),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_color_manual(values = c("grey40", "red1", "blue"), labels = c("Control", "Drought", "Watered")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,.85)) +
  labs(y = "Mortality") +
  facet_wrap(~ Species) 

ggsave(Fig.MortA.SI, filename = "Figures/Fig-SI1-a.tiff", width = 6, height = 2.5, units = "in", dpi = 600)

Fig.MortT.SI <- ggplot(dem.sum3[dem.sum3$strat == "Drought tolerator",], aes(y = p.mort, x = Subplot, col = Treat.Code, group = Treat.Code)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = p.mort - se, ymax = p.mort + se, width = 0.06)) +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(size = 10), 
        #plot.title = element_text(size=30, face="bold", vjust = 2),
        axis.title.y = element_text(size = 12), 
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = "none",
        legend.key.size = unit(1.5, 'lines'),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_color_manual(values = c("grey40", "red1", "blue"), labels = c("Control", "Drought", "Watered")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,.85)) +
  labs(y = "Mortality") +
  facet_wrap(~ Species)

ggsave(Fig.MortT.SI, filename = "Figures/Fig-SI1-b.tiff", width = 6, height = 2.5, units = "in", dpi = 600)

#### Fig. S5: Spp seed set ####
flo.seed.sum3 <- summarySE(flo.seed, measurevar = "n.seed.ind", groupvars = c("Treat.Code", "Subplot", "Species", "strat"), na.rm = T)
flo.seed.sum3$Species <- factor(flo.seed.sum3$Species, levels = c("Agoseris heterophylla", "Lasthenia californica", "Plantago erecta", "Clarkia purpurea", "Hemizonia congesta", "Calycadenia pauciflora"))
flo.seed.sum3$strat <- revalue(flo.seed.sum3$strat, c("SA" = "Drought avoider", "ST" = "Drought tolerator"))
flo.seed.sum3$Subplot <- revalue(flo.seed.sum3$Subplot, c("No Grass" = "No grass"))

Fig.SI2.a <- ggplot(flo.seed.sum3[flo.seed.sum3$strat == "Drought avoider",], aes(y = n.seed.ind, x = Subplot, col = Treat.Code, group = Treat.Code)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = n.seed.ind - se, ymax = n.seed.ind + se, width = 0.06)) +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(size = 10), 
        #plot.title = element_text(size=30, face="bold", vjust = 2),
        axis.title.y = element_text(size = 12), 
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = c(0.9, 0.68),
        legend.key.size = unit(1.5, 'lines'),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_color_manual(values = c("grey40", "red1", "blue"), labels = c("Control", "Drought", "Watered")) +
  scale_y_continuous(limits = c(0,75)) +
  labs(y = "Seed Set") +
  facet_wrap(~ Species) 

ggsave(Fig.SI2.a, filename = "Figures/Fig-SI2-a.tiff", width = 6, height = 2.5, units = "in", dpi = 600)

Fig.SI2.b <- ggplot(flo.seed.sum3[flo.seed.sum3$strat == "Drought tolerator",], aes(y = n.seed.ind, x = Subplot, col = Treat.Code, group = Treat.Code)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = n.seed.ind - se, ymax = n.seed.ind + se, width = 0.06)) +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(size = 10), 
        #plot.title = element_text(size=30, face="bold", vjust = 2),
        axis.title.y = element_text(size = 12), 
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = "none",
        legend.key.size = unit(1.5, 'lines'),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_color_manual(values = c("grey40", "red1", "blue"), labels = c("Control", "Drought", "Watered")) +
  labs(y = "Seed Set") +
  scale_y_continuous(limits = c(0,75)) +
  facet_wrap(~ Species) 

ggsave(Fig.SI2.b, filename = "Figures/Fig-SI2-b.tiff", width = 6, height = 2.5, units = "in", dpi = 600)

#### Fig. S6: Spp germination ####
germ.sum3 <- summarySE(filter(dem, !(Species == "Lasthenia californica" & Year == 2016), !(Species == "Plantago erecta" & Year == 2016)), measurevar = "p.germ", groupvars = c("strat", "Species"), na.rm = T)
#germ.sum3 <- summarySE(dem, measurevar = "p.germ", groupvars = c("strat", "Species"), na.rm = T)
germ.sum3$strat <- revalue(germ.sum3$strat, c("SA" = "Acquisitive", "ST" = "Conservative"))
germ.sum3$Species <- factor(germ.sum3$Species, levels = c("Agoseris heterophylla", "Lasthenia californica", "Plantago erecta", "Clarkia purpurea", "Hemizonia congesta", "Calycadenia pauciflora"))

addline_format <- function(x,...){
    gsub('\\s','\n',x)
}

Fig.SI3 <- ggplot(germ.sum3, aes(y = p.germ, x = Species, fill = strat)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = p.germ - se, ymax = p.germ + se, width = 0.06)) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = c(.8, .8),
        legend.key.size = unit(1.5, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
  scale_x_discrete(labels = addline_format(c("Agoseris heterophylla", "Lasthenia californica", "Plantago erecta", "Clarkia purpurea", "Hemizonia congesta", "Calycadenia pauciflora"))) +
  labs(y = "Germination")

ggsave(Fig.SI3, filename = "Figures/Fig-SI3.tiff", width = 6, height = 3, units = "in", dpi = 600)

#### Fig. S7: Spp seed survival ####
surv.sum2 <- summarySE(sb, measurevar = "p.surv", groupvars = c("strat", "Species", "Treat.Code"), na.rm = T)
surv.sum2$Treat.Code <- factor(surv.sum2$Treat.Code, c("D", "C", "W"))
surv.sum2$strat <- revalue(surv.sum2$strat, c("SA" = "Drought avoider", "ST" = "Drought tolerator"))
surv.sum2$Treat.Code <- revalue(surv.sum2$Treat.Code, c("D" = "Drought", "C" = "Control", "W" = "Watered"))
surv.sum2$Species <- factor(surv.sum2$Species, levels = c("Agoseris heterophylla", "Lasthenia californica", "Plantago erecta", "Clarkia purpurea", "Hemizonia congesta", "Calycadenia pauciflora"))

Fig.SI5.a <- ggplot(surv.sum2[surv.sum2$strat == "Drought avoider",], aes(y = p.surv, x = Species, fill = Treat.Code, group = Treat.Code)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  geom_errorbar(position = position_dodge(0.9), aes(ymin = p.surv - se, ymax = p.surv + se, width = 0.1)) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size=30, face="bold", vjust = 2),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 10),
        legend.position = c(.1, .74),
        legend.key.size = unit(1.5, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
  labs(y = "Seed Survival") +
  #scale_x_discrete(labels = addline_format(c("Lasthenia californica", "Plantago erecta", "Agoseris heterophylla", "Clarkia purpurea", "Hemizonia congesta", "Calycadenia pauciflora"))) +
  scale_fill_manual(values = c("red1", "grey40", "blue"))

ggsave(Fig.SI5.a, filename = "Figures/Fig-SI5-a.tiff", width = 6, height = 2.5, units = "in", dpi = 600)

Fig.SI5.b <- ggplot(surv.sum2[surv.sum2$strat == "Drought tolerator",], aes(y = p.surv, x = Species, fill = Treat.Code, group = Treat.Code)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  geom_errorbar(position = position_dodge(0.9), aes(ymin = p.surv - se, ymax = p.surv + se, width = 0.1)) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size=30, face="bold", vjust = 2),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 10),
        legend.position = "none",
        legend.key.size = unit(1.5, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  labs(y = "Seed Survival") +
  #scale_x_discrete(labels = addline_format(c("Lasthenia californica", "Plantago erecta", "Agoseris heterophylla", "Clarkia purpurea", "Hemizonia congesta", "Calycadenia pauciflora"))) +
  scale_fill_manual(values = c("red1", "grey40", "blue"))

ggsave(Fig.SI5.b, filename = "Figures/Fig-SI5-b.tiff", width = 6, height = 2.5, units = "in", dpi = 600)

#### Fig. S8: Greenhouse based traits ####
plot.trait2 <- ggplot(trait.w, aes(x = D13C.GH, y = RGRt.GH)) +
  theme_classic() +
  geom_smooth(method = "lm", se = F, col = "black", size = 1) +
  geom_text(aes(label = Species.short), hjust = .5, vjust = .5) +
  theme(axis.text = element_text(size = 10), 
        plot.title = element_text(size = 30, face="bold", vjust = 2),
        axis.title = element_text(size = 13), 
        strip.text = element_text(size = 15),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_x_reverse(lim = c(25.3, 22.3)) +
  labs(x = "Carbon isotope discrimination", y = "Relative growth rate")

# plot.trait3 <- ggplot(trait.w, aes(x = D13C.GH, y = RGR.la.GH)) +
#   theme_classic() +
#   geom_smooth(method = "lm", se = F, col = "black", size = 1) +
#   geom_text(aes(label = Species.short), hjust = .5, vjust = .5) +
#   theme(axis.text = element_text(size = 10), 
#         plot.title = element_text(size = 30, face="bold", vjust = 2),
#         axis.title = element_text(size = 13), 
#         strip.text = element_text(size = 15),
#         axis.line = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
#   scale_x_reverse(lim = c(25.2, 22.5)) +
#   labs(x = "Carbon isotope discrimination", y = "Relative growth rate")

ggsave(plot.trait2, filename = "Figures/trait-gh.tiff", width = 3, height = 2.5, units = "in", dpi = 600)

#### Table S1: observations per analysis ####
dem.sum3 <- dem.sum3[,c(1:3,5)]
names(dem.sum3)[4] <- "n.mort"

flo.seed.sum3 <- flo.seed.sum3[,c(1:3,5)]
names(flo.seed.sum3)[4] <- "n.set"

full.sum2  <- full.sum2[,c(1:3,5)]
names(full.sum2)[4] <- "n.lam"

surv.sum2 <- surv.sum2 [,c(1:4)]
names(surv.sum2)[4] <- "n.sdsur"

# pop.sum <- summarySE(dem, measurevar = "germ.proj", groupvars = c("Treat.Code", "Subplot", "Species", "strat"), na.rm = T)
# pop.sum <- pop.sum[,c(1:3,6)]
# names(pop.sum)[4] <- "n.pop"

#SI.table1 <- merge(pop.sum, full.sum2, by = c("Species", "Treat.Code", "Subplot"))
SI.table1 <- merge(full.sum2, dem.sum3, by = c("Species", "Treat.Code", "Subplot"))
SI.table1 <- merge(SI.table1, flo.seed.sum3, by = c("Species", "Treat.Code", "Subplot"))

write.csv(SI.table1, "Figures/SI-table1.csv", row.names = F)
write.csv(surv.sum2, "Figures/SI-table1ext.csv", row.names = F)
