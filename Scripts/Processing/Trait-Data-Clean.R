### Prepping all my trait data! ###

rm(list=ls())

#### Load Libraries ####
library(ggplot2)
library(plyr)
library(dplyr)
library(growthrates)
library(lubridate)
library(cowplot)
library(Rmisc)
library(nlme)
library(mvtnorm)
library(lme4)

#### Load Datasets ####

## Carbon Isotope Data
cn <- read.csv("Data/Traits/Carbon-Isotope/C-N-Results.csv")

## SLA Data
sla <- read.csv("Data/Traits/Trait-Data/SLA_2017.csv")

# Relative growth rate (biomass)
la.dim.f <- read.csv("Data/Traits/Trait-Data/dim-sla.csv")
la.dim.f$Date <- as.Date(la.dim.f$Date, "%m/%d/%y")

biomass.f <- read.csv("Data/Traits/Trait-Data/Field_FinalHarvest.csv")
traits.f <- read.csv("Data/Traits/Trait-Data/Field-traits-2017.csv")
traits.f$Date <- as.Date(traits.f$Date, "%m/%d/%y")

traits.gh <- read.csv("Data/Traits/Trait-Data/Greenhouse-traits-2017.csv")
traits.gh$Date <- as.Date(traits.gh$Date, "%m/%d/%y")

# Other trait data
Species.full <- c("Agoseris heterophylla", "Clarkia purpurea", "Plantago erecta", "Calycadenia pauciflora", "Hemizonia congesta", "Lasthenia californica")
trait <- read.csv("/Users/Marina/Google Drive/02_McLaughlin_80Sites_Organized/Modified_CombinedFiles/McL_80SitesSpeciesTraits_012615.csv")[,c(2,3)] #trait data
trait <- filter(trait, Species_Name %in% Species.full)
trait$Species_Name <- revalue(trait$Species_Name, c("Agoseris heterophylla" = "AGHE", "Clarkia purpurea" = "CLPU", "Lasthenia californica" = "LACA","Plantago erecta" = "PLER", "Hemizonia congesta" = "HECO", "Calycadenia pauciflora" = "CAPA"))

## Seed mass data
# seed.mass <- read.csv("~Marina/Documents/UC-Davis/Research/Phenology/Original/Traits/Kew database.csv")
# seed.mass <- seed.mass[seed.mass$id %in% Species.full,]
# seed.mass <- ddply(seed.mass, .(id), summarize, seed_mass_g = mean(seed_mass_g))
#
# # OR ##
#
# seed.mass <- read.csv("Data/Seed_weights.csv")
# seed.mass <- seed.mass[seed.mass$Species %in% Species.full,]

#### PREP DATASETS ####

### F: Final Biomass ####
biomass.f <- biomass.f[,c(1,2,5:7)]
biomass.f <- biomass.f[complete.cases(biomass.f),]
biomass.f <- filter(biomass.f, Species != "CETR")
biomass.f$total.g <- biomass.f$stem.g + biomass.f$leaf.g + biomass.f$repr.g
biomass.f <- ddply(biomass.f, .(Species), summarize, stem.g = mean(stem.g), leaf.g = mean(leaf.g), repr.g = mean(repr.g), total.g = mean(total.g))
biomass.f$gh_f <- "F"
biomass.f$root.g <- NA
biomass.f$root.shoot <- NA
biomass.f <- biomass.f[,c(1:4,7,5,6,8)]

### GH: Final Biomass ####
traits.gh <- filter(traits.gh, Species != "CETR")
biomass.gh <- traits.gh[,c(1:3,8:11)]
biomass.gh[is.na(biomass.gh)] <- 0
biomass.gh$total.g <- biomass.gh$root.g + biomass.gh$stem.g + biomass.gh$leaf.g + biomass.gh$repr.g
max.dates <- ddply(biomass.gh, .(Species), summarize, Date = max(Date))
biomass.gh <- merge(biomass.gh, max.dates, by = c("Species","Date"), all.x = F)
biomass.gh <- ddply(biomass.gh, .(Species), summarize, stem.g = mean(stem.g), leaf.g = mean(leaf.g), repr.g = mean(repr.g), root.g = mean(root.g), total.g = mean(total.g))
biomass.gh$gh_f <- "GH"
biomass.gh$root.shoot <- biomass.gh$root.g/(biomass.gh$stem.g + biomass.gh$leaf.g)
biomass <- rbind(biomass.f, biomass.gh)
rm(biomass.f, biomass.gh, max.dates)

### GH: Seedling root:shoot ratio ####
r.s <- filter(traits.gh, Date == "2017-02-21")
r.s <- r.s[,c(2,3,8:10)]
r.s[is.na(r.s)] <- 0
r.s$root.shoot <- r.s$root.g/(r.s$stem.g +r.s$leaf.g)
r.s <- ddply(r.s, .(Species), summarize, root.shoot = mean(root.shoot))
r.s$gh_f <- "GH"

#### GH: RGR (Mass) ####
simp.rgr <- traits.gh[,c(1:3,7:11)]
simp.rgr[is.na(simp.rgr)] <- 0

# get rid of CETR
simp.rgr <- filter(simp.rgr, Species != "CETR")

# add together all biomass
simp.rgr$total.g <- simp.rgr$root.g + simp.rgr$stem.g + simp.rgr$leaf.g + simp.rgr$repr.g

# simplify dataset for easier analysis
simp.rgr <- simp.rgr[,c(1,2,9)]
simp.rgr$log.m <- log(simp.rgr$total.g)
simp.rgr <- ddply(simp.rgr, .(Species), transform, day_number = as.numeric(difftime(Date, min(Date), units = "days")))

# ggplot(simp.rgr, aes(x = Date, y = total.g)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = F) +
#   facet_wrap(~ Species)
# 
# ggplot(simp.rgr, aes(x = Date, y = log.m)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = F) +
#   facet_wrap(~ Species) # even taking the log, they are clearly still asymptotic

# Calculate RGR using Paine-2012 suggestions, 4-parameter logistic fit best so continue with that
n.preds <- 100   
# before i wasnt adjusting for different number of days between species
Length <- list(PLER = data.frame(X = seq(min(simp.rgr$day_number), max(simp.rgr[simp.rgr$Species == "PLER",]$day_number), length = n.preds)), HECO = data.frame(X = seq(min(simp.rgr$day_number), max(simp.rgr[simp.rgr$Species == "HECO",]$day_number), length = n.preds)), CAPA = data.frame(X = seq(min(simp.rgr$day_number), max(simp.rgr[simp.rgr$Species == "CAPA",]$day_number), length = n.preds)), LACA = data.frame(X = seq(min(simp.rgr$day_number), max(simp.rgr[simp.rgr$Species == "LACA",]$day_number), length = n.preds)), AGHE = data.frame(X = seq(min(simp.rgr$day_number), max(simp.rgr[simp.rgr$Species == "AGHE",]$day_number), length = n.preds)), CLPU = data.frame(X = seq(min(simp.rgr$day_number), max(simp.rgr[simp.rgr$Species == "CLPU",]$day_number), length = n.preds)))
               
SPP <- unique(simp.rgr$Species)

rgr.avg <- data.frame(Species = SPP, RGRt = rep(NA, 6), RGRt.hi = rep(NA, 6), RGRt.lo = rep(NA, 6))

summarizer <- function(dat, alpha){ # this function returns confidence envelopes around growth trajectories, and growth rates. 
	n <- length(dat)
	quantiles <- c(alpha/2, 1-(alpha/2))
	CIs <- data.frame(matrix(NA, ncol(dat[[1]]), n*2))
	names(CIs) <- paste(rep(names(dat), each = 2), c("lo", "hi"), sep = ".")
	for(i in 1:n){
		CIs[,(2*i-1):(2*i)] <- t(apply(dat[[i]],    2, quantile, quantiles, na.rm = T))
		}
	return(CIs)
}

output.fpl.gnls <- function(fit, times, CI = F, LOG = F, alpha = 0.05){
	params <- coef <- coef(fit)
	M0 = params[1]; K = params[2]; xmid = params[3]; r = params[4]
	fitted <- fit$fitted
	resid  <- fit$residuals
	data <- data.frame(fitted = fitted, resid = resid)
	eq   <- bquote(paste(.(round(M0, 4))+ (.(round(K-M0, 4)) /1+e^{(.(round(xmid, 3))*t)/.(round(r, 3))})))
	mss <- sum((fitted - mean(fitted))^2)
	rss <- sum(resid^2)
	R2  <- mss/(mss + rss)
	rmse <- sqrt(rss)
	AIC <- AIC(fit)
	summary <- c(R2 = R2, AIC = AIC, RMSE = rmse)
	exp <- exp((xmid-times)/r)
	rates = data.frame(
		times = times,
		M    = M0+(K-M0)/(1+exp),
		AGR  = ((K-M0)*exp)/(r*(1+exp)^2)
		)
	rates$RGRt <- rates$AGR/rates$M
	rates$RGRm <- (M0-rates$M)*(K-rates$M)/(r*rates$M*(M0-K))
	if(LOG == T){
		rates$RGRt <- rates$AGR
		rates$RGRm <- (M0-rates$M)*(K-rates$M)/(r*(M0-K))
		rates$AGR  <- rates$AGR*exp(rates$M)
		}
	if(CI == T){
		cov   <- fit$varBeta
		x <- data.frame(rmvnorm(n=1000, mean=coef, sigma=cov))
		names(x) <- c("M0", "K", "xmid", "r")
		M <- AGR <- RGRt <- RGRm <- matrix(NA, ncol = length(times), nrow = nrow(x))
		for(i in 1:nrow(x)){
			x.i  <- x[i,]
			exp <- exp((x.i$xmid-times)/x.i$r)
			M[i,]     <- x.i$M0 + (x.i$K-x.i$M0)/(1+exp)
			AGR[i,]   <- ((x.i$K-x.i$M0)*exp)/(x.i$r*(1+exp)^2)
			RGRt[i,] <- AGR[i,]/M[i,]
			RGRm[i,]  <- x.i$r*(1 - (M[i,]/x.i$K)^(1/x.i$xmid))
#			RGRm[i,]  <- ((x.i$M0-M[i,])*(x.i$K-M[i,]))/(M[i,]*x.i$r*(x.i$M0-x.i$K))
			if(LOG == T){
				RGRt[i,] <- AGR[i,]
				RGRm[i,] <- (x.i$M0-M[i,])*(x.i$K-M[i,])/(x.i$r*(x.i$M0-x.i$K))
				AGR[i,]  <- AGR[i,]*exp(M[i,]) 
				}
		}
		CIs <- summarizer(list(M = M, AGR = AGR, RGRt = RGRt, RGRm = RGRm), alpha)
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = cbind(rates, CIs))
		} else {
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = rates)
		}
	return(out)
	}

for(i in SPP){
  
# ## Logistic fit ##
# tmp.logis  <- getInitial(log.m ~ SSlogis(day_number, Asym, xmid, scal), data = simp.rgr[simp.rgr$Species == i,])
# fit.logis <- gnls(log.m ~ SSlogis(day_number, Asym, xmid, scal), data = simp.rgr[simp.rgr$Species == i,], weights= varExp(form = ~ fitted(.)))
# out.logis <- output.logis.gnls(fit.logis, Xes_asymp$X, CI = T)

## 4-p Logistic fit ##
#tmp.fpl  <- getInitial(log.m ~ SSfpl(day_number, A, B, xmid, scal), data = simp.rgr[simp.rgr$Species == i,])
fit.fpl <- gnls(log.m ~ SSfpl(day_number, A, B, xmid, scal), data = simp.rgr[simp.rgr$Species == i,], weights= varExp(form = ~ fitted(.)))
out.fpl <- output.fpl.gnls(fit.fpl, Length[[i]]$X, CI = T, LOG = T)

# Compare all models for each species
par(mfrow = c(1,3))

### Summary  
plot(1:10, 31:40, type = "n", ann = F, axes = F)
text(0.5, 39, i, pos = 4)
text(0.5, 38, expression(paste("Model             ", R^2, "    ", "AIC", sep = "")), pos = 4)
text(0.5,  37, sprintf("4-p Logistic  %1.3f %4.1f", out.fpl$summary[1],   out.fpl$summary[2]),   pos = 4, col = "#1B9E77", font = 2)
#text(0.5, 34, sprintf("Logistic      %1.3f %4.1f", out.logis$summary[1], out.logis$summary[2]), pos = 4, col = COL[7])
par(xpd = F, family = "Helvetica", pty = "s")

### Biomass
plot(simp.rgr[simp.rgr$Species == i,]$day_number, simp.rgr[simp.rgr$Species == i,]$log.m, xlim = c(0, 150), xlab = "Days since germination", ylab = "Biomass (g)", col = "black", axes = F)
axis(1); axis(2, at = log(c(0.001, 0.01, 0.1, 1, 5, 15)), labels = c(0.001, 0.01, 0.1, 1, 5, 15), las =1)
#lines(out.logis$rates$times, out.logis$rates$M, col = COL[7])
lines(out.fpl$rates$times,   out.fpl$rates$M,   col = "#1B9E77")

##RGRt
plot(Length[[i]]$X,  out.fpl$rates$RGRt, xlab = "Days since germination", ylab = expression(paste("RGR ", (g%.%g^-1%.%day^-1))), type = "n", xlim = c(0, 150), ylim = c(0, 0.35))
#lines(Length[[i]]$X, out.logis$rates$RGRt, col = COL[7])  # Logistic
lines(Length[[i]]$X, out.fpl$rates$RGRt,   col = "#1B9E77")  # 4-param Logistic

# average RGR
  rgr.avg[rgr.avg$Species == i, 2] <- mean(out.fpl$rates$RGRt)  
  rgr.avg[rgr.avg$Species == i, 3] <- mean(out.fpl$rates$RGRt.hi)
  rgr.avg[rgr.avg$Species == i, 4] <- mean(out.fpl$rates$RGRt.lo)

}

# output: rgr.avg
rgr.avg$gh_f <- "GH"

rm(fit.fpl, out.fpl, Length, i, n.preds, simp.rgr)

##GH: RGR (Leaf Area) ####
la.rgr <- traits.gh[,c(1:3,12)]
la.rgr <- la.rgr[complete.cases(la.rgr),]

# get rid of CETR
la.rgr <- filter(la.rgr, Species != "CETR")

la.rgr$log.la <- log(la.rgr$leaf.area.cm)
la.rgr <- ddply(la.rgr, .(Species), transform, day_number = as.numeric(difftime(Date, min(Date), units = "days")))

n.preds <- 100   
# before i wasnt adjusting for different number of days between species
Length <- list(PLER = data.frame(X = seq(min(la.rgr$day_number), max(la.rgr[la.rgr$Species == "PLER",]$day_number), length = n.preds)), HECO = data.frame(X = seq(min(la.rgr$day_number), max(la.rgr[la.rgr$Species == "HECO",]$day_number), length = n.preds)), CAPA = data.frame(X = seq(min(la.rgr$day_number), max(la.rgr[la.rgr$Species == "CAPA",]$day_number), length = n.preds)), LACA = data.frame(X = seq(min(la.rgr$day_number), max(la.rgr[la.rgr$Species == "LACA",]$day_number), length = n.preds)), AGHE = data.frame(X = seq(min(la.rgr$day_number), max(la.rgr[la.rgr$Species == "AGHE",]$day_number), length = n.preds)), CLPU = data.frame(X = seq(min(la.rgr$day_number), max(la.rgr[la.rgr$Species == "CLPU",]$day_number), length = n.preds)))

rgr.la.gh <- data.frame(Species = SPP, RGRt = rep(NA, 6), RGRt.hi = rep(NA, 6), RGRt.lo = rep(NA, 6))


for(i in SPP){
  
## 4-p Logistic fit ##
#tmp.fpl  <- getInitial(log.la ~ SSfpl(day_number, A, B, xmid, scal), data = la.rgr[la.rgr$Species == i,])
fit.fpl <- gnls(log.la ~ SSfpl(day_number, A, B, xmid, scal), data = la.rgr[la.rgr$Species == i,], weights= varExp(form = ~ fitted(.)))
out.fpl <- output.fpl.gnls(fit.fpl, Length[[i]]$X, CI = T, LOG = T)

# Compare all models for each species
par(mfrow = c(1,3))

### Summary  
plot(1:10, 31:40, type = "n", ann = F, axes = F)
text(0, 39, i, pos = 4)
text(0, 38, expression(paste("Model             ", R^2, "    ", "AIC", sep = "")), pos = 4)
text(0,  37, sprintf("4-p Logis  %1.3f %4.1f", out.fpl$summary[1],   out.fpl$summary[2]),   pos = 4, col = "#1B9E77", font = 2)
par(xpd = F, family = "Helvetica", pty = "s")

### Biomass
plot(la.rgr[la.rgr$Species == i,]$day_number, la.rgr[la.rgr$Species == i,]$log.la, xlim = c(0, max(la.rgr[la.rgr$Species == i,]$day_number)), xlab = "Days since germination", ylab = "Total leaf area (cm2)", col = "black", axes = F)
axis(1); axis(2, at = log(c(0.001, 0.01, 0.1, 1, 5, 15)), labels = c(0.001, 0.01, 0.1, 1, 5, 15), las =1)
lines(out.fpl$rates$times,   out.fpl$rates$M,   col = "#1B9E77")

##RGRt
plot(Length[[i]]$X,  out.fpl$rates$RGRt, xlab = "Days since germination", ylab = expression(paste("RGR ", (g%.%g^-1%.%day^-1))), type = "n", xlim = c(0, max(la.rgr[la.rgr$Species == i,]$day_number)), ylim = c(0, 0.3))
lines(Length[[i]]$X, out.fpl$rates$RGRt,   col = "#1B9E77")  # 4-param Logistic

# average RGR
  rgr.la.gh[rgr.la.gh$Species == i, 2] <- mean(out.fpl$rates$RGRt)
  rgr.la.gh[rgr.la.gh$Species == i, 3] <- mean(out.fpl$rates$RGRt.hi)
  rgr.la.gh[rgr.la.gh$Species == i, 4] <- mean(out.fpl$rates$RGRt.lo)

}

## F: RGR (Leaf Area) ####

# First look for outliers
la.dim.f.long <- melt(la.dim.f[,c(1:4,10,13,16,19,22)], id.vars = c("Date","Species", "Individual", "leaf"), measure.vars = c("leafarea1","leafarea2","leafarea3","leafarea4","leafarea5"), variable_name = "ID")

la.dim.f.long$Date <- as.factor(la.dim.f.long$Date)
la.dim.f.long$Individual <- as.factor(la.dim.f.long$Individual)

ggplot(filter(la.dim.f.long, Species == "HECO"), aes(x = Date, y = value, col = Individual, by = Individual)) +
  geom_boxplot()

# merge with other leaf area data
la.dim.f <- ddply(la.dim.f, .(Date, Species, Individual, leaf), transform, avg.la = sum(leafarea1, leafarea2, leafarea3, leafarea4, leafarea5)/5)

la.dim.f <- ddply(la.dim.f, .(Date, Species, Individual), summarize, total.la = sum(n.leaves*avg.la))

la.dim.f <- la.dim.f[complete.cases(la.dim.f),]
traits.f <- traits.f[,c(1,2,3,8)]
names(traits.f)[4] <- "total.la"
traits.f <- traits.f[complete.cases(traits.f),]
traits.f <- rbind(traits.f, la.dim.f)

traits.f <- filter(traits.f, Species != "CETR")

ggplot(traits.f, aes(x = Date, y = total.la)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~ Species)

ggplot(traits.f, aes(x = Date, y = log(total.la))) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~ Species) 

traits.f$log.la <- log(traits.f$total.la)

# Prep for loop
traits.f <- ddply(traits.f, .(Species), transform, day_number = as.numeric(difftime(Date, min(Date), units = "days")))

n.preds <- 100   
# before i wasnt adjusting for different number of days between species
Length <- list(PLER = data.frame(X = seq(min(traits.f$day_number), max(traits.f[traits.f$Species == "PLER",]$day_number), length = n.preds)), HECO = data.frame(X = seq(min(traits.f$day_number), max(traits.f[traits.f$Species == "HECO",]$day_number), length = n.preds)), CAPA = data.frame(X = seq(min(traits.f$day_number), max(traits.f[traits.f$Species == "CAPA",]$day_number), length = n.preds)), LACA = data.frame(X = seq(min(traits.f$day_number), max(traits.f[traits.f$Species == "LACA",]$day_number), length = n.preds)), AGHE = data.frame(X = seq(min(traits.f$day_number), max(traits.f[traits.f$Species == "AGHE",]$day_number), length = n.preds)), CLPU = data.frame(X = seq(min(traits.f$day_number), max(traits.f[traits.f$Species == "CLPU",]$day_number), length = n.preds)))

rgr.la.f <- data.frame(Species = SPP, RGRt = rep(NA, 6), RGRt.hi = rep(NA, 6), RGRt.lo = rep(NA, 6))

SPP1 <- c("AGHE", "CLPU", "HECO", "CAPA")

for(i in SPP1){
  
## 4-p Logistic fit ##
#tmp.fpl  <- getInitial(log.la ~ SSfpl(day_number, A, B, xmid, scal), data = traits.f[traits.f$Species == i,])
fit.fpl <- gnls(log.la ~ SSfpl(day_number, A, B, xmid, scal), data = traits.f[traits.f$Species == i,], weights= varExp(form = ~ fitted(.)))
out.fpl <- output.fpl.gnls(fit.fpl, Length[[i]]$X, CI = T, LOG = T)

# Compare all models for each species
par(mfrow = c(1,3))

### Summary  
plot(1:10, 31:40, type = "n", ann = F, axes = F)
text(0, 39, i, pos = 4)
text(0, 38, expression(paste("Model             ", R^2, "    ", "AIC", sep = "")), pos = 4)
text(0,  37, sprintf("4-p Logis  %1.3f %4.1f", out.fpl$summary[1],   out.fpl$summary[2]),   pos = 4, col = "#1B9E77", font = 2)
par(xpd = F, family = "Helvetica", pty = "s")

### Biomass
plot(traits.f[traits.f$Species == i,]$day_number, traits.f[traits.f$Species == i,]$log.la, xlim = c(0, max(traits.f[traits.f$Species == i,]$day_number)), xlab = "Days since germination", ylab = "Total leaf area (cm2)", col = "black", axes = F)
axis(1); axis(2, at = log(c(0.001, 0.01, 0.1, 1, 5, 15)), labels = c(0.001, 0.01, 0.1, 1, 5, 15), las =1)
lines(out.fpl$rates$times,   out.fpl$rates$M,   col = "#1B9E77")

##RGRt
plot(Length[[i]]$X,  out.fpl$rates$RGRt, xlab = "Days since germination", ylab = expression(paste("RGR ", (g%.%g^-1%.%day^-1))), type = "n", xlim = c(0, max(traits.f[traits.f$Species == i,]$day_number)), ylim = c(0, 0.1))
lines(Length[[i]]$X, out.fpl$rates$RGRt,   col = "#1B9E77")  # 4-param Logistic

# average RGR
  rgr.la.f[rgr.la.f$Species == i, 2] <- mean(out.fpl$rates$RGRt)
  rgr.la.f[rgr.la.f$Species == i, 3] <- mean(out.fpl$rates$RGRt.hi)
  rgr.la.f[rgr.la.f$Species == i, 4] <- mean(out.fpl$rates$RGRt.lo)

}

# LACA and PLER have too few data points to fit a 4 parameter logistic, and for some reason its not letting me even fit a 3 parameter logistic, exp will have to work
Init.exp <- function(mCall,LHS,data) {
 	xy    <- sortedXyData(mCall[["X"]],LHS,data)
	r     <- (coef(lm(log(y) ~ x, xy))[2]) # Use the slope from a linear fit to the data as an initial guess for r
	M0    <- min(xy$y)		             # Use the minimum y value as an initial guess for A
 	value <- c(M0, r)
 	names(value) <- mCall[c("M0", "r")]
 	return(value)
	}
fmla.exp <- as.formula("~M0*exp(r*X)")

SS.exp   <- selfStart(fmla.exp, initial = Init.exp, parameters = c("M0", "r"))
output.exp.nls <- function(fit, times, CI = F, alpha = 0.05){
	params <- coef(fit)
	names(params) <- NULL
	M0 <- params[1]; r <- params[2]
	fitted <- fit$m$fitted()
	resid  <- fit$m$resid()
	data <- data.frame(fitted = fitted, resid = resid)
	eq   <- bquote(.(round(M0, 3)) * e^{.(round(r, 3))*t})
	mss <- sum((fitted - mean(fitted))^2)
	rss <- sum(resid^2)
	R2  <- mss/(mss + rss)
	rmse <- sqrt(rss)
	AIC <- AIC(fit)
	summary <- c(R2 = R2, AIC = AIC, RMSE = rmse)
	rates = data.frame(
		times = times,
		M    =   M0*exp(r*times),
		AGR  = r*M0*exp(r*times),
		RGRt = r,
		RGRm = r
		)
	if(CI == T){
		cov  <- summary(fit)$cov
		x    <- data.frame(rmvnorm(n=1000, mean=coef(fit), sigma=cov))
		M <- AGR <- RGRt <- RGRm <- matrix(NA, ncol = length(times), nrow = nrow(x))
		for(i in 1:nrow(x)){
			x.i  <- x[i,]
			M[i,]     <-       x.i$M0*exp(x.i$r*times)
			AGR[i,]   <- x.i$r*x.i$M0*exp(x.i$r*times)
			RGRt[i,] <- AGR[i,]/M[i,]
			RGRm[i,]  <- x.i$r
		}
		CIs <- summarizer(list(M = M, AGR = AGR, RGRt = RGRt, RGRm = RGRm), alpha)
		names(params) <- c("M0", "r")
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = cbind(rates, CIs))
		} else {
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = rates)
		}
	return(out)
}

SPP2 <- c("LACA", "PLER")

for(i in SPP2){
  
## Exp fit
fit.exp  <-  nls(total.la ~ SS.exp(day_number, r, M0), data = traits.f[traits.f$Species == i,]) 
out.exp  <- output.exp.nls(fit.exp, Length[[i]]$X, T)

# Compare all models for each species
par(mfrow = c(1,3))

### Summary  
plot(1:10, 31:40, type = "n", ann = F, axes = F)
text(0.5, 39, i, pos = 4)
text(0.5, 38, expression(paste("Model             ", R^2, "    ", "AIC", sep = "")), pos = 4)
text(0.5, 36, sprintf("Exponential      %1.3f %4.1f", out.exp$summary[1], out.exp$summary[2]), pos = 4, col = "#7570B3")
par(xpd = F, family = "Helvetica", pty = "s")

### Biomass
plot(traits.f[traits.f$Species == i,]$day_number, traits.f[traits.f$Species == i,]$total.la, xlim = c(0, 115), xlab = "Days since germination", ylab = "Total leaf area (cm2)", col = "black", axes = F)
axis(1); axis(2, at = log(c(0.001, 0.01, 0.1, 1, 5, 15)), labels = c(0.001, 0.01, 0.1, 1, 5, 15), las =1)
lines(out.exp$rates$times, out.exp$rates$M, col = "#7570B3")

##RGRt
plot(Length[[i]]$X,  out.fpl$rates$RGRt, xlab = "Days since germination", ylab = expression(paste("RGR ", (g%.%g^-1%.%day^-1))), type = "n", xlim = c(0, 115), ylim = c(0, 0.05))
lines(Length[[i]]$X, out.exp$rates$RGRt, col = "#7570B3")  # Exponential

# average RGR
  rgr.la.f[rgr.la.f$Species == i, 2] <- mean(out.exp$rates$RGRt)
  rgr.la.f[rgr.la.f$Species == i, 3] <- mean(out.exp$rates$RGRt.hi)
  rgr.la.f[rgr.la.f$Species == i, 4] <- mean(out.exp$rates$RGRt.lo)

}

rm(fit.fpl, out.fpl, fit.exp, out.exp, Length, i, n.preds,summarizer, output.exp.nls, output.fpl.gnls, Init.exp, fmla.exp, SS.exp, SPP1, SPP2)

# bind together field and greenhouse ones
rgr.la.f$gh_f <- "F"
rgr.la.gh$gh_f <- "GH"

rgr.la <- rbind(rgr.la.f, rgr.la.gh)

names(rgr.la)[2] <- "RGR.la"
rm(rgr.la.f, rgr.la.gh, la.dim.f, la.dim.f.long, la.rgr, traits.f, traits.gh)

## GH & F: SLA ####
sla$SLA <- sla$leaf.area/sla$mass

# Final SLA
sla <- ddply(sla, .(Species, gh_f), summarize, SLA = mean(SLA))

sla.short <- sla[sla$Species %in% SPP,]

#### GH & F: Isotopes ####
names(cn)[3] <- "gh_f"
cn$gh_f <- revalue(cn$gh_f, c("G" = "GH"))
cn$CN <- cn$C.ug/cn$N.ug

ggplot(cn, aes(x = Species, y = d13C, by = gh_f, col = gh_f)) +
  geom_boxplot()

ggplot(cn, aes(x = Species, y = CN, by = gh_f, col = gh_f)) +
  geom_boxplot()

d13C <- ddply(cn, .(Species, gh_f), summarize, d13C = mean(d13C))
CN <- ddply(cn, .(Species, gh_f), summarize, CN = mean(CN))

d13C$D13C <- (-0.008 - (d13C$d13C * 0.001))/(1 + (d13C$d13C * 0.001))*1000 # convert relative delta to Discrimination

CN <- filter(CN, Species %in% SPP)

d13C.short <- filter(d13C, Species %in% SPP)

#### Final Traits ####
final.trait <- merge(biomass, CN, by = c("Species", "gh_f"))
final.trait <- merge(final.trait, rgr.avg[,c(1,2,5)], by = c("Species", "gh_f"), all.x = T)
final.trait <- merge(final.trait, rgr.la[,c(1,2,5)])
final.trait <- merge(final.trait, d13C.short, by = c("Species", "gh_f"), all.x = T)
final.trait <- merge(final.trait, sla.short, by = c("Species", "gh_f"), all.x = T)

rm(biomass, cn, CN.short, d13C.short, rgr.avg, rgr.la, sla.short, SPP, r.s, seed.mass)

# Wide traits #
final.trait.w <- reshape(final.trait, timevar = "gh_f", idvar = "Species", direction = "wide")

# Full set of 10 species: SLA vs D13C #
all.sp <- merge(sla[sla$gh_f == "F",], d13C, by = c("Species", "gh_f"))
all.sp <- filter(all.sp, Species != "LOWR")

rm(sla, d13C)

### Graphs ####
plots <- "Figures/"

###
# Greenhouse
###
gh1 <- ggplot(filter(final.trait, gh_f == "GH"), aes(x = d13C, y = RGRt)) +
  geom_text(aes(label = Species), hjust = 0.5, vjust = 0.5) +
  theme_bw() +
  xlim(-32.5, -29.5) +
  labs(title = "GH: WUE v RGR (mass)", y = "RGR (mass)", x = "dC13") +
  theme(plot.title = element_text(hjust = 0.5))

gh2 <- ggplot(filter(final.trait, gh_f == "GH"), aes(x = SLA, y = RGRt)) +
  geom_text(aes(label = Species), hjust = 0.5, vjust = 0.5) +
  theme_bw() +
  xlim(50, 375) +
  labs(title = "GH: SLA v RGR (mass)", y = "RGR (mass)", x = "SLA") +
  theme(plot.title = element_text(hjust = 0.5))

gh3 <- ggplot(filter(final.trait, gh_f == "GH"), aes(x = d13C, y = SLA)) +
  geom_text(aes(label = Species), hjust = 0.5, vjust = 0.5) +
  theme_bw() +
  xlim(-32.5, -29.5) +
  labs(title = "GH: WUE v SLA", y = "SLA", x = "dC13") +
  theme(plot.title = element_text(hjust = 0.5))

gh4 <- ggplot(filter(final.trait, gh_f == "GH"), aes(x = d13C, y = RGR.la)) +
  geom_text(aes(label = Species), hjust = 0.5, vjust = 0.5) +
  theme_bw() +
  xlim(-32.5, -29.5) +
  labs(title = "GH: WUE v RGR (leaf area)", y = "RGR.la", x = "dC13") +
  theme(plot.title = element_text(hjust = 0.5))

# comparing RGR.la field to RGR.la GH
gh5 <- ggplot(filter(final.trait, gh_f == "GH"), aes(x = RGR.la, y = RGRt)) +
  geom_text(aes(label = Species), hjust = 0.5, vjust = 0.5) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1) +
  xlim(.035, .11) +
  labs(title = "GH: RGR (mass) v RGR (la)", y = "RGR (biomass)", x = "RGR (leaf area)") +
  theme(plot.title = element_text(hjust = 0.5))

greenhouse <- grid.arrange(gh1, gh4, gh3, gh2, gh5, nrow = 2)
ggsave(greenhouse, filename = "Greenhouse-traits.pdf", path = plots, width = 300, height = 169, units = "mm", dpi = 600)


###
# Comparisons
###

c1 <- ggplot(tmp, aes(x = RGR.la.F, y = RGR.la.GH)) +
  geom_text(aes(label = Species), hjust = 0.5, vjust = 0.5) +
  theme_bw() +
   xlim(.0035, .03) +
  labs(title = "Relative Growth Rate (leaf area)", y = "Greenhouse", x = "Field") +
  theme(plot.title = element_text(hjust = 0.5))

c2 <- ggplot(tmp, aes(x = SLA.F, y = SLA.GH)) +
  geom_text(aes(label = Species), hjust = 0.5, vjust = 0.5) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1) +
  xlim(40, 450) +
  labs(title = "SLA", y = "Greenhouse", x = "Field") +
  theme(plot.title = element_text(hjust = 0.5))

c3 <- ggplot(tmp, aes(x = d13C.F, y = d13C.GH)) +
  geom_text(aes(label = Species), hjust = 0.5, vjust = 0.5) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1) +
  xlim(-32.2,-27.8) +
  labs(title = "dC13 (WUE)", y = "Greenhouse", x = "Field") +
  theme(plot.title = element_text(hjust = 0.5))

compare <- grid.arrange(c1, c2, c3, nrow = 1)
ggsave(compare, filename = "compare-RGR.pdf", path = plots, width = 240, height = 80, units = "mm", dpi = 600)

###
# Field
###

ggplot(all.sp, aes(x = d13C, y = log(SLA))) +
  geom_text(aes(label = Species), hjust = 0.5, vjust = 0.5) +
  theme_classic() +
  xlim(-32, -27.5) +
  labs(title = "Field WUE v Field SLA") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_smooth(method = "lm", se = F, col = "red")


plot.trait.sla13 <- ggplot(all.sp, aes(x = D13C, y = SLA)) +
  theme_classic() +
  geom_smooth(method = "lm", se = F) +
  #geom_text(aes(label = Species), size = 3, hjust = .5, vjust = .5) +
  geom_text(all.sp[all.sp$Species %in% SPP,], aes(label = Species), col = "red", size = 3, hjust = .5, vjust = .5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.line = element_blank(),
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 13)) +
  scale_x_reverse(limits = c(24.7,20.3)) +
  labs(x = "Carbon isotope discrimination (∆, \u2030)", y = expression("SLA ("*cm^"2"*" leaf area / g dry mass)"))

ggsave(plot.trait.sla13, filename = "Figures/ESA-2018/Trait-sla13.pdf", width = 90, height = 80, units = "mm", dpi = 600)


f1 <- ggplot(filter(final.trait, gh_f == "F"), aes(x = d13C, y = SLA)) +
  geom_text(aes(label = Species), hjust = 0.5, vjust = 0.5) +
  theme_bw() +
  xlim(-32, -27.5) +
  labs(title = "Field WUE v Field SLA") +
  theme(plot.title = element_text(hjust = 0.5))

f2 <- ggplot(filter(final.trait, gh_f == "F"), aes(x = d13C, y = RGR.la)) +
  geom_text(aes(label = Species), hjust = 0.5, vjust = 0.5) +
  theme_bw() +
  xlim(-32, -27.5) +
  labs(title = "Field WUE v Field RGR (leaf area)") +
  theme(plot.title = element_text(hjust = 0.5))
  

field <- grid.arrange(f1, f2, nrow = 1)
ggsave(field, filename = "Field-traits.pdf", path = plots, width = 169, height = 80, units = "mm", dpi = 600)

#### Write Tables ####
final.trait$Species <- revalue(final.trait$Species, c("AGHE" = "Agoseris heterophylla", "CLPU" = "Clarkia purpurea", "LACA" = "Lasthenia californica", "PLER" = "Plantago erecta", "HECO" = "Hemizonia congesta", "CAPA" = "Calycadenia pauciflora"))

final.trait.w$Species <- revalue(final.trait.w$Species, c("AGHE" = "Agoseris heterophylla", "CLPU" = "Clarkia purpurea", "LACA" = "Lasthenia californica", "PLER" = "Plantago erecta", "HECO" = "Hemizonia congesta", "CAPA" = "Calycadenia pauciflora"))


write.table(final.trait, "Data/Post-Processing/final-traits.csv", sep = ",", row.names = F)

write.table(final.trait.w, "Data/Post-Processing/final-traits-w.csv", sep = ",", row.names = F)

write.table(all.sp, "Data/Post-Processing/final-sla-13c.csv", sep = ",", row.names = F)

#### PCA Traits ####

# Field traits
trait <- final.trait
trait$gh_f <- revalue(trait$gh_f, c("F" = "Fi"))
trait.pca <- prcomp(trait[trait$gh_f == "Fi", c(9:10)], scale = T)
biplot(trait.pca)
summary(trait.pca)

new <- as.data.frame(trait.pca$x)
new$rows <- row.names(trait.pca$x)
trait$rows <- row.names(trait)
new <- merge(trait[trait$gh_f == "Fi",c(1,12)], new[,c(1,3)], by = "rows")

## Greenhouse traits
trait.pca <- prcomp(trait[trait$gh_f == "GH", c(8,10)], scale = T)
biplot(trait.pca)
summary(trait.pca)

new <- as.data.frame(trait.pca$x)
new$rows <- row.names(trait.pca$x)
trait$rows <- row.names(trait)
new <- merge(trait[trait$gh_f == "GH",c(1,12)], new[,c(1,3)], by = "rows")
