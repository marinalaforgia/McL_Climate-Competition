## Script to clean up demography (germ and mortality) data 2017; output: dem-data-17.csv and flo-data-17.csv ##

rm(list = ls())

#### Load libraries and data ####
library(plyr)
library(dplyr)
library(ggplot2)

d <- read.csv("Data/Census-Data/2016-2017/Total-Census-2017.csv") # demography data
treat <- read.csv("Data/Setup/Marina-Treatment-30.csv") # treatment data
added <- read.csv("Data/Setup/Exp-Seed-Species.csv") # seeds sown
viab <- read.csv("Data/Viability/Viability-Overall.csv") # seed viability data
  
#### Prep Data Sets ####

###
# Clean up demography data
###

d$Subplot <- revalue(d$Subplot, c("G" = "Grass", "N" = "No Grass", "T" = "Thatch"))
d$Species <- revalue(d$Species, c("AGHE" = "Agoseris heterophylla", "CLPU" = "Clarkia purpurea", "LACA" = "Lasthenia californica", "PLER" = "Plantago erecta", "HECO" = "Hemizonia congesta", "CAPA" = "Calycadenia pauciflora"))

# Remove CETR since this species barely germinated and none survived to flower
d <- filter(d, Species != "CETR")
d[is.na(d)] <- 0
d$Species <- factor(d$Species)

###
# Separate flowering and mortality data
###

# flowering
flo <- d[,c(1:4,10:18)]
flo <- flo[!rowSums(flo[ ,5:13]) == 0,]
flo$Year <- 2017
flo <- flo[,c(14,1:13)]

# mortality and germination
d <- d[,c(1:9,18,19)]


###
# Adjust seed added by viability data
###
added$Species <- revalue(added$Species, c("AGHE" = "Agoseris heterophylla", "CLPU" = "Clarkia purpurea", "LACA" = "Lasthenia californica", "PLER" = "Plantago erecta", "HECO" = "Hemizonia congesta", "CAPA" = "Calycadenia pauciflora"))
added <- filter(added, Species %in% unique(d$Species), exp.year == 2017)[,c(2,4)] # avg number of seed added per subplot in 2017
viab <- viab[c(1,6:9,12), c(1,9)] 
added <- merge(added, viab, by = "Species")
added$viable <- added$avg.num.per.sub*added$p.viable
added$viable <- as.integer(added$viable)


###
# Check Mortality and germ data for weird things
###

# add together germ.in and germ.out as well as dead.in and dead.out
d$germ.proj <- d$N.Germ.in + d$N.Germ.out # This is for mortality data
d$germ.tot <- d$N.Germ.in + d$N.Thin.in # This is germination data
d$mort <- d$N.dead.in + d$N.dead.out

d <- d[,c(1:4,13,12,14,10,11)]

# Summarize data across dates to get per plot data
d$Date <- as.Date(d$Date, "%m/%d/%y")
d$Month <- as.numeric(format(as.Date(d$Date), "%m"))
#d.plot <- ddply(d, c("Plot", "Subplot", "Species"), summarize, germ.tot = sum(germ.tot), germ.proj = sum(germ.proj), mort = sum(mort), done = sum(N.done), .drop = F)
d.plot <- ddply(d, c("Plot", "Subplot", "Species"), summarize, germ.tot = sum(germ.tot), germ.proj = sum(germ.proj), trt.mort = sum(mort[which(Month %in% c(12, 1, 2))]), pst.mort = sum(mort[which(Month %in% c(3:9))]), tot.mort = sum(mort), .drop = F)

#### Viab Opt 1####
added$viable <- as.integer(added$viable)
d.plot <- merge(d.plot, added[,c(1,4)], by = "Species") 

# With this option, how many plots have negative numbers?
d.plot$total <- d.plot$viable - d.plot$germ.tot
d.plot.neg <- filter(d.plot, total < 0) # AGHE seed so low in mass, makes sense many are over esp in plots with high seed dispersal last year, one weird plantago plot, the rest are in range

# Adjust plots with too many, cap at avg.
d.plot$germ.tot <- ifelse(d.plot$germ.tot > d.plot$viable, d.plot$viable, d.plot$germ.tot)

d.plot <- d.plot[,c(1:9)]
rm(d.plot.neg)

# Check mortality numbers to make sure they dont exceed germination
d.plot$dead.check <- d.plot$germ.proj - d.plot$tot.mort
d.plot.neg <- filter(d.plot, dead.check < 0) # only one had mortality exceed germ and only by one
d.plot$tot.mort <- ifelse(d.plot$dead.check < 0, d.plot$germ.proj, d.plot$tot.mort)
d.plot <- d.plot[,1:9]
rm(d.plot.neg)

#merge with treatment data
d.plot <- merge(d.plot, unique(treat[,c(1,3)], by = "Plot"))

rm(viab, added, treat)

d.plot$Year <- 2017
d.plot <- d.plot[,c(11,10,1,3,2,4:9)]

write.table(d.plot, "Data/Post-Processing/dem-data-17.csv", row.names = F, sep = ",")
write.table(flo, "Data/Census-Data/flo-data-17.csv", row.names = F, sep = ",")

