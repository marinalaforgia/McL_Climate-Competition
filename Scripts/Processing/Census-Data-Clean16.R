## Script to clean up demography (germ and mortality) data 2016; output: dem-data-16.csv and flo-data-16.csv ##

rm(list = ls())

#### Load libraries and data ####
library(plyr)
library(dplyr)
library(ggplot2)

d <- read.csv("Data/Census-Data/2015-2016/Total-Census-2016.csv") # demography data
treat <- read.csv("Data/Setup/Marina-Treatment-30.csv") # treatment data
added <- read.csv("Data/Setup/Exp-Seed-Species.csv") # seeds sown
viab <- read.csv("Data/Viability/Viability-Overall.csv") # seed viability data

#### Prep data to clean ####
treat$Subplot <- revalue(treat$Subplot, c("A" = "Grass", "B" = "No Grass"))
d$Subplot <- revalue(d$Subplot, c("A" = "Grass", "B" = "No Grass"))
d$Species <- revalue(d$Species, c("AGHE" = "Agoseris heterophylla", "CLGR" = "Clarkia purpurea", "LACA" = "Lasthenia californica", "PLER" = "Plantago erecta", "HECO" = "Hemizonia congesta"))

# Remove LUBI since this species barely germinated and none survived to flower
d <- filter(d, Species != "LUBI")
d[is.na(d)] <- 0
d$Species <- factor(d$Species)

# removed plots 85, 86, and 87 due to burrow damage early in the season
# grass plots in 82 and 89 because of less than 10& grass cover
d <- filter(d, Plot != 87, Plot != 85, Plot != 86, !(Plot == 82 & Subplot == "Grass"), !(Plot == 89 & Subplot == "Grass"))

###
# Separate flowering and germ/mortality data
###

# flowering
flo <- d[,c(1:4,8:15)]
flo <- flo[!rowSums(flo[,5:12]) == 0,]
flo$Year <- 2016
flo <- flo[,c(13,1:12)]

# mortality and germination
d <- d[,c(1:7,19)]

###
# Remove unused plots, plots with mortality due to insects
###

# get rid of plots I stopped using halfway through the year
d <- filter(d, d$Plot %in% treat$Plot)

# get rid of plots with mortality due to insects/burrows
d$n.died <- d$n.died - d$Deaths.bur.ins
d <- d[,c(1:7)]

###
# Adjust seed added by viability data
###
added$Species <- revalue(added$Species, c("AGHE" = "Agoseris heterophylla", "CLPU" = "Clarkia purpurea", "LACA" = "Lasthenia californica", "PLER" = "Plantago erecta", "HECO" = "Hemizonia congesta", "CAPA" = "Calycadenia pauciflora"))
added <- filter(added, Species %in% unique(d$Species), exp.year == 2016)[,c(2,4)] # avg number of seed added per subplot in 2016
viab <- viab[c(1:5), c(1,9)] # viability data
added <- merge(added, viab, by = "Species")
added$viable <- added$avg.num.per.sub*added$p.viable
added$viable <- as.integer(added$viable)

#### Check Mortality and germ data ####

# Summarize data across dates to get per plot data
d$Date <- as.Date(d$Date, "%m/%d/%y")
d$Month <- as.numeric(format(as.Date(d$Date), "%m"))


d.plot <- ddply(d, c("Plot", "Subplot", "Species"), summarize, germ.tot = sum(n.germ), germ.proj = sum(n.germ) - sum(n.thin), trt.mort = sum(n.died[which(Month %in% c(12, 1, 2))]), pst.mort = sum(n.died[which(Month %in% c(3:9))]), tot.mort = sum(n.died), .drop = F)
d.plot <- merge(d.plot, added[,c(1,4)], by = "Species") 

# With this option, how many plots have negative numbers?
d.plot$total <- d.plot$viable - d.plot$germ.tot

d.plot.neg <- filter(d.plot, total < 0) # 8 HECO... only 2 or 3 really weird ones with 20-40 over, not outside of range of possibility...

# Adjust plots with too many, cap at avg.
d.plot$germ.tot <- ifelse(d.plot$germ.tot > d.plot$viable, d.plot$viable, d.plot$germ.tot)

# check mortality doesnt exceed germ
d.plot$total <- d.plot$germ.proj - d.plot$tot.mort
d.plot.neg <- filter(d.plot, total < 0) # none, we good

#merge with treatment data
d.plot <- merge(d.plot, unique(treat[,c(1,3)], by = "Plot"))

d.plot$Year <- 2016
d.plot <- d.plot[,c(12,11,1,3,2,4:9)]

write.table(d.plot, "Data/Post-Processing/dem-data-16.csv", row.names = F, sep = ",")
write.table(flo, "Data/Census-Data/flo-data-16.csv", row.names = F, sep = ",")

