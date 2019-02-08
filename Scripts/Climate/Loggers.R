# Processing Climate Data for Watering Experiment #
rm(list=ls())

library(stringr)
library(plyr)
library(ggplot2)
library(Rmisc)
library(dplyr)
library(reshape)
library(lme4)
library(lsmeans)
library(multcomp)

#### Logger data ####

# 2015 - 2016
climate16.17 <- read.csv("Data/Climate/Loggers/2015-2016/17-W 9Jun16-0856.csv")[,1:3]
climate16.18 <- read.csv("Data/Climate/Loggers/2015-2016/18-WC 9Jun16-0859.csv")[,1:3]
climate16.29 <- read.csv("Data/Climate/Loggers/2015-2016/29-W 9Jun16-0811.csv")[,1:3]
climate16.30 <- read.csv("Data/Climate/Loggers/2015-2016/30-WC 9Jun16-0817.csv")[,1:3]
climate16.85 <- read.csv("Data/Climate/Loggers/2015-2016/85 9Jun16-0906.csv")[,1:3]
climate16.86 <- read.csv("Data/Climate/Loggers/2015-2016/86 9Jun16-0903.csv")[,1:3]
climate16.97 <- read.csv("Data/Climate/Loggers/2015-2016/97 9Jun16-0823.csv")[,1:3]
climate16.98 <- read.csv("Data/Climate/Loggers/2015-2016/98 9Jun16-0820.csv")[,1:3]

# 2016 - 2017
climate17.17 <- read.csv("Data/Climate/Loggers/2016-2017/17-W 3Nov17-1533.csv")[,1:3]
climate17.18 <- read.csv("Data/Climate/Loggers/2016-2017/18-WC 3Nov17-1320.csv")[,1:3]
climate17.29 <- read.csv("Data/Climate/Loggers/2016-2017/29-W 3Nov17-1226.csv")[,1:3]
climate17.30 <- read.csv("Data/Climate/Loggers/2016-2017/30-WC 3Nov17-1232.csv")[,1:3]
climate17.85 <- read.csv("Data/Climate/Loggers/2016-2017/85 3Nov17-1529.csv")[,1:3]
climate17.86 <- read.csv("Data/Climate/Loggers/2016-2017/86 3Nov17-1526.csv")[,1:3]
climate17.97 <- read.csv("Data/Climate/Loggers/2016-2017/97 3Nov17-1357.csv")[,1:3]
climate17.98 <- read.csv("Data/Climate/Loggers/2016-2017/98 3Nov17-1524.csv")[,1:3]

# Treatment Data
treat <- read.csv("Data/Marina-Treatment-30.csv") # treatment data

# Process logger data
climate16.17 <- climate16.17[3:nrow(climate16.17),]
climate16.18 <- climate16.18[3:nrow(climate16.18),]
climate16.29 <- climate16.29[3:nrow(climate16.29),]
climate16.30 <- climate16.30[3:nrow(climate16.30),]
climate16.85 <- climate16.85[3:nrow(climate16.85),]
climate16.86 <- climate16.86[3:nrow(climate16.86),]
climate16.97 <- climate16.97[3:nrow(climate16.97),]
climate16.98 <- climate16.98[3:nrow(climate16.98),]

colnames(climate16.17) <- c("Time", "Moisture", "Temp")
colnames(climate16.18) <- c("Time", "Moisture", "Temp")
colnames(climate16.29) <- c("Time", "Moisture", "Temp")
colnames(climate16.30) <- c("Time", "Moisture", "Temp")
colnames(climate16.85) <- c("Time", "Moisture", "Temp")
colnames(climate16.86) <- c("Time", "Moisture", "Temp")
colnames(climate16.97) <- c("Time", "Moisture", "Temp")
colnames(climate16.98) <- c("Time", "Moisture", "Temp")

climate16.17$Plot <- 17
climate16.18$Plot <- 18
climate16.29$Plot <- 29
climate16.30$Plot <- 30
climate16.85$Plot <- 85
climate16.86$Plot <- 86
climate16.97$Plot <- 97
climate16.98$Plot <- 98

climate17.17 <- climate17.17[3:nrow(climate17.17),]
climate17.18 <- climate17.18[3:nrow(climate17.18),]
climate17.29 <- climate17.29[3:nrow(climate17.29),]
climate17.30 <- climate17.30[3:nrow(climate17.30),]
climate17.85 <- climate17.85[3:nrow(climate17.85),]
climate17.86 <- climate17.86[3:nrow(climate17.86),]
climate17.97 <- climate17.97[3:nrow(climate17.97),]
climate17.98 <- climate17.98[3:nrow(climate17.98),]

colnames(climate17.17) <- c("Time", "Moisture", "Temp")
colnames(climate17.18) <- c("Time", "Moisture", "Temp")
colnames(climate17.29) <- c("Time", "Moisture", "Temp")
colnames(climate17.30) <- c("Time", "Moisture", "Temp")
colnames(climate17.85) <- c("Time", "Moisture", "Temp")
colnames(climate17.86) <- c("Time", "Moisture", "Temp")
colnames(climate17.97) <- c("Time", "Moisture", "Temp")
colnames(climate17.98) <- c("Time", "Moisture", "Temp")

climate17.17$Plot <- 17
climate17.18$Plot <- 18
climate17.29$Plot <- 29
climate17.30$Plot <- 30
climate17.85$Plot <- 85
climate17.86$Plot <- 86
climate17.97$Plot <- 97
climate17.98$Plot <- 98

climate16 <- rbind(climate16.17, climate16.18, climate16.29, climate16.30, climate16.85, climate16.86, climate16.97, climate16.98)
climate17 <- rbind(climate17.17, climate17.18, climate17.29, climate17.30, climate17.85, climate17.86, climate17.97, climate17.98)
climate16$Year <- 2016
climate17$Year <- 2017
climate <- rbind(climate16, climate17)
rm(climate16.17, climate16.18, climate16.29, climate16.30, climate16.85, climate16.86, climate16.97, climate16.98, climate16.17, climate16.18, climate16.29, climate16.30, climate16.85, climate16.86, climate16.97, climate16.98, climate16, climate17)

test <- as.data.frame(cbind(as.character(climate$Time), str_split_fixed(climate$Time, " ", 2), climate$Plot))
colnames(test) <- c("old", "Date", "Time", "Plot")
colnames(climate)[1] <- "old"
climate <- merge(climate, test, by = c("old", "Plot"))
climate <- climate[,2:ncol(climate)]
rm(test)
climate$Moisture <- as.numeric(as.character(climate$Moisture))
climate$Temp <- as.numeric(as.character(climate$Temp))
climate$Date <- as.Date(climate$Date, "%m/%d/%Y")
climate.d <- ddply(climate, .(Plot, Date), summarize, Moisture = mean(Moisture), Temp = mean(Temp))

# Graph
ggplot(climate.d, aes(x = Date, y = Moisture)) +
  geom_line() +
  facet_wrap(~ Plot)

climate.d <- merge(climate.d, unique(treat[,c(1,3)]), by = "Plot")

climate.d.sum <- summarySE(climate.d, measurevar = "Moisture", groupvars = c("Date", "Treat.Code"))

ggplot(climate.d.sum, aes(x = Date, y = Moisture)) +
  geom_line() +
  facet_wrap(~ Treat.Code) # why are my watered plots drier than my drought plots?

climate.win16 <- filter(climate.d.sum, Date > as.Date("2015-12-01"), Date < as.Date("2016-03-01"))

ggplot(climate.win16, aes(x = Date, y = Moisture)) +
  geom_line() +
  facet_wrap(~ Treat.Code)

climate.win17 <- filter(climate.d.sum, Date > as.Date("2016-12-01"), Date < as.Date("2017-03-01"))

ggplot(climate.win17, aes(x = Date, y = Moisture)) +
  geom_line() +
  facet_wrap(~ Treat.Code)

# conclusion = loggers worthless; treatments = even more worthless
rm(climate.d.sum, climate.d, climate.win16, climate.win17, climate)

#### TDR Data ####

# 2015-2016
tdr <- read.csv("Data/Climate/TDR/soil water data winter 2015-2016.csv")

colnames(tdr) <- c("Plot","Treat.Code","12/11/15","1/2/16","3/2/16","3/2/16.d")

# Deep soil moisture
tdr.d <- tdr[-1, c(1, 2, 6)]

# Shallow soil moisture
tdr<- tdr[-1,1:5]
tdr$Treat.Code <- revalue(tdr$Treat.Code, c("WC" = "C", "SC" = "C", "S" = "D"))

tdr <- melt(tdr, id.vars = c("Plot", "Treat.Code"))
colnames(tdr)[3:4] <- c("Date", "Moisture")
tdr$Moisture <- as.numeric(as.character(tdr$Moisture))
tdr.mine <- filter(tdr, Plot %in% unique(treat$Plot))
summarySE(tdr.mine, measurevar = "Moisture", groupvars = "Treat.Code")

m.moist <- lmer(Moisture ~ Treat.Code + (1|Date), data = tdr.mine)
summary(m.moist)

K <- rbind("D - C" = c(0, -1, 0),
           "W - C" = c(0, -1, 1))

summary(glht(m.moist, linfct = K), test = adjusted("BH")) # Shallow water significant

# Deep soil moisture
tdr.d$Treat.Code <- revalue(tdr.d$Treat.Code, c("WC" = "C", "SC" = "C", "S" = "D"))
colnames(tdr.d)[3] <- "Moisture"
tdr.d$Moisture <- as.numeric(as.character(tdr.d$Moisture))
tdr.mine <- filter(tdr.d, Plot %in% unique(treat$Plot))
summarySE(tdr.mine, measurevar = "Moisture", groupvars = "Treat.Code")

m.moist <- lm(Moisture ~ Treat.Code, data = tdr.mine)
summary(m.moist)

K <- rbind("D - C" = c(0, -1, 0),
           "W - C" = c(0, -1, 1))

summary(glht(m.moist, linfct = K), test = adjusted("BH")) # deep water not sig diff in drought plots, but watered plots sig wetter


