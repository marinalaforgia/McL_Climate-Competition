# Script to update Precip Data each year from Knoxville

rm(list=ls())

#### READ ME ####
# See Read Me file in annual update folder for more information on downloading precip data

# Update working directory to match your computer's file path
setwd("/Users/Marina/Google Drive/02_McLaughlin_80Sites_Organized") 

# Update file name for daily data below with current date 
# Update current year 
currentyr <- 2019

#### Load libraries ####
#install.packages("plyr")
library(plyr)
#install.packages("dplyr")
library(dplyr)

#### Read in data ####
daily <- read.csv("/Users/Marina/Documents/UC-Davis/Projects/McL_Climate-Competition/Data/Climate/Weather-station/knoxville_daily_data120319.csv", header = F, sep = ",", skip = 6, strip.white = T, blank.lines.skip = T, fileEncoding="UTF-8") # update with current date

# there are 6 lines of text above data, check to make sure this isn't cutting off actual data
head(daily)

# column names
colnames(daily) = c("date","year", "day.of.year","day.of.run", "solar.radiation.kWhr.m2", "ave.wind.speed.mps", "avg.wind.direction.deg", "max.wind.gust.mps", "avg.temp.C", "max.avg.temp.C", "min.avg.temp.C", "avg.rel.humid", "max.avg.rel.humid", "min.avg.rel.humid", "total.precip.mm")

#### Check daily data for outliers ####
# Record and remove outliers:
# Three days with over 250 mm of rainfall: 03/18/2004, 2/25/2003, 2/12/2004
daily[daily$total.precip.mm == 99999,]$total.precip.mm <- NA
boxplot(as.numeric(daily$total.precip.mm))
daily <- filter(daily, total.precip.mm < 250) # also gets rid of NAs

#### Monthly summaries ####
daily$date <- as.Date(daily$date, "%m/%d/%Y")
daily$month <- as.numeric(format(daily$date, "%m"))

# Water Year rainfall
month <- ddply(daily, c("year", "month"), summarize, precip = sum(total.precip.mm))

sepdec.annual <- ddply(month, "year", summarize, precip = sum(precip[which(month %in% 9:12)]))
janaug.annual <- ddply(month, "year", summarize, precip = sum(precip[which(month %in% 1:8)]))

k <- data.frame(Year = c(month$year[1]+1):currentyr)
k$ppt.yr <- sepdec.annual$precip[1:nrow(sepdec.annual)-1] + janaug.annual$precip[2:nrow(sepdec.annual)]

# Winter Rainfall
D.annual <- ddply(month, "year", summarize, precip = sum(precip[month == 12]))
JF.annual <- ddply(month, "year", summarize, precip = sum(precip[which(month %in% 1:2)]))

k$ppt.winter <- D.annual$precip[1:nrow(sepdec.annual)-1] + JF.annual$precip[2:nrow(sepdec.annual)]

colMeans(k[-31,])


# Temperature summer
month <- ddply(daily, c("year", "month"), summarize, temp = mean(max.avg.temp.C))

junaug.annual <- ddply(month, "year", summarize, temp = mean(temp[which(month %in% 6:8)]))
colMeans(junaug.annual[-31,], na.rm = T)


# Winter temp
D.annual <- ddply(month, "year", summarize, temp = mean(temp[month == 12]))
JF.annual <- ddply(month, "year", summarize, temp = mean(temp[which(month %in% 1:2)]))

k$win.t <- D.annual$temp[1:nrow(sepdec.annual)-1] + JF.annual$temp[2:nrow(sepdec.annual)]

colMeans(k[-31,])

daily.16.17 <- filter(daily, year == 2015 | year == 2016 | year == 2017)
daily.16.17$date <- as.Date(daily.16.17$date, "%m/%d/%Y")
daily.16.17$month <- as.numeric(format(daily.16.17$date, "%m"))
daily.16.17 <- filter(daily.16.17, month == 9 | month == 10 | month == 11, total.precip.mm > 0)
daily.16.17 <- daily.16.17[,c(1,2,15,16)]
