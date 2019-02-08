#### Grass Cover Data ####
rm(list=ls())

treat <- read.csv("Data/Setup/Marina-Treatment-30.csv")

## Grass 2016 ##
grass.16 <- read.csv("Data/Census-Data/2015-2016/Grass-Cover-2016.csv")
grass.16$Year <- 2016
grass.16$Subplot <- revalue(grass.16$Subplot, c("A" = "Grass"))

## Grass 2017 ##
grass.17 <- read.csv("Data/Census-Data/2016-2017/Grass-Cover-2017.csv")
grass.17 <- ddply(grass.17, .(Plot, Subplot), summarize, Cover = sum(Cover))
grass.17$Subplot <- revalue(grass.17$Subplot, c("G" = "Grass", "T" = "Thatch"))
grass.17$Year <- 2017
grass.17$Cover <- ifelse(grass.17$Cover > 100, 100, grass.17$Cover)


grass <- rbind(grass.16, grass.17)
grass <- merge(grass, unique(treat[,c(1,3)]), by = "Plot")

write.table(grass, "Data/Post-Processing/grass-cover.csv", sep = ",", row.names = F)
