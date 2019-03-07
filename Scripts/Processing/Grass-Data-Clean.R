#### Grass Cover Data ####
rm(list=ls())

treat <- read.csv("Data/Setup/Marina-Treatment-30.csv")

#### Summed by Plot ####
## Grass 2016 ##
grass.16 <- read.csv("Data/Census-Data/2015-2016/Grass-Cover-2016.csv")
grass.16$Year <- 2016
grass.16$Subplot <- revalue(grass.16$Subplot, c("A" = "Grass"))

## Grass 2017 ##
grass.17.sp <- read.csv("Data/Census-Data/2016-2017/Grass-Cover-2017.csv")
grass.17 <- ddply(grass.17.sp, .(Plot, Subplot), summarize, Cover = sum(Cover))
grass.17$Subplot <- revalue(grass.17$Subplot, c("G" = "Grass", "T" = "Thatch"))
grass.17$Year <- 2017
grass.17$Cover <- ifelse(grass.17$Cover > 100, 100, grass.17$Cover)


grass <- rbind(grass.16, grass.17)
grass <- merge(grass, unique(treat[,c(1,3)]), by = "Plot")

write.table(grass, "Data/Post-Processing/grass-cover.csv", sep = ",", row.names = F)

#### Summed by Species ####
grass.16.sp <- read.csv("Data/Census-Data/2015-2016/grass-sp-cover-2016.csv")
grass.16.sp$Year <- 2016
grass.16.sp <- filter(grass.16.sp, Plot != 87, Plot != 85, Plot != 86) # removed plots 85, 86, and 87 due to burrow damage early in the season

grass.17.sp <- grass.17.sp[,c(1:4)]
grass.17.sp$Year <- 2017
grass.17.sp$Cover <- grass.17.sp$Cover/100
grass.sp <- rbind(grass.16.sp, grass.17.sp)

#### Density of individuals ####
grass.dens <- read.csv("Data/Census-Data/2015-2016/grass-density-2016.csv")
grass.sp <- merge(grass.sp, grass.dens[,c(1,8)], by.x = "Species", by.y = "Grass")
grass.sp$dens.est <- grass.sp$Cover*grass.sp$avg.sub.100p

grass.sp$Species <- revalue(grass.sp$Species, c("AVFA" = "Avena fatua", "BRHO" = "Bromus hordeaceus", "TACA" = "Taeniatherum caput-medusae", "VUMI" = "Vulpia microstachys", "LOMU" = "Lolium multiflorum"))

#### Combining with Seet Set ####
grass.seed <- read.csv("Data/Post-Processing/McL_Grass_Seed-set.csv")[,-c(7,8)]
grass.seed$Subplot <- revalue(grass.seed$Subplot, c("Thatch" = "T", "Grass" = "G"))

grass.sp <- merge(grass.sp, grass.seed, by = c("Year", "Plot", "Subplot", "Species"), all = T)

# fix missing treatments
treat <- read.csv("Data/Setup/Marina-Treatment-30.csv")
grass.sp <- merge(grass.sp[,-8], unique(treat[,c(1,3)]), by = "Plot")

# remove drought plots 2017
grass.sp <- filter(grass.sp, !(Year == 2017 & Treat.Code == "D"))

write.table(grass.sp, "Data/Post-Processing/McL_Grass_density_seed.csv", sep = ",", row.names = F)
