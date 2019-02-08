old <- read.csv("Raw.Data/Census.Datasheets/April-3-2016-Census.csv")
head(old)
levels(old$Species)
Plot <- unique(old$Plot)
Plot <- rep(Plot, each = 12)
Subplot <- rep(c("A","B"), each = 6)
Subplot <- rep(Subplot, times = 40)
Species <- rep(levels(old$Species), times = 80)
new <- cbind(Plot,Subplot,Species)
write.csv(new, "datasheet.csv")

# now with new plots
old.2 <- read.csv("Raw.Data/Census.Datasheets/May-8-2016-Census.csv")
head(old.2)
levels(old$Species)
Plot <- unique(old.2$Plot)
Plot <- rep(Plot, each = 12)
Subplot <- rep(c("A","B"), each = 6)
Subplot <- rep(Subplot, times = 30)
Species <- rep(levels(old$Species), times = 60)
new <- cbind(Plot,Subplot,Species)
write.csv(new, "datasheet2.csv")


#d.sum.30 <- filter(d.sum, n.total > 30, Species != "LUBI") # check these numbers
#length(unique(d.sum.30$Plot))

#d.sum.20 <- filter(d.sum, n.total < 15, Species != "LUBI")
#length(unique(d.sum.20$Plot))

#d.sum$n.thin <- ifelse(d.sum$n.total<30, 0, d.sum$n.total - 30)
#d.sum$n.need <- ifelse(d.sum$n.total>20, 0, 20 - d.sum$n.total)

#d.sum.noL <- filter(d.sum, Species != "LUBI")


# create census sheets for 2017 season
library(plyr)
dem <- read.csv("Data/Census Data/2016-2017/Total-Census-2017.csv")
dem[is.na(dem)] <- 0
new <- ddply(dem, .(Plot, Subplot, Species), summarize, N.Germ.in = sum(N.Germ.in), N.Germ.out = sum(N.Germ.out), N.dead.in = sum(N.dead.in), N.dead.out = sum(N.dead.out))
new$tot.in <- new$N.Germ.in - new$N.dead.in
new$tot.out <- new$N.Germ.out - new$N.dead.out
new <- new[,c(1,2,3,8,9)]

write.table(new, "4-Census-Apr9.csv", sep = ",", row.names = F)

dem <- read.csv("Data/Census Data/2016-2017/5-Census-May19.csv")
dem[is.na(dem)] <- 0
new <- ddply(dem, .(Plot, Sub, Species), summarize, t.in = t.in - d.in - N.done, t.out = t.out - d.out)
write.table(new, "6-Census-Jun26.csv", sep = ",", row.names = F)
