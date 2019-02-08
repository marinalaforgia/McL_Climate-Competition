library(dplyr)
s <- read.csv("Anu_Cover2015.csv")
a <- read.csv("/Users/Marina/Documents/UC-Davis/02_McLaughlin_80Sites_Organized/Modified_CombinedFiles/McL_80SitesSpeciesTraits_012615.csv")
s.2 <- reshape(s, varying = 2:100, v.names = "Cover", timevar = "Species", times=names(s)[2:100], direction = "long")

s.2 <- filter(s.2, Cover != 0)

s.2 <- s.2[,1:3]
s.3 <- merge(s.2, a[,c(1,16)], by.x = "Species", by.y = "Species_Name")

s.2$Species[gsub("[.]", " ", s.2$Species) 
write.table(s.2, "Anu_Species2015.csv", sep = ",", row.names = F)
