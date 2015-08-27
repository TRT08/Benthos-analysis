library(doBy)
library(reshape)
library("ggplot2")
library(vegan)
library(reshape2)
library(ggplot2)
library(scales)
library(lubridate)
library(extrafont)
library(ggplot2)
library(plyr)
library(outliers)
library(forams)
library(FSA)
library(vegan)
library(stringr)
library(RODBC)

# Get PCR data
db <- "F:/DATA/SLBE/AVBOT Database.accdb"
con2 <- odbcConnectAccess2007(db)
PCR <- sqlFetch(con2, "Updated PCR Results")
close(con2)
rm(con2,db)
colnames(PCR) <- make.names(colnames(PCR), unique = TRUE)

#Restructure PCR data  to rollup with other data
PCRw <- dcast(PCR, YearSER + Site + Year + Event ~ Medium_Abbrev, value.var="Qpcr_Score", mean, na.rm = TRUE)
PCRw[is.na(PCRw )] <- NA
nam <- names(PCRw)[5:10]
names(PCRw)[5:10] <-  str_c(nam, '.PCR_Score')

PCRw$RenameSite[PCRw$Site %in% c("GH2","GHA2")] <- "GH2/GHA2"
PCRw$RenameSite[PCRw$Site %in% c("GH1","GHA1","GHHN")] <- "GH1/GHA1"
PCRw$RenameSite[PCRw$Site %in% c("GHC1","GHD1")] <- "GHC1/GHD1"
PCRw$RenameSite[PCRw$Site %in% c("SMC2","SMD3")] <- "SMC2/SMD3"
PCRw$RenameSite[PCRw$YearSER=="2011-106"] <- "SMA1"
PCRw$RenameSite[PCRw$YearSER=="2012-115"] <- "SMC2"
PCRw$RenameSite2 <- as.factor(ifelse(is.na(PCRw$RenameSite), as.character(PCRw$Site), as.character(PCRw$RenameSite)))
PCRw <- PCRw[ , -which(names(PCRw) %in% c("RenameSite","Site"))]
names(PCRw)[names(PCRw)=="RenameSite2"] <- "Site"
names(PCRw)[names(PCRw)=="YearSER"] <- "PCR_YearSER"

PCRw$Event <- as.factor(PCRw$Event)
PCRw$Year <- as.factor(PCRw$Year)


##########################################################  
d <- CompiledBenthosLogandFreq
d$Month <- month(d$Date)
d$Year <- as.factor(d$Year)
d<- join(d, PCRw, by=c("Year","Event","Site"), type="left", match="first")
d$DayNum <- yday(d$Date)

c <- CompiledBenthosLog
#levels(c$Status)
selected<-c("Processed")
c <- c[c$Status %in% selected,]
c$Year <- as.factor(c$Year)
c<- join(c, PCRw, by=c("Year","Event","Site"), type="left", match="first")
c$Month <- month(c$Date)
c$DayNum <- yday(c$Date)

rm(PCR,PCRw)
#########################################################
##############################################################################
###Get diversity by YEARSER

freq <- cast(d, YearSER ~ Family, value='Total.Organisms.sum', sum)
#freq$Sample <- paste(freq[,1], freq[,2], freq[,3]) 
#row.names(freq)<-freq$Sample
row.names(freq)<-freq$YearSER

freq <- freq[ ,-which(names(freq) %in% c("YearSER"))] ###Choose cols/animals to remove
freq <- freq[, colSums(freq == 0) != nrow(freq)] #Remove any columns that are all zeroes
freq_mat<- data.matrix(freq) ######CONVERT DATAFRAME TO MATRIX!!!!!#######
#shannon diversity
di <- diversity(freq_mat, index = "shannon")
di2<- data.frame(di)
names(di2)[names(di2)=="di"]<- "FamilyShannonDI"
di2$YearSER <- row.names(di2)
c <- join(c, di2, by="YearSER", type="left") ####ADD THIS BACK TO THE COMPILED LOG

rm(di2, di, freq, freq_mat)
###############################################################################
#####REPLACE MISSING WEIGHTS WITH AVERAGES ############

#Get the perc. weight with some estimations...
#First, how many biomass estimations are missing from the benthos data?
#length(which(is.na(d$Av.biomass.mg))) # How many missing
#length(which(!is.na(d$Av.biomass.mg))) # Not missing
#Proportion missing:
#length(which(is.na(d$Av.biomass.mg)))/(length(which(is.na(d$Av.biomass.mg)))+length(which(!is.na(d$Av.biomass.mg))))
#let's try to fill in some of these missing biomasses with benthos data!
ben1<-join(NA.omit.biomass,CompiledBenthosLog,by="YearSER")
avbioben<- summaryBy(biomass.mg ~ Taxon,data=ben1, FUN=c(mean))
avbiobenOrd<- summaryBy(biomass.mg ~ Order,data=ben1, FUN=c(mean))
names(avbioben)[names(avbioben)=="biomass.mg.mean"] <-"biomass.mg.mean.general"
names(avbiobenOrd)[names(avbiobenOrd)=="biomass.mg.mean"] <-"biomass.mg.mean.general"
d<-join(d,avbioben,by=c("Taxon"))

d$est.final.biomass.mg<-ifelse(is.na(d$Sum.biomass.mg), d$biomass.mg.mean.general*d$Total.Organisms.sum, d$Sum.biomass.mg)

#Fill in estimates if none have been measured in benthos yet
#Bytho - Bilkovic and Lehman 1997 biomass for size1-4 Lake Michigan in ug
d$est.final.biomass.mg[d$Taxon =="Bythotrephes"] <- ifelse(
  is.na(d$est.final.biomass.mg[d$Taxon =="Bythotrephes"]),
  mean(c(133.9*0.001,621.6*0.001)), 
  d$est.final.biomass.mg[d$Taxon =="Bythotrephes"])

d$est.final.biomass.mg[d$Order=="Copepoda (sub class)"] <- ifelse(
  is.na(d$est.final.biomass.mg[d$Order=="Copepoda (sub class)"]),
  avbiobenOrd$biomass.mg.mean.general[avbiobenOrd$Order == "Copepoda (sub class)"],
  d$est.final.biomass.mg[d$Order=="Copepoda (sub class)"])

rm(avbiobenOrd, avbioben)
#Now how many are missing?
#length(which(is.na(d$est.final.biomass.mg))) # How many missing
#length(which(!is.na(d$est.final.biomass.mg))) # Not missing
#length(which(is.na(d$est.final.biomass.mg)))/(length(which(is.na(d$est.final.biomass.mg)))+length(which(!is.na(d$est.final.biomass.mg))))

#missingvals<- subset(d, is.na(d$est.final.biomass.mg))

#### Forams Abundance-Biomass Calculation Curves ######
#See Chap 8 Primer-E Manual for explanation

#Need data in this format:
#row.names (Species)  N (abundance)  Biomass (Total biomass)

agg2<- summaryBy(Total.Organisms.sum  + est.final.biomass.mg ~ YearSER + Taxon, data=d, FUN=c(sum),na.rm=TRUE)
agg2<- agg2[complete.cases(agg2),]
agg2<- agg2[-which(agg2$Taxon == "Quagga"), ] ###Remove mussels because very heavy compared to little guys
agg2<- agg2[-which(agg2$Taxon == "Zebra"), ] 

colnames(agg2)[3:4] <- c("N", "Biomass")

cols <- c("YearSER","Taxon")
agg2$combo <- do.call(paste, c(agg2[cols], sep="-"))
row.names(agg2) <- agg2$combo
agg2<- agg2[ , -which(names(agg2) %in% c("combo","Taxon"))]

Wstat <- list()
WstatL95 <- list()
WstatR95 <- list()

for (cat in unique(agg2$YearSER)){
  e <- subset(agg2, agg2$YearSER == cat)
  x <- abc(e[,2:3])
  WstatL95[[cat]] <- x@W.Stat[[1]] #Get the left 95% W stat
  Wstat[[cat]] <- x@W.Stat[[2]] #Get the W stat
  WstatR95[[cat]] <- x@W.Stat[[3]] #Get the right 95% W stat
}

Wstat <- do.call(rbind.data.frame, Wstat)
colnames(Wstat)[1]<- "ABC.W.Stat"
WstatL95 <- do.call(rbind.data.frame, WstatL95)
colnames(WstatL95)[1]<- "ABC.W.stat.L95" #left 95% W stat
WstatR95  <- do.call(rbind.data.frame, WstatR95)
colnames(WstatR95)[1]<- "ABC.W.stat.R95" #right 95% W stat

W.Stats<- do.call(cbind, list(Wstat, WstatL95, WstatR95))

W.Stats$YearSER <- row.names(W.Stats)

c <- join(c, W.Stats, by="YearSER", type="left")

rm(agg2, Wstat,WstatR95,WstatL95, e, W.Stats)

#########Add the critter counts to the compiled log ###########
library(reshape2)
critter.reshape <-dcast(d,  YearSER ~ Order, value.var="Total.Organisms.sum",fun.aggregate=sum ) 
colnames(critter.reshape)[-1] <- paste("Ord", colnames(critter.reshape)[-1], sep = ".")
c <- join(c, critter.reshape, by="YearSER", type="left")

critter.reshapeFAM <-dcast(d,  YearSER ~ Family, value.var="Total.Organisms.sum",fun.aggregate=sum ) 
colnames(critter.reshapeFAM)[-1] <- paste("Fam", colnames(critter.reshapeFAM)[-1], sep = ".")
c <- join(c, critter.reshapeFAM, by="YearSER", type="left")

#######Add richness#########
library(vegan)

c <- c[with(c, order(YearSER)), ]

row.names(critter.reshape)<-critter.reshape$YearSER
critter.reshape <- critter.reshape[,-1]   
c$Order.Richness <- specnumber(critter.reshape)

row.names(critter.reshapeFAM)<-critter.reshapeFAM$YearSER
critter.reshapeFAM <- critter.reshapeFAM[,-1]   
c$Family.Richness <- specnumber(critter.reshapeFAM)

rm(critter.reshape, critter.reshapeFAM)
