
####################################
####GETTING STARTED - RUN FIRST ####
####################################

setwd("F:/DATA/SLBE/R scripts/Benthos biomass and frequency/")

require(reshape) || install.packages("reshape") 
require(plyr) || install.packages("plyr") 
require(RODBC) || install.packages("RODBC") 
require(splitstackshape) || install.packages("splitstackshape") 
require(data.table) || install.packages("data.table") 
require(car) || install.packages("car") 
require(doBy) || install.packages("doBy")
require(chron) || install.packages("chron")
require(SDMTools) || install.packages("SDMTools")

options(scipen=999) #Keeps scientific notation away

##########################################
#######GET DATA FROM DATABASE       ######
##########################################

# Get benthos/scope data from Database
db <- "F:/DATA/SLBE/AVBOT Database.accdb"
con2 <- odbcConnectAccess2007(db)
#sqlTables(con2, tableType = "TABLE")$TABLE_NAME   #list table names
RawBenthosLog <- sqlFetch(con2, "Inventory Control Log - Benthos")
RawBenthosData <- sqlFetch(con2, "Compiled Benthic Analysis")
RawBenthosTaxa <- sqlFetch(con2, "Potential Benthos Taxon List")
Rawmussel.lengths <- sqlFetch(con2, "Mussel Lengths")
Rawscopes <- sqlFetch(con2, "Microscope Magnifications")
SerialLogBenthos <- sqlFetch(con2, "Serial Log - Benthos")
Video <- sqlFetch(con2, "Serial Log - Video")
BT <- sqlFetch(con2, "BT Environmental Data")
HOBO <- sqlFetch(con2, "Combined")
Coords<- sqlFetch(con2, "Site Coordinates - Fixed and Random")
Subst <- sqlFetch(con2, "Benthos Substrate Composition")
BenGIS <- sqlFetch(con2, "Benthos Buffer GIS Zonal Stats")
close(con2)

#Make database colnames appropriate for R to use
colnames(RawBenthosLog) <- make.names(colnames(RawBenthosLog), unique = TRUE)
colnames(RawBenthosData) <- make.names(colnames(RawBenthosData), unique = TRUE)
colnames(RawBenthosTaxa) <- make.names(colnames(RawBenthosTaxa), unique = TRUE)
colnames(Rawmussel.lengths) <- make.names(colnames(Rawmussel.lengths), unique = TRUE)
colnames(Rawscopes) <- make.names(colnames(Rawscopes), unique = TRUE)
colnames(SerialLogBenthos) <- make.names(colnames(SerialLogBenthos), unique = TRUE)
colnames(Video) <- make.names(colnames(Video), unique = TRUE)
colnames(BT) <- make.names(colnames(BT), unique = TRUE)
colnames(HOBO) <- make.names(colnames(HOBO), unique = TRUE)
colnames(Coords) <- make.names(colnames(Coords), unique = TRUE)
colnames(Subst) <- make.names(colnames(Subst), unique = TRUE)
colnames(BenGIS) <- make.names(colnames(BenGIS), unique = TRUE)
DownCasts <- BT[which(BT$Cast.Direction=="Down"),]

##########################################
#######   PREPARE THE BENTHOS LOG   ######
##########################################                                              


#Add latitude/longitude, water temp and sample time
#First get the more accurate lat/longs from the benthos log
RawBenthosLog <- join(RawBenthosLog, SerialLogBenthos[ , c("YearSER","Water.Temp","Time","Notes","Lat","Long")], by ="YearSER", type = "left", match = "all")
#Then add in the more general site coordinates
RawBenthosLog <- join(RawBenthosLog, Coords[,c("Site","Year","Latitude","Longitude")], by=c("Site","Year"), type="left", match="first")
#Merge the two lat/longs together, adding more general ones to the more specific ones
RawBenthosLog$LatitudeCor <- ifelse(is.na(RawBenthosLog$Lat), RawBenthosLog$Latitude, RawBenthosLog$Lat)
RawBenthosLog$LongitudeCor <- ifelse(is.na(RawBenthosLog$Long), RawBenthosLog$Longitude, RawBenthosLog$Long)
#Remove unwanted lat long columns
RawBenthosLog <- RawBenthosLog[ , -which(names(RawBenthosLog) %in% c("Lat","Long","Latitude","Longitude"))] 
#Rename the good ones
setnames(RawBenthosLog, old = c('LatitudeCor','LongitudeCor'), new = c('Latitude','Longitude'))

#convert feet to meters and add standard depths if true depths not available!
RawBenthosLog$Depth.m.standard[RawBenthosLog$Site %in% c("GH1","GHA1","GHB1","GHC1","GHD1","GHN1","SBD1","SBN1","SM1","SMA1","SMB1","SMC1","SMD1","SMN1","C2","A1","B1","C1","D1","E1")] <- 10
RawBenthosLog$Depth.m.standard[RawBenthosLog$Site %in% c("GH2","GHA2","GHB2","GHC2","GHD2","GHN2","SBD2","SBN2","SM2","SMA2","SMB2","SMC2","SMD2","SMN2","A2","B2","C3","C3 (569)","D2","E2")] <- 20
RawBenthosLog$Depth.m.standard[RawBenthosLog$Site %in% c("GHD3","GHN3","SBD3","SBN3","SMD3","SMN3")]<- 30
RawBenthosLog$Depth.m <- RawBenthosLog$Depth.ft * 0.3048 #convert feet to meters
RawBenthosLog$Depth.m[is.na(RawBenthosLog$Depth.m)] <- RawBenthosLog$Depth.m.standard[is.na(RawBenthosLog$Depth.m)] #combine the standard depths (10,20,30) with actual depths into one column
RawBenthosLog <- RawBenthosLog[ , -which(names(RawBenthosLog) %in% "Depth.ft")] #remove intermediate step columns

#Add general location names (e.g. South Manitou)
RawBenthosLog$GeneralLoc[RawBenthosLog$Site %in% c("GH","GH1","GH2","GHA1" ,"GHA2","GHB1","GHB2","GHC1","GHC2","GHD1","GHD2","GHD3","GHHD","GHHN","GHN1","GHN2","GHN3","GH-S1","GH-S4","GH-SWLeg")] <- "GoodHarbor"
RawBenthosLog$GeneralLoc[RawBenthosLog$Site %in% c("SM1","SM2","SMA1","SMA2","SMB1","SMB2","SMC1","SMC2","SMD1","SMD2","SMD3","SMN1","SMN2","SMN3","D1","D2","CG","SM","SM-S2","SM-S3")] <- "SouthManitou"
RawBenthosLog$GeneralLoc[RawBenthosLog$Site %in% c("SBD1","SBD2","SBD3","SBN1","SBN2","SBN3","B1","B2")] <- "SleepingBearShoals"
RawBenthosLog$GeneralLoc[RawBenthosLog$Site %in% c("E1","E2")] <- "NorthManitou"
RawBenthosLog$GeneralLoc[RawBenthosLog$Site %in% c("C1","C2","C3","C3 (569)")] <- "PlatteBay"
RawBenthosLog$GeneralLoc[RawBenthosLog$Site %in% c("A1","A2")] <- "PyramidPoint"

#####Die-off year
RawBenthosLog$DieoffYear[RawBenthosLog$Year==2010]<- TRUE
RawBenthosLog$DieoffYear[RawBenthosLog$Year==2011]<- FALSE
RawBenthosLog$DieoffYear[RawBenthosLog$Year==2012]<- TRUE
RawBenthosLog$DieoffYear[RawBenthosLog$Year==2013]<- FALSE

#####Find the difference between the standard depth (10,20,30m) and the actual depth of ponar grab
RawBenthosLog$DepthDifference.m <- abs(RawBenthosLog$Depth.m.standard-RawBenthosLog$Depth.m)

##############ADD VIDEO DATA#############
#Figure out site designation from video log by rolling up video SERs with the closest site and date
VideoGrab<- Video[ , c("YearSER","Site","Date","SedRate","MussRate", "CladRate","Site.Condition.For.Event")]
names(VideoGrab)[names(VideoGrab)=="YearSER"] <- "Video.YearSER"
DT <- data.table(RawBenthosLog, key = c("Site","Date"))
tm <- data.table(VideoGrab, key = key(DT))
test <- tm[DT, roll='nearest', allow.cartesian=TRUE]
test <- test[!duplicated(test$YearSER),] 
RawBenthosLog <- data.frame(test)

##Add the date of the VideoSER to get difference in videograb
Videodateyear <- Video[,c("Date","YearSER")]
names(Videodateyear)[names(Videodateyear)=="Date"] <- "Video.Date"
names(Videodateyear)[names(Videodateyear)=="YearSER"] <- "Video.YearSER"
Videodateyear<- unique(Videodateyear)
test2 <- join(RawBenthosLog, Videodateyear, by ="Video.YearSER", type = "left", match = "all")
withVideo <- data.frame(test2)
withVideo$Video.Time.Diff.Days<- abs(difftime(withVideo$Video.Date, withVideo$Date, units="days"))

##############ADD BT DATA FOR 2013##################
dc <- data.table(DownCasts, key=c("Site","Date"))
dc <- data.table(DownCasts)
maxdepth<- dc[, list(Depth..m.=max(Depth..m.)), by=c("SER","Site","Date")]
maxdepthBT <- join(maxdepth, DownCasts, by =c("SER","Site","Date","Depth..m."), type = "left")
maxdepthBT <- data.frame(maxdepthBT)
names(maxdepthBT)[names(maxdepthBT)=="Depth..m."] <- "Depth.m"

#Roll up BT data with diets using the closest site,date, and depth compared to the average trap/net depth
DT <- data.table(withVideo, key = c("Site","Date"))
tm <- data.table(maxdepthBT, key = key(DT))
test <- tm[DT, roll='nearest', allow.cartesian=TRUE]
test <- test[!duplicated(test$YearSER),]

##Add the date of the BT drop
BTdateyear <- BT[,c("Date","YearSer")]
names(BTdateyear)[names(BTdateyear)=="Date"] <- "BT.Date"
BTdateyear<- unique(BTdateyear)
test2 <- join(test, BTdateyear, by ="YearSer", type = "left", match = "all")
withBT <- data.frame(test2)
names(withBT)[names(withBT)=="YearSer"] <- "BT.YearSer"
names(withBT)[names(withBT)=="Depth.m"] <- "Max.BT.depth.m"
withBT$BT.Time.Diff.Days<- abs(difftime(withBT$BT.Date, withBT$Date, units="days"))

#remove extra columns created by merge with video and BT
withBT <- withBT[ , -which(names(withBT) %in% c("ID1","ID","Year","SER","Pressure..psi.","Flag","Scan.Count","Original.Fluorescence.ECO.AFL.FL..mg.m.3.","Cast.Direction","Elapsed.Time..s."))] 

###GET RID OF PESKY COLUMN NAME CHANGES FROM DATATABLE
colnames(withBT) <- sub("i\\.", "", colnames(withBT))
colnames(withBT) <- sub("\\.1", "", colnames(withBT))

#Add difference in max.bt.depth and the average trap/net depth
withBT$BT.depth.and.ponar.diff.m <- abs(withBT$Max.BT.depth.m-withBT$Depth.m)

#####DEP NONDEP DESIGNATIONS
withBT$DepNonDep[withBT$Site %in% c("GHD1","GHD2","GHD3","GHHD","SMD1","SMD2","SMD3","SBD1","SBD2","SBD3")]<- "DEP"
withBT$DepNonDep[withBT$Site %in% c("GHHN","GHN1","GHN2","GHN3","SMN1","SMN2","SMN3","SBN1","SBN2","SBN3")]<-"NONDEP"
withBT$DepNonDep[withBT$Site.Condition.For.Event=="SLOUGHED" & withBT$Year=="2012"]<-"DEP"
withBT$DepNonDep[withBT$Site.Condition.For.Event %in% c("BARE","LIVE") & withBT$Year=="2012"]<-"NONDEP"

#######HOBO DATA ROLLUP###########
HOBO$Time2 <- as.POSIXlt(HOBO$Time) 
HOBO$Time2 <- times(format(HOBO$Time2, "%H:%M:%S"))
HOBO <- subset(HOBO, select = -Time )
names(HOBO)[names(HOBO)=="Time2"] <- "Time" 
HOBO <-  within(HOBO, { timestamp=format(as.POSIXct(paste(HOBO$Date, HOBO$Time)), "%m/%d/%Y %H:%M:%S") }) #Make datetime

###Remove first and last HOBO readings to account for pick up and drop off times)
HOBO2 <- HOBO[(order(HOBO$YearSER, HOBO$ID)), ]

highest<-by(HOBO2, HOBO2$YearSER, tail, n=1)
lowest<-by(HOBO2, HOBO2$YearSER, head, n=1)

highestd<-do.call("rbind", as.list(highest))
lowestd<-do.call("rbind", as.list(lowest))

HOBO3 <- HOBO2[!HOBO2$ID %in% highestd$ID, ]
HOBO4 <- HOBO3[!HOBO3$ID %in% lowestd$ID, ]

HOBO4 <- HOBO4[ , -which(names(HOBO4) %in% c("ID","Year","Event","SER","Gear","Depth","Time","timestamp"))] 

names(HOBO4)[names(HOBO4)=="Temp_C"] <- "HOBO_Temp_C" 
names(HOBO4)[names(HOBO4)=="Intensity_Lux"] <- "HOBO_Intensity_Lux" 
names(HOBO4)[names(HOBO4)=="YearSER"] <- "HOBO_YearSER" 

###Now get temp and light summaries###

HOBOTempMean <- summaryBy(HOBO_Temp_C ~ Date + Site, data=HOBO4, FUN=mean, na.rm=TRUE)
Lux <- HOBO4[ -which(HOBO4$HOBO_Intensity_Lux ==0),]
Lux <- Lux[!(is.na(Lux$HOBO_Intensity_Lux)),] 
HOBOLuxMean <- summaryBy(HOBO_Intensity_Lux  ~ Date + Site, data=Lux, FUN=mean, na.rm=TRUE)

HOBO5 <- join(HOBOTempMean,HOBOLuxMean, by = "Date", type = 'full')
HOBOdateyear <- HOBO4[,c("Date","HOBO_YearSER","Site")]
HOBOdateyear<- unique(HOBOdateyear)
HOBO5 <- join(HOBO5,HOBOdateyear, by = c("Site","Date"), type = 'full')

#Roll up HOBO data with diets using the closest site,date
DT <- data.table(withBT, key = c("Site","Date"))
tm <- data.table(HOBO5, key = key(DT))
withBT <- tm[DT, roll='nearest', allow.cartesian=TRUE]
withBT <- data.frame(withBT)
withBT <- withBT[!duplicated(withBT$YearSER),]

names(HOBOdateyear)[names(HOBOdateyear)=="Date"] <- "HOBO.Date"
HOBOdateyear <- HOBOdateyear[,c("HOBO.Date","HOBO_YearSER")]

##Add the date of the HOBO
withBT<- join(withBT, HOBOdateyear, by ="HOBO_YearSER", type = "left", match = "first")
withBT$HOBO.Time.Diff.Days<- abs(difftime(withBT$HOBO.Date, withBT$Date, units="days"))

###########################################
###ROLL UP NOWCAST WEATHER DATA###########
#IMPORTANT NOTE: many lat/long points lay outside the model parameters. so instead some locations were changed:
#SMN3 data was used for SMN1 and SMN2
#For GHN1-3: 44.9845,-85.7913
#For SBN1-3: 44.9119, -86.0189
library(lubridate)

# get the nowcast weather data
NC <- read.csv("F:/DATA/SLBE/NowCast Data/CombinedFixedSites.csv")
NC$DateTime <- strptime(as.character(NC$Date.Time.GMT.0500.), "%m/%d/%Y %H:%M")
withBT$Time <- format(withBT$Time, format="%H:%M:%S") 

NAs <- is.na(withBT$Time)
withBT$Time[NAs] <- "12:30:00" #replace all missing times with standard time

withBT$DateTime <- as.POSIXct(paste(withBT$Date, withBT$Time), format="%Y-%m-%d %H:%M:%S")

withBT$mins24DateTime <- withBT$DateTime - as.difftime(24, unit="hours")

tempwithBT <- withBT[withBT$Status=="Processed",]

Now <- list()
for (i in 1:nrow(tempwithBT)){
  print(as.character(tempwithBT$YearSER[i]))
  e <-  NC[NC$Site %in% tempwithBT$Site[i], ]
  e$YearSER <- tempwithBT$YearSER[i]
  int <- new_interval(tempwithBT$mins24DateTime[i], tempwithBT$DateTime[i])
  Now[[i]] <- e[e$DateTime %within% int,]
}

NowW<- do.call(rbind.data.frame, Now)
Nowcast <- NowW[!(is.na(NowW$Date.Time.GMT.0500.)), ]
rm(tempwithBT, NC, NowW,e,int,NAs)

Nowcast[,c(3:26)] <- as.numeric(as.character(unlist(Nowcast[,c(3:26)])))

###Get circular averages
AirVel <- Nowcast[,c("Site","YearSER","Air.Velocity.Direction.Degrees.0.from.North.")]
AirVel <- AirVel[complete.cases(AirVel),]
circ.AirVelocityDir<- tapply(AirVel$Air.Velocity.Direction.Degrees.0.from.North., list(AirVel$YearSER), function(x) circular.averaging(x, deg = TRUE))
circ.AirVelocityDir <- as.data.frame(circ.AirVelocityDir)
circ.AirVelocityDir$YearSER <- rownames(circ.AirVelocityDir) 
colnames(circ.AirVelocityDir)[which(names(circ.AirVelocityDir) == "circ.AirVelocityDir")] <- "circ.AirVelocityDir.mean"
circ.AirVelocityDir$circ.AirVelocityDir.mean<- as.numeric(circ.AirVelocityDir$circ.AirVelocityDir.mean)
circ.AirVelocityDir$YearSER<- as.factor(circ.AirVelocityDir$YearSER)
circ.AirVelocityDir <- circ.AirVelocityDir[!(is.na(circ.AirVelocityDir$circ.AirVelocityDir.mean)), ]

WaterVel <- Nowcast[,c("Site","YearSER","Water.Velocity.at.Surface.Direction.Degrees.0.toward.North.")]
WaterVel <- WaterVel[complete.cases(WaterVel),]
circ.WaterVelocityDir<- tapply(WaterVel$Water.Velocity.at.Surface.Direction.Degrees.0.toward.North., list(WaterVel$YearSER), function(x) circular.averaging(x, deg = TRUE))
circ.WaterVelocityDir <- as.data.frame(circ.WaterVelocityDir)
circ.WaterVelocityDir$YearSER <- rownames(circ.WaterVelocityDir) 
colnames(circ.WaterVelocityDir)[which(names(circ.WaterVelocityDir) == "circ.WaterVelocityDir")] <- "circ.WaterVelocityDir.mean"
circ.WaterVelocityDir$circ.WaterVelocityDir.mean<- as.numeric(circ.WaterVelocityDir$circ.WaterVelocityDir.mean)
circ.WaterVelocityDir$YearSER<- as.factor(circ.WaterVelocityDir$YearSER)
circ.WaterVelocityDir <- circ.WaterVelocityDir[!(is.na(circ.WaterVelocityDir$circ.WaterVelocityDir.mean)), ]

DAvWaterVel <- Nowcast[,c("Site","YearSER","Depth.Averaged.Water.Velocity.Direction.Degrees.0.toward.North.")]
DAvWaterVel <- DAvWaterVel[complete.cases(DAvWaterVel),]
circ.DAvWaterVelocityDir<- tapply(DAvWaterVel$Depth.Averaged.Water.Velocity.Direction.Degrees.0.toward.North., list(DAvWaterVel$YearSER), function(x) circular.averaging(x, deg = TRUE))
circ.DAvWaterVelocityDir <- as.data.frame(circ.DAvWaterVelocityDir)
circ.DAvWaterVelocityDir$YearSER <- rownames(circ.DAvWaterVelocityDir) 
colnames(circ.DAvWaterVelocityDir)[which(names(circ.DAvWaterVelocityDir) == "circ.DAvWaterVelocityDir")] <- "circ.DAvWaterVelocityDir.mean"
circ.DAvWaterVelocityDir$circ.DAvWaterVelocityDir.mean<- as.numeric(circ.DAvWaterVelocityDir$circ.DAvWaterVelocityDir.mean)
circ.DAvWaterVelocityDir$YearSER<- as.factor(circ.DAvWaterVelocityDir$YearSER)
circ.DAvWaterVelocityDir <- circ.DAvWaterVelocityDir[!(is.na(circ.DAvWaterVelocityDir$circ.DAvWaterVelocityDir.mean)), ]

WavesVel <- Nowcast[,c("Site","YearSER","Wave.Direction.Degrees.0.toward.North.")]
WavesVel <- WavesVel[complete.cases(WavesVel),]
circ.WavesVelocityDir<- tapply(WavesVel$Wave.Direction.Degrees.0.toward.North., list(WavesVel$YearSER), function(x) circular.averaging(x, deg = TRUE))
circ.WavesVelocityDir <- as.data.frame(circ.WavesVelocityDir)
circ.WavesVelocityDir$YearSER <- rownames(circ.WavesVelocityDir) 
colnames(circ.WavesVelocityDir)[which(names(circ.WavesVelocityDir) == "circ.WavesVelocityDir")] <- "circ.WavesVelocityDir.mean"
circ.WavesVelocityDir$circ.WavesVelocityDir.mean<- as.numeric(circ.WavesVelocityDir$circ.WavesVelocityDir.mean)
circ.WavesVelocityDir$YearSER<- as.factor(circ.WavesVelocityDir$YearSER)
circ.WavesVelocityDir <- circ.WavesVelocityDir[!(is.na(circ.WavesVelocityDir$circ.WavesVelocityDir.mean)), ]

sumNow<- summaryBy(. ~ YearSER, data=Nowcast, FUN=mean, na.rm=TRUE)

sumNow<- sumNow[ ,-which(names(sumNow) %in% c("Day.of.Year.mean", 
                    "Water.Velocity.at.Surface.Direction.Degrees.0.toward.North..mean",
                    "Depth.Averaged.Water.Velocity.Direction.Degrees.0.toward.North..mean",
                    "Wave.Direction.Degrees.0.toward.North..mean",
                    "Air.Velocity.Direction.Degrees.0.from.North..mean"))]

sumAll<- join_all(list(sumNow, circ.WavesVelocityDir, circ.DAvWaterVelocityDir, circ.WaterVelocityDir, circ.AirVelocityDir), by="YearSER", type="left")
rm(sumNow, circ.WavesVelocityDir, circ.DAvWaterVelocityDir, circ.WaterVelocityDir, circ.AirVelocityDir)

withNowcast<- join(withBT, sumAll, by="YearSER", type="left")

rm(sumAll, withBT, WaterVel,WavesVel,Nowcast,DAvWaterVel, AirVel)


############################
#######RENAME FINAL benthos LOG
#Join with the benthos GIS data
CompiledBenthosLog <- join(withNowcast, BenGIS, by="Site", type="left")

SimpClasses <- read.csv("F:/DATA/SLBE/Arc/Extracted Arc Data/BenthosData/10mSimpclass.csv")
SimpClasses2013 <- read.csv("F:/DATA/SLBE/Arc/Extracted Arc Data/BenthosData/10mSimpclass2013.csv")

SimpClasses <- SimpClasses[ , c("Site","Predom.Simp.Ben.Class")]
SimpClasses2013 <- SimpClasses2013[ , c("YearSER","Predom.Simp.Ben.Class")]
names(SimpClasses2013)[2]<- "Predom.Simp.Ben.Class.2013"

CompiledBenthosLog <- join(CompiledBenthosLog, SimpClasses, by="Site", type="left")
CompiledBenthosLog <- join(CompiledBenthosLog, SimpClasses2013, by="YearSER", type="left")

CompiledBenthosLog$Predom.Simp.Ben.Class2 <-as.factor(ifelse(CompiledBenthosLog$Year== "2013", as.character(CompiledBenthosLog$Predom.Simp.Ben.Class.2013), as.character(CompiledBenthosLog$Predom.Simp.Ben.Class)))
CompiledBenthosLog <- CompiledBenthosLog[ , -which(names(CompiledBenthosLog) %in% c("Predom.Simp.Ben.Class","Predom.Simp.Ben.Class.2013"))]
names(CompiledBenthosLog)[names(CompiledBenthosLog)=="Predom.Simp.Ben.Class2"] <- "Predom.Simp.Ben.Class"

#Join with the substrate classifications

cols<- c("DepNonDep","SiteCondition", "CladMap", "Sand", "Silt","Clay","Shells","LiveMussels",                              
"DOMdecomp","Slimebacteria","Muckmud","Diatoms","Clad","Chara",                                  
"Rocky","GravelPebbles","Black","Badsmells" ,"Bare","DecompAlgae")

colnames(Subst)[colnames(Subst) %in% cols] <- paste("Sub", colnames(Subst)[colnames(Subst) %in% cols], sep = ".")

Subst <- Subst[ , c("YearSER", "Sub.DepNonDep","Sub.SiteCondition", "Sub.CladMap", "Sub.Sand","Sub.Shells","Sub.LiveMussels",                              
                    "Sub.DOMdecomp","Sub.Slimebacteria","Sub.Muckmud","Sub.Diatoms","Sub.Clad","Sub.Chara",                                  
                    "Sub.Rocky","Sub.GravelPebbles","Sub.Black","Sub.Badsmells" ,"Sub.Bare","Sub.DecompAlgae")]

CompiledBenthosLog <- join(CompiledBenthosLog, Subst, by="YearSER", type="left", match="first")

rm(Video, cols, BenGIS, Subst, DT,tm,test,test2,BTdateyear,BT,dc,Videodateyear,maxdepth,maxdepthBT,withVideo,VideoGrab,DownCasts,
  SerialLogBenthos,HOBOdateyear,Lux,HOBOTempMean,HOBOLuxMean,HOBO2,HOBO3,HOBO4,HOBO5,HOBO,highestd,lowestd,withNowcast)

#####RENAME COMMON SITES (same throughout years) BY PROxIMITY
CompiledBenthosLog$RenameSite[CompiledBenthosLog$Site %in% c("GH2","GHA2")] <- "GH2/GHA2"
CompiledBenthosLog$RenameSite[CompiledBenthosLog$Site %in% c("GH1","GHA1","GHHN")] <- "GH1/GHA1"
CompiledBenthosLog$RenameSite[CompiledBenthosLog$Site %in% c("GHC1","GHD1")] <- "GHC1/GHD1"
CompiledBenthosLog$RenameSite[CompiledBenthosLog$Site %in% c("SMC2","SMD3")] <- "SMC2/SMD3"
CompiledBenthosLog$RenameSite2 <- as.factor(ifelse(is.na(CompiledBenthosLog$RenameSite), as.character(CompiledBenthosLog$Site), as.character(CompiledBenthosLog$RenameSite)))
CompiledBenthosLog <- CompiledBenthosLog[ , -which(names(CompiledBenthosLog) %in% c("RenameSite","Site"))]
names(CompiledBenthosLog)[names(CompiledBenthosLog)=="RenameSite2"] <- "Site"

###COMPILE THE COMBINED DEPNONDEPS and SITE CONDITIONS##
CompiledBenthosLog$Sub.DepNonDep <- as.character(CompiledBenthosLog$Sub.DepNonDep)
CompiledBenthosLog$Site.Condition.For.Event <- as.character(CompiledBenthosLog$Site.Condition.For.Event)
CompiledBenthosLog$Sub.SiteCondition <- as.character(CompiledBenthosLog$Sub.SiteCondition)
CompiledBenthosLog$All.DepNonDep <- ifelse(CompiledBenthosLog$Year<=2011, CompiledBenthosLog$Sub.DepNonDep, CompiledBenthosLog$DepNonDep) 
CompiledBenthosLog$All.SiteCondition <- ifelse(CompiledBenthosLog$Year<=2011, CompiledBenthosLog$Sub.SiteCondition, CompiledBenthosLog$Site.Condition.For.Event) 

#Make backups for later use
BenthosLog <- RawBenthosLog 
BenthosData <- RawBenthosData
BenthosTaxa <- RawBenthosTaxa
mussel.lengths <- Rawmussel.lengths
Scopes <- Rawscopes

##########################################
#######START BIOMASS / LENGTH CALCS ######
##########################################

ben <- subset(BenthosData, select = -1) #Remove "ID" column
names(ben)[names(ben)=="Taxon"] <- "ID"  #Prep for joining
ben <- join(ben, BenthosTaxa, by ="ID", type = "left", match = "all")

#Remove taxon with ONLY counts - keep only the Total Organism counts (e.g. remove all chironomid heads if they don't count as total organisms)
ben <- ben[!(is.na(ben$Total.Organisms)),] #Remove NA total orgs
ben <-ben[(ben$Total.Organisms>0),] #Remove total org=0

#Melt diet data so each measured animal has its own row--make sure to include all columns except the "X1","X2",etc.

keep.cols <- names(ben)[!names(ben) %in% c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10")]

meltedben<- melt(ben, id=c(keep.cols))

#Change names to melt with microscope magnifiers
names(meltedben)[names(meltedben)=="variable"] <- "Animal.measured"
names(meltedben)[names(meltedben)=="value"] <- "Length"

#Merge scopes multiplier and data
ben2 <- join(meltedben, Scopes, by = c("Scope.Name","Mag"), type = "left", match = "all")
ben2 <- ben2[ , -which(names(ben2) %in% c("ID"))]

#Find true length of organisms (in mm)
ben2$Length <- as.numeric(ben2$Length, na.rm = TRUE)
ben2$mm.length <- ben2$Length * ben2$Multiplier

##############MUSSEL DATA##########################
#Fix up mussel length data to reshape and bind
#mussel.lengths <- subset(mussel.lengths, select = -1)
names(mussel.lengths)[names(mussel.lengths)=="Date.Processed"] <- "Processed.Date..End."
names(mussel.lengths)[names(mussel.lengths)=="Size.Class"] <- "mm.length"

#Reshape mussels to bind with the benthos data
#Function:
expandrows <- function(dataset, count, count.is.col = TRUE) {
  if (!isTRUE(count.is.col)) {
    if (length(count) == 1) {
      dataset[rep(rownames(dataset), each = count), ]
    } else {
      if (length(count) != nrow(dataset)) {
        stop("Expand vector does not match number of rows in data.frame")
      }
      dataset[rep(rownames(dataset), count), ]
    }
  } else {
    dataset[rep(rownames(dataset), dataset[[count]]), 
            splitstackshape:::othernames(dataset, count)]
  }
}

muss <- expandrows(mussel.lengths, "Count", count.is.col = TRUE)
row.names(muss) <- seq(nrow(muss)) #Remove weird row names created by previous script

#Remove bad columns
muss <- muss[ , -which(names(muss) %in% c("Year",	"SER",	"Date",	"Processed.Date..End.",	"Site", "ID"))]
  
#merge the mussels with the serial info
muss<- join(muss, BenthosLog, by ="YearSER", type = "left", match = "all")
muss <- muss[ , -which(names(muss) %in% c("Sample.Type","Status","Location","Specific.Location","Storage.Notes","Priority"))]

#Merge mussels with BenthosTaxa list
muss$ID <- rep(NA, nrow(muss)) #make a blank column to fill with ID numbers
muss[muss$Taxon == "Quagga", ][, "ID"] <- 47
muss[muss$Taxon == "Zebra", ][, "ID"] <- 53
muss <- join(muss, BenthosTaxa, by ="ID", type = "left", match = "all")

nms <- colnames(ben2)  
Missing <- setdiff(nms, names(muss))  # Find names of missing columns
muss[Missing] <- "NA"                    # Add them, filled with '0's
muss <- muss[nms]  

#Bind the mussel lengths and the benthic data!!!!!!
ben.all.L <- rbind(muss, ben2)

######Biomass determination by animal
##Isolate complete cases
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}
NA.omit.biomass <- completeFun(ben.all.L, c("a", "mm.length"))

names(NA.omit.biomass)[names(NA.omit.biomass)=="a"] <- "aX"
names(NA.omit.biomass)[names(NA.omit.biomass)=="b"] <- "bX"
names(NA.omit.biomass)[names(NA.omit.biomass)=="c"] <- "cX"

##GET BIOMASSES#####
#Note: make sure all mm.length = mm and all mass=mg!
#Note2: log in r = ln, and log10 in r=regular log

NA.omit.biomass$biomass.mg<- 
ifelse(NA.omit.biomass$Mass.Unit=="mg", with(NA.omit.biomass, (aX*(mm.length^bX)*factor)),
ifelse(NA.omit.biomass$Mass.Unit=="ug", with(NA.omit.biomass, (aX*(mm.length^bX)*factor)/1000), 
ifelse(NA.omit.biomass$Mass.Unit=="ln_mg", with(NA.omit.biomass, (exp(aX+log(mm.length)*bX))),
ifelse(NA.omit.biomass$Mass.Unit=="ln_ug", with(NA.omit.biomass, (exp(aX+log(mm.length)*bX))/1000), 
ifelse(NA.omit.biomass$Mass.Unit=="log_mg", with(NA.omit.biomass, 10^(aX+(bX*(log10(mm.length))))),
ifelse(NA.omit.biomass$Mass.Unit=="ng", with(NA.omit.biomass, ((aX-(bX*(mm.length*1000)) + cX*(mm.length*1000)^2)/1000000)), 
ifelse(NA.omit.biomass$Mass.Unit=="Cerco_mg", with(NA.omit.biomass, exp(aX*log(mm.length)-bX)), 
ifelse(NA.omit.biomass$Mass.Unit=="Diff_ug", with(NA.omit.biomass, (aX+mm.length^bX)/1000), 
ifelse(NA.omit.biomass$Mass.Unit=="Alon_ug", with(NA.omit.biomass, (aX*(mm.length/1000)^bX)/1000),  
ifelse(NA.omit.biomass$Mass.Unit=="Illyo_g", with(NA.omit.biomass, ((aX*mm.length)-bX)/1000), 
ifelse(NA.omit.biomass$Mass.Unit=="Dia_ug", with(NA.omit.biomass, (aX*(cX*mm.length)^bX)/1000),
ifelse(NA.omit.biomass$Mass.Unit=="Actual_ug", with(NA.omit.biomass,  aX/1000), 
ifelse(NA.omit.biomass$Mass.Unit=="Limno_ug", with(NA.omit.biomass, (10^(aX*mm.length-bX))/1000), 
ifelse(NA.omit.biomass$Mass.Unit=="Ostra_ug", with(NA.omit.biomass, (10^(aX*log10(mm.length)))/1000),       
ifelse(NA.omit.biomass$Mass.Unit=="log_ug", with(NA.omit.biomass, (10^(aX+(bX*(log10(mm.length)))))/1000),NA)))))))))))))))

#####Renaming Freq/biomass critters with heads
NA.omit.biomass$Taxon[which(NA.omit.biomass$Taxon =="Chironomid pupa HW")] <- "Chironomid pupa TL"
NA.omit.biomass$Taxon[which(NA.omit.biomass$Taxon =="Chironomidae HW")] <- "Chironomidae TL"
NA.omit.biomass$Taxon[which(NA.omit.biomass$Taxon =="Chironominae HW")] <- "Chironominae TL"
NA.omit.biomass$Taxon[which(NA.omit.biomass$Taxon =="Tanypodinae HW")] <- "Tanypodinae TL"
NA.omit.biomass$Taxon[which(NA.omit.biomass$Taxon =="Orthocladiinae HW")] <- "Orthocladiinae TL" 
NA.omit.biomass$Taxon[which(NA.omit.biomass$Taxon =="Chironomidae HW - Head Only")] <- "Chironomidae TL" 
NA.omit.biomass$Taxon[which(NA.omit.biomass$Taxon =="Chironominae HW - Head Only")] <- "Chironominae TL" 
NA.omit.biomass$Taxon[which(NA.omit.biomass$Taxon =="Tanypodinae HW - Head Only")] <- "Tanypodinae TL" 
NA.omit.biomass$Taxon[which(NA.omit.biomass$Taxon =="Orthocladiinae HW - Head Only")] <- "Orthocladiinae TL" 
levels(NA.omit.biomass$Taxon)[levels(NA.omit.biomass$Taxon) %in%  c("Lirceus BW", "Lirceus HW")] <- "Lirceus"
levels(NA.omit.biomass$Taxon)[levels(NA.omit.biomass$Taxon) %in%  c("Caecidotea BW", "Caecidotea HW")] <- "Caecidotea"
levels(NA.omit.biomass$Taxon)[levels(NA.omit.biomass$Taxon) %in%  c("Bythotrephes Body", "Bythotrephes Body + Spine", "Bythotrephes Spine")] <- "Bythotrephes"

NA.omit.biomass$Taxon <- droplevels(NA.omit.biomass$Taxon)

####Get averages of biomass for each taxon per SER
mean.biomass.mg <- cast(NA.omit.biomass,  Taxon ~ YearSER, mean, value='biomass.mg')
mean.biomass.mg  <- melt(mean.biomass.mg, id=c("YearSER"))
row.names(mean.biomass.mg) <- seq(nrow(mean.biomass.mg)) 
names(mean.biomass.mg)[names(mean.biomass.mg)=="value"] <- "Av.biomass.mg"

####Get number of animals measured for each taxon per SER
length.biomass.mg <- cast(NA.omit.biomass,  Taxon ~ YearSER, length, value='biomass.mg')
length.biomass.mg  <- melt(length.biomass.mg, id=c("YearSER"))
row.names(length.biomass.mg) <- seq(nrow(length.biomass.mg)) 
names(length.biomass.mg)[names(length.biomass.mg)=="value"] <- "Num.critters.averaged.4.biomass"
length.biomass.mg[length.biomass.mg == 0] <- NA

####Get biomass sum for each taxon per SER
sum.biomass.mg <- cast(NA.omit.biomass,  Taxon ~ YearSER, sum, value='biomass.mg')
sum.biomass.mg  <- melt(sum.biomass.mg, id=c("YearSER"))
row.names(sum.biomass.mg) <- seq(nrow(sum.biomass.mg)) 
names(sum.biomass.mg)[names(sum.biomass.mg)=="value"] <- "Sum.biomass.mg"
sum.biomass.mg[sum.biomass.mg == 0] <- NA


##########################################
#######     FREQUENCY DATA         ######
##########################################
BenthosFreq1 <- subset(BenthosData, select = -1)
names(BenthosFreq1)[names(BenthosFreq1)=="Taxon"] <- "ID"  #Prep for joining
BenthosFreq1 <- join(BenthosFreq1, BenthosTaxa, by ="ID", type = "left", match = "all")

#Remove unneeded columns
Freq.All.benthos <- BenthosFreq1[ , -which(names(BenthosFreq1) %in% c("X1","X2","X3","X4", "X5","X6","X7","X8","X9","X10","Scope.Name", "Mag",  "ID","a","b","c","factor","Mass.Unit","Source","Equation" ))]

#####Renaming Freq/biomass critters with heads
Freq.All.benthos$Taxon[which(Freq.All.benthos$Taxon =="Chironomid pupa HW")] <- "Chironomid pupa TL"
Freq.All.benthos$Taxon[which(Freq.All.benthos$Taxon =="Chironomidae HW")] <- "Chironomidae TL"
Freq.All.benthos$Taxon[which(Freq.All.benthos$Taxon =="Chironominae HW")] <- "Chironominae TL"
Freq.All.benthos$Taxon[which(Freq.All.benthos$Taxon =="Tanypodinae HW")] <- "Tanypodinae TL"
Freq.All.benthos$Taxon[which(Freq.All.benthos$Taxon =="Orthocladiinae HW")] <- "Orthocladiinae TL" 
Freq.All.benthos$Taxon[which(Freq.All.benthos$Taxon =="Chironomidae HW - Head Only")] <- "Chironomidae TL" 
Freq.All.benthos$Taxon[which(Freq.All.benthos$Taxon =="Chironominae HW - Head Only")] <- "Chironominae TL" 
Freq.All.benthos$Taxon[which(Freq.All.benthos$Taxon =="Tanypodinae HW - Head Only")] <- "Tanypodinae TL" 
Freq.All.benthos$Taxon[which(Freq.All.benthos$Taxon =="Orthocladiinae HW - Head Only")] <- "Orthocladiinae TL" 
levels(Freq.All.benthos$Taxon)[levels(Freq.All.benthos$Taxon) %in%  c("Lirceus BW", "Lirceus HW")] <- "Lirceus"
levels(Freq.All.benthos$Taxon)[levels(Freq.All.benthos$Taxon) %in%  c("Caecidotea BW", "Caecidotea HW")] <- "Caecidotea"
levels(Freq.All.benthos$Taxon)[levels(Freq.All.benthos$Taxon) %in%  c("Bythotrephes Body", "Bythotrephes Body + Spine", "Bythotrephes Spine")] <- "Bythotrephes"

Freq.All.benthos$Taxon <- droplevels(Freq.All.benthos$Taxon)

Freq.sum <- summaryBy(Total.Organisms + Count ~ YearSER + Taxon, data=Freq.All.benthos, FUN=sum, na.rm=TRUE)

Dudeinfo <- Freq.All.benthos[ ,-which(names(Freq.All.benthos) %in% c("Total.Organisms","Count"))]

Freq.sum.plus.info <- join(Freq.sum, Dudeinfo, by =c("YearSER","Taxon"), type = "left", match="first")

#Bind the frequency and the biomass data
combo<- join(length.biomass.mg, mean.biomass.mg, by =c("YearSER","Taxon"), type = "left", match = "all")
combo<- join(combo, sum.biomass.mg, by =c("YearSER","Taxon"), type = "left", match = "all")
Freq.All.benthos<- join(Freq.sum.plus.info, combo, by =c("YearSER","Taxon"), type="left")

#Add column of proportion of total critters used to get biomass estimates
Freq.All.benthos$Prpn.biomass.calc <- Freq.All.benthos$Num.critters.averaged.4.biomass/Freq.All.benthos$Total.Organisms.sum

#Make a total biomass column 
Freq.All.benthos$Final.biomass.mg <- Freq.All.benthos$Av.biomass.mg*Freq.All.benthos$Total.Organisms.sum
Freq.All.benthos[is.na(Freq.All.benthos)] <- NA #Turn NANs into NAs

#remove intermediate data frames
rm(combo,length.biomass.mg,mean.biomass.mg, sum.biomass.mg, Scopes,Freq.sum,Freq.sum.plus.info,Rawscopes,
   BenthosData,BenthosFreq1,BenthosLog,BenthosTaxa,Dudeinfo,RawBenthosData,RawBenthosLog,RawBenthosTaxa,
   Rawmussel.lengths,ben,ben2,meltedben,muss,mussel.lengths)

##########################################
####### COMPILE Benthos LOG WITH taxon ######
##########################################

CompiledBenthosLogandFreq<- join(Freq.All.benthos, CompiledBenthosLog, by =c("YearSER"), type="left")

#write.csv(CompiledBenthosLogandFreq,"Benthos.csv")

rm(Missing, Now, con2,db, highest, i, lowest,nms,completeFun,expandrows)