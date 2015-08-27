library(vegan)

#####PERMANOVA--example from vingette
data(dune)
data(dune.env)

adonis(dune ~ Management*A1, data=dune.env, permutations=99)

#####lets try it with benthos data...then probably fail...because i dont care about science
##########################################################  
d <- CompiledBenthosLogandFreq
d$Month <- month(d$Date)
d$Year <- as.factor(d$Year)
d<- join(d, PCRw, by=c("Year","Event","Site"), type="left", match="first")
d$DayNum <- yday(d$Date)

c <- CompiledBenthosLog
levels(c$Status)
selected<-c("Processed")
c <- c[c$Status %in% selected,]
names(c)[names(c)=="Year.1"]<- "Year"
c$Year <- as.factor(c$Year)
c<- join(c, PCRw, by=c("Year","Event","Site"), type="left", match="first")
c$Month <- month(c$Date)
c$DayNum <- yday(c$Date)
##########################################################
#Get your species data
freq <- cast(d, YearSER ~ Family, value='Total.Organisms.sum', sum)

#Figure out which environmental vars you want
names(c)
environ <- c[,c("YearSER","DayNum","Depth.m.standard","Year","Month", "Latitude","Longitude","DieoffYear","Clad.")]

environSERs <- environ$YearSER[is.na(environ$Clad.)]
environ<- environ[!is.na(environ$Clad.),]
freq <- freq[!freq$YearSER %in% environSERs, ]

row.names(environ)<- environ$YearSER
environ$YearSER<-NULL
mode(environ$DieoffYear) <- "integer"

row.names(freq)<- freq$YearSER
freq$YearSER<-NULL

environ$NumMussels <- freq$Dreissenidae

freq <- subset(freq, select = -c(Dreissenidae) )

environ<-data.frame(lapply(environ,as.numeric))

###adonis(freq ~ DayNum*Depth.m.standard*Year*Month*Latitude*Longitude*DieoffYear*Clad.*NumMussels, data=environ, method="bray",permutations=999)
###lat/long may need to be removed bc it makes the permanova blow up

adonis(freq ~ DayNum*Depth.m.standard*Year*Month*
         DieoffYear*Clad.*NumMussels, data=environ, method="bray",permutations=999)


