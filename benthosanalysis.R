library(doBy)
library(reshape)
library(ggplot2)
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

##########################################################  
#  RUN THESE FIRST TO GET THE NECESSARY FILES FOR THIS SCRIPT TO RUN ########
source("F:/DATA/SLBE/R scripts/General-Scripts/GobyThemes.R") #Load GGPLOT themes
source("biomassfreqscript.R",echo=TRUE, max.deparse.length=10000)
source("biomass_additional_functions.R")
setwd("F:/DATA/SLBE/Manuscripts/Benthos-analysis/") #set directory 

#d = CompiledBenthosLogandFreq
#c = CompiledBenthosLog

###Remove all pelagic animals from the benthos######
#table(d$Habitat)
d <- d[d$Habitat %in% c("Benthic", "Epibenthic", "Migrator"),]

###Convert ponar area to meters squared
#Ponar width: 9inx9in, 22.86cmx22.86c
#ponar area: 522.5796 cm2
#meter area: 10000 cm2

ponar <- function(x){
  (x * 10000)/522.5796
}

##########################################################  
#####Check benthos for outliers#####
outties<-NA.omit.biomass[,c("Taxon","YearSER","mm.length")]
names(outties)[names(outties)=="Length"]<- "Microscope.Length"
source("FindOutliers.R",echo=TRUE, max.deparse.length=10000)
outord["checked"] <- NA
#write.csv(outord, "outord.csv") #run this the first time

#Add a column in excel called "checked" and add a 1 if checked and okay.
benchecked<-read.csv("outord.csv")
benout<- join(outord,benchecked[,c("Taxon","YearSER","checked")], by=c("YearSER","Taxon"),type="left")
benout <- benout[!benout$checked==1, ]

##Check benthos for incorrect measurements
NA.omit.biomass$Mag <- as.numeric(NA.omit.biomass$Mag)
DiffMags <-  NA.omit.biomass[!is.na(NA.omit.biomass$Mag),]
DiffMags <- DiffMags[!DiffMags$Measurement.Standard == DiffMags$Mag,]
DiffMags <- DiffMags[!(DiffMags$Measurement.Standard == 3 & DiffMags$Mag == 3.2),] #leica equiv.
DiffMags <- DiffMags[!(DiffMags$Measurement.Standard == 3 & DiffMags$Mag == 2.5),] #unitron
DiffMags <- DiffMags[!(DiffMags$Measurement.Standard == 3 & DiffMags$Mag == 25),] #margi's scope
DiffMags <- DiffMags[!(DiffMags$Measurement.Standard == 3 & DiffMags$Mag > 3),] #measured higher than 3 for more precision

#now compare this to the benthos outliers, 
#this will tell if there are any issues with possible wrong magnifications
wrongMags <- join(benout, DiffMags, by=c("Taxon", "YearSER"), type="inner")

####################################################
#Which samples should we double check the ID?

table(CompiledBenthosLogandFreq$Taxon)

#Double-check these yearSERs
DoubleCheckID <- CompiledBenthosLogandFreq[,c(1:2,5)][CompiledBenthosLogandFreq$Taxon %in% c("Amphipod","Calanoid","Gammarid sp",
            "Chydorid sp", "Copepod", "Cyclopoid","Daphnia galeata mendotae","Gastropod TL",
            "Water Mite TL"),]
write.csv(DoubleCheckID, "C:/Users/trtucker/Desktop/DoubleCheckID.csv")

#ones to sort chironomids out of for pat
ChiroSERS <- unique(CompiledBenthosLogandFreq[,1][CompiledBenthosLogandFreq$Taxon %in% c("Chironomid pupa TL",
              "Chironomidae TL","Chironominae TL", "Orthocladiinae TL","Tanypodinae TL")])
write.csv(ChiroSERS, "C:/Users/trtucker/Desktop/ChiroSERS.csv")

#Priority ID chironomids to family
ChiroFamilySERS <- CompiledBenthosLogandFreq[,c(1:2,5)][CompiledBenthosLogandFreq$Taxon %in% c("Chironomidae TL"),]
write.csv(ChiroFamilySERS, "C:/Users/trtucker/Desktop/ChiroFamilySERS.csv")

rm(ChiroFamilySERS, ChiroSERS, DoubleCheckID )

#################################################################
#Make basic density / biomass per square meter table#
Density <- cast(d, TaxonomicGroup + Family ~ Year, value='Total.Organisms.sum', sum)
Density[3:6] <- ponar(Density[3:6]) #convert to sq. m.
Biomass <- cast(d, TaxonomicGroup + Family ~ Year, value='est.final.biomass.mg', sum)
Biomass[3:6] <- ponar(Biomass[3:6]) #convert to sq. m.
cols <- c("2010", "2011", "2012", "2013")
colnames(Biomass)[colnames(Biomass) %in% cols] <- paste("Mass", colnames(Biomass)[colnames(Biomass) %in% cols], sep = ".")
DensTable <- join(Density, Biomass, by=c("TaxonomicGroup", "Family"), type="left")
DensTable[,3:10]<-round(DensTable[,3:10], 1)
rm(Density, Biomass, cols)

#################################################################
#make nice graphs with error bars, ala Nalepa by year
#Mean number of orgs v. month v. site type v. taxon

source("F:/DATA/SLBE/R scripts/General-Scripts/summary.R") #gets mean, se, sd, ci

summ<- summarySE(d, measurevar="Total.Organisms.sum", groupvars=c("TaxonomicGroup", "Month","All.SiteCondition","Depth.m.standard"))  #"All.DepNonDep"
#summ[is.na(summ)] <- 0
summ$Total.Organisms.sum <- ponar(summ$Total.Organisms.sum) #convert to sq. m.
summ$Year <- as.numeric(as.character(summ$Year))
summ$Event <- as.numeric(as.character(summ$Event))
summ <- summ[summ$TaxonomicGroup %in% names(which(table(summ$TaxonomicGroup) > 3)), ] #remove ones with few obs.
summ$Depth.m.standard <- as.factor(summ$Depth.m.standard)
summ <- summ[!summ$Depth.m.standard == 30,]

ggplot(summ, aes(x=Month, y=Total.Organisms.sum, group=Depth.m.standard, color=Depth.m.standard)) + 
  geom_errorbar(aes(ymin=Total.Organisms.sum-se, ymax=Total.Organisms.sum+se), width=.5) + geom_point() + Goby_theme +
  geom_line(linetype="dotted") + facet_grid(TaxonomicGroup ~ All.SiteCondition + ., scales="free")+
  labs(x="Month", y="Mean number of organisms", title="Benthos over time") +  
  theme(strip.text.y = element_text(size = 9, angle = 0)) 

ggsave(filename = "F:/DATA/SLBE/Manuscripts/Benthos-analysis/Figs/NumOverTime.png")                                                                        "cm", "mm"), dpi = 300, ...)

######################################################
#####Rarefaction curves#####
###Are there enough samples for each year? 
freq <- cast(d, Year+ YearSER ~ Family, value='Total.Organisms.sum', sum)
row.names(freq)<-freq$YearSER
freq <- freq[ ,-which(names(freq) %in% c("YearSER"))] ###Choose cols/animals to remove
freq <- freq[, colSums(freq == 0) != nrow(freq)] 

opar <- par() 
op <- par(mfrow=c(2,2),col.lab="black",col.main="black")

for (i in levels(freq$Year)){
  freqYEAR <- freq[which(freq$Year==i), ]
  freqYEAR <- freqYEAR[,-1]
  freqYEAR<- data.matrix(freqYEAR)
  rareYEAR<-specaccum(freqYEAR, method="rarefaction")
  plot(rareYEAR, xlab="Number of benthos samples", ylab="Number of taxa", 
     main=paste(i), ci.type = "polygon", ci.col="Gray", ci.lty=0)
}

par(opar)

#######################################################################
####Bar Plot taxa over Time by Year####
Grab<- d[ , c("Month","Order","Year","Total.Organisms.sum","Final.biomass.mg","Family", "GeneralLoc","Depth.m.standard","All.DepNonDep")]

#By year, month
agg <- aggregate(data=Grab, Total.Organisms.sum~ Year+Month+Family, function(x) sum(x))
ggplot(agg, aes(x=Month, y=Total.Organisms.sum, fill=Family)) +
  geom_bar(stat="identity", position = "fill") +
  Goby_theme + scale_x_continuous(breaks=pretty_breaks(n=10))+
  labs(x="Month", y="Proportion of sample taxa", title="Benthos over time")+
  facet_grid(Year ~ .)

#By depth, year, month
agg <- aggregate(data=Grab, Total.Organisms.sum~ Year+Month+Order+Depth.m.standard, function(x) sum(x))
ggplot(agg, aes(x=Month, y=Total.Organisms.sum, fill=Order)) +
  geom_bar(stat="identity", position = "fill") +
  Goby_theme + scale_x_continuous(breaks=pretty_breaks(n=10))+
  labs(x="Month", y="Proportion of sample taxa", title="Benthos over time")+
  facet_grid(Year ~ Depth.m.standard + .)

#By All.DepNonDep, year, month
agg <- aggregate(data=Grab, Total.Organisms.sum~ Year+Month+Order+All.DepNonDep, function(x) sum(x))
ggplot(agg, aes(x=Month, y=Total.Organisms.sum, fill=Order)) +
  geom_bar(stat="identity", position = "fill") +
  Goby_theme + scale_x_continuous(breaks=pretty_breaks(n=10))+
  labs(x="Month", y="Proportion of sample taxa", title="Benthos over time")+
  facet_grid(Year ~ All.DepNonDep + .)

#By general loc, year, month
agg <- aggregate(data=Grab, Total.Organisms.sum~ Year+Month+Order+GeneralLoc, function(x) sum(x))
agg <- agg[which(agg$GeneralLoc==c('GoodHarbor','SouthManitou')), ]
ggplot(agg, aes(x=Month, y=Total.Organisms.sum, fill=Order)) +
  geom_bar(stat="identity", position = "fill") +
  Goby_theme + scale_x_continuous(breaks=pretty_breaks(n=10))+
  labs(x="Month", y="Proportion of sample taxa", title="Benthos over time")+
  facet_grid(Year ~ GeneralLoc + .)

#####PLOT SHANNON #####

ggplot(c,aes(x=DayNum,y=FamilyShannonDI, color=Year)) + 
  geom_point() +
  facet_grid(All.SiteCondition ~ Depth.m.standard + .)+
  xlab("DayNum") + ylab("Shannon Diversity Index") + Goby_theme

####Basic differences between site types - alive, bare, and sloughed####
#Are more animals captured in the clad. samples?

#First, compare total orgs. by substrate alone
sum.d <- summaryBy(Total.Organisms.sum ~ All.SiteCondition + YearSER, data = d, FUN=sum)
sum.d$All.SiteCondition <- as.factor(sum.d$All.SiteCondition)
kruskal.test(Total.Organisms.sum.sum ~ All.SiteCondition, data = sum.d) 
boxplot(Total.Organisms.sum.sum~All.SiteCondition,data=sum.d, main="Total orgs. by site condition", xlab="Site condition", ylab="Total number of organisms")

#Next, compare them by taxonomic group and substrate
sum.d3 <- summaryBy(Total.Organisms.sum ~ All.SiteCondition + TaxonomicGroup + YearSER, data = d, FUN=sum)
ggplot(sum.d3,aes(x=TaxonomicGroup,y=Total.Organisms.sum.sum,fill=TaxonomicGroup)) + geom_boxplot() + 
  facet_grid(All.SiteCondition~.,scales="free",space="free") + coord_cartesian(ylim = c(0, 2000)) +
   Goby_theme + theme(axis.text.x=element_text(angle = 90, hjust = 0))

ggsave(filename = "F:/DATA/SLBE/Manuscripts/Benthos-analysis/Figs/NumbySub.png")

#finally, do another test to look at which organisms are very different and which sites are different from eachother
library(mgcv)
allyears <- gam(Total.Organisms.sum.sum ~ All.SiteCondition + TaxonomicGroup, data=sum.d3)
gam.check(allyears)
summary(allyears)

########Change in benthos over time? by loc? ##########
d3<- d[ , which(names(d) %in% c("Month","Order","Total.Organisms.sum","Year","GeneralLoc"))]
#d3$Month <- as.factor(d3$Month)
#levels(d3$Month)[levels(d3$Month) %in%  c("5","6","7","8")] <- "Summer"
#levels(d3$Month)[levels(d3$Month) %in%  c("9","10","11")] <- "Fall"
test<- summaryBy(Total.Organisms.sum ~ Order + Month + GeneralLoc + Year,data=d3, FUN=c(sum))
test<-test[!is.na(test$Order),]
test[is.na(test)] <- 0 


glm.model = glm(Total.Organisms.sum.sum ~ Order * GeneralLoc * Year, data=test, family = poisson)
anova(glm.model, test="Chisq")


###############################################################################
#### NMDS  ####
#Prep sample for Bray-curtis######
#freq <- cast(d, All.DepNonDep + YearSER ~ Family, value='Total.Organisms.sum', sum)

freq <- cast(d, YearSER ~ Family, value='Total.Organisms.sum', sum)

#names(freq)

#Remove NA values you don't want
freq<- freq[!(is.na(freq$YearSER)), ]

#Combine year and site to make sample name
#freq$Sample <- paste(freq[,1], freq[,2]) #if two vars
freq$Sample <- freq[,1] #if one var

row.names(freq)<-freq$Sample

####if mussels is a treatment
#hist(freq$Dreissenidae, breaks=30, xaxp  = c(0, 1200, 30))
freq$NumMussels <- ifelse(freq$Dreissenidae>50, "LotsofMuss", "FewMuss")
freq <- subset(freq, select = -c(Dreissenidae) )

###MAKE TREATMENTS FOR LATER
treat <- freq$NumMussels

###Choose cols/animals to remove
freq <- freq[ , -which(names(freq) %in% c("NumMussels","??UPDATE","Fish Egg","Sample","YearSER"))]

#Remove any columns that are all zeroes
freq <- freq[, colSums(freq == 0) != nrow(freq)] 

freq_mat<- data.matrix(freq) #######VERY IMPORTANT###CONVERT DATAFRAME TO MATRIX!!!!!#######

#shannon diversity
diversity(freq_mat, index = "shannon")

####BRAY
#binary=FALSE means you look at the number of individuals.  
#TRUE would give the result for presence-absence (Sorenson's index)
bc<-vegdist(freq_mat, method="bray", binary=FALSE)

####NMDS using Bray-Curtis ordination
freq.mds<-metaMDS(freq_mat, distance = "bray", k = 2) 
stressplot(freq.mds)

plot(freq.mds) #plots the ordination axes
points(freq.mds, display = c("sites", "species"))#displays both sites and species on the same plot.  
text(freq.mds, display = c("sites", "species"))

ordiplot(freq.mds,type="n")
ordihull(freq.mds,groups=treat,draw="polygon",col="grey90", label=T)
orditorp(freq.mds,display="species",col="red",air=0.01)

#add the sample labels#
orditorp(freq.mds,display="sites",col=c(rep("green",5),rep("blue",5)),
         air=0.01,cex=1.25)

####NMDS using mussels as the treatment

####################   PERMANOVA    ####################

#Better than ANOSIM http://www.esajournals.org/doi/abs/10.1890/12-2010.1

############## TEST FOR MULTIVARIATE SPREAD####
#http://thebiobucket.blogspot.com/2011/04/assumptions-for-permanova-with-adonis.html
#An important assumtption for PERMANOVA is same "multivariate spread" among groups, 
#which is similar to variance homogeneity in univariate ANOVA.

library(vegan)

# two similar populations:
dat1a<-matrix(sample(c(0,1,1,1),200,replace=T),10,20)
dat1b<-matrix(sample(c(0,1,1,1),200,replace=T),10,20)

# generating a third sample from the same population, but with reduced
# number of species occurrences. this set will have higher
# beta-diversity (or "multivariate spread"):
dat2<-matrix(sample(c(0,0,0,1),200,replace=T),10,20)

# distance matrices:
fac<-gl(2,10)
dist11<-vegdist(rbind(dat1a,dat1b))
dist12<-vegdist(rbind(dat1a,dat2)) 
  
# when computing sets with same beta-dispersion we get a
# correct adonis result with no sign. group differences
# in species composition:
anova(betadisper(dist11,fac))
adonis(rbind(dat1a,dat1b)~fac)

# when using sets with different beta-diversity you may
# get false significant adonis results - the location/composition
# is actually the same (!) this result is due to different
# multivariate spread in dat1 and dat2:
anova(betadisper(dist12,fac))
adonis(rbind(dat1a,dat2)~fac)

# see ordination diagram where location (centroids) between dat1 and dat2
# is not shifted more than for dat1a dat1b, still you yield a (false)
# sign. adonis result
# plot:

windows(10,5)

opar<-par()
par(mfrow=c(1,2))
plot(meta11<-metaMDS(dist11,zerodist=ignore),type="n",
     main="same beta-disp\nsame location")
points(meta11,select=which(fac==1),col="red")
points(meta11,select=which(fac==2),col="blue")
ordispider(meta11,group=fac)

plot(meta12<-metaMDS(dist12,zerodist=ignore),type="n",
     main="diff beta-disp\nsame location")
points(meta12,select=which(fac==1),col="red")
points(meta12,select=which(fac==2),col="blue")
ordispider(meta12,group=fac)

par(opar)


#################################################

#Get your species data
freq <- cast(d, YearSER ~ Family, value='Total.Organisms.sum', sum)

#Figure out which environmental vars you want
names(c)
environ <- c[,c("YearSER","DayNum","Depth.m.standard","Year","Month", "DieoffYear","All.SiteCondition")] #"Latitude","Longitude",

#environSERs <- environ$YearSER[is.na(environ$All.SiteCondition)]
#environ<- environ[!is.na(environ$All.SiteCondition),]
#freq <- freq[!freq$YearSER %in% environSERs, ]

row.names(environ)<- environ$YearSER
environ$YearSER<-NULL
mode(environ$DieoffYear) <- "integer"

#Make site condition into a dummy
dummy <- as.data.frame(model.matrix( ~ All.SiteCondition - 1, data=environ))
environ<- cbind(dummy,environ)
environ$All.SiteCondition <-NULL

row.names(freq)<- freq$YearSER
freq$YearSER<-NULL

environ$NumMussels <- freq$Dreissenidae
freq <- subset(freq, select = -c(Dreissenidae) )
environ<-data.frame(lapply(environ,as.numeric))

##########
#Test for multivariate normality
library(MVN)
#https://cran.r-project.org/web/packages/MVN/vignettes/MVN.pdf
freqN <- mardiaTest(freq, qqplot = TRUE)
environN<- mardiaTest(environ, qqplot = TRUE)
result

#Multivariate homogeneity of groups dispersions (variances)
betadisper(d, group, type = c("centroid", "median"))
#http://r-sig-ecology.471788.n2.nabble.com/multiple-factors-in-vegan-betadisper-td5556506.html

##########

###adonis(freq ~ DayNum*Depth.m.standard*Year*Month*Latitude*Longitude*DieoffYear*Clad.*NumMussels, data=environ, method="bray",permutations=999)
###lat/long may need to be removed bc it makes the permanova blow up

perma <- adonis(freq ~ DayNum*Depth.m.standard*Year*All.SiteConditionBARE*All.SiteConditionLIVE*All.SiteConditionSLOUGHED*NumMussels, data=environ, method="bray",permutations=999)


########################################################################################
#### Forams Abundance-Biomass Calculation Curves ######
#See Chap 8 Primer-E Manual for explanation

#Need data in this format:
#row.names (Species)  N (abundance)  Biomass (Total biomass)

#REMEMBER that the mussels are all removed from the W-stat calculations

####PRINT A PLOT FOR EACH YEARSER TO A PDF ########
somePDFPath = "F:/DATA/SLBE/Manuscripts/Benthos-analysis/Figs/ABCplots.pdf"
pdf(file=somePDFPath)
par(mfrow=c(3,2));

for (cat in unique(agg2$YearSER)){
  e <- subset(agg2, agg2$YearSER == cat)
  plot(abc(e[,2:3]),
       main=paste("YearSER", cat))
}

dev.off() 

###plotting W-stats etc. ###

c$Year <- as.numeric(c$Year)

boxplot(ABC.W.Stat~GeneralLoc,data=c, main="ABC W-Stat", xlab="General Location", ylab="ABC W-Stat")
boxplot(ABC.W.Stat~Year,data=c, main="ABC W-Stat", xlab="Year", ylab="ABC W-Stat")
boxplot(ABC.W.Stat~All.DepNonDep,data=c, main="ABC W-Stat", xlab="All.DepNonDep", ylab="ABC W-Stat")
boxplot(ABC.W.Stat~Depth.m.standard,data=c, main="ABC W-Stat", xlab="Depth", ylab="ABC W-Stat")

#wilcox tests for W stats
wilcox.test(ABC.W.Stat~Depth.m.standard, data = c)
wilcox.test(ABC.W.Stat~All.DepNonDep, data = c)

#2-way anova on the ABC stat
int <- aov(ABC.W.Stat ~ Depth.m.standard*All.DepNonDep, data=c)
summary(int)


######Come up with a method to determine site condition######
#####classification tree ###

CompiledBenthosLog$DayNum <- yday(CompiledBenthosLog$Date)

rpartanaly <- CompiledBenthosLog[,c("All.All.DepNonDep","Predominant.Benthic.Class","BroadBPIMEAN","CurveMEAN",
                                    "DepthMEAN","SlopemMEAN","CosAspMEAN", "SinAspMEAN","AspectMEAN" ,
                                    "DayNum","All.SiteCondition","Predom.Simp.Ben.Class", "Sub.CladMap", "YearSER")]

rpartanaly <- rpartanaly[!is.na(rpartanaly$All.SiteCondition),]

rpartanaly <- rpartanaly[!rpartanaly$All.SiteCondition %in% c("ALIVE?", "ALIVE?SLOUGHED?", "BARE or ALIVE", "DEAD"),]
rpartanaly <- rpartanaly[!rpartanaly$All.All.DepNonDep == "DEP?",]

library(rpart)
library(rpart.plot)	

# grow tree 
fit <- rpart(All.SiteCondition ~ All.All.DepNonDep + DayNum  + Predominant.Benthic.Class + SlopemMEAN + 
               AspectMEAN + DepthMEAN + BroadBPIMEAN + Sub.CladMap,
             method="class", data=rpartanaly)

printcp(fit) # display the results 
plotcp(fit) # visualize cross-validation results 
summary(fit) # detailed summary of splits

#Variables actually used in tree construction:
#[1] All.All.DepNonDep             AspectMEAN                BroadBPIMEAN              DayNum                   
#[5] DepthMEAN                 Predominant.Benthic.Class Sub.CladMap              

#Root node error: 322/621 = 0.51852

#n=621 (30 observations deleted due to missingness)

#CP nsplit rel error  xerror     xstd
#1 0.391304      0   1.00000 1.00000 0.038669
#2 0.208075      1   0.60870 0.60870 0.035968
#3 0.071429      2   0.40062 0.40062 0.031396
#4 0.024845      4   0.25776 0.25776 0.026335
#5 0.020186      5   0.23292 0.25155 0.026064
#6 0.012422      8   0.16770 0.25466 0.026200
#7 0.010000     10   0.14286 0.23602 0.025363

#Node number 1: 621 observations,    complexity param=0.3913043
#predicted class=BARE      expected loss=0.5185185  P(node) =1
#class counts:   299   146   176
#probabilities: 0.481 0.235 0.283 
#left son=2 (395 obs) right son=3 (226 obs)

prp(fit, uniform=TRUE, main="Classification Tree for Site Condition")

####PREDICT###
classthese <- CompiledBenthosLog[,c("YearSER","All.All.DepNonDep","Predominant.Benthic.Class","BroadBPIMEAN","CurveMEAN",
                                    "DepthMEAN","SlopemMEAN","CosAspMEAN", "SinAspMEAN","AspectMEAN" ,
                                    "DayNum","All.SiteCondition","Predom.Simp.Ben.Class", "Sub.CladMap")]

classthese <- classthese[!classthese$All.SiteCondition %in% c("LIVE", "SLOUGHED", "BARE"),]

row.names(classthese) <- classthese$YearSER

preds <- predict(fit, new=classthese, type = c("prob"))

preds <- as.data.frame(preds)

preds$YearSER <- row.names(preds)

combo <- CompiledBenthosLog[,c("YearSER","Date","Site","Year")]

preds <- join(preds, combo, by="YearSER", type="left")
rm(combo)

preds <- preds[order(preds$Year,preds$Site,preds$Date), ]
preds$Pred.Site.Cond <- colnames(preds)[apply(preds[,c("BARE","LIVE","SLOUGHED")],1,which.max)]

#write.csv(preds, "predicted.site.conditions.csv")

##########################################################################
######Run cluster analysis#########

hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) as.dist((1-cor(t(x)))/2)

d <- distfunc(rog_diet2011)
fit <- hclustfunc(d)
plot(fit, hang=-1)

####### Anosim analysis#####################################################
#WORSE THAN PERMANOVA http://www.esajournals.org/doi/abs/10.1890/12-2010.1
library(sinkr)

#Get your species data
freq <- cast(d, YearSER ~ Family, value='Total.Organisms.sum', sum)

#Figure out which environmental vars you want
names(c)
environ <- c[,c("YearSER","DayNum","Depth.m.standard","Year","CladRate")]

#"CloudCoverfraction.mean", "ModelWaterLevelmeter.mean",     
#"WaterVelocityatSurfacems.mean", "DepthAveragedWaterVelocityms.mean","SignificantWaveHeightmeter.mean", "WavePeriodsecond.mean","circ.WavesVelocityDir.mean",                
#"circ.DAvWaterVelocityDir.mean","circ.WaterVelocityDir.mean"

environSERs <- environ$YearSER[is.na(environ$CladRate)]
environ<- environ[!is.na(environ$CladRate),]
freq <- freq[!freq$YearSER %in% environSERs, ]

row.names(environ)<- environ$YearSER
environ$YearSER<-NULL

row.names(freq)<- freq$YearSER
freq$YearSER<-NULL

#hist(freq$Dreissenidae, breaks=30, xaxp  = c(0, 1200, 30))
#environ$HiMussels<- ifelse(freq$Dreissenidae > 50, 1, 0)
environ$NumMussels <- freq$Dreissenidae

freq <- subset(freq, select = -c(Dreissenidae) )

environ<-data.frame(lapply(environ,as.numeric))

res <- bioEnv(freq, environ,
              fix.dist.method="bray", var.dist.method="euclidean",
              scale.fix=FALSE, scale.var=TRUE)
res

summary(res)