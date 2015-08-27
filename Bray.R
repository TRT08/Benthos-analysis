#Prep sample for Bray-curtis

library(vegan)

freq <- freq.year2

names(freq)

#Remove NA values you don't want
#freq<- freq[!(is.na(freq$DepNonDep)), ]

###MAKE TREATMENTS FOR LATER
Year<- freq$Year
GeneralLoc <- freq$GeneralLoc
DepNonDep<-freq$DepNonDep

#Combine year and site to make sample name
#freq$Sample <- paste(freq[,1], freq[,2]) 

row.names(freq)<-freq$YearSER

###Choose cols/animals to remove
freq <- freq[ ,-which(names(freq) %in% c("YearSER","DepNonDep" ))]

#Remove any columns that are all zeroes
freq <- freq[, colSums(freq == 0) != nrow(freq)] 

#######VERY IMPORTANT#########
######CONVERT DATAFRAME TO MATRIX!!!!!#######
freq_mat<- data.matrix(freq)

#shannon diversity
diversity(freq_mat, index = "shannon")


####Rarefaction curve
spi<-specaccum(freq, method="rarefaction")
plot(spi)

####BRAY
#binary=FALSE means you look at the number of individuals.  
#TRUE would give the result for presence-absence (Sorenson's index)
bc<-vegdist(freq_mat, method="bray", binary=FALSE)

####NMDS####
#using Bray-Curtis ordination
freq.mds<-metaMDS(freq_mat, distance = "bray", k = 2) 
stressplot(freq.mds)

plot(freq.mds) #plots the ordination axes
points(freq.mds, display = c("sites", "species"))#displays both sites and species on the same plot.  
text(freq.mds, display = c("sites", "species"))

treat <- DepNonDep
ordiplot(freq.mds,type="n")
ordihull(freq.mds,groups=treat,draw="polygon",col="grey90", label=T)
orditorp(freq.mds,display="species",col="red",air=0.01)
orditorp(freq.mds,display="sites",col=c(rep("green",5),rep("blue",5)),
         air=0.01,cex=1.25)


########OTHER CRAP##################
###Select a year to run...or run all of them together (proceed to next step)
rog_diet2011 <-rog_diet3[311:660,]
rog_diet2012 <-rog_diet3[10:310,]
rog_diet2013 <-rog_diet3[1:9,]


###Run NMDS
rog_diet.nmds <-monoMDS(rog_diet.dis)   
plot(rog_diet.nmds)


###Run cluster analysis

hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) as.dist((1-cor(t(x)))/2)

d <- distfunc(rog_diet2011)
fit <- hclustfunc(d)
plot(fit, hang=-1)