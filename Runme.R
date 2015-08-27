#########################################
##########BENTHOS MEGA SCRIPT############
#### for biomass and frequency calcs#####
#########################################

###Code by Taaja Tucker
###Last update 9/25/14

###RUN FIRST###
setwd("F:/DATA/SLBE/Manuscripts/Benthos-analysis/") #set directory
logfile <- file("output.log") #creates a backlog of the script activity (look in folder)
sink(logfile, append=TRUE) 
sink(logfile, append=TRUE, type="message")
source("biomassfreqscript.R",echo=TRUE, max.deparse.length=10000)
source("biomass_additional_functions.R",echo=TRUE, max.deparse.length=10000)

sink()  # Restore output to console
sink(type="message")
rm(NA.omit.biomass)
#Now manually open the output.log file in the folder and check for any warnings
###############

##Files available to work with:
#1. CompiledBenthosLogandFreq - MOST IMPORTANT FINAL PRODUCT - contains frequency/count data for mussel sheet + 
#benthos data (Total Orgs and Count) but also the total biomass for each taxon. 
#You can see the intermediate steps in these calculation processes in the following columns:
#Av.biomass.mg- the average biomass (mg) per measured critter
#Sum.biomass.mg - The sume of the biomass of each critter per taxon per SER
#Prpn.biomass.calc - The ratio of critters measured / number counted
#Final.biomass.mg -  Av.biomass.mg*Total.Organisms (this is inherently the sum of biomass if Prpn.biomass.calc=1)

#2.Freq.All.Benthos - AN INTERMEDIATE FILE - contains all individual critters from mussel sheet and benthos data and lengths (mm) calculated 
#for all non-missing critters. Takes only "Total Organisms" and not "Count" into account to prevent double biomass
#counts (e.g. chironomid heads vs bodies).

#3. CompiledBenthosLog - This file has information for each ponar drop (whether or not it was processed) and shows the associated BT info, the video ratings, weather, and more!

###GUIDE TO SOME OF THE THINGS YOU'LL FIND IN THE FILE: CompiledBenthosLogandFreq
# "Total.Organisms.sum" - Sum of all organisms found per SER / taxon. Heads are included if they were counted as an organism and rolled up into one category
# e.g. A chironomidae head was counted as one organism and 5 chironomidae TL were counted. This value would = 6, and the taxon would be "Chironomidae"                       
# "Num.critters.averaged.4.biomass"
#"Av.biomass.mg"  - the average biomass of the animals measured                           
#"Sum.biomass.mg" - the total biomass of ONLY the animals measured                           
#"Prpn.biomass.calc" - The proportion of the animals that were measured compared to the total counted     
#"Final.biomass.mg"   -  Av.biomass.mg  *  Total.Organisms.sum                                           
#"Max.BT.depth.m" - The maximum depth the BT was dropped to by SER, which was used roll up BT data with benthos data
#"BT.YearSer"  - the YearSER of the BT drop    
#"BT.Time.Diff.Days" - the difference in the number of days between the BT drop and the ponar drop
#"BT.depth.and.ponar.diff.m"  - the difference in the depth between the maximum BT drop depth and the ponar drop depth
#"DepthDifference.m" - the difference between the standard depth (e.g. 10, 20, 30 m) and the actual depht of the ponar drop                                                          
#"Video.Time.Diff.Days"  - the difference in the number of days between the video recording and the ponar drop 
#"GeneralLoc" - the general location of the sample (e.g. South Manitou)                                
#"DepNonDep"  - Based on the site designations (2013) or the video ratings (from 2012).                                
#"DieoffYear"  - True or False. 2010&2012=dieoff years, 2011+2013 not die off years.                                                                                                    
                    
#Below are examples of calculations to use. You can make your own by changing the variables in the equations!
######BIOMASS CALCULATIONS########
##Calculate total biomass
names(CompiledBenthosLogandFreq)
sum(CompiledBenthosLogandFreq$Final.biomass.mg, na.rm = TRUE)
summary(CompiledBenthosLogandFreq$Final.biomass.mg)

######Biomass determination by SER
#mean biomass
mean.biomass.mg <- cast(CompiledBenthosLogandFreq, Year + Site.Condition.For.Event + Depth.m.standard ~ TaxonomicGroup, mean, value='Final.biomass.mg') 
#Total biomass
sum.biomass.mg <- cast(CompiledBenthosLogandFreq, TaxonomicGroup ~  Year + Site.Condition.For.Event + Depth.m.standard, sum, value='Final.biomass.mg')

######Freq determination by Year and general location
freq.year2 <- cast(CompiledBenthosLogandFreq,  YearSER + DepNonDep~ Taxon, value='Total.Organisms.sum', sum)



