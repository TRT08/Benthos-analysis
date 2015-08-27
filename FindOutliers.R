
# create a factor version of Taxon a factor, ordered by the median length (for plotting)
medianlength <- tapply(outties$mm.length, outties$Taxon, median, na.rm=TRUE)
itemfactor <- factor(outties$Taxon, levels=names(sort(medianlength)))

# add median length to data frame
index <- tapply(outties$mm.length, outties$Taxon)
outties$medianlen <- medianlength[index]

# boxplots of all the diet item sizes
windows(h=9, w=6.5)
par(mar=c(4, 12, 1, 1), cex=0.6)
boxplot(outties$mm.length ~ itemfactor, horizontal=TRUE, las=1, xlab="Length  (mm)", ylab="", main="")

rm(medianlength, itemfactor)

# outlier measure: the number of interquartile ranges beyond the lower and upper quartiles

# calculate quartiles and interquartile range
q1 <- tapply(outties$mm.length, outties$Taxon, quantile, 0.25, na.rm=TRUE)[index]
q3 <- tapply(outties$mm.length, outties$Taxon, quantile, 0.75, na.rm=TRUE)[index]
iqr <- q3 - q1

# set outlier measure to zero
outmeas <- rep(0, length(iqr))
# change it to NA if the interquartile range is missing
outmeas[is.na(iqr)] <- NA
# calculate the number of interquartile ranges below the first quartile
sel <- !is.na(outties$mm.length) & outties$mm.length < q1
outmeas[sel] <- ((outties$mm.length - q1)/iqr)[sel]
# calculate the number of interquartile ranges above the third quartile
sel <- !is.na(outties$mm.length) & outties$mm.length > q3
outmeas[sel] <- ((outties$mm.length - q3)/iqr)[sel]
# save the result to the outties data frame
outties$outmeas <- outmeas

rm(index, q1, q3, iqr, outmeas, sel)

# sort the data by outmeas (in absolute value) to look at potential outliers
outord <- outties[order(-abs(outties$outmeas), outties$medianlen, outties$mm.length), ]

# a common cutoff is 1.5, which is a good place to start,
# but you may find that you want to look at more or fewer potential outliers than are suggested by this cutoff
sel <- !is.na(outord$outmeas) & abs(outord$outmeas) > 1.5
outord<- outord[sel, ]

rm(sel)

