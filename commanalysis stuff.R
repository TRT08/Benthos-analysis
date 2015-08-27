require(vegan) 
data(varespec) 
data(varechem) 
alt.varechem <- varechem  
alt.varechem$N <- log(alt.varechem$N) 
res<-bioenv(wisconsin(varespec), alt.varechem) 
res

library(sinkr) #better version
res <- bioEnv(wisconsin(varespec), varechem,
              fix.dist.method="bray", var.dist.method="euclidean",
              scale.fix=FALSE, scale.var=TRUE)
res

data(bryceveg) # returns a vegetation data.frame
dis.bc <- dsvdis(bryceveg,'bray/curtis') # returns a dissimilarity matrix
clust <- sample(1:5,nrow(bryceveg),replace=TRUE)
s <- indval(bryceveg,clust)

library(vegan)
data(dune)
data(dune.env)
dune.dist <- vegdist(dune)
attach(dune.env)
dune.ano <- anosim(dune.dist, Management)
summary(dune.ano)
plot(dune.ano)