library("tidyr")
library("dplyr")
library("poppr")
library("pegas")
library("ape")
library("adegenet")
library("ade4")


s = 16
run = 1
gen = 500

##############################################################
## Convert the data
##############################################################
Mydata = read.csv(paste(getwd(),'/para_set_',s,'/model_run_',run,'/paraset_',s,'_offspring_map_',gen,'.csv',sep=""))
locus = Mydata %>% select(.,genotypeA,genotypeB,genotypeC)
locus$genotype1 = paste(Mydata$neut1a,Mydata$neut1b,sep="")
locus$genotype2 = paste(Mydata$neut2a,Mydata$neut2b,sep="")
locus$genotype3 = paste(Mydata$neut3a,Mydata$neut3b,sep="")
locus$genotype4 = paste(Mydata$neut4a,Mydata$neut4b,sep="")
locus$genotype5 = paste(Mydata$neut5a,Mydata$neut5b,sep="")
locus$genotype6 = paste(Mydata$neut6a,Mydata$neut6b,sep="")
locus$genotype7 = paste(Mydata$neut7a,Mydata$neut7b,sep="")
locus$genotype8 = paste(Mydata$neut8a,Mydata$neut8b,sep="")
locus$genotype9 = paste(Mydata$neut9a,Mydata$neut9b,sep="")
locus$genotype10 = paste(Mydata$neut10a,Mydata$neut10b,sep="")
locus$genotype11 = paste(Mydata$neut11a,Mydata$neut11b,sep="")
locus$genotype12 = paste(Mydata$neut12a,Mydata$neut12b,sep="")
locus$genotype13 = paste(Mydata$neut13a,Mydata$neut13b,sep="")
locus$genotype14 = paste(Mydata$neut14a,Mydata$neut14b,sep="")
locus$genotype15 = paste(Mydata$neut15a,Mydata$neut15b,sep="")
locus$genotype16 = paste(Mydata$neut16a,Mydata$neut16b,sep="")
locus$genotype17 = paste(Mydata$neut17a,Mydata$neut17b,sep="")
locus$genotype18 = paste(Mydata$neut18a,Mydata$neut18b,sep="")
locus$genotype19 = paste(Mydata$neut19a,Mydata$neut19b,sep="")
locus$genotype20 = paste(Mydata$neut20a,Mydata$neut20b,sep="")
locus$genotype21 = paste(Mydata$neut21a,Mydata$neut21b,sep="")
locus$genotype22 = paste(Mydata$neut22a,Mydata$neut22b,sep="")
locus$genotype23 = paste(Mydata$neut23a,Mydata$neut23b,sep="")
locus$genotype24 = paste(Mydata$neut24a,Mydata$neut24b,sep="")

#Convert matrix to a  genind object (for the package adegenet). #The genind object can then easily be converted into loci objects (package pegas) (i.e. Mydata2)
Mydata1 = df2genind(locus, ploidy = 2, sep="")
Mydata2 = genind2loci(Mydata1)

##############################################################
## Number allelic differences between two individuals
##############################################################
distgenDISS = diss.dist(Mydata1, percent = FALSE, mat = TRUE)

library(lattice) 
levelplot(distgenDISS[1:100,1:100])
heatmap(distgenDISS[1:100,1:100])
heatmap(distgenDISS[101:200,101:200])
