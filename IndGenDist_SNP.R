library(dplyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)

library("poppr")
library("pegas")
library("ape")
library("adegenet")
library("ade4")

##############################################################
## Convert the data
##############################################################

Mydata <- read.table("Master_Pinus_data_genotype.txt", header = TRUE, check.names = FALSE)   
dim(Mydata) 

ind <- as.character(Mydata$tree_id) # individual labels 
population <- as.character(Mydata$state) # population labels
county <- Mydata$county 

#Create matrix with only genotypes (keeping only first 100 SNPs)
#Convert matrix to a  genind object (for the package adegenet). #The genind object can then easily be converted into loci objects (package pegas) (i.e. Mydata2)
locus <- Mydata[, -c(1, 2, 3, 4, 105:ncol(Mydata))] 
Mydata1 <- df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep="")
Mydata2 <-genind2loci(Mydata1)

##############################################################
## Individual genetic distance: euclidean distance
##############################################################

#Use the function dist() from adegenet 
#Use euclidean distance among vector of allele frequencies
distgenEUCL <- dist(Mydata1, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
hist(distgenEUCL)

##############################################################
## Individual genetic distance: num. loci at which ind. differ
##############################################################

#Option pairwise.deletion = FALSE in the command dist.gene() removes all loci with one missing values
distgenDIFF <- dist.gene(Mydata2, method="pairwise", pairwise.deletion = FALSE, variance = FALSE)
hist(distgenDIFF)

#Get percent missing data per population
missing_data <- info_table(Mydata1, type = "missing")
sum(missing_data["Total", 1:100] > 0)
barplot(missing_data["Total", 1:100], xlab = "Locus", ylab = "Percent Missing")

#Keep loci with missing values
distgenDIFF <- dist.gene(Mydata2, method="pairwise", pairwise.deletion = TRUE, variance = FALSE)
hist(distgenDIFF)

##############################################################
## Number allelic differences between two individuals
##############################################################
distgenDISS <- diss.dist(Mydata1, percent = FALSE, mat = FALSE)
hist(distgenDISS)
