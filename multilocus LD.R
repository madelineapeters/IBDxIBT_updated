library(pegas)

#####as.loci function#####
# S3 method for class 'data.frame'
as.loci(x, allele.sep = "/|", col.pop = NULL, col.loci = NULL, ...)

#Arguments
  #x - an object of class "loci" or "genind", a data frame, a factor, or a vector of mode character.
  #allele.sep - the character(s) separating the alleles for each locus in the data file (a forward slash by default)
  #col.pop - specifies whether one of the column of the data file identifies the population; default NULL, otherwise an integer or a character giving the number or the name of the column.
  #col.loci - a vector of integers or of characters specifying the indices or the names of the columns that are loci. By default, all columns are taken as loci except the one labelled "population", if present or specified.
  #ploidy - the ploidy level (see details).
  #na.alleles - a vector of charater strings giving the alleles to be treated as missing data.

#####LD functions#####
#These two functions analyse linkage disequilibrium in the case of phased (LD) or unphased (LD2) genotypes
LD(x, locus = c(1, 2), details = TRUE) #based on the observed frequencies of haplotypes
LD2(x, locus = c(1, 2), details = TRUE) #based on the observed frequencies of different genotypes

#Arguments
  #x - an object of class "loci"
  #locus - a vector of two integers giving the loci to analyse
  #details - a logical value indicating whether to print the correlation matrix among alleles

#Values
  #For LD: if details = TRUE, a named list with the following elements:
      #Observed frequencies - the counts of haplotypes in the data.
      #Expected frequencies - the expected frequencies of haplotypes computed from the observed proportions of alleles under the assumption of no linkage disequilibrium.
      #Correlations among alleles - the observed correlations among alleles from both loci.
      #LRT (G-squared) - the likelihood-ratio test of the null hypothesis of no linkage disequilibrium.
      #Pearson’s test (chi-squared) - the chi-squared test based on haplotypes counts.
      #T2 the T2 test with its number of degrees of freedom (df) 
  #For LD2: if details = TRUE, a named list with two elements:
      #Delta the correlations among alleles (denoted Delta in Schaid 2004)
      #T2 the T2 test with its number of degrees of freedom (df).

#####Run linkage analyses#####
#load in required packages
library(pegas)

#set base directory
base_wd<-"~/498/final sets"

#parameters
pop_size<-400
no_neutral<-10
no_FL<-5
no_runs<-5

#start loop over parameter sets
for (x in c(54,55,56,57)) {

  #start loop over runs
  for (r in 1:no_runs) {
    
    #set working directory
    start_wd<-(paste(base_wd, paste("para_set", x, sep="_"), paste("model_run_", r, sep=""), sep="/"))
    setwd(start_wd)
    
    #create dataframes to hold stats outputs (WHAT STATS DO I SAVE JFC)
    LD.r<-as.data.frame(matrix(nrow=11, ncol=5))
    mantel.sig<-as.data.frame(matrix(nrow=11, ncol=5))
    
    #start loop over generations
    for (i in c(1, seq(from = 50, to = 250, by = 50))) {
      #read in generation dataframe
      offspring <- read.csv(paste(getwd(), paste("offspring_map", i, sep="_"), sep="/"))
      
      #create dataframe to store numerical allele scores
      offspring_map<-as.data.frame(matrix(nrow=pop_size, ncol=2*5+5))
      
      #loop over population and number of neutral alleles to fill in allele score dataframe
      for (k in 1:pop_size){
        for (j in 1:5)
        {
          offspring_map[k,1+(2*(j-1))]<-if (offspring[k,(2*5)+(2*(j-1))]=="1") {
            "a"} else if (offspring[k,(2*5)+(2*(j-1))]=="0.5") {"a"} else {"A"}
          offspring_map[k,2+(2*(j-1))]<-if (offspring[k,(2*5)+1+(2*(j-1))]=="1") {
            "a"} else if (offspring[k,(2*5)+1+(2*(j-1))]=="0.5") {"a"} else {"A"}
          
          #fill in columns that store loci scores
          offspring_map[k,2*5+(j)]<-paste(offspring_map[k,1+(2*(j-1))],offspring_map[k,2+(2*(j-1))],sep="/")
        }
      }
      
      #after filling in columns with loci scores, remove allele score columns 
      offspring_map<-subset.data.frame(offspring_map, select= (2*5+1):(2*5+5))
      
      ####################################################################################
      
      #transform df to LD compatible object
      LD.offspring_map<-as.loci(offspring_map, allele.sep = "/|", col.pop = NULL, col.loci = NULL, ...)
      
      LD.out<-LD(LD.offspring_map, locus = c(1, 2, 3), details = TRUE)
      
      
      
    } #end for loop over generations
   
    #Remove text from Mantel output dataframes
    mantel.r <- as.data.frame(sapply(mantel.r,gsub,pattern="Mantel statistic r: ", replacement=""))
    mantel.sig <- as.data.frame(sapply(mantel.sig,gsub,pattern="Significance: ", replacement=""))
    
    write.csv(mantel.sig, paste("mantel.FL.loci.sig.", r, ".csv", sep=""))
    write.csv(mantel.r, paste("mantel.FL.loci.r.", r, ".csv", sep=""))
    
  } #end for loop over runs
} #end loop over parameter sets

