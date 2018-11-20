#load required packages
library(vegan)
library(tidyr)
library(dplyr)

for (e in seq(from = 54, to = 65, by = 1)) {
  
  library(vegan)
  library(tidyr)
  library(dplyr)
  library(doParallel)
  
  #set parameters
  no_runs<-5
  FLmean<-100
  no_neutral<-10
  pop_size<-400
  no_FL<-5
  
  #start loop over runs
  for (o in 1:no_runs) {
    
    for(f in 1:no_neutral){
    
    LD.FLN.D.avg<-as.data.frame(matrix(nrow=11, ncol=no_neutral))
    LD.FLN.r.avg<-as.data.frame(matrix(nrow=11, ncol=no_neutral))
    
    #create output dataframe
    LD.df.D<-as.data.frame(matrix(nrow=11, ncol=100))
    LD.df.r<-as.data.frame(matrix(nrow=11, ncol=100))
    
    #set working directory
    start_wd<-(paste("C:/Users/Madeline/Desktop/Weis lab/EEB498/Updated_sets", paste("para_set", e, sep="_"), paste("model_run_", o, sep=""), sep="/"))
    setwd(start_wd)

    #start loop over generations
    for (i in c(1, seq(from = 50, to = 500, by = 50))) {
      #read in generation dataframe
      offspring <- read.csv(paste(getwd(), paste("offspring_map", i, sep="_"), sep="/"))
      
      #create dataframe to store numerical allele scores
      offspring_map<-as.data.frame(matrix(nrow=pop_size, ncol=no_neutral))
      names(offspring_map)<-c("1a","2a","3a","4a","5a", "6a", "7a", "8a", "9a", "10a")
      #loop over population and number of neutral alleles to fill in allele score dataframe
      for (k in 1:pop_size){
        for (j in 1:no_neutral)
        {
          offspring_map[k,j]<-if (offspring[k,10+(2*no_FL)+(2*(j-1))] == "D") {
            1} else {0}
        } #end loop over j
      } #end loop over k

      
      ####################################################################################

        A<-f #for focal locus
          for (B in 1:(no_neutral)){ #for comparative locus
            if (A != B) {
              p.a<-sum(offspring_map[,A])/(pop_size) #probability of having A allele
              p.b<-sum(offspring_map[,B])/(pop_size)  #probability of having B allele
              p.ab<-0
              for (x in 1:pop_size) {
                if (offspring_map[x,A]+offspring_map[x,B] == 2) {p.ab<-p.ab+1}
              }
              p.ab<-p.ab/pop_size #probability of having A and B allele
              
              D<-(p.ab-(p.a*p.b)) #disequilibrium coefficient
              r<-D/sqrt(p.a*(1-p.a)*p.b*(1-p.b)) #correlation coefficient
              if (i == 1) {
                LD.df.D[1,(B)]<-D
                LD.df.r[1,(B)]<-r
              } else if (i > 1) {
                LD.df.D[((i/50)+1),(B)]<-D
                LD.df.r[((i/50)+1),(B)]<-r
              } #end if/else
            }
          } #end b
    } #end for loop over generations
    
    #remove columns with NA
    LD.df.D<-LD.df.D[ , colSums(is.na(LD.df.D)) == 0]
    LD.df.D<-abs(LD.df.D)
    
    x<-1
    for(c in 2:100) {
      if (is.na(unique(LD.df.r[,c]))) {x <- c(x,c)}
    }
    LD.df.r<-select(LD.df.r, -x)
    LD.df.r<-abs(LD.df.r)
    
    LD.FLN.D.avg[f]<-rowMeans(LD.df.D)
    LD.FLN.r.avg[f]<-rowMeans(LD.df.r, na.rm=TRUE)
  } #end loop over neutral loci
    
    write.csv(LD.FLN.D.avg, paste("LD.N.D.", o, ".csv", sep=""))
    write.csv(LD.FLN.r.avg, paste("LD.N.r.", o, ".csv", sep=""))
  } #end for loop over runs
  
} #end for loop over parameter sets
