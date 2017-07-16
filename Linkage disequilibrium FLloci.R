#load required packages
library(vegan)
library(tidyr)
library(dplyr)
library(doParallel)

cl<- makeCluster(3)
registerDoParallel(cl)

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
    
    #create output dataframe
    LD.df.D<-as.data.frame(matrix(nrow=11, ncol=88))
    names(LD.df.D)[1:20]<-c("1.2","1.3","1.4","1.5","2.")
    LD.df.r<-as.data.frame(matrix(nrow=11, ncol=88))
    
    #set working directory
    start_wd<-(paste("C:/Users/Madeline/Desktop/Weis lab/EEB498/Updated_sets", paste("para_set", e, sep="_"), paste("model_run_", o, sep=""), sep="/"))
    setwd(start_wd)

    #start loop over generations
    for (i in c(1, seq(from = 50, to = 500, by = 50))) {
      #read in generation dataframe
      offspring <- read.csv(paste(getwd(), paste("offspring_map", i, sep="_"), sep="/"))
      
      #create dataframe to store numerical allele scores
      offspring_map<-subset(offspring, select=10:(9+2*no_FL))
      
      #loop over population and number of neutral alleles to fill in allele score dataframe
      for (k in 1:pop_size){
        for (j in 1:no_FL)
        {
          offspring_map[k,1+(2*(j-1))]<-if (offspring[k,10+(2*(j-1))] >= 0) {
            1} else {0}
          offspring_map[k,2+(2*(j-1))]<-if (offspring[k,11+(2*(j-1))] >= 0) {
            1} else {0}
          
          #fill in columns that store loci scores
          offspring_map[k,2*no_FL+(j)]<-sum(offspring_map[k,1+(2*(j-1))],offspring_map[k,2+(2*(j-1))])
        }
      }
      
      hapseq.1<-seq(from=1, to=2*no_FL, by=2)
      hapseq.2<-seq(from=2, to=2*no_FL, by=2)
      
      offspring.hap1<-select(offspring_map, hapseq.1)
      offspring.hap2<-select(offspring_map, hapseq.2)
      
      ####################################################################################

        for (A in 1:(no_FL-1)){ #for focal locus
          for (B in (A+1):no_FL){ #for comparative locus
            p.a<-sum(offspring.hap1[,A])/(pop_size) #probability of having A allele
            p.b<-sum(offspring.hap1[,B])/(pop_size)  #probability of having B allele
            p.ab<-0
            for (x in 1:pop_size) {
              if (offspring.hap1[x,A]+offspring.hap1[x,B] == 2) {p.ab<-p.ab+1}
            }
            p.ab<-p.ab/pop_size #probability of having A and B allele
            
            D<-(p.ab-(p.a*p.b)) #disequilibrium coefficient
            r<-D/sqrt(p.a*(1-p.a)*p.b*(1-p.b)) #correlation coefficient
            if (is.na(r)) {r <- "N"}
            if (i == 1) {
              LD.df.D[1,(A*(B))]<-D
              LD.df.r[1,(A*(B))]<-r
            } else if (i > 1) {
              LD.df.D[((i/50)+1),(A*(B))]<-D
              LD.df.r[((i/50)+1),(A*(B))]<-r
            } #end if/else
          } #end b
        } #end a
    } #end for loop over generations
    
    #remove columns with NA
    LD.df.D<-LD.df.D[ , colSums(is.na(LD.df.D)) == 0]
    names(LD.df.D)<-c("1.2","1.3","1.4","1.5","2.3","2.4","2.5","3.4","3.5","4.5")
    
    LD.df.r<-LD.df.r[ , colSums(is.na(LD.df.r)) == 0]
    LD.df.r[is.na(LD.df.r)] <- 0
    names(LD.df.r)<-c("1.2","1.3","1.4","1.5","2.3","2.4","2.5","3.4","3.5","4.5")
    
    write.csv(LD.df.D, paste("LD.df.D.", o, ".csv", sep=""))
    write.csv(LD.df.r, paste("LD.df.r.", o, ".csv", sep=""))
    
  } #end for loop over runs
  
} #end for loop over parameter sets
