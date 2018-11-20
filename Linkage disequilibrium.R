#load required packages
library(vegan)
library(tidyr)
library(dplyr)
library(doParallel)

cl<- makeCluster(3)
registerDoParallel(cl)

foreach (e = seq(from = 54, to = 65, by = 1)) %dopar% {
  #set parameters
  no_runs<-10
  FLmean<-100
  no_neutral<-10
  pop_size<-400
  no_FL<-5
  
  #start loop over runs
  for (r in 1:no_runs) {
    
    #create output dataframe
    LD.df.D<-as.data.frame(matrix(nrow=11, ncol=88))
    LD.df.r<-as.data.frame(matrix(nrow=11, ncol=88))
    
    #set working directory
    start_wd<-(paste("C:/Users/Madeline/Desktop/Weis lab/EEB498/Updated_sets", paste("para_set", e, sep="_"), paste("model_run_", r, sep=""), sep="/"))
    setwd(start_wd)

    #start loop over generations
    for (i in c(1, seq(from = 50, to = 500, by = 50))) {
      #read in generation dataframe
      offspring <- read.csv(paste(getwd(), paste("offspring_map", i, sep="_"), sep="/"))
      
      #start loop to set up simplified neutral locus matrix
      for (j in 1:no_neutral)
      {
        names(offspring)[10+(2*5)+(2*(j-1))]<-paste("neut",j,"a",sep="")
        names(offspring)[11+(2*5)+(2*(j-1))]<-paste("neut",j,"b",sep="")
      }
      
      #create dataframe to store numerical allele scores
      offspring_map<-as.data.frame(matrix(nrow=pop_size, ncol=2*no_neutral+no_neutral))
      
      #loop over population and number of neutral alleles to fill in allele score dataframe
      for (k in 1:pop_size){
        for (j in 1:no_neutral)
        {
          offspring_map[k,1+(2*(j-1))]<-if (offspring[k,10+(2*5)+(2*(j-1))]=="d") {
            0} else {1}
          offspring_map[k,2+(2*(j-1))]<-if (offspring[k,11+(2*5)+(2*(j-1))]=="d") {
            0} else {1}
          
          #fill in columns that store loci scores
          offspring_map[k,2*no_neutral+(j)]<-sum(offspring_map[k,1+(2*(j-1))],offspring_map[k,2+(2*(j-1))])
        }
      }
      
      hapseq.1<-seq(from=1, to=2*no_neutral, by=2)
      hapseq.2<-seq(from=2, to=2*no_neutral, by=2)
      
      offspring.hap1<-select(offspring_map, hapseq.1)
      offspring.hap2<-select(offspring_map, hapseq.2)
      ####################################################################################
      
      r.fun<-function(a, b) {
        p.a<-sum(offspring.hap1[,a])/(pop_size) #probability of having A allele
        p.b<-sum(offspring.hap2[,b])/(pop_size)  #probability of having B allele
        p.ab<-0
        for (x in 1:pop_size) {
          if (offspring.hap1[x,a]+offspring.hap1[x,b] == 4) {p.ab<-p.ab+1}
        }
        p.ab<-p.ab/pop_size #probability of having A and B allele
        
        D<-(p.ab-(p.a*p.b)) #disequilibrium coefficient
        r<-D/sqrt(p.a*(1-p.a)*p.b*(1-p.b)) #correlation coefficient
      } #end r.fun function

        for (A in 1:(no_neutral-1)){ #for focal locus
          for (B in (A+1):no_neutral){ #for comparative locus
            p.a<-sum(offspring.hap1[,A])/(pop_size) #probability of having A allele
            p.b<-sum(offspring.hap2[,B])/(pop_size)  #probability of having B allele
            p.ab<-0
            for (x in 1:pop_size) {
              if (offspring.hap1[x,A]+offspring.hap1[x,B] == 2) {p.ab<-p.ab+1}
            }
            p.ab<-p.ab/pop_size #probability of having A and B allele
            
            D<-(p.ab-(p.a*p.b)) #disequilibrium coefficient
            r<-D/sqrt(p.a*(1-p.a)*p.b*(1-p.b)) #correlation coefficient
            if (i == 1) {
              LD.df.D[1,(A*(B-1))]<-D
              LD.df.r[1,(A*(B-1))]<-r
            } else if (i > 1) {
              LD.df.D[((i/50)+1),(A*(B-1))]<-D
              LD.df.r[((i/50)+1),(A*(B-1))]<-r
            } #end if/else
          } #end b
        } #end a
    } #end for loop over generations
    
    write.csv(LD.df.D, paste("LD.df.D.", r, ".csv", sep=""))
    write.csv(LD.df.r, paste("LD.df.r.", r, ".csv", sep=""))
    
  } #end for loop over runs
  
} #end for loop over parameter sets

library(ggplot2)
r.plot<-ggplot()+geom_line(aes(y=LD.df.r$V1))
plot(r.plot)

plot(LD.df.r$V1)
