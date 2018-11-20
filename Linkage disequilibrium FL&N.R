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
    
    LD.FLN.D.avg<-as.data.frame(matrix(nrow=11, ncol=no_FL))
    LD.FLN.r.avg<-as.data.frame(matrix(nrow=11, ncol=no_FL))
    
    for (f in 1:no_FL){
    #create output dataframe
    LD.df.D<-as.data.frame(matrix(nrow=11, ncol=15))
    LD.df.r<-as.data.frame(matrix(nrow=11, ncol=15))
    
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
      
      
      #create dataframe to store numerical allele scores
      offspring_FL<-as.data.frame(matrix(nrow=pop_size), col=1)
      names(offspring_FL)<-"FL.1a"
      
      #create dataframe for FL allele scores
      j<-f
      for (k in 1:pop_size){
          offspring_FL[k,1]<-if (offspring[k,10+(2*(j-1))] >= 0) {
            1} else {0}
      }
      
      N.FL.df<-bind_cols(offspring_FL, offspring_map)
      
      ####################################################################################

        A<-1 #for focal locus
          for (B in 2:(no_neutral+1)){ #for comparative locus
            p.a<-sum(N.FL.df[,A])/(pop_size) #probability of having A allele
            p.b<-sum(N.FL.df[,B])/(pop_size)  #probability of having B allele
            p.ab<-0
            for (x in 1:pop_size) {
              if (N.FL.df[x,A]+N.FL.df[x,B] == 2) {p.ab<-p.ab+1}
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
          } #end b
    } #end for loop over generations
    
    #remove columns with NA
    LD.df.D<-LD.df.D[ , colSums(is.na(LD.df.D)) == 0]
    LD.df.D<-abs(LD.df.D)
    
    x<-1
    for(c in 2:15) {
      if (is.na(unique(LD.df.r[,c]))) {x <- c(x,c)}
    }
    LD.df.r<-select(LD.df.r, -x)
    LD.df.r<-abs(LD.df.r)
    
    LD.FLN.D.avg[f]<-rowMeans(LD.df.D)
    LD.FLN.r.avg[f]<-rowMeans(LD.df.r, na.rm=TRUE)
    
    } #end loop over FL loci 
    write.csv(LD.FLN.D.avg, paste("LD.FLN.D.", o, ".csv", sep=""))
    write.csv(LD.FLN.r.avg, paste("LD.FLN.r.", o, ".csv", sep=""))
  } #end for loop over runs
  
} #end for loop over parameter sets
