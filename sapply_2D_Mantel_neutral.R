#load required packages
library(tidyr)
library(dplyr)
library(doParallel)

cl<- makeCluster(2)
registerDoParallel(cl)

foreach (e = starting.set1:final.set1) %dopar% {

#start loop over runs
for (r in 1:no_runs) {

  #set working directory
  start_wd<-(paste(base_wd, sub_wd, paste("para_set", e, sep="_"), paste("model_run_", r, sep=""), sep="/"))
  setwd(start_wd)
  
  #create dataframes to hold stats outputs
  mantel.r<-as.data.frame(matrix(nrow=11, ncol=no_neutral))
  mantel.sig<-as.data.frame(matrix(nrow=11, ncol=no_neutral))
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
    
    #after filling in columns with loci scores, remove allele score columns 
    offspring_map<-subset.data.frame(offspring_map, select= (2*no_neutral+1):(2*no_neutral+no_neutral))
    
    ####################################################################################
    #set up geographic distance matrix
    distance_matrix<-matrix(nrow=pop_size, ncol=pop_size)
    
    #fill in geographic distance matrix
    for (m in 1:pop_size) {
      for (p in 1:pop_size) {
        
        distance_matrix[m,p]<-sqrt((offspring$X_pos[m]-offspring$X_pos[p])^2 + (offspring$Y_pos[m] - offspring$Y_pos[p])^2)	
      }}
  
    #start loop over neutral loci to run and record Mantel test  
    for (n in 1:no_neutral){
      
      #set up matrix to store per-locus genetic distance 
      neut_dist_matrix<-matrix(nrow=pop_size, ncol=pop_size)
      
      #fill in genetic distance matrix
      for (m in 1:pop_size) {
        for (p in 1:pop_size) {
          neut_dist_matrix[m,p]<-abs(offspring_map[m,(1+(n-1))]-offspring_map[p,(1+(n-1))])
        }
      }
      
      #run Mantel test
      if (sum(neut_dist_matrix) == 0) {
        if (i == 1) {
          mantel.r[i,n]<-paste("NA")
          mantel.sig[i,n]<-paste("NA")
        } else {
          mantel.r[((i/50)+1),n]<-paste("NA")
          mantel.sig[((i/50)+1),n]<-paste("NA")
        }
      } else {
        library(vegan)
        neut_mantel<-mantel(distance_matrix, neut_dist_matrix, permutations=1000)
        
        #store output from Mantel test
        if (i == 1) {
          mantel.r[i,n]<-paste((capture.output(neut_mantel))[7])
          mantel.sig[i,n]<-paste((capture.output(neut_mantel))[8])
        } else {
          mantel.r[((i/50)+1),n]<-paste((capture.output(neut_mantel))[7])
          mantel.sig[((i/50)+1),n]<-paste((capture.output(neut_mantel))[8])
        } #end store output if/else
      } #end run Mantel if/else
      
    } #end loop over n

  } #end for loop over generations
  
  
  #Remove text from Mantel output dataframes
  mantel.r <- as.data.frame(sapply(mantel.r,gsub,pattern="Mantel statistic r: ", replacement=""))
  mantel.sig <- as.data.frame(sapply(mantel.sig,gsub,pattern="Significance: ", replacement=""))
  
  write.csv(mantel.sig, paste("2D.mantel.sig.", r, ".csv", sep=""))
  write.csv(mantel.r, paste("2D.mantel.r.", r, ".csv", sep=""))
  
} #end for loop over runs

} #end for loop over parameter sets

