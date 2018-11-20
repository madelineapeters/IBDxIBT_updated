#load required packages
library(vegan)
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
  mantel.FL.r<-as.data.frame(matrix(nrow=11, ncol=1))
  mantel.FL.sig<-as.data.frame(matrix(nrow=11, ncol=1))
  
  #start loop over generations
  for (i in c(1, seq(from = 50, to = 500, by = 50))) {
    #read in generation dataframe
    offspring_map <- read.csv(paste(getwd(), paste("offspring_map", i, sep="_"), sep="/"))
    
    ####################################################################################
    #set up geographic distance matrix
    distance_matrix<-matrix(nrow=pop_size, ncol=pop_size)
    
    #fill in geographic distance matrix
    for (m in 1:pop_size) {
      for (p in 1:pop_size) {
        distance_matrix[m,p]<-sqrt((offspring_map$X_pos[m]-offspring_map$X_pos[p])^2 + (offspring_map$Y_pos[m] - offspring_map$Y_pos[p])^2)	
      }}
  
    #set up and fill in FL day matrix  
    FL_matrix<-matrix(nrow=pop_size, ncol=pop_size)
    
    for (m in 1:pop_size) {
      for (p in 1:pop_size) {
        
        FL_matrix[m,p]<-abs((offspring_map$FLday[m]-offspring_map$FLday[p])/FLmean)
      }}
      
    #run Mantel test
    library(vegan)
    FL_mantel<-mantel(distance_matrix, FL_matrix, permutations=1000)
      
      #store output from Mantel test
      if (i == 1) {
        mantel.FL.r[i,1]<-paste((capture.output(FL_mantel))[7])
        mantel.FL.sig[i,1]<-paste((capture.output(FL_mantel))[8])
      } else {
        mantel.FL.r[((i/50)+1),1]<-paste((capture.output(FL_mantel))[7])
        mantel.FL.sig[((i/50)+1),1]<-paste((capture.output(FL_mantel))[8])
      } #end store output if/else

  } #end for loop over generations
  
  
  #Remove text from Mantel output dataframes
  mantel.r <- as.data.frame(sapply(mantel.FL.r,gsub,pattern="Mantel statistic r: ", replacement=""))
  mantel.sig <- as.data.frame(sapply(mantel.FL.sig,gsub,pattern="Significance: ", replacement=""))
  
  write.csv(mantel.FL.sig, paste("2D.mantel.FL.sig.", r, ".csv", sep=""))
  write.csv(mantel.FL.r, paste("2D.mantel.FL.r.", r, ".csv", sep=""))
  
} #end for loop over runs

} #end for loop over parameter sets

