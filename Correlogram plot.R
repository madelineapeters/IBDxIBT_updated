
#load required packages
library(tidyr)
library(dplyr)
library(doParallel)
library(spatial)

#set up for DoParallel
cl<- makeCluster(3)
registerDoParallel(cl)

#read in stats parameters
source(file=paste("~/498/stats parameters.R"))

#set up dataframe to store corr values
corr.df<-as.data.frame(matrix(nrow=20, ncol=4))
names(corr.df)<-c("Dist",1:9)

#start parallel loop
for (e in c(55, 57, 59)) {
  
  ############Set-up############ 
  #re-read in stats parameters for each core
  source(file=paste("~/498/stats parameters.R"))
  
  #select run
  r<-1
  
  #select generation
  i<-500
    
  #set working directory
  start_wd<-(paste(base_wd, sub_wd, paste("para_set", e, sep="_"), paste("model_run_", r, sep=""), sep="/"))
  setwd(start_wd)
  
  #read in generation dataframe
  offspring <- read.csv(paste(start_wd, paste("offspring_map", i, sep="_"), sep="/"))
  
  #set up and fill in FL day matrix  
  FL_matrix<-matrix(nrow=pop_size, ncol=pop_size)
  
  for (m in 1:pop_size) {
    for (p in 1:pop_size) {
      FL_matrix[m,p]<-abs((offspring$FLday[m]-offspring$FLday[p])/FLmean)
    } #end p loop
  } #end m loop
  
  #set up and fill in geographic distance matrix
  distance1<-matrix(nrow=pop_size, ncol=pop_size)
  distance2<-matrix(nrow=pop_size, ncol=pop_size)

  for (m in 1:pop_size) {
    for (p in 1:pop_size) {
      distance1[m,p]<-sqrt((offspring$X_pos[m]-offspring$X_pos[p])^2 + (offspring$Y_pos[m] - offspring$Y_pos[p])^2)
      distance2[m,p]<-sqrt(((sqrt(pop_size)+1)-abs(offspring$X_pos[m]-offspring$X_pos[p]))^2 + ((sqrt(pop_size)+1)-abs(offspring$Y_pos[m]-offspring$Y_pos[p]))^2)
    }
  }
  
  distance1<-distance1/sum(distance1, na.rm=TRUE)
  distance2<-distance2/sum(distance2, na.rm=TRUE)
  distance_matrix<-distance1+distance2
  distance_matrix<-distance_matrix/sum(distance_matrix, na.rm=TRUE)
  
  #run correlogram function
  FL.mantel<-mantel.correlog(D.eco=FL_matrix, D.geo=distance_matrix, nperm=1000)

  #fill in correlogram dataframe)
  if (e == 54 || 55) {
    for (s in 9:17){
      FL.string<-capture.output(FL.mantel)[s]
      FL.string<-substr(FL.string, 10,19)
      corr.df[(s-8),1]<-FL.string
    } #end for s loop
  } #end if loop
  for (s in 9:17){
    FL.string<-capture.output(FL.mantel)[s]
    FL.string<-substr(FL.string, 33, 44)
    if (e == 54 || 56 || 58){
      corr.df[(s-8),((e/2)-25)]<-FL.string
    } else {
      corr.df[(s-8),(((e-1)/2)-25)]<-FL.string
    } #end else
  } #end for s loop 

} #end loop over parameter sets   

write.csv(corr.df, "corr.df.csv")
