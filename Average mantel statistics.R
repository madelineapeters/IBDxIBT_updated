#load required packages
library(tidyr)
library(dplyr)
library(data.table)

for (e in c(52,53)) {
  #set parameters
  no_runs<-3
  FLmean<-100
  no_neutral<-10
  pop_size<-400
  
  for (r in 1:no_runs) {
    #Set the working directory
    start_wd<-(paste("C:/Users/Madeline/Desktop/Weis lab/EEB498/", paste("para_set", e, sep="_"), paste("/model_run_", r, sep=""), sep=""))
    setwd(start_wd)
    
    #ensure all dataframes are loaded into the environment
    for (df in list.files(pattern = 'mantel.sig.*')){
      assign(df, read.csv(df, header = TRUE))
    } 
    for (df in list.files(pattern = 'mantel.r.*')){
      assign(df, read.csv(df, header = TRUE))
    }
  } #end loop over runs
  
  #put dataframes into lists
  dflist.sig<-lapply(ls(pattern = "mantel.sig.*"), get)
  dflist.r<-lapply(ls(pattern = "mantel.r.*"), get)
  
  #bind data frames and rename columns
  total.sig<-bind_cols(dflist.sig)
  for (g in 1:(3*no_neutral+no_runs)){
    names(total.sig)[g]<-paste(g)
  }
  total.r<-bind_cols(dflist.r)
  for (g in 1:(3*no_neutral+no_runs)){
    names(total.r)[g]<-paste(g)
  }
  
  #remove X columns (labelling generations)
  for( r in 1:no_runs) {
    total.sig<-select(total.sig, -(10*(r-1)+1))
  }
  for( r in 1:no_runs) {
    total.r<-select(total.r, -(10*(r-1)+1))
  }
  
  #sum over rows to get average statistics
  total.sig<-rowMeans(total.sig, na.rm = TRUE, dims = 1)
  total.r<-rowMeans(total.r, na.rm = TRUE, dims = 1)

  #reset the working directory
  start_wd<-(paste("C:/Users/Madeline/Desktop/Weis lab/EEB498/", paste("para_set", e, sep="_"), sep=""))
  setwd(start_wd)
  
  #save averaged dataframes
  write.csv(total.sig, paste("mantel.sig.avg.csv"), row.names=F)
  write.csv(total.r, paste("mantel.r.avg.csv"), row.names=F)
  
  #clear environment before next parameter set
  rm(list=ls())
  
} #end loop over parameter sets