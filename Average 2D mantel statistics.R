#load required packages
library(tidyr)
library(dplyr)
library(data.table)

for (e in 64) {
  #set parameters
  no_runs<-10
  FLmean<-100
  no_neutral<-10
  pop_size<-400
  
  for (r in 1:no_runs) {
    #Set the working directory
    start_wd<-(paste("G:/final sets", paste("para_set", e, sep="_"), paste("model_run_", r, sep=""), sep="/"))
    setwd(start_wd)
    
    #ensure all dataframes are loaded into the environment
    df.sig<-paste("2D.mantel.sig.", r, sep="")
    assign(df.sig, read.csv(paste(df.sig, ".csv", sep=""), header = TRUE))
    df.r<-paste("2D.mantel.r.", r, sep="")
    assign(df.r, read.csv(paste(df.r, ".csv", sep=""), header = TRUE))
    
  } #end loop over runs
  
  rm(df.sig)
  rm(df.r)
  
  #put dataframes into lists
  dflist.sig<-lapply(ls(pattern = "2D.mantel.sig.*"), get)
  dflist.r<-lapply(ls(pattern = "2D.mantel.r.*"), get)
  
  #bind data frames and rename columns
  total.sig<-bind_cols(dflist.sig)
  for (g in 1:(no_runs*no_neutral+no_runs)){
    names(total.sig)[g]<-paste(g)
  }
  total.r<-bind_cols(dflist.r)
  for (g in 1:(no_runs*no_neutral+no_runs)){
    names(total.r)[g]<-paste(g)
  }
  
  #remove X columns (labelling generations)
  for(r in 1:no_runs) {
    if(r == 1){
      total.sig<-select(total.sig, -(1))
      total.r<-select(total.r, -(1))
    } else{
      total.sig<-select(total.sig, -((r-1)*no_neutral+1))
      total.r<-select(total.r, -((r-1)*no_neutral+1))
    }
  }
  
  #sum over rows to get average statistics
  total.sig<-rowMeans(total.sig, na.rm = TRUE, dims = 1)
  total.r<-rowMeans(total.r, na.rm = TRUE, dims = 1)

  #reset the working directory
  setwd('../')
  
  #save averaged dataframes
  write.csv(total.sig, paste("2D.mantel.sig.avg.csv"), row.names=F)
  write.csv(total.r, paste("2D.mantel.r.avg.csv"), row.names=F)
  
  #clear environment before next parameter set
  rm(list=ls())
  
} #end loop over parameter sets