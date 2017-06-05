#load required packages
library(tidyr)
library(dplyr)

for (e in 42:63) {
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
    df.FL.sig<-paste("2D.mantel.FL.sig.", r, ".csv", sep="")
    df.FL.sig<-as.data.frame(assign(df.FL.sig, read.csv(df.FL.sig, header = TRUE)))
    df.FL.r<-paste("2D.mantel.FL.r.", r, ".csv", sep="")
    df.FL.r<-as.data.frame(assign(df.FL.r, read.csv(df.FL.r, header = TRUE)))
    
  } #end loop over runs
  
  #put dataframes into lists
  dflist.sig<-lapply(ls(pattern = "2D.mantel.FL.sig.*"), get)
  dflist.r<-lapply(ls(pattern = "2D.mantel.FL.r.*"), get)
  
  
  #bind data frames and rename columns
  total.sig<-bind_cols(dflist.sig)
  for (g in 1:(2*no_runs)){
    names(total.sig)[g]<-paste(g)
  }
  total.r<-bind_cols(dflist.r)
  for (g in 1:(2*no_runs)){
    names(total.r)[g]<-paste(g)
  }
  
  #remove X columns (labelling generations)
  for(r in 1:no_runs) {
    if (r == 1){
      total.sig<-select(total.sig, -(1))
      total.r<-select(total.r, -(1))
    } else {
      total.sig<-select(total.sig, -(r))
      total.r<-select(total.r, -(r))
    }
  }
  
  #sum over rows to get average statistics
  total.sig<-rowMeans(total.sig, na.rm = TRUE, dims = 1)
  total.r<-rowMeans(total.r, na.rm = TRUE, dims = 1)

  #reset the working directory
  setwd('../')
  
  #save averaged dataframes
  write.csv(total.sig, paste("FL.mantel.sig.avg.csv"), row.names=F)
  write.csv(total.r, paste("FL.mantel.r.avg.csv"), row.names=F)
  
  #clear environment before next parameter set
  rm(list=ls())
  
} #end loop over parameter sets