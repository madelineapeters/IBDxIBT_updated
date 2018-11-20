#load required packages
library(tidyr)
library(dplyr)
library(data.table)

for(e in starting.set1:ending.set2) {
  
  for (r in 1:no_runs) {
    #Set the working directory
    start_wd<-(paste(base_wd, sub_wd, paste("para_set", e, sep="_"), paste("model_run_", r, sep=""), sep="/"))
    setwd(start_wd)
    
    #ensure all dataframes are loaded into the environment
    df.var<-paste("neutral.dev.p.", r, sep="")
    assign(df.var, read.csv(paste(df.var, ".csv", sep=""), header = TRUE))
  } #end loop over runs
  
    rm(df.var)
    
  #put dataframes into lists
  dflist.var<-lapply(ls(pattern = "neutral.dev.p.*"), get)
 
  #bind data frames and rename columns
  total.var<-bind_cols(dflist.var)
  for (g in 1:(no_runs*no_neutral+no_runs)){
    names(total.var)[g]<-paste(g)
  }
  
  #remove X columns (labelling generations)
  for (r in 1:no_runs) {
    if (r == 1){
      total.var<-select(total.var, -(1))
    } else {
      total.var<-select(total.var, -(10*(r-1)+1))
    }
  }
  
  #sum over rows to get average statistics
  total.var<-rowMeans(total.var, na.rm = TRUE, dims = 1)

  #reset the working directory
  setwd('../')
  
  #save averaged dataframes
  write.csv(total.var, paste("neutral.dev.p.avg.csv"), row.names=F)
  
  #clear environment before next parameter set
  rm(dflist.var)
  
} #end loop over parameter sets