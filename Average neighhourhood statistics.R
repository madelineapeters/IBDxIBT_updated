#load required packages
library(tidyr)
library(dplyr)
library(data.table)

#set parameters
no_runs<-5
FLmean<-100
no_neutral<-10
pop_size<-100
no_gen<-100

#create generation label dataframe
plot.df<-as.data.frame(matrix(nrow=no_gen, ncol=1))
plot.df[1:no_gen,1]<-(1:no_gen)
names(plot.df)[1]<-"gen"

#start loop over parameter sets to make average dataframes
for (e in seq(from=82, to=85, by=1)) {
  
  #loop over runs to call in run-specific dataframes
  for (r in 1:no_runs) {
    #Set the working directory
    start_wd<-(paste("C:/Users/Madeline/Desktop/Weis lab/EEB498", paste("para_set", e, sep="_"), paste("model_run_", r, sep=""), sep="/"))
    setwd(start_wd)
    
    #ensure all dataframes are loaded into the environment
    df.N<-paste("sum.N.", r, sep="")
    assign(df.N, read.csv("sum.N.csv", header = TRUE))
  } #end loop over runs 
  
  #remove dummy variable
  rm(df.N)
  
  #put dataframes into lists
  dflist.N<-lapply(ls(pattern = "sum.N.*"), get) #to produce dataframe for neighbourhood size
  dflist.d<-dflist.N #to produce dataframe for density
  
  #bind data frames
  total.N<-bind_cols(dflist.N)
  names(total.N)[1:(no_runs*3)]<-c(1:(no_runs*3))
  total.d<-bind_cols(dflist.d)
  names(total.d)[1:(no_runs*3)]<-c(1:(no_runs*3))
  
  #remove unwanted columns
  for(r in 1:no_runs) {
      total.N<-select(total.N, -((r+1):(r+2))) #removes variance and density
      total.d<-select(total.d, -(r:(r+1))) #removes neighbourhood size and variance
  }
  
  #sum over rows to get average statistics
  total.N<-mutate(total.N, avg=rowMeans(total.N, na.rm = TRUE, dims = 1))
  total.d<-mutate(total.d, avg=rowMeans(total.d, na.rm = TRUE, dims = 1))

  #reset the working directory
  setwd('../')
  
  #save averaged dataframes
  write.csv(total.N, paste("avg.d.csv"), row.names=F)
  write.csv(total.d, paste("avg.d.csv"), row.names=F)
  
  #clear environment before next parameter set
  rm(total.N, total.d, sum.N.1, sum.N.2, sum.N.3, sum.N.4, sum.N.5)
  
} #end loop over parameter sets

#start loop to call in average dataframes for parameter sets to compare
for (e in c(84,85)) {
  #Set the working directory
  start_wd<-(paste("C:/Users/Madeline/Desktop/Weis lab/EEB498", paste("para_set", e, sep="_"), sep="/"))
  setwd(start_wd)
  
  #ensure all dataframes are loaded into the environment
  df.N<-paste("avg.d.", e, sep="")
  assign(df.N, read.csv("avg.d.csv", header = TRUE))
  
  df.d<-paste("avg.d.", e, sep="")
  assign(df.d, read.csv("avg.d.csv", header = TRUE))
  
} #end loop over runs 

rm(df.N, df.d)

avg.d.84<-bind_cols(avg.d.84, plot.df)
avg.d.84<-bind_cols(avg.d.84, plot.df)
avg.d.85<-bind_cols(avg.d.85, plot.df)
avg.d.85<-bind_cols(avg.d.85, plot.df)
avg.d.84[1:6]<-avg.d.84[1:6]/4
avg.d.84[1:6]<-avg.d.84[1:6]/4
avg.d.85[1:6]<-avg.d.85[1:6]/4

library(ggplot2)

d.plot<-ggplot()+geom_line(data=avg.d.84, aes(x=gen, y=avg), color="blue", size=2)+
  geom_line(data=avg.d.85, aes(x=gen, y=avg), color="darkgreen", size=2)+
  geom_line(data=avg.d.82, aes(x=gen, y=avg), color="blue", size=2)+
  geom_line(data=avg.d.83, aes(x=gen, y=avg), color="darkgreen", size=2)+
  xlab("Generation")+ylab("Effective density")+theme_bw()+
  ylim(0,1)

plot(d.plot)

geom_line(data=avg.d.84, aes(x=gen, y=`X3`))+
  geom_line(data=avg.d.85, aes(x=gen, y=`X3`))+
  geom_line(data=avg.d.84, aes(x=gen, y=`X6`))+
  geom_line(data=avg.d.85, aes(x=gen, y=`X6`))+
  geom_line(data=avg.d.84, aes(x=gen, y=`X9`))+
  geom_line(data=avg.d.85, aes(x=gen, y=`X9`))+
  geom_line(data=avg.d.84, aes(x=gen, y=`X12`))+
  geom_line(data=avg.d.85, aes(x=gen, y=`X12`))+
  geom_line(data=avg.d.84, aes(x=gen, y=`X15`))+
  geom_line(data=avg.d.85, aes(x=gen, y=`X15`))+

