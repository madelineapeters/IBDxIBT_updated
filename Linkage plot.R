#load required packages
library(vegan)
library(tidyr)
library(dplyr)

#set parameters
no_runs<-5
FLmean<-100
no_neutral<-10
pop_size<-400
no_FL<-5

#make generation dataframe
plot.df<-as.data.frame(matrix(nrow=11, ncol=1))
names(plot.df)<-"gen"
plot.df[1:11,1]<-c(1, seq(from=50, to=500, by=50))

for (e in c(54,56,58)) {
  #set parameters
  no_runs<-5
  FLmean<-100
  no_FL<-10
  pop_size<-400
  no_FL<-10
  

  start_wd<-(paste("C:/Users/Madeline/Desktop/Weis lab/EEB498/Updated_sets", paste("para_set", e, sep="_"), sep="/"))
  setwd(start_wd)
    
  #ensure all dataframes are loaded into the environment
  df.sig<-paste("LD.N.D.avg.", e, sep="")
  assign(df.sig, read.csv("LD.D.avg.csv", header = TRUE))
  
  rm(df.sig)
  
} #next e
    
dflist.sig<-lapply(ls(pattern = "LD.N.D.avg.*"), get)

LD.D<-bind_cols(dflist.sig)
names(LD.D)<-c("P55","P57","P59")

LD.D<-bind_cols(plot.df, LD.D)

library(ggplot2)
LD.D.plot<-ggplot(data=LD.D, aes(x=gen))+
  geom_line(aes(y=P55), colour="darkgreen", size=2)+
  geom_line(aes(y=P57), colour="yellow", size=2)+
  geom_line(aes(y=P59), colour="darkblue", size=2)+
  xlab("Generation (0-500)")+ylab("Disequilibrium coefficient")+
  theme_bw()+theme(axis.text.x = element_blank())+ylim(-0.01,0.15)
plot(LD.D.plot)
