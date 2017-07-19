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
plot.df<-as.data.frame(matrix(nrow=100, ncol=1))
names(plot.df)<-"gen"
plot.df[1:100,1]<-c(1, seq(from=2, to=100, by=1))

for (e in c(86, 87, 88, 89)) {
  #set parameters
  no_runs<-3
  FLmean<-100
  no_neutral<-10
  pop_size<-400
  no_FL<-5
  

  for (r in 1:no_runs) {
    start_wd<-(paste("C:/Users/Madeline/Desktop/Weis lab/EEB498", paste("para_set", e, sep="_"), paste("model_run", r, sep="_"), sep="/"))
    setwd(start_wd)
      
    #ensure all dataframes are loaded into the environment
    df.sig<-paste("self.", r, sep="")
    assign(df.sig, read.csv("self.df.csv", header = TRUE))
    
    rm(df.sig)
    
  } #end r loop
  
  dflist.sig<-lapply(ls(pattern = "self.*"), get)
  
  LD.d<-bind_cols(dflist.sig)
  LD.d<-as.data.frame(rowMeans(LD.d))
  plot.df<-bind_cols(plot.df, LD.d)
  
  rm(self.1, self.2, self.3)
  
} #next e

names(plot.df)[2:5]<-c("P86","P87","P88","P89")

library(ggplot2)
LD.D.plot<-ggplot(data=plot.df, aes(x=gen))+
  geom_line(aes(y=P86), colour="darkblue", size=2)+
  geom_line(aes(y=P87), colour="darkblue", size=2)+
  geom_line(aes(y=P88), colour="darkgreen", size=2)+
  geom_line(aes(y=P89), colour="darkgreen", size=2)+
  xlab("Generation (0-500)")+ylab("Selfing rate")+
  theme_bw()+theme(axis.text.x = element_blank())+ylim(0,1)
plot(LD.D.plot)
