library(dplyr)
library(stats)

seq<-c(1, seq(from=50, to=500, by=50))

for (p in 111){
  for (r in 1:3){
    start_wd<-(paste("C:/Users/Madeline/Desktop/Weis lab/EEB498", paste("para_set", p, sep="_"), paste("model_run", r, sep="_"), sep="/"))
    setwd(start_wd)
    
    Poisson.df<-as.data.frame(matrix(nrow=11, ncol=6))
    names(Poisson.df)<-c("Mom.mean","Mom.var", "Mom.p", "Dad.mean", "Dad.var", "Dad.p")
    
    for (z in seq) {
      NH <- read.csv(paste(getwd(), paste("NH", z, "csv", sep="."), sep="/"))
      pop_size<-dim(NH)[1]
      
      Mother.df<-as.data.frame(table(NH$Mother))
      Zero.mom<-as.data.frame(matrix(nrow=(pop_size-dim(Mother.df)[1]), ncol=2))
      names(Zero.mom)<-c("Var1", "Freq")
      Zero.mom[,2]<-0
      Mother.comb<-bind_rows(Mother.df, Zero.mom)
      
      Father.df<-as.data.frame(table(NH$Father))
      Zero.dad<-as.data.frame(matrix(nrow=(pop_size-dim(Father.df)[1]), ncol=2))
      names(Zero.dad)<-c("Var1", "Freq")
      Zero.dad[,2]<-0
      Father.comb<-bind_rows(Father.df, Zero.dad)

      
      poisson.mom<-poisson.test(sum(Mother.comb$Freq), length(Mother.comb$Freq), mean(Mother.comb$Freq))
      poisson.dad<-poisson.test(sum(Father.comb$Freq), length(Father.comb$Freq), mean(Father.comb$Freq))

      Poisson.df$Mom.mean[match(z, seq)]<-mean(Mother.comb$Freq)
      Poisson.df$Mom.var[match(z, seq)]<-var(Mother.comb$Freq)
      Poisson.df$Dad.mean[match(z, seq)]<-mean(Father.comb$Freq)
      Poisson.df$Dad.var[match(z, seq)]<-var(Father.comb$Freq)
      Poisson.df$Mom.p[match(z, seq)]<-poisson.mom$p.value
      Poisson.df$Dad.p[match(z, seq)]<-poisson.dad$p.value
    } #end z
    
    write.csv(Poisson.df, "Poisson.df.csv")
    
  } #end r
} #end p
    