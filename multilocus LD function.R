#load required packages
library(dplyr)
      
#set base directory
base_wd<-"~/498/final sets"

#parameters
pop_size<-400
no_neutral<-10
no_FL<-5
no_runs<-5

#make generation label dataframe
plot.df<-as.data.frame(matrix(nrow=11, ncol=1))
plot.df[1:11,1]<-c(1, seq(from=50, to=500, by=50))
names(plot.df)[1]<-"gen"

#start loop over parameter sets
for (x in c(54,55,56,57)) {
  
  #create dataframes to hold stats output
  LD.5<-as.data.frame(matrix(nrow=11, ncol=5))
  LD.4<-as.data.frame(matrix(nrow=11, ncol=5))
  LD.3<-as.data.frame(matrix(nrow=11, ncol=5))
  LD.avg<-as.data.frame(matrix(nrow=11, ncol=3))
  
  #start loop over runs
  for (o in 1:no_runs) {
    
    #set working directory
    start_wd<-(paste(base_wd, paste("para_set", x, sep="_"), paste("model_run_", o, sep=""), sep="/"))
    setwd(start_wd)
    
    #start loop over generations
    for (i in c(1, seq(from = 50, to = 500, by = 50))) {
      #read in generation dataframe
      offspring <- read.csv(paste(getwd(), paste("offspring_map", i, sep="_"), sep="/"))
      
      #create dataframe to store numerical allele scores
      offspring_map<-as.data.frame(matrix(nrow=pop_size, ncol=2*5))
      
      #loop over population and number of neutral alleles to fill in allele score dataframe
      for (k in 1:pop_size){
        for (j in 1:5)
        {
          offspring_map[k,1+(2*(j-1))]<-if (offspring[k,(2*5)+(2*(j-1))]=="1") {
            1} else if (offspring[k,(2*5)+(2*(j-1))]=="0.5") {1} else {0}
          #offspring_map[k,2+(2*(j-1))]<-if (offspring[k,(2*5)+1+(2*(j-1))]=="1") {
          #1} else if (offspring[k,(2*5)+1+(2*(j-1))]=="0.5") {1} else {0}

        }
      }
      
      offspring_map<-offspring_map[ , colSums(is.na(offspring_map)) == 0]
      
      ####################################################################################
      
      multiloc.LD<-function(df, var1=c(D), A, B, C=NA, D=NA, E=NA){
        if ((is.na(C))&(is.na(D))&(is.na(E))) {
          
          offspring_sub<-select(df, A, B)
          p.a<-sum(offspring_sub[,1])/(pop_size) #probability of having A allele
          p.b<-sum(offspring_sub[,2])/(pop_size)  #probability of having B allele
          c.ab<-0
          for (r in 1:pop_size){
            if (sum(offspring_sub[r,]) == 2) {c.ab<-c.ab+1}
          }
          p.ab<-c.ab/pop_size
          
          D<-(p.ab-(p.a*p.b)) #disequilibrium coefficient
          r<-D/sqrt(p.a*(1-p.a)*p.b*(1-p.b)) #correlation coefficient
          
        } else if ((is.na(D))&(is.na(E))) {
          
          offspring_sub<-select(df, A, B, C)
          p.a<-sum(offspring_sub[,1])/(pop_size) #probability of having A allele
          p.b<-sum(offspring_sub[,2])/(pop_size)  #probability of having B allele
          p.c<-sum(offspring_sub[,3])/(pop_size)  #probability of having C allele
          c.abc<-0
          for (r in 1:pop_size){
            if (sum(offspring_sub[r,]) == 3) {c.abc<-c.abc+1}
          }
          p.abc<-c.abc/pop_size
          
          D<-(p.abc-(p.a*p.b*p.c)) #disequilibrium coefficient
          r<-D/sqrt(p.a*(1-p.a)*p.b*(1-p.b)*p.c*(1-p.c)) #correlation coefficient
        } else if (is.na(E)) {
          offspring_sub<-select(df, A, B, C, D)
          p.a<-sum(offspring_sub[,1])/(pop_size) #probability of having A allele
          p.b<-sum(offspring_sub[,2])/(pop_size) #probability of having B allele
          p.c<-sum(offspring_sub[,3])/(pop_size) #probability of having C allele
          p.d<-sum(offspring_sub[,4])/(pop_size) #probability of having D allele
          c.abcd<-0
          for (r in 1:pop_size){
            if (sum(offspring_sub[r,]) == 4) {c.abcd<-c.abcd+1}
          }
          p.abcd<-c.abcd/pop_size
          
          D<-(p.abcd-(p.a*p.b*p.c*p.d)) #disequilibrium coefficient
          r<-D/sqrt(p.a*(1-p.a)*p.b*(1-p.b)*p.c*(1-p.c)*p.d*(1-p.d)) #correlation coefficient
        } else {
          offspring_sub<-select(df, A, B, C, D, E)
          p.a<-sum(offspring_sub[,1])/(pop_size) #probability of having A allele
          p.b<-sum(offspring_sub[,2])/(pop_size) #probability of having B allele
          p.c<-sum(offspring_sub[,3])/(pop_size) #probability of having C allele
          p.d<-sum(offspring_sub[,4])/(pop_size) #probability of having D allele
          p.e<-sum(offspring_sub[,5])/(pop_size) #probability of having E allele
          c.abcde<-0
          for (r in 1:pop_size){
            if (sum(offspring_sub[r,]) == 5) {c.abcde<-c.abcde+1}
          }
          p.abcde<-c.abcde/pop_size
          
          D<-(p.abcde-(p.a*p.b*p.c*p.d*p.e)) #disequilibrium coefficient
          r<-D/sqrt(p.a*(1-p.a)*p.b*(1-p.b)*p.c*(1-p.c)*p.d*(1-p.d)*p.e*(1-p.e)) #correlation coefficient
        }
      
        out<-var1
        return(out)
      } #end function
      
      #calculate frequency of full haplotype
      c.abcde<-0
      for (r in 1:pop_size){
        if (sum(offspring_map[r,]) == 5) {c.abcde<-c.abcde+1}
      }
      p.ABCDE<-c.abcde/pop_size
      
      #calculate frequency for 4/5 loci
      c.abcd<-0
      for (r in 1:pop_size){
        if (sum(offspring_map[r,1:4]) == 4) {c.abcd<-c.abcd+1}
      }
      p.ABCD<-c.abcd/pop_size
      
      #calculate frequency for 3/5 loci
      c.abc<-0
      for (r in 1:pop_size){
        if (sum(offspring_map[r,1:3]) == 3) {c.abc<-c.abc+1}
      }
      p.ABC<-c.abc/pop_size
      
      p.A<-sum(offspring_map[,1])/(pop_size)
      p.B<-sum(offspring_map[,2])/(pop_size)
      p.C<-sum(offspring_map[,3])/(pop_size)
      p.D<-sum(offspring_map[,4])/(pop_size)
      p.E<-sum(offspring_map[,5])/(pop_size)
      
      D.ABCD<-multiloc.LD(df=offspring_map, A=1, B=2, C=3, D=4)
      D.ABCE<-multiloc.LD(df=offspring_map, A=1, B=2, C=3, D=5)
      D.ACDE<-multiloc.LD(df=offspring_map, A=1, B=3, C=4, D=5)
      D.ABCE<-multiloc.LD(df=offspring_map, A=1, B=2, C=3, D=5)
      D.BCDE<-multiloc.LD(df=offspring_map, A=2, B=3, C=4, D=5)
      D.ABDE<-multiloc.LD(df=offspring_map, A=1, B=2, C=4, D=5)
      
      D.ABC<-multiloc.LD(df=offspring_map, A=1, B=2, C=3)
      D.ABD<-multiloc.LD(df=offspring_map, A=1, B=2, C=4)
      D.ABE<-multiloc.LD(df=offspring_map, A=1, B=2, C=5)
      D.ACE<-multiloc.LD(df=offspring_map, A=1, B=3, C=5)
      D.ACD<-multiloc.LD(df=offspring_map, A=1, B=3, C=4)
      D.ADE<-multiloc.LD(df=offspring_map, A=1, B=4, C=5)
      D.BCD<-multiloc.LD(df=offspring_map, A=2, B=3, C=4)
      D.BCE<-multiloc.LD(df=offspring_map, A=2, B=3, C=5)
      D.BDE<-multiloc.LD(df=offspring_map, A=2, B=4, C=5)
      D.CDE<-multiloc.LD(df=offspring_map, A=3, B=4, C=5)
      
      D.AB<-multiloc.LD(df=offspring_map, A=1, B=2)
      D.AC<-multiloc.LD(df=offspring_map, A=1, B=3)
      D.AD<-multiloc.LD(df=offspring_map, A=1, B=4)
      D.AE<-multiloc.LD(df=offspring_map, A=1, B=5)
      D.BC<-multiloc.LD(df=offspring_map, A=2, B=3)
      D.BD<-multiloc.LD(df=offspring_map, A=2, B=4)
      D.BE<-multiloc.LD(df=offspring_map, A=2, B=5)
      D.CD<-multiloc.LD(df=offspring_map, A=3, B=4)
      D.CE<-multiloc.LD(df=offspring_map, A=3, B=5)
      D.DE<-multiloc.LD(df=offspring_map, A=4, B=5)
      
      D.3<-p.ABC - p.A*D.BC - p.B*D.AC - p.C*D.AB
      
      D.4<-
      p.ABCD - p.A*D.BCD - p.B*D.ACD - p.C*D.ABD - p.D*D.ABC -
      p.A*p.B*D.CD - p.A*p.C*D.BD - p.A*p.D*D.BC - 
        p.B*p.C*D.AD - p.B*p.D*D.AC -
        p.C*p.D*D.AB - 
      D.AB*(D.CD)-D.AC*(D.BD)-D.AD*(D.BC)
      
      D.5<-
        p.ABCDE - p.A*D.BCDE - p.B*D.ACDE - p.C*D.ABDE - p.D*D.ABCE - p.E*D.ABCD
      p.A*p.B*D.CDE - p.A*p.C*D.BDE - p.A*p.D*D.BCE - p.A*p.E*D.BCD -
        p.B*p.C*D.ADE - p.B*p.D*D.ACE - p.B*p.E*D.ACD -
        p.C*p.D*D.ABE - p.C*p.E*D.ABD -
        p.D*p.E*D.ABC -
      p.A*p.B*p.C*D.DE - p.A*p.B*p.D*D.CE - p.A*p.B*p.E*D.CD - p.A*p.C*p.D*D.BE - p.A*p.C*p.E*D.BD - p.A*p.D*p.E*D.BC -
        p.B*p.C*p.D*D.AE - p.B*p.D*p.E*D.AC - p.B*p.C*p.E*D.AD -
        p.C*p.D*p.E*D.AB -
      D.AB*(D.CD+D.CE+D.DE)-D.AC*(D.BD+D.BE+D.DE)-D.AD*(D.BC+D.BE+D.CE)-D.AE*(D.BC+D.BD+D.CD)-D.BC*(D.DE)-D.BD*(D.CE)-D.BE*(D.CD)
      
      if (i == 1){LD.5[1,o]<-D.5} else {LD.5[(1+(i/50)),o]<-D.5}
      if (i == 1){LD.4[1,o]<-D.4} else {LD.4[(1+(i/50)),o]<-D.4}
      if (i == 1){LD.3[1,o]<-D.3} else {LD.3[(1+(i/50)),o]<-D.3}

    } #end for loop over generations
    
  } #end for loop over runs
  
  LD.5<- mutate(LD.5, avg=rowMeans(LD.5[2:5])) %>%
    abs() %>% 
    bind_cols(plot.df, .)  
  LD.4<- mutate(LD.4, avg=rowMeans(LD.4[2:5])) %>%
    abs() %>% 
    bind_cols(plot.df, .)  
  LD.3<- mutate(LD.3, avg=rowMeans(LD.3[2:5])) %>%
    abs() %>% 
    bind_cols(plot.df, .)  

 LD.avg[,1]<-LD.5$avg
 LD.avg[,2]<-LD.4$avg
 LD.avg[,3]<-LD.3$avg
 LD.avg<-bind_cols(plot.df, LD.avg)

} #end loop over parameter sets

library(ggplot2)
LDavg.plot<-ggplot(data=LD.avg, aes(x=gen))+geom_line(aes(y=V1), colour="red")+
  geom_line(aes(y=V2), colour="blue")+
  geom_line(aes(y=V3))+
  #geom_line(aes(y=V4))+
  #geom_line(aes(y=V5))+
  #geom_line(aes(y=avg), size=2, colour="red")+
  xlab("Generation")+ylab("Disequilibrium coefficient D")+ggtitle(paste("Average - Parameter set", x, sep=" "))+theme_bw()+ylim(0,1)
plot(LDavg.plot)

LD.plot<-ggplot(data=LD.5, aes(x=gen))+geom_line(aes(y=V1))+
  geom_line(aes(y=V2))+
  geom_line(aes(y=V3))+
  geom_line(aes(y=V4))+
  geom_line(aes(y=V5))+
  geom_line(aes(y=avg), size=2, colour="blue")+
  xlab("Generation")+ylab("Disequilibrium coefficient D")+ggtitle(paste("Parameter set", x, sep=" "))+theme_bw()+ylim(0,1)
plot(LD.plot)
