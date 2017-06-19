#load required packages
library(tidyr)
library(dplyr)
library(doParallel)

cl<- makeCluster(3)
registerDoParallel(cl)

source(file=paste("~/498/stats parameters.R"))

foreach (e = starting.set2:ending.set2) %dopar% {
  
  source(file=paste("~/498/stats parameters.R"))
  
  #make data frame for average values
  run.avg.df<-as.data.frame(matrix(nrow=11, ncol=1))
  avg.df<-as.data.frame(matrix(nrow=11, ncol=no_runs))
  for (r in 1:no_runs) {
    
    #make allele data frame
    allele.df<-as.data.frame(matrix(nrow=11, ncol=no_neutral))
    
    #set working directory
    start_wd<-paste(base_wd,sub_wd, paste("para_set", e, sep="_"), paste("model_run", r, sep="_"), sep="/")
    setwd(start_wd)
    
    #start loop over generations
    for (i in c(1, seq(from = 50, to = 500, by = 50))) {
      #read in generation dataframe
      offspring <- read.csv(paste(getwd(), paste("offspring_map", i, sep="_"), sep="/"))
      
      #start loop to set up simplified neutral locus matrix
     # for (j in 1:no_neutral)
      #{
       # names(offspring)[10+(2*5)+(2*(j-1))]<-paste("neut",j,"a",sep="")
        #names(offspring)[11+(2*5)+(2*(j-1))]<-paste("neut",j,"b",sep="")
     # }
      
      #create dataframe to store numerical allele scores
      offspring_map<-as.data.frame(matrix(nrow=pop_size, ncol=2*no_neutral+no_neutral))
      
      #loop over population and number of neutral alleles to fill in allele score dataframe
      for (k in 1:pop_size){
        for (j in 1:no_neutral)
        {
          offspring_map[k,1+(2*(j-1))]<-if (offspring[k,10+(2*(j-1))]=="1") {
            0} else if (offspring[k,10+(2*(j-1))]=="0.5") {0} else {1}
          offspring_map[k,2+(2*(j-1))]<-if (offspring[k,11+(2*(j-1))]=="1") {
            0} else if (offspring[k,11+(2*(j-1))]=="0.5") {0} else {1}
          
          #fill in columns that store loci scores
          offspring_map[k,2*no_neutral+(j)]<-sum(offspring_map[k,1+(2*(j-1))],offspring_map[k,2+(2*(j-1))])
        }
      }
      
      #after filling in columns with loci scores, remove allele score columns 
      offspring_map<-subset.data.frame(offspring_map, select= (2*no_neutral+1):(2*no_neutral+no_neutral))
      
      #start loop over neutral loci to run and record Mantel test  
      
      for (n in 1:no_neutral){
        
        p<-sum(offspring_map[,n])/(2*pop_size)
        var_p<-p*(1-p)
        dev.p<-abs(p-0.5)
        
        if (i == 1){
          allele.df[1,n]<-var_p
        } else {
          allele.df[((i/50)+1),n]<-var_p
        }
        
      } #end loop over n
      
    } #end for loop over generations
    write.csv(allele.df, paste("FL.var.", r, ".csv", sep=""))
    
    #make average data frame and plot
    allele.df$avg<-rowMeans(allele.df[1:no_neutral])
    avg.df[,r]<-allele.df$avg

  } #end for loop over runs
  
  start_wd<-paste(base_wd, sub_wd, paste("para_set", e, sep="_"), sep="/")
  setwd(start_wd)
  
  runavg.df<-rowMeans(avg.df)
  run.avg.df<-as.data.frame(run.avg.df)
  write.csv(avg.df, paste("Sep.FL.var.", e, ".csv", sep=""))
  write.csv(run.avg.df, paste("Avg.FL.var.", e, ".csv", sep=""))

}