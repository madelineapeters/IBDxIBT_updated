
#load required packages
library(vegan)
library(tidyr)
library(dplyr)

#set up dataframe to store corr values
corr.df<-as.data.frame(matrix(nrow=20, ncol=6))
names(corr.df)<-c("Dist",1:9)
corr.df[2:6]<-0

#start loop
for (e in c(54, 56, 58)) {
  
############Set-up############ 
  #re-read in stats parameters for each core
  source(file=paste("~/498/stats parameters.R"))
  
  #select run
  r<-1
  
  #select generation
  i<-500
  
  #select trait
  trait<-"neutral"
  
  #set working directory
  start_wd<-(paste(base_wd, sub_wd, paste("para_set", e, sep="_"), paste("model_run_", r, sep=""), sep="/"))
  setwd(start_wd)
  
  #read in generation dataframe
  offspring <- read.csv(paste(start_wd, paste("offspring_map", i, sep="_"), sep="/"))
  
  #set up and fill in geographic distance matrix
  distance1<-matrix(nrow=pop_size, ncol=pop_size)
  distance2<-matrix(nrow=pop_size, ncol=pop_size)
  
  for (m in 1:pop_size) {
    for (p in 1:pop_size) {
      distance1[m,p]<-sqrt((offspring$X_pos[m]-offspring$X_pos[p])^2 + (offspring$Y_pos[m] - offspring$Y_pos[p])^2)
      distance2[m,p]<-sqrt(((sqrt(pop_size)+1)-abs(offspring$X_pos[m]-offspring$X_pos[p]))^2 + ((sqrt(pop_size)+1)-abs(offspring$Y_pos[m]-offspring$Y_pos[p]))^2)
    }
  }
  
  distance1<-distance1/sum(distance1, na.rm=TRUE)
  distance2<-distance2/sum(distance2, na.rm=TRUE)
  distance_matrix<-distance1+distance2
  distance_matrix<-distance_matrix/sum(distance_matrix, na.rm=TRUE)

############Trait matrix############  
  if (trait == "FLday") {
    #set up and fill in FL day matrix  
    trait_matrix<-matrix(nrow=pop_size, ncol=pop_size)
    
    for (m in 1:pop_size) {
      for (p in 1:pop_size) {
        trait_matrix[m,p]<-abs((offspring$FLday[m]-offspring$FLday[p])/FLmean)
      } #end p loop
    } #end m loop
    
    #run correlogram function
    trait.mantel<-mantel.correlog(D.eco=trait_matrix, D.geo=distance_matrix, nperm=1000)
    
    #fill in correlogram dataframe)
    if (e == 54 || 55) {
      for (s in 9:17){
        trait.string<-capture.output(trait.mantel)[s]
        trait.string<-substr(trait.string, 10,19)
        corr.df[(s-8),1]<-trait.string
      } #end for s loop
    } #end if loop
    for (s in 9:17){
      trait.string<-capture.output(trait.mantel)[s]
      trait.string<-substr(trait.string, 33, 44)
      if (e == 54 || 56 || 58){
        corr.df[(s-8),((e/2)-25)]<-trait.string
      } else {
        corr.df[(s-8),(((e-1)/2)-25)]<-trait.string
      } #end else
    } #end for s loop 
    
  } else if (trait == "neutral") {
    
    #start loop to set up simplified neutral locus matrix
    for (j in 1:no_neutral){
      names(offspring)[10+(2*5)+(2*(j-1))]<-paste("neut",j,"a",sep="")
      names(offspring)[11+(2*5)+(2*(j-1))]<-paste("neut",j,"b",sep="")
    }
    
    #create dataframe to store numerical allele scores
    offspring_map<-as.data.frame(matrix(nrow=pop_size, ncol=2*no_neutral+no_neutral))
    
    #loop over population and number of neutral alleles to fill in allele score dataframe
    for (k in 1:pop_size) {
      for (j in 1:no_neutral) {
        offspring_map[k,1+(2*(j-1))]<-if (offspring[k,10+(2*5)+(2*(j-1))]=="d") {
          0} else {1}
        offspring_map[k,2+(2*(j-1))]<-if (offspring[k,11+(2*5)+(2*(j-1))]=="d") {
          0} else {1}
        
        #fill in columns that store loci scores
        offspring_map[k,2*no_neutral+(j)]<-sum(offspring_map[k,1+(2*(j-1))],offspring_map[k,2+(2*(j-1))])
      } #end j loop
    } #end k loop
    
    #after filling in columns with loci scores, remove allele score columns 
    offspring_map<-subset.data.frame(offspring_map, select= (2*no_neutral+1):(2*no_neutral+no_neutral))
    
    #start loop over neutral loci to set-up genetic distance matrix and record output 
    for (n in 1:no_neutral){
      
      #set up matrix to store per-locus genetic distance 
      trait_matrix<-matrix(nrow=pop_size, ncol=pop_size)
      
      #fill in genetic distance matrix
      for (m in 1:pop_size) {
        for (p in 1:pop_size) {
          trait_matrix[m,p]<-abs(offspring_map[m,(1+(n-1))]-offspring_map[p,(1+(n-1))])
        } #end p loop
      } #end m loop
      
      #run correlogram function
      if (sum(trait_matrix) > 0) {
        
        #run correlog
        trait.mantel<-mantel.correlog(D.eco=trait_matrix, D.geo=distance_matrix, nperm=1000)
        
        #fill in correlogram dataframe)
        if (e == 54 || 55) {
          for (s in 9:17){
            trait.string<-capture.output(trait.mantel)[s]
            trait.string<-substr(trait.string, 10,19)
            corr.df[(s-8),1]<-trait.string
          } #end for s loop
        } #end if loop
        for (s in 9:17){
          trait.string<-capture.output(trait.mantel)[s]
          trait.string<-as.numeric(substr(trait.string, 33, 44))
          if (e == 54 || 56 || 58){
            corr.df[(s-8),((e/2)-25)]<-corr.df[(s-8),((e/2)-25)]+trait.string
            if (n == no_neutral) {corr.df[(s-8),(((e)/2)-25)]<-corr.df[(s-8),(((e)/2)-25)]/no_neutral}
          } else {
            corr.df[(s-8),(((e-1)/2)-25)]<-corr.df[(s-8),(((e-1)/2)-25)]+trait.string
            if (n == no_neutral) {corr.df[(s-8),(((e-1)/2)-25)]<-corr.df[(s-8),(((e-1)/2)-25)]/no_neutral}
          } #end else
        } #end for s loop 
      } #end if statement
    } #end n loop
  } else if (trait == "FLloci") {
    
      #create dataframe to store numerical allele scores
      offspring_map<-as.data.frame(matrix(nrow=pop_size, ncol=2*no_FL+no_FL))
    
      #loop over population and number of FL alleles to fill in allele score dataframe
      for (k in 1:pop_size){
        for (j in 1:5)
        {
          offspring_map[k,1+(2*(j-1))]<-if (offspring[k,(2*5)+(2*(j-1))]=="1") {
            0} else if (offspring[k,(2*5)+(2*(j-1))]=="0.5") {0} else {1}
          offspring_map[k,2+(2*(j-1))]<-if (offspring[k,(2*5)+1+(2*(j-1))]=="1") {
            0} else if (offspring[k,(2*5)+1+(2*(j-1))]=="0.5") {0} else {1}
          
          #fill in columns that store loci scores
          offspring_map[k,2*5+(j)]<-sum(offspring_map[k,1+(2*(j-1))],offspring_map[k,2+(2*(j-1))])
        } #end j loop
      } #end k loop
    
      #after filling in columns with loci scores, remove allele score columns 
      offspring_map<-subset.data.frame(offspring_map, select= (2*no_FL+1):(2*no_FL+no_FL))
      
      #start loop over FL loci to set-up genetic distance matrix and record output 
      for (n in 1:no_FL){
        
        #set up matrix to store per-locus genetic distance 
        trait_matrix<-matrix(nrow=pop_size, ncol=pop_size)
        
        #fill in genetic distance matrix
        for (m in 1:pop_size) {
          for (p in 1:pop_size) {
            trait_matrix[m,p]<-abs(offspring_map[m,(1+(n-1))]-offspring_map[p,(1+(n-1))])
          } #end p loop
        } #end m loop
        
        #run correlogram function
        if (sum(trait_matrix) > 0) {
          
          #run correlog
          trait.mantel<-mantel.correlog(D.eco=trait_matrix, D.geo=distance_matrix, nperm=1000)
          
          #fill in correlogram dataframe)
          if (e == 54 || 55) {
            for (s in 9:17){
              trait.string<-capture.output(trait.mantel)[s]
              trait.string<-substr(trait.string, 10,19)
              corr.df[(s-8),1]<-trait.string
            } #end for s loop
          } #end if loop
          for (s in 9:17){
            trait.string<-capture.output(trait.mantel)[s]
            trait.string<-as.numeric(substr(trait.string, 33, 44))
            if (e == 54 || 56 || 58){
              corr.df[(s-8),((e/2)-25)]<-corr.df[(s-8),((e/2)-25)]+trait.string
              if (n == no_FL) {corr.df[(s-8),(((e)/2)-25)]<-corr.df[(s-8),(((e)/2)-25)]/no_FL}
            } else {
              corr.df[(s-8),(((e-1)/2)-25)]<-corr.df[(s-8),(((e-1)/2)-25)]+trait.string
              if (n == no_FL) {corr.df[(s-8),(((e-1)/2)-25)]<-corr.df[(s-8),(((e-1)/2)-25)]/no_FL}
            } #end else
          } #end for s loop 
        } #end if statement
      } #end n loop
    }

} #end loop over parameter sets   

write.csv(corr.df, paste("corr.df", e, trait, "csv", sep="."))
