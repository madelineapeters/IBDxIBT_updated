
for (e in starting.set1:ending.set2) {

  #start loop over runs
  for (r in 1:no_runs){
    
    #set working directory
    start_wd<-(paste(base_wd, sub_wd, paste("para_set", e, sep="_"), paste("model_run_", r, sep=""), sep="/"))
    setwd(start_wd)
    
    #create dataframes to hold stats outputs
    Fst_tab<-matrix(nrow=11, ncol=no_runs)
    for (g in 1:no_runs){
      names(Fst_tab)[g]<-paste("run",g, sep="_")
    }
    
    #start loop over generations
    for (i in c(1, seq(from = 50, to = 500, by = 50))) {
      #read in generation dataframe
      offspring <- read.csv(paste(getwd(), paste("offspring_map", i, sep="_"), sep="/"))
      
      #start loop to set up simplified neutral locus matrix
      for (j in 1:no_neutral)
      {
        names(offspring)[10+(2*5)+(2*(j-1))]<-paste("neut",j,"a",sep="")
        names(offspring)[11+(2*5)+(2*(j-1))]<-paste("neut",j,"b",sep="")
      }
      
      #create dataframe to store numerical allele scores
      offspring_map<-as.data.frame(matrix(nrow=pop_size, ncol=2*no_neutral+no_neutral))
      
      #loop over population and number of neutral alleles to fill in allele score dataframe
      for (k in 1:pop_size){
        for (j in 1:no_neutral)
        {
          offspring_map[k,1+(2*(j-1))]<-if (offspring[k,10+(2*5)+(2*(j-1))]=="d") {
            0} else {1}
          offspring_map[k,2+(2*(j-1))]<-if (offspring[k,11+(2*5)+(2*(j-1))]=="d") {
            0} else {1}
          
          #fill in columns that store loci scores
          offspring_map[k,2*no_neutral+(j)]<-sum(offspring_map[k,1+(2*(j-1))],offspring_map[k,2+(2*(j-1))])
        }
      }
      
      #after filling in columns with loci scores, remove allele score columns 
      offspring_map<-subset.data.frame(offspring_map, select= (2*no_neutral+1):(2*no_neutral+no_neutral))
      
      ####################################################################################
      
      #set up pairwise Fst matrix
      pair_matrix<-matrix(nrow=pop_size, ncol=pop_size)
      
      #fill in total distance matrix
      for (m in 1:pop_size) {
        for (p in 1:pop_size){
          #To estimate FST, take the following steps:
          
          #Find the allele frequencies for each subpopulation.
          Pm<-sum(offspring_map[m,])/(2*no_neutral)
          Pp<-sum(offspring_map[p,])/(2*no_neutral)
          Qm<-1-Pm
          Qp<-1-Pp
          
          #Find the average allele frequencies for the total population.
          Ptotal<-sum(Pm,Pp)/2
          Qtotal<-1-Ptotal
          
          #Calculate the heterozygosity (2pq) for each subpopulation.
          Hm<-2*Pm*Qm
          Hp<-2*Pp*Qp
          
          #Calculate the average of these subpopulation heterozygosities. This is HS.
          HS<-sum(Hm, Hp)/2
          
          #Calculate the heterozygosity based on the total population allele frequencies. This is HT.
          HT<-(2*(Ptotal)*Qtotal)
          
          #Finally, calculate FST=(HT-HS)/HT.
          Fst<-(HT-HS)/HT
          
          #Store Fst value
          pair_matrix[m,p]<-Fst
        } #next p
      } #next m
      
      #find average pairwise Fst for across whole matix
      Fst_avg<-mean(pair_matrix, na.rm = TRUE)
      
      #Store average pairwise Fst in data frame
      if (i == 1){
        Fst_tab[i,r]<-Fst_avg
      } else {
        Fst_tab[((i/50)+1),r]<-Fst_avg
      }
      
    } #end for loop over generations
    
  } #end for loop over runs
  
  #Calculate average generational Fst values across runs
  Fst_tab<-as.data.frame(Fst_tab)
  Fst_run_avg<-matrix(nrow=11, ncol=1)
  Fst_run_avg<-rowMeans(Fst_tab, na.rm = TRUE, dims = 1)
  
  #Reset working directory
  start_wd<-(paste(base_wd, "final sets", paste("para_set", e, sep="_"), sep="/"))
  setwd(start_wd)
  
  #Write .csv file for Fst_tab	
  write.csv(Fst_run_avg, paste("Fst.", e, ".csv", sep=""))
  
} #end for loop over parameter sets
