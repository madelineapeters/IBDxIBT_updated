#Set parameters
fl_mean<-100 #mean flowering day in generation 1
fl_seq<-c(2, 3, 5, 8, 11, 13, 15, 17, 15, 13, 11, 8, 5, 3, 2) #flowers produced per day
duration<-length(fl_seq) #length of flowering period
no_run<-3 #number of runs
seq<-c(1, seq(from=10, to=100, by=10)) #generations to test
seed_lambda<-2

#Start loop over parameter sets
for (e in 88:89) {
  
  #set metres per unit based on parameter set
  m_per_unit<-if ((e == 86) | (e == 87)) {0.5} else {1}
  
  #Start loop over runs
  for (r in 1:no_run) {
    
    S.N.df<-as.data.frame(matrix(nrow=length(seq), ncol=3))
    names(S.N.df)<-c("gen", "S", "num.pool")
    S.N.df[1]<-seq
    
    #Start loop over generations
    for (g in seq) {
      
      #Set working directory
      start_wd<-(paste("~/498/Neighbourhood genetics", paste("para_set", e, sep="_"), paste("model_run", r, sep="_"), sep="/"))
      setwd(start_wd)
      
      #read in generation dataframe
      offspring_map <- read.csv(paste(getwd(), paste("NH", g, "csv", sep="."), sep="/"))

      #fill in flowering schedule
      days<-c(min(offspring_map$FLday):(max(offspring_map$FLday)+(duration-1)))

      #set up matrix to fill in flowers produced per individual per day over the entire flowering season      
      flowers<-as.data.frame(matrix(nrow=pop_size, ncol=length(days)))
      names(flowers)<-paste("d", days, sep="")

      #start loop over individuals to fill in flowers matrix
      for (i in 1:pop_size){
        flowers[i, which(days==offspring_map$FLday[i]):(which(days==offspring_map$FLday[i])+duration-1)]<-fl_seq
      } #end loop over individuals

      flowers[is.na(flowers)]<-0 #replaces NAs with 0

      #set up M matrix: proportion of all flowers over season on mom 'n' on day 'd'
      mat<-as.matrix(flowers/sum(flowers))

      #set up P matrix: proportion of flowers open on dad 'n' relative to all flowers open on day 'd' (not across the entire season)
      pat<-matrix(nrow=nrow(flowers), ncol=ncol(flowers))
      
      for (t in 1:ncol(pat))	{	
        pat[,t]<-flowers[,t]/sum(flowers[,t])
      }

      #fix P matrix to deal with days where NO flowers are produced in the population (this removes NaN)
      if (any(colSums(flowers)==0)){pat[,which(colSums(flowers)==0)]<-0}


      #create mating opportunity matrix (places moms as rows and dads as columns)
      mat_opp<-mat %*% t(pat)

      #create dispersal mating matrix (for a torus)
      distance1<-matrix(nrow=pop_size, ncol=pop_size)
      distance2<-matrix(nrow=pop_size, ncol=pop_size)
      for (m in 1:pop_size) {
        for (p in 1:pop_size) {
          distance1[m,p]<-sqrt((offspring_map$X_pos[m]-offspring_map$X_pos[p])^2 + (offspring_map$Y_pos[m] - offspring_map$Y_pos[p])^2)
          distance2[m,p]<-sqrt(((sqrt(pop_size)+1)-abs(offspring_map$X_pos[m]-offspring_map$X_pos[p]))^2 + ((sqrt(pop_size)+1)-abs(offspring_map$Y_pos[m]-offspring_map$Y_pos[p]))^2)
        }
      }
      
      distance1<-distance1*m_per_unit
      distance2<-distance2*m_per_unit

      #run through pollen dispersal function
      distance_for_mating1<-distance1
      distance_for_mating2<-distance2
      pollen_dispersal1<-(exp(-(seed_lambda*distance_for_mating1)))
      pollen_dispersal2<-(exp(-(seed_lambda*distance_for_mating2)))
    
      #recalibrate mating opportunity matrix to account for space
      mat_opp_adj<-mat_opp * (pollen_dispersal1 + pollen_dispersal2)
      mat_opp_adj<-mat_opp_adj/sum(mat_opp_adj)

      #eigenvalue extration and calculation of S statistic
      EV<-eigen(mat_opp_adj, sym=FALSE)
      ev.vector<-EV$values
      ev.1<-ev.vector[1]
      ev.sum<-sum(ev.vector)
      S<-(ev.1/sum(ev.vector))
      num_pools<-ev.sum/ev.1
      
      #store output
      if (g == 1) {
        S.N.df[g,2]<-S
        S.N.df[g,3]<-num_pools
      } else {
        S.N.df[(g/10)+1,2]<-S
        S.N.df[(g/10)+1,3]<-num_pools
      }

    } #end loop over g
    
    write.csv(S.N.df, paste("S.N.df", r, "csv", sep="."))
    
  } #end loop over runs
} #end loop over parameter sets
