library(tidyr)
library(dplyr)

###Create and set directory for saved files###
dir.create(path=paste(getwd(), paste("para_set", para_set, sep="_"), paste("model_run", r, sep="_"), sep="/"), recursive=T)
setwd(paste(getwd(), paste("para_set", para_set, sep="_"), paste("model_run", r, sep="_"), sep="/"))

#create neighbourhood size dataframe
sum.N<-as.data.frame(matrix(nrow=no_gen,ncol=3))
names(sum.N)<-c("N.size","var","d")

#create selfing rate dataframe
self.df<-as.data.frame(matrix(nrow=no_gen, ncol=1))
names(self.df)<-"self.rate"

#Step 1) Create table to hold parental information
parents<-as.data.frame(matrix(nrow=pop_size, ncol=(3+(2*FT_loci)+3)))

names(parents)[1:3]<-c("Mother","Father","FLday")
for (n in 1:FT_loci)
{
  names(parents)[4+2*(n-1)]<-paste("loc",n,"a",sep="")
  names(parents)[5+2*(n-1)]<-paste("loc",n,"b",sep="")
}
names(parents)[(ncol(parents)-3): ncol(parents)]<-c("FTnoise","position","X_pos","Y_pos")

#Step 2) Generate flowering time onsets
#2a Give each individual its flowering time genotype (2 haploid complements)
for (n in 1:FT_loci)
{
  parents[,4+2*(n-1)]<-sample(FT_alleles, size=pop_size, replace=T)
  parents[,5+2*(n-1)]<-sample(FT_alleles, size=pop_size, replace=T)
}

#2b) Add environmental variance to flowering time
parents$FTnoise<-rnorm(n=pop_size, mean=0, sd=fl_noise)

#2c) Add genertic and environmetal to some positive baseline to get each individual's flowering time
parents$FLday<-floor(fl_mean + rowSums(parents[,4:(3+2*FT_loci)])+parents$FTnoise)


#Step 3) Give each individual a position in space
#3a) From the matrix of available positions in space, develop a list of positions from which to sample
for (y in 1:ncol(space_20by20_sq_grid)) {
  if (y ==1 ) {
    pos_vector<-c(as.character(space_20by20_sq_grid[,y])) }
  else {
    pos_vector<-c(pos_vector, as.character(space_20by20_sq_grid[,y])) 
  }
}

#3b) Give everyone a position from the set of possible positions, and then extract their x and y coordinates
parents$position<-(pos_vector)	
parents$X_pos<-sapply(X=parents$position, FUN=function(X) {
  as.numeric(strsplit(X, split=",")[[1]][1])
})
parents$Y_pos<-sapply(X=parents$position, FUN=function(X) {
  as.numeric(strsplit(X, split=",")[[1]][2])
})	

###Run offspring generations(generation 2 through no_gen)###
for (g in 1:no_gen){
  
  
  #Step 2) Fill in flowering schedule
  #2a) Set days of flowering schedule
  days<-c(min(parents$FLday):(max(parents$FLday)+(duration-1)))
  
  #2b) Fill in flowering schedule using parental generation abundances
  flowers<-as.data.frame(matrix(nrow=pop_size, ncol=length(days)))
  names(flowers)<-paste("d", days, sep="")
  for (i in 1:pop_size){
    flowers[i, which(days==parents$FLday[i]):(which(days==parents$FLday[i])+duration-1)]<-fl_seq
  }
  flowers[is.na(flowers)]<-0
  
  #2c) Account for inbreeding depression if applicable
  if (LA_inbreeding_flowers == "yes") {
    for (i in 1:nrow(flowers)){
      flowers[i,]<-flowers[i,]*(fl_ID_weight)*(1-parents$ID_coef[i])
    }
  }
  
  #Step 3) Create mating opportunity matrix based on phenology
  if (phenology== "yes") {
    #3a) Set maternal matrix as the proportion of flowers made by a single individual relative to those made by the total population (both over entire season)
    mat<-as.matrix(flowers/sum(flowers))
    #3b) Set paternal matrix as the proportion of flowers made by a single individual relative to the entire population
    pat<-matrix(nrow=nrow(flowers), ncol=ncol(flowers))
    for (t in 1:ncol(pat))	{	
      pat[,t]<-flowers[,t]/sum(flowers[,t])
    }
    
    #3c) Fix pat matrix to deal with days where NO flowers are produced in the population (this removes NaN)
    if (any(colSums(flowers)==0)){
      pat[,which(colSums(flowers)==0)]<-0}
    
    #3d) Place moms as rows and dads as columns
    mat_opp<-mat %*% t(pat)
    
  } else if (phenology=="no") {
    mat_opp<-matrix(nrow=pop_size, ncol=pop_size)
    mat_opp[,]<-1/(pop_size^2)
    
  }
  
  #Step 4) Develop pairwise space (Euclidean distance) matrix and run through the pollen dispersal function
  #4a) Calculate probability based on distance
  if (torus == "no") {
    distance1<-matrix(nrow=pop_size, ncol=pop_size)
    distance2<-matrix(nrow=pop_size, ncol=pop_size)
    for (m in 1:pop_size) {
      for (p in 1:pop_size) {
        distance1[m,p]<-sqrt((parents$X_pos[m]-parents$X_pos[p])^2 + (parents$Y_pos[m] - parents$Y_pos[p])^2)	
      }
    }
    distance1<-distance1*m_per_unit
    distance2[,]<-0
    
    #4a1) Run through pollen dispersal function
    if (space_pollen=="yes"){
      distance_for_mating1<-distance1
      distance_for_mating2<-distance2
    } else if (space_pollen=="no") {
      distance_for_mating1<-matrix(nrow=pop_size, ncol=pop_size)
      distance_for_mating1[,]<-1
      distance_for_mating2[,]<-0
    }
    
    pollen_dispersal1<-(exp(-(seed_lambda*distance_for_mating1)))
    pollen_dispersal2<-(exp(-(seed_lambda*distance_for_mating2)))
    
  } else if (torus == "yes") {
    distance1<-matrix(nrow=pop_size, ncol=pop_size)
    distance2<-matrix(nrow=pop_size, ncol=pop_size)
    for (m in 1:pop_size) {
      for (p in 1:pop_size) {
        distance1[m,p]<-sqrt((parents$X_pos[m]-parents$X_pos[p])^2 + (parents$Y_pos[m] - parents$Y_pos[p])^2)
        distance2[m,p]<-sqrt(((sqrt(pop_size)+1)-abs(parents$X_pos[m]-parents$X_pos[p]))^2 + ((sqrt(pop_size)+1)-abs(parents$Y_pos[m]-parents$Y_pos[p]))^2)
        
      }
    }
    distance1<-distance1*m_per_unit
    distance2<-distance2*m_per_unit
    
    #4a1) Run through pollen dispersal function
    if (space_pollen=="yes"){
      distance_for_mating1<-distance1
      distance_for_mating2<-distance2
    } else if (space_pollen=="no") {
      distance_for_mating1<-matrix(nrow=pop_size, ncol=pop_size)
      distance_for_mating1[,]<-1
      distance_for_mating2<-matrix(nrow=pop_size, ncol=pop_size)
      distance_for_mating2[,]<-1
    }
    pollen_dispersal1<-(exp(-(seed_lambda*distance_for_mating1)))
    pollen_dispersal2<-(exp(-(seed_lambda*distance_for_mating2)))
  } #end of if torus "yes"
  
  #Step 5) Recalibrate mating opportunity matrix to account for space
  mat_opp_adj<-mat_opp * (pollen_dispersal1 + pollen_dispersal2)
  mat_opp_adj<-mat_opp_adj/sum(mat_opp_adj)
  
  #Step 9) Generate offspring
  offspring<-as.data.frame(matrix(nrow=pop_size, ncol=ncol(parents)))
  names(offspring)<-names(parents)
  
  #9a) Designate spatial locations, assuming each space will be occupied by a single plant in the next generation
  offspring$position<-parents$position
  offspring$X_pos<-parents$X_pos
  offspring$Y_pos<-parents$Y_pos	
  
  #9b) Start loop over offspring
  for (o in 1:pop_size) {
    
    #9b1: Calculate dispersal of seeds to location of offspring[o] as a weight on relative fitness
    if (space_seed == "yes"){
      seed_dispersal1<-(exp(-(seed_lambda*distance1[o,])))} else if (space_seed == "no") {
        seed_dispersal1<-rep(1, times=pop_size)}
    if (space_seed == "yes"){	  
      seed_dispersal2<-(exp(-(seed_lambda*distance2[o,])))} else if (space_seed == "no") {
        seed_dispersal2<-rep(1, times=pop_size)}
    seed_dispersal<-seed_dispersal1+seed_dispersal2
    
    #9b2) Choose a mom based on relative fitness and chance of getting a seed to the specified location
    ##NOTE- Edge plants at an inherent disadvantage if a torus is not used.
    mom<-sample(c(1:pop_size), size=1, prob=(1*seed_dispersal)/sum(1*seed_dispersal))
    
    #9b3) If there are no sires for chosen mom, choose another mom
    while (rowSums(mat_opp_adj)[mom]==0)	{mom<-sample(c(1:pop_size), size=1, prob=1)}
    
    #9b4) Choose sire for selected mom
    pat_prob<-mat_opp_adj[mom,]/rowSums(mat_opp_adj)[mom]
    dad<-sample(x=c(1:pop_size), size=1, prob=pat_prob)
    
    #9b5) Record mom and dad
    offspring$Mother[o]<-mom
    offspring$Father[o]<-dad		
    
    #9b10) Mendellian inheritance of flowering time loci 
    for (n in 1:(FT_loci)) {
      offspring[o,4 + 2*(n-1)]<-sample(x=parents[mom,(4+2*(n-1)):(5 + 2*(n-1))], size=1)
      offspring[o,5 + 2*(n-1)]<-sample(x=parents[dad,(4+2*(n-1)):(5 + 2*(n-1))], size=1)
    }	#close loop over the flowering time loci
    
    
  }	#close loop over offspring
  
  #9c) Add environmental variance to flowering time, and sum flwg time columns to calculate flowering time	
  offspring$FTnoise<-rnorm(n=pop_size, mean=0, sd=fl_noise)
  offspring$FLday<-floor(fl_mean + rowSums(offspring[,4:(3+2*FT_loci)]) + offspring$FTnoise) 
  
  NH<-subset.data.frame(offspring, select=c(FLday, Mother, Father, X_pos, Y_pos))
  
  #set up starting flowering distribution
  fl_df<-as.data.frame(fl_seq)
  fl_df<-mutate(fl_df, day=c(-7:7))
  
  #fill in mom coordinates
  for (m in 1:pop_size) {
    mom.x<-NH$Mother[m]
    NH$M.x[m]<-sapply(X=pos_vector[mom.x], FUN=function(X) {
      as.numeric(strsplit(X, split=",")[[1]][1])
    })
    NH$M.y[m]<-sapply(X=pos_vector[mom.x], FUN=function(X) {
      as.numeric(strsplit(X, split=",")[[1]][2])
    })
    
  }
  
  #fill in dad coorindates
  for (p in 1:pop_size) {
    dad.x<-NH$Father[p]
    NH$P.x[p]<-sapply(X=pos_vector[dad.x], FUN=function(X) {
      as.numeric(strsplit(X, split=",")[[1]][1])
    })
    NH$P.y[p]<-sapply(X=pos_vector[dad.x], FUN=function(X) {
      as.numeric(strsplit(X, split=",")[[1]][2])
    })
  }
  
  #score selfing
  for (k in 1:pop_size) {
    NH$self[k]<-if (NH$Mother[k] == NH$Father[k]) {1} else if (NH$Mother[k] != NH$Father[k]) {0}
  }
  
  #calculate Euclidean distance for each parent
  for (m in 1:pop_size){
    NH$mom.dist[m]<-sqrt((NH$M.x[m]-NH$X_pos[m])^2+(NH$M.y[m]-NH$Y_pos[m])^2)
  }
  for (m in 1:pop_size){
    NH$dad.dist[m]<-sqrt((NH$P.x[m]-NH$X_pos[m])^2+(NH$P.y[m]-NH$Y_pos[m])^2)
  }
  
  #calculate the mean distance
  mean.dist<-sum(mean(NH$mom.dist), mean(NH$dad.dist))/2
  
  #calculate variance for each offspring
  NH<-mutate(NH, ind.var=(mom.dist-mean.dist)^2+(dad.dist-mean.dist)^2)
  
  #calculate population variance
  var<-sum(NH$ind.var)/(2*pop_size)
  
  #sample offspring
  off.samp<-sample(c(1:pop_size), size=100, replace=FALSE)
  off.df<-as.data.frame(matrix(nrow=100,ncol=1))
  #start loop over sampled offspring
  for (w in off.samp){
    
    if (phenology == "no") {off.df[match(c(w),off.samp),]<- 1/(m_per_unit*m_per_unit)} else if (phenology == "yes") {
      
      ind.df<-as.data.frame(matrix(nrow=1, ncol=1))
      ind.df.t<-as.data.frame(matrix(nrow=1, ncol=1))
      a<-NH$X_pos[w]
      b<-NH$Y_pos[w]
      FLday.o<-NH$FLday[w]
      p.x<-round(sqrt(var))+a #round so can translate to cartesian coordinates
      if(p.x > sqrt(pop_size)) {p.x <- sqrt(pop_size)}
      n.x<-a-round(sqrt(var))
      if(n.x < 1) {n.x <- 1}
      for (x in n.x:p.x){
        p.y<-if (var-(x-a)^2 < 0) {b} else {round(sqrt(var-(x-a)^2))+b}
        if(p.y > sqrt(pop_size)) {p.y <- sqrt(pop_size)}
        n.y <- if (var-(x-a)^2 < 0) {b} else {b-round(sqrt(var-(x-a)^2))}
        if(n.y < 1) {n.y <- 1}
        for (z in n.y:p.y) {
          NH.temp<-filter(NH, X_pos==x, Y_pos==z)
          #set up integral for overlap in flowering schedule
          FLday.t<-NH.temp$FLday[1]-FLday.o #paired indivdual relative flowering day
          fl_df.t<-fl_df
          fl_df.t[2]<-fl_df.t[2]+FLday.t
          fl_df.n <- as.data.frame((min(fl_df[2],fl_df.t[2]):max(fl_df[2],fl_df.t[2])))
          names(fl_df.n)<-"day"
          fl_df.n<-full_join(fl_df.n,fl_df,by="day")
          fl_df.n<-full_join(fl_df.n,fl_df.t,by="day")
          fl_df.n[is.na(fl_df.n)] <- 0
          fl_df.n<-mutate(fl_df.n, AUC=min(fl_seq.x,fl_seq.y))
          for (i in 1:length(fl_df.n$day)){
            fl_df.n$AUC[i]<-sapply(X=fl_df.n$fl_seq.x[i], Y=fl_df.n$fl_seq.y[i], FUN=function(X,Y) {
              min(X,Y)
            })
          } #end i loop
          dif.t<-sum(fl_df.n$AUC)/sum(fl_df.n$fl_seq.x)
          ind.df.t[1,1]<-dif.t
          ind.df<-bind_rows(ind.df, ind.df.t)
        } #next y position
      } #next x position
      off.df[match(c(w),off.samp),]<-sum(ind.df[,1], na.rm=TRUE)/length(ind.df[,1])
    } #end if/else
  } #next offspring
  
  write.csv(off.df, paste("off.df.", g, ".csv", sep=""))
  write.csv(NH, paste("NH.", g, ".csv", sep=""))
  
  sum.N[g,1]<-(4*3.14*var*mean(off.df$V1, na.rm=TRUE))
  sum.N[g,2]<-var
  sum.N[g,3]<-mean(off.df$V1, na.rm=TRUE)
  
  self.df[g,1]<-mean(NH$self)
  
  write.csv(sum.N,
            "sum.N.csv", row.names=F)
  write.csv(self.df, 
            "self.df.csv", row.names=F)
  
  if (g<no_gen){
    offspring->parents}
  
}