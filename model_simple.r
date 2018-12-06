library(dplyr)

dir.create(path=paste(getwd(), paste("para_set", x, sep="_"), sep="/"), recursive=T)

setwd(paste(getwd(), paste("para_set", para_set, sep="_"), sep="/"))

r=1
dir.create(path=paste(getwd(), paste("model_run", r, sep="_"), sep="/"), recursive=T)

#dir.create(path=paste(getwd(), paste("para_set", para_set, sep="_"), paste("model_run", r, sep="_"), sep="/"), recursive=T)
#directory is created using the para_set number read in from the parameter script and parameters csv and the model_run number from the commander script

setwd(paste(getwd(), paste("model_run", r, sep="_"), sep="/"))
#setwd(paste(getwd(), paste("para_set", para_set, sep="_"), paste("model_run", r, sep="_"), sep="/"))
#directory is set relative to the current working directory, which is set in the commander script

output<-as.data.frame(matrix(nrow=no_gen, ncol=11))
#no_gen is taken from the parameters script
names(output)<-c("gen","FTmean","FTvar","FTh2","IBmean","pA","pB","pC","FisA","FisB","FisC")

##TABLES HOLDING DATA WHERE THERE ARE SEVERAL VALUES PER GENERATION

#TABLE TO HOLD FLOWERING TIMES BY GENERATION
onset_tab<-(matrix(nrow=pop_size, ncol=no_gen))

#TABLE FOR A-LOCUS GENOTYPES BY GENERATION
Ageno_tab<-matrix(nrow=pop_size, ncol=no_gen)

#TABLE FOR B-LOCUS GENOTYPES BY GENERATION
Bgeno_tab<-matrix(nrow=pop_size, ncol=no_gen)

#TABLE FOR C-LOCUS GENOTYPES BY GENERATION
Cgeno_tab<-matrix(nrow=pop_size, ncol=no_gen)

sum.N<-as.data.frame(matrix(nrow=no_gen, ncol=3))
self.df<-as.data.frame(matrix(nrow=no_gen, ncol=1))

parents<-as.data.frame(matrix(nrow=pop_size, ncol=(10+2*FT_loci)+4+2*neut_loci))
#2 columns are needed for each locus so that each allele per locus can be stored

names(parents)[1:10]<-c("Mother","Father","IB_coef","FLday","A1","A2","B1","B2","C1","C2")
#A, B and C are three additional neutral loci used primarily for heatmap purposes
#A can also be set to be under selection using the parameter .csv

for (n in 1:FT_loci) {
  names(parents)[11+2*(n-1)]<-paste("loc",n,"a",sep="")
  names(parents)[12+2*(n-1)]<-paste("loc",n,"b",sep="")
}

names(parents)[(ncol(parents)-3): ncol(parents)]<-c("FTnoise","position","X_pos","Y_pos")

for (j in 1:neut_loci) {
  names(parents)[11+(2*FT_loci)+(2*(j-1))]<-paste("neut",j,"a",sep="")
  names(parents)[12+(2*FT_loci)+(2*(j-1))]<-paste("neut",j,"b",sep="")
}

#name assignments are set up so column positions are relative to the number of neutral loci specified

for (n in 1:FT_loci) {
  parents[,11+2*(n-1)]<-sample(FT_alleles, size=pop_size, replace=T)
  parents[,12+2*(n-1)]<-sample(FT_alleles, size=pop_size, replace=T)
}

#potential alleles are drawn from FT_alleles, which are specified in the parameters script
#to compare to neutral loci, best to make each FT locus biallelic

parents$FTnoise<-rnorm(n=pop_size, mean=0, sd=fl_noise)

parents$FLday<-floor(fl_mean + rowSums(parents[,11:(10+2*FT_loci)])+parents$FTnoise)

parents$A1<-sample(c("A","a"), size=pop_size, prob=c(freqA, 1-freqA), replace=T)
parents$A2<-sample(c("A","a"), size=pop_size, prob=c(freqA, 1-freqA), replace=T)

parents$B1<-sample(c("B","b"), size=pop_size, prob=c(freqB, 1-freqB), replace=T)
parents$B2<-sample(c("B","b"), size=pop_size, prob=c(freqB, 1-freqB), replace=T)

parents$C1<-sample(c("C","c"), size=pop_size, prob=c(freqC, 1-freqC), replace=T)
parents$C2<-sample(c("C","c"), size=pop_size, prob=c(freqC, 1-freqC), replace=T)

#starts population at HWE
#allele frequencies are typically starting at 0.5 (can be adjusted to alter starting conditions)

freqD<-0.5

for (j in 1:neut_loci)
{
  parents[,11+(2*FT_loci)+(2*(j-1))]<-sample(c("D","d"), size=pop_size, prob=c(freqD, 1-freqD), replace=T)
  parents[,12+(2*FT_loci)+(2*(j-1))]<-sample(c("D","d"), size=pop_size, prob=c(freqD, 1-freqD), replace=T)
}

#assign parents a space
comment.out = function(){
  for (y in 1:ncol(space_20by20_sq_grid)) {
    if (y == 1) {
      pos_vector<-c(as.character(space_20by20_sq_grid[,y]))
    } else {
      pos_vector<-c(pos_vector, as.character(space_20by20_sq_grid[,y]))
    }
  }
  
  parents$position<-(pos_vector)
  parents$X_pos<-sapply(X=parents$position, FUN=function(X) {
    as.numeric(strsplit(X, split=",")[[1]][1])
  })
  parents$Y_pos<-sapply(X=parents$position, FUN=function(X) {
    as.numeric(strsplit(X, split=",")[[1]][2])
  })
}

parents$X_pos = xy.coor[,1]
parents$Y_pos = xy.coor[,2]

#coancestry
coancestry_p<-matrix(nrow=pop_size, ncol=pop_size)

coancestry_p[,]<-0
for (i in 1:pop_size) {
  coancestry_p[i,i]<-0.5
}

parents$IB_coef<-0

for (g in 1:no_gen) {
  
  output$gen[g]<-g
  output$FTmean[g]<-mean(parents$FLday)
  output$FTvar[g]<-var(parents$FLday)
  output$FTh2[g]<-var(rowSums(parents[,11:(10+2*FT_loci)]))/var(parents$FLday)
  output$IBmean[g]<-mean(parents$IB_coef)
  output$pA[g]<-length(which(c(parents$A1, parents$A2)=="A"))/(2*pop_size)
  output$pB[g]<-length(which(c(parents$B1, parents$B2)=="B"))/(2*pop_size)
  output$pC[g]<-length(which(c(parents$C1, parents$C2)=="C"))/(2*pop_size)
  
  output$FisA[g]<-((2*output$pA[g]*(1-output$pA[g])) - (length(which(Ageno_tab[,g]=="aA"))/pop_size))/(2*output$pA[g]*(1-output$pA[g]))
  output$FisB[g]<-((2*output$pB[g]*(1-output$pB[g])) - (length(which(Bgeno_tab[,g]=="bB"))/pop_size))/(2*output$pB[g]*(1-output$pB[g]))
  output$FisC[g]<-((2*output$pC[g]*(1-output$pC[g])) - (length(which(Cgeno_tab[,g]=="cC"))/pop_size))/(2*output$pC[g]*(1-output$pC[g]))
  
  days<-c(min(parents$FLday):(max(parents$FLday)+(duration-1)))
  
  flowers<-as.data.frame(matrix(nrow=pop_size, ncol=length(days)))
  #flowers matrix is as wide as there are days in the flowering season, so an individual's flowers on each day are recorded
  names(flowers)<-paste("d", days, sep="")
  
  for (i in 1:pop_size){
    flowers[i, which(days==parents$FLday[i]):(which(days==parents$FLday[i])+duration-1)]<-fl_seq
  } #applies flowering schedule from parameters script to this range of days for an individual
  
  flowers[is.na(flowers)]<-0 #removes NA days
  
  if (LA_inbreeding_flowers == "yes") {
    for (i in 1:nrow(flowers)){
      flowers[i,]<-flowers[i,]*(fl_ID_weight)*(1-parents$ID_coef[i])
    }
  }
  #lowers fitness of individuals if they are inbred by lowering their flower output
  
  if (phenology== "yes") {
    #Set maternal matrix as the proportion of flowers made by a single individual relative to those made by the total population (both over entire season)
    mat<-as.matrix(flowers/sum(flowers))
    #dimensions are pop_size x length(days)
    #mat flower adjusted for total flowers made by the whole population over the entire flowering season
    
    #Set paternal matrix as the proportion of flowers made by a single individual relative to the entire population
    pat<-matrix(nrow=pop_size, ncol=length(days))
    
    #pat flower are adjusted for other flowers open on each day of an individual's flowering peiod
    for (t in 1:ncol(pat))	{	
      pat[,t]<-flowers[,t]/sum(flowers[,t])
    }
    #t is a day in the flowering season, because the columns in flower are days
    
    #Fix pat matrix to deal with days where NO flowers are produced in the population (this removes NaN)
    if (any(colSums(flowers)==0)){
      pat[,which(colSums(flowers)==0)]<-0
    }
    
    #Place moms as rows and dads as columns
    mat_opp<-mat %*% t(pat)
    
  } else if (phenology=="no") {
    mat_opp<-matrix(nrow=pop_size, ncol=pop_size)
    mat_opp[,]<-1/(pop_size^2)
  } #if phenology is not affecting mating probabilities, then everyone has equal opportunity so prob for individuals i and j is just 1/n
  
  if (compatibility == "yes") {
    S_allele_matching<-matrix(nrow=pop_size, ncol=pop_size)
    for (m in 1:pop_size) {
      for (p in 1:pop_size) {
        S_allele_matching[m, p]<- if (m == p) {0} else if (m != p) {1}
      } #assuming an infinite number of SI alleles, such that selfing cannot occur but individual i can mate with any other individual in the population (technically)
    }
  } else if (compatibility=="no") {
    S_allele_matching<-matrix(nrow=pop_size, ncol=pop_size)
    S_allele_matching[,]<-1
  } #if SI is not included, then everyone just gets a 1
  
  mat_opp_adj<-mat_opp * S_allele_matching
  
  #Re-adjust everything so that matrix sums to 1 (these need to be probability matrices)
  mat_opp_adj<-mat_opp_adj/sum(mat_opp_adj)

  distance1<-matrix(nrow=pop_size, ncol=pop_size)
  distance2<-matrix(nrow=pop_size, ncol=pop_size)
  
  for (m in 1:pop_size) {
    for (p in 1:pop_size) {
      distance1[m,p]<-sqrt((parents$X_pos[m]-parents$X_pos[p])^2 + (parents$Y_pos[m] - parents$Y_pos[p])^2)	
    } #distance is just the Euclidean distance between two points
  }
  distance1<-distance1*m_per_unit
  distance2[,]<-0 #when not on a torus, distance is not considered in the second direction
  
  #6a1) Run through pollen dispersal function
  if (space_pollen=="yes"){ #does pollen have limited dispersal?
    distance_for_mating1<-distance1
    distance_for_mating2<-distance2
  } else if (space_pollen=="no") {
    distance_for_mating1<-matrix(nrow=pop_size, ncol=pop_size)
    distance_for_mating2<-matrix(nrow=pop_size, ncol=pop_size)
    distance_for_mating1[,]<-1
    distance_for_mating2[,]<-0
  }
  
  pollen_dispersal1<-(((0.375)^2)/(2*pi))*exp(-0.375*distance_for_mating1) #if space_pollen is no, everyone gets the same value
    pollen_dispersal2<-0 #everyone gets the same value here as well
  
  mat_opp_adj2<-mat_opp_adj * (pollen_dispersal1 + pollen_dispersal2) #because pollen can be dispersed in either direction 1 OR direction 2, the probabilities are additive
  mat_opp_adj2<-mat_opp_adj2/sum(mat_opp_adj2) #readjusts the matrix to a probability
  
  parents$zFLday<-(parents$FLday-mean(parents$FLday))/sd(parents$FLday) #deviation from mean flowering day
  parents$fitness<-NA
  for (i in 1:pop_size){
    genotype<-paste(parents$A1[i], parents$A2[i], sep="") 
    parents$fitness[i]<-if (genotype=="aa") {
      1 + selection_aa*parents$zFLday[i] } else if (genotype=="AA"){1 + selection_AA*parents$zFLday[i]} else {1 + selection_Aa*parents$zFLday[i]}
  } #sets fitness as values around 1 based on A locus and in what direction of the mean flowering day the individual sits (will be relative because z is calulated relative to population)
  
  parents$fitness[which(parents$fitness<0)]<-0
  
  if (LA_inbreeding_fecundity == "yes") {
    parents$fitness<-parents$fitness * (1-parents$IB_coef)
  }
  
  offspring<-as.data.frame(matrix(nrow=pop_size, ncol=ncol(parents)))
  names(offspring)<-names(parents)
  
  offspring$position<-parents$position
  offspring$X_pos<-parents$X_pos
  offspring$Y_pos<-parents$Y_pos
  
  for (o in 1:pop_size) {
    
    #Calculate dispersal of seeds to location of offspring[o] as a weight on relative fitness and distance
    if (space_seed == "yes"){
      seed_dispersal1<-(exp(-(seed_lambda*distance1[o,])))
      seed_dispersal<-seed_dispersal1
    } else if (space_seed == "no") {
      seed_dispersal1<-rep(1, times=pop_size)
      seed_dispersal<-seed_dispersal1
    } 
    seed_dispersal<-seed_dispersal1
    
    if ((space_seed == "yes") & (torus == "yes")) { 
      seed_dispersal2<-(exp(-(seed_lambda*distance2[o,])))
      seed_dispersal<-seed_dispersal1+seed_dispersal2
    } else if ((space_seed == "no") & (torus == "yes")) {
      seed_dispersal2<-rep(1, times=pop_size)
      seed_dispersal<-seed_dispersal1+seed_dispersal2
    }
    
    #Choose a mom based on 1) distance to seed location and 2) relative amount of pollen received
    dist.prob = (parents$fitness*seed_dispersal)/sum(parents$fitness*seed_dispersal)
    pollen.prob = rowSums(mat_opp_adj2)/sum(rowSums(mat_opp_adj2))
    joint.prob = (dist.prob*pollen.prob)/sum(dist.prob*pollen.prob)
    #Choose a mom based on relative fitness and chance of getting a seed to the specified location
    ##NOTE- Edge plants at an inherent disadvantage if a torus is not used.
    mom<-sample(c(1:pop_size), size=1, prob=joint.prob)
    
    #If there are no sires for chosen mom, choose another mom
    while (rowSums(mat_opp_adj2)[mom]==0) {
      mom<-sample(c(1:pop_size), size=1, prob=(parents$fitness)/sum(parents$fitness))
    }
    
    #Choose sire for selected mom
    pat_prob<-mat_opp_adj2[mom,]/rowSums(mat_opp_adj2)[mom]
    dad<-sample(x=c(1:pop_size), size=1, prob=pat_prob)
    
    #9b5) Record mom, dad, and inbreeding coefficient
    offspring$Mother[o]<-mom
    offspring$Father[o]<-dad
    offspring$IB_coef[o]<-coancestry_p[mom,dad]
    
    #If early acting inbreeding-depression is in place, check to see if IB_coef is below threshold
    #If it's too high, then discard the seed and choose again
    if (EA_inbreeding == "yes") {
      while (offspring$IB_coef[o] > EA_id_threshold) {
        #Choose a new mom
        mom<-sample(c(1:pop_size), size=1, prob=(parents$fitness*seed_dispersal)/sum(parents$fitness*seed_dispersal))
        #If there are no sires for chosen mom, choose another.
        while (rowSums(mat_opp_adj2)[mom]==0) {
          mom<-sample(c(1:pop_size), size=1, prob=(parents$fitness)/sum(parents$fitness))
        }
      } #close while statement regarding threshold
    } #close if over whether early-acting ID in is place
    
    #A-locus inhertiance (one from mom, one from dad)
    offspring$A1[o]<-as.character(sample(x=parents[mom, 5:6], size=1))
    offspring$A2[o]<-as.character(sample(x=parents[dad, 5:6], size=1))
    
    #B-locus inhertiance (one from mom, one from dad)
    offspring$B1[o]<-as.character(sample(x=parents[mom, 7:8], size=1))
    offspring$B2[o]<-as.character(sample(x=parents[dad, 7:8], size=1))
    
    #C-locus inhertiance (one from mom, one from dad)
    offspring$C1[o]<-as.character(sample(x=parents[mom, 9:10], size=1))
    offspring$C2[o]<-as.character(sample(x=parents[dad, 9:10], size=1))
    
    #Mendellian inheritance of flowering time loci 
    for (n in 1:(FT_loci)) {
      offspring[o,11 + 2*(n-1)]<-sample(x=parents[mom,(11+2*(n-1)):(12 + 2*(n-1))], size=1)
      offspring[o,12 + 2*(n-1)]<-sample(x=parents[dad,(11+2*(n-1)):(12 + 2*(n-1))], size=1)
    } #close loop over the flowering time loci
    
    #Mutation at flowering time loci
    if (mut == "yes") {
      for (n in 1:(FT_loci)) {
        mutation<-sample(1:(1/mut_rate), size=1) #assuming mutation rate is for flowering time loci
        if (mutation == 1){
          effect<-sample((-1:1),size=1)
          parent_mut<-sample(1:2,size=1)
          if (parent_mut == 1) {
            offspring[o,11+2*(n-1)]<-offspring[o,11+2*(n-1)]+effect
          } else {
            offspring[o,12+2*(n-1)]<-offspring[o,12+2*(n-1)]+effect
          }
        }
      }
    }
    
    #Mendellian inheritance of remaining neutral loci
    for (j in 1:(neut_loci)) {
      offspring[o,11+(2*FT_loci)+(2*(j-1))]<-as.character(sample(x=parents[mom,(11+(2*FT_loci)+(2*(j-1))):(12+(2*FT_loci)+(2*(j-1)))], size=1))
      offspring[o,12+(2*FT_loci)+(2*(j-1))]<-as.character(sample(x=parents[dad,(11+(2*FT_loci)+(2*(j-1))):(12+(2*FT_loci)+(2*(j-1)))], size=1))
    } #close loop over remaining neutral loci
    
    #Mutation at neutral loci
    if (mut == "yes") {
      for (n in 1:(neut_loci)) {
        mutation<-sample(1:(1/mut_rate), size=1) #assuming mutation rate is for flowering time loci
        if (mutation == 1){
          effect<-sample(c("D","d"),size=1)
          parent_mut<-sample(1:2,size=1)
          if (parent_mut == 1) {
            offspring[o,11+(2*FT_loci)+(2*(n-1))]<-effect
          } else {
            offspring[o,12+(2*FT_loci)+(2*(n-1))]<-effect
          }
        }
      }
    }
  } #close loop over offspring
  
  offspring$FTnoise<-rnorm(n=pop_size, mean=0, sd=fl_noise)
  offspring$FLday<-floor(fl_mean + rowSums(offspring[,11:(10+2*FT_loci)]) + offspring$FTnoise) 
  
  coancestry_o<-matrix(nrow=pop_size, ncol=pop_size)
  
  for (i in 2:pop_size){
    for (j in 1:(i-1)){	
      Mi<-offspring$Mother[i]
      Mj<-offspring$Mother[j]
      Fi<-offspring$Father[i]
      Fj<-offspring$Father[j]
      
      pathA<-(1/4)*coancestry_p[Mi, Mj]
      pathB<-(1/4)*coancestry_p[Fi, Mj]
      pathC<-(1/4)*coancestry_p[Mi, Fj]
      pathD<-(1/4)*coancestry_p[Fi, Fj]
      
      coancestry_o[i, j]<-(pathA + pathB + pathC + pathD)
      coancestry_o[j, i]<-(pathA + pathB + pathC + pathD)
    }
  }
  
  for (i in 1:pop_size){
    coancestry_o[i,i]<-(1/2)*(1 + offspring$IB_coef[i])
  }
  
  if ((g == 1)||(g == 10)||(g == 20)||(g == 30)||(g == 40)||(g == 50)||(g == 60)||(g == 70)||(g == 80)||(g == 90)||(g == 100)||(g == 150)||(g == 200)||(g == 250)||(g == 300)||(g == 350)||(g == 400)||(g == 450)||(g == 500)||(g == 550)||(g == 600)||(g == 650)||(g == 700)||(g == 750)||(g == 800)){
    
    offspring_map<-as.data.frame(matrix(nrow=pop_size, ncol=(11+(2*FT_loci)+(2*neut_loci))))
    names(offspring_map)[1:11]<-c("FLday","X_pos","Y_pos","genotypeA","mapA","genotypeB","mapB","genotypeC","mapC","Mother","Father")
    
    offspring_map$FLday<-offspring$FLday
    offspring_map$X_pos<-offspring$X_pos
    offspring_map$Y_pos<-offspring$Y_pos
    offspring_map$Mother<-offspring$Mother
    offspring_map$Father<-offspring$Father
    
    for (n in 1:FT_loci)
    {
      names(offspring_map)[12+2*(n-1)]<-paste("loc",n,"a",sep="")
      names(offspring_map)[13+2*(n-1)]<-paste("loc",n,"b",sep="")
    }
    
    for (i in 1:pop_size){
      for (n in 1:FT_loci)
      {
        offspring_map[i,(12+2*(n-1))]<-offspring[i,(11+2*(n-1))]
        offspring_map[i,(13+2*(n-1))]<-offspring[i,(12+2*(n-1))]
      }
    }
    
    for (j in 1:neut_loci)
    {
      names(offspring_map)[12+(2*FT_loci)+(2*(j-1))]<-paste("neut",j,"a",sep="")
      names(offspring_map)[13+(2*FT_loci)+(2*(j-1))]<-paste("neut",j,"b",sep="")
    }
    
    for (i in 1:pop_size){
      for (j in 1:neut_loci)
      {
        offspring_map[i,12+(2*FT_loci)+(2*(j-1))]<-offspring[i,12+(2*FT_loci)+(2*(j-1))]
        offspring_map[i,13+(2*FT_loci)+(2*(j-1))]<-offspring[i,13+(2*FT_loci)+(2*(j-1))]
      }
    }
    
    for (i in 1:pop_size){
      offspring_map$genotypeA[i]<-paste(offspring$A1[i], offspring$A2[i], sep="")	
      offspring_map$mapA[i]<-if (offspring_map$genotypeA[i]=="aa") {
        0} else if (offspring_map$genotypeA[i]=="AA"){2} else {1}
    }
    
    for (i in 1:pop_size){
      offspring_map$genotypeB[i]<-paste(offspring$B1[i], offspring$B2[i], sep="")	
      offspring_map$mapB[i]<-if (offspring_map$genotypeB[i]=="bb") {
        0} else if (offspring_map$genotypeB[i]=="BB"){2} else {1}
    }
    
    for (i in 1:pop_size){
      offspring_map$genotypeC[i]<-paste(offspring$C1[i], offspring$C2[i], sep="")	
      offspring_map$mapC[i]<-if (offspring_map$genotypeC[i]=="cc") {
        0} else if (offspring_map$genotypeC[i]=="CC"){2} else {1}
    }
    
    write.csv(offspring_map, paste("paraset_",x,"_offspring_map_", g, ".csv",sep="")
              , row.names=F)
    
  }
  #end loop over g if
  
  comment.out = function(){
    if ((g == 1)||(g == 50)||(g == 100)||(g == 150)||(g == 200)||(g == 250)||(g == 300)||(g == 350)||(g == 400)||(g == 450)||(g == 500)){
      
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
    } #end loop over g if
  }
  
  write.csv(output, "data_by_generation.csv")
  
  write.csv(onset_tab,
            "Flowering_time_onsets_over_generations.csv", row.names=F)
  
  write.csv(Ageno_tab,
            "A_locus_genotypes_over_generations.csv", row.names=F)
  
  write.csv(Bgeno_tab,
            "B_locus_genotypes_over_generations.csv", row.names=F)
  
  write.csv(Cgeno_tab,
            "C_locus_genotypes_over_generations.csv", row.names=F)
  
  if (g<no_gen){
    offspring->parents
    coancestry_o->coancestry_p}
}
