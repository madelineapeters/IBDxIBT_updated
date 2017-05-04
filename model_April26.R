############## MODEL #################
##Cleaned 26/04/2017
##Includes torus spacing
##Includes mutation
##Includes customisable number of neutral loci

###Create and set directory for saved files###
dir.create(path=paste(getwd(), paste("para_set", para_set, sep="_"), paste("model_run", no_runs, sep="_"), sep="/"), recursive=T)
setwd(paste(getwd(), paste("para_set", para_set, sep="_"), paste("model_run", no_runs, sep="_"), sep="/"))


###Create table holding data with one value per generation###
output<-as.data.frame(matrix(nrow=no_gen, ncol=11))
names(output)<-c("gen","FTmean","FTvar","FTh2","IBmean","pA","pB","pC","FisA","FisB","FisC")


###Create tables holding data with several values per generation###

#Create table for flowering time by generation
onset_tab<-(matrix(nrow=pop_size, ncol=no_gen))

#Create table for A-locus genotypes by generation
Ageno_tab<-matrix(nrow=pop_size, ncol=no_gen)

#Create table for B-locus genotypes by generation
Bgeno_tab<-matrix(nrow=pop_size, ncol=no_gen)

#Create table for C-locus genotypes by generation
Cgeno_tab<-matrix(nrow=pop_size, ncol=no_gen)


###Create parental generation (generation 1)###

#Step 1) Create table to hold parental information
	parents<-as.data.frame(matrix(nrow=pop_size, ncol=(10+(2*FT_loci)+4+(2*neut_loci))))

	names(parents)[1:10]<-c("Mother","Father","IB_coef","FLday","A1","A2","B1","B2","C1","C2")
	for (n in 1:FT_loci)
	{
	  names(parents)[11+2*(n-1)]<-paste("loc",n,"a",sep="")
	  names(parents)[12+2*(n-1)]<-paste("loc",n,"b",sep="")
	}
	names(parents)[(ncol(parents)-3): ncol(parents)]<-c("FTnoise","position","X_pos","Y_pos")
	for (j in 1:neut_loci)
	{
	  names(parents)[11+(2*FT_loci)+(2*(j-1))]<-paste("neut",j,"a",sep="")
	  names(parents)[12+(2*FT_loci)+(2*(j-1))]<-paste("neut",j,"b",sep="")
	}
#Step 2) Generate flowering time onsets
	#2a Give each individual its flowering time genotype (2 haploid complements)
	for (n in 1:FT_loci)
	{
	  parents[,11+2*(n-1)]<-sample(FT_alleles, size=pop_size, replace=T)
	  parents[,12+2*(n-1)]<-sample(FT_alleles, size=pop_size, replace=T)
	}

	#2b) Add environmental variance to flowering time
	parents$FTnoise<-rnorm(n=pop_size, mean=0, sd=fl_noise)

	#2c) Add genertic and environmetal to some positive baseline to get each individual's flowering time
	parents$FLday<-floor(fl_mean + rowSums(parents[,11:(10+2*FT_loci)])+parents$FTnoise)

#Step 3) Assign genotypes at the neutral loci
	#3a) Give each individual a genotype at the A-locus
	parents$A1<-sample(c("A","a"), size=pop_size, prob=c(freqA, 1-freqA), replace=T)
	parents$A2<-sample(c("A","a"), size=pop_size, prob=c(freqA, 1-freqA), replace=T)

	#3b) Give each individual a genotype at the B-locus
	parents$B1<-sample(c("B","b"), size=pop_size, prob=c(freqB, 1-freqB), replace=T)
	parents$B2<-sample(c("B","b"), size=pop_size, prob=c(freqB, 1-freqB), replace=T)

	#3c) Give each individual a genotype at the C-locus
	parents$C1<-sample(c("C","c"), size=pop_size, prob=c(freqC, 1-freqC), replace=T)
	parents$C2<-sample(c("C","c"), size=pop_size, prob=c(freqC, 1-freqC), replace=T)
	
	#3d) Give each individual a genotype at its remaining neutral loci
	for (j in 1:neut_loci)
	{
	  parents[,11+(2*FT_loci)+(2*(j-1))]<-sample(c("D","d"), size=pop_size, prob=c(freqD, 1-freqD), replace=T)
	  parents[,12+(2*FT_loci)+(2*(j-1))]<-sample(c("D","d"), size=pop_size, prob=c(freqD, 1-freqD), replace=T)
	}

#Step 4) Give each individual a position in space
	#4a) From the matrix of available positions in space, develop a list of positions from which to sample
	for (y in 1:ncol(space_50by50_sq_grid)) {
	  if (y ==1 ) {
		pos_vector<-c(as.character(space_50by50_sq_grid[,y])) }
	  else {
		pos_vector<-c(pos_vector, as.character(space_50by50_sq_grid[,y])) 
	  }
	}
	
	#4b) Give everyone a position from the set of possible positions, and then extract their x and y coordinates
	parents$position<-(pos_vector)	
	parents$X_pos<-sapply(X=parents$position, FUN=function(X) {
	  as.numeric(strsplit(X, split=",")[[1]][1])
	})
	parents$Y_pos<-sapply(X=parents$position, FUN=function(X) {
	  as.numeric(strsplit(X, split=",")[[1]][2])
	})	


#Step 5) Calculate coancestry and inbreeding coefficients of parents in generation 1
	#5a) Create coancestry matrix
	coancestry_p<-matrix(nrow=pop_size, ncol=pop_size)
	
	#5b) Set coancestry values (all equal 0, except along diagonal, at which it equals 0.5)
	coancestry_p[,]<-0
	for (i in 1:pop_size) {
	  coancestry_p[i,i]<-0.5
	}

	#5c) Set inbreeding coefficient for all pairings not on the diagonal
	parents$IB_coef<-0

###Run offspring generations(generation 2 through no_gen)###
for (g in 1:no_gen){
  
#Step 1) Store generation g conditions into the various tables

	#1a) Fill table holding data with one value per generation
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
  
  	#1b) Fill table holding data for flowering time onset over generations
	onset_tab[,g]<-parents$FLday
  
   	#1c) Fill table holding data for neutral loci over generations
	Ageno_tab[,g]<-apply(X=parents[,5:6], MARGIN=1, FUN=function(X){
	paste(sort(X)[1], sort(X)[2], sep="")
	})

	Bgeno_tab[,g]<-apply(X=parents[,7:8], MARGIN=1, FUN=function(X){
	paste(sort(X)[1], sort(X)[2], sep="")
	})

	Cgeno_tab[,g]<-apply(X=parents[,9:10], MARGIN=1, FUN=function(X){
	paste(sort(X)[1], sort(X)[2], sep="")
	})
  
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
  
#Step 4) Develop a compatibility matrix based on self-incompatibility alleles
  if (compatibility == "yes") {
    S_allele_matching<-matrix(nrow=pop_size, ncol=pop_size)
    for (m in 1:pop_size) {
		for (p in 1:pop_size) {
			mat_expr1<-parents$S1[m]
			mat_expr2<-parents$S2[m]
			pat_expr1<-parents$S1[p]
			pat_expr2<-parents$S2[p]
			
			S_allele_matching[m, p]<- if (identical(mat_expr1, pat_expr1)==T) {0} else if (identical(mat_expr1, pat_expr2)==T) {0} else if (identical(mat_expr2, pat_expr1)==T) {0} else if (identical(mat_expr2, pat_expr2)==T) {0} else {1} 
		}
	}
  } else if (compatibility=="no") {
    S_allele_matching<-matrix(nrow=pop_size, ncol=pop_size)
    S_allele_matching[,]<-1
  }
  
#Step 5) Recalibrate the mating opportunity matrix to account for self-incompatibility
  mat_opp_adj<-mat_opp * S_allele_matching
  
  #5a) Re-adjust everything so that matrix sums to 1
  mat_opp_adj<-mat_opp_adj/sum(mat_opp_adj)

#Step 6) Develop pairwise space (Euclidean distance) matrix and run through the pollen dispersal function
	#6a) Calculate probability based on distance
	if (torus == "no") {
	  distance1<-matrix(nrow=pop_size, ncol=pop_size)
	  distance2<-matrix(nrow=pop_size, ncol=pop_size)
	  for (m in 1:pop_size) {
	    for (p in 1:pop_size) {
	      distance[m,p]<-sqrt((parents$X_pos[m]-parents$X_pos[p])^2 + (parents$Y_pos[m] - parents$Y_pos[p])^2)	
	    }
	  }
	  distance1<-distance*m_per_unit
	  distance2<-0
	  
	  #6a1) Run through pollen dispersal function
	  if (space_pollen=="yes"){
	    distance_for_mating<-distance
	  } else if (space_pollen=="no") {
	    distance_for_mating<-matrix(nrow=pop_size, ncol=pop_size)
	    distance_for_mating[,]<-1
	  }
	  pollen_dispersal1<-(((0.375)^2)/(2*pi))*exp(-0.375*distance_for_mating)	#function taken from Lavigne et al. 1998, technically only appropriate for distances >3 (45% of pollen), but we don't have anything else to describe what happens within 3m...
	  pollen_dispersal2<-0
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
		
		#6a1) Run through pollen dispersal function
		if (space_pollen=="yes"){
			distance_for_mating1<-distance1
			distance_for_mating2<-distance2
		} else if (space_pollen=="no") {
			distance_for_mating1<-matrix(nrow=pop_size, ncol=pop_size)
			distance_for_mating1[,]<-1
			distance_for_mating2<-matrix(nrow=pop_size, ncol=pop_size)
			distance_for_mating2[,]<-1
		}
		pollen_dispersal1<-(((0.375)^2)/(2*pi))*exp(-0.375*distance_for_mating1)
		pollen_dispersal2<-(((0.375)^2)/(2*pi))*exp(-0.375*distance_for_mating2)
	} #end of if torus "yes"
 
#Step 7) Recalibrate mating opportunity matrix to account for space
	mat_opp_adj2<-mat_opp_adj * (pollen_dispersal1 + pollen_dispersal2)
	mat_opp_adj2<-mat_opp_adj2/sum(mat_opp_adj2)
 
#Step 8) Calculate relative fitness based on A locus and flowering time		
	parents$zFLday<-(parents$FLday-mean(parents$FLday))/sd(parents$FLday)
	parents$fitness<-NA
	for (i in 1:pop_size){
		genotype<-paste(parents$A1[i], parents$A2[i], sep="")	
		parents$fitness[i]<-if (genotype=="aa") {
			1 + selection_aa*parents$zFLday[i] } else if (genotype=="AA"){1 + selection_AA*parents$zFLday[i]} else {1 + selection_Aa*parents$zFLday[i]}
	}			
  
  #8a) Correct for fitness values
  parents$fitness[which(parents$fitness<0)]<-0	
  
  #8b) Adjust fitness according to (1-IB_coef)	
	if (LA_inbreeding_fecundity == "yes") {
		parents$fitness<-parents$fitness * (1-parents$IB_coef)
	} 	
  
  
#Step 9) Generate offspring - the probability of choosing a given mom depends on her fitness as defined by the A-locus and flowering time.
	#Dad's identity depends on adjusted mating matrix.
	#If a chosen mom has no compatible crosses available (i.e., if rowSums(mat_opp_adj2)[m]==0, then we have to choose another mom.
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
			seed_dispersal<-(exp(-(seed_lambda*(distance1[o,]+distance2[o,]))))} else if (space_seed == "no") {
			seed_dispersal<-rep(1, times=pop_size)}
		
		#9b2) Choose a mom based on relative fitness and chance of getting a seed to the specified location
			##NOTE- Edge plants at an inherent disadvantage if a torus is not used.
		mom<-sample(c(1:pop_size), size=1, prob=(parents$fitness*seed_dispersal)/sum(parents$fitness*seed_dispersal))
    
		#9b3) If there are no sires for chosen mom, choose another mom
		while (rowSums(mat_opp_adj2)[mom]==0)	{mom<-sample(c(1:pop_size), size=1, prob=(parents$fitness)/sum(parents$fitness))}
    
    
		#9b4) Choose sire for selected mom
		pat_prob<-mat_opp_adj2[mom,]/rowSums(mat_opp_adj2)[mom]
		dad<-sample(x=c(1:pop_size), size=1, prob=pat_prob)
    
		#9b5) Record mom, dad, and inbreeding coefficient
		offspring$Mother[o]<-mom
		offspring$Father[o]<-dad		
		offspring$IB_coef[o]<-coancestry_p[mom,dad]
    
		#9b6) #If early acting inbreeding-depression is in place, check to see if IB_coef is below threshold
		#If it's too high, then discard the seed and choose again
		if (EA_inbreeding == "yes") {
			while (offspring$IB_coef[o] > EA_id_threshold) {
				#12a: Choose a mom
				mom<-sample(c(1:pop_size), size=1, prob=(parents$fitness*seed_dispersal)/sum(parents$fitness*seed_dispersal))
				#If there are no sires for chosen mom, choose another.
				while (rowSums(mat_opp_adj2)[mom]==0) {
					mom<-sample(c(1:pop_size), size=1, prob=(parents$fitness)/sum(parents$fitness))
				}
			} #close while statement regarding threshold	
		} #close if over whether early-acting ID in is place.
    
    
		#9b7) A-locus inhertiance (one from mom, one from dad)
		offspring$A1[o]<-as.character(sample(x=parents[mom, 5:6], size=1))
		offspring$A2[o]<-as.character(sample(x=parents[dad, 5:6], size=1))
		
		#9b8) B-locus inhertiance (one from mom, one from dad)
		offspring$B1[o]<-as.character(sample(x=parents[mom, 7:8], size=1))
		offspring$B2[o]<-as.character(sample(x=parents[dad, 7:8], size=1))
		
		#9b9) C-locus inhertiance (one from mom, one from dad)
		offspring$C1[o]<-as.character(sample(x=parents[mom, 9:10], size=1))
		offspring$C2[o]<-as.character(sample(x=parents[dad, 9:10], size=1))
		
		#9b10) Mendellian inheritance of flowering time loci 
		for (n in 1:(FT_loci)) {
		  offspring[o,11 + 2*(n-1)]<-sample(x=parents[mom,(11+2*(n-1)):(12 + 2*(n-1))], size=1)
		  offspring[o,12 + 2*(n-1)]<-sample(x=parents[dad,(11+2*(n-1)):(12 + 2*(n-1))], size=1)
		}	#close loop over the flowering time loci
		
		#9b11) Mutation at flowering time loci (locus wide mutation rate u = 0.001)
		  if (mut == "yes") {
  		  for (n in 1:(FT_loci)) {
  			mutation<-sample(1:((1/mut_rate)*n), size=1)
  			if (mutation == 1){
  			  effect<-sample((-2:2),size=1)
  			  parent_mut<-sample(1:2,size=1)
  			  if (parent_mut == 1) {
  				  offspring[o,11+2*(n-1)]<-offspring[o,11+2*(n-1)]+effect*0.5
  			  } else {
  				  offspring[o,12+2*(n-1)]<-offspring[o,12+2*(n-1)]+effect*0.5
  			  }
  			}
  		  }
		  }
		#9b12) Mendellian inheritance of remaining neutral loci
		for (j in 1:(neut_loci)) {
		  offspring[o,11+(2*FT_loci)+(2*(j-1))]<-as.character(sample(x=parents[mom,(11+(2*FT_loci)+(2*(j-1))):(12+(2*FT_loci)+(2*(j-1)))], size=1))
		  offspring[o,12+(2*FT_loci)+(2*(j-1))]<-as.character(sample(x=parents[dad,(11+(2*FT_loci)+(2*(j-1))):(12+(2*FT_loci)+(2*(j-1)))], size=1))
		}	#close loop over remaining neutral loci
  }	#close loop over offspring
  
  #9c) Add environmental variance to flowering time, and sum flwg time columns to calculate flowering time	
  offspring$FTnoise<-rnorm(n=pop_size, mean=0, sd=fl_noise)
  offspring$FLday<-floor(fl_mean + rowSums(offspring[,11:(10+2*FT_loci)]) + offspring$FTnoise) 
	
#Step 10) Calculate coancestry in the new generation
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

  
#Step 11) Repeat model using offspring as parents
  #place within 'if' statement so that both the parent and offspring files of final generation are preserved for debugging
  #11a) Save information about loci and space every 50 generations
  if ((g == 1)||(g == 50)||(g == 100)||(g == 150)||(g == 200)||(g == 250)||(g == 300)||(g == 350)||(g == 400)||(g == 450)||(g == 500)){
    
    offspring_map<-as.data.frame(matrix(nrow=pop_size, ncol=(9+(2*FT_loci)+(2*neut_loci))))
    names(offspring_map)[1:9]<-c("FLday","X_pos","Y_pos","genotypeA","mapA","genotypeB","mapB","genotypeC","mapC")
    
    offspring_map$FLday<-offspring$FLday
    offspring_map$X_pos<-offspring$X_pos
    offspring_map$Y_pos<-offspring$Y_pos
    
    for (n in 1:FT_loci)
    {
      names(offspring_map)[10+2*(n-1)]<-paste("loc",n,"a",sep="")
      names(offspring_map)[11+2*(n-1)]<-paste("loc",n,"b",sep="")
    }
    
    for (i in 1:pop_size){
      for (n in 1:FT_loci)
      {
        offspring_map[i,(10+2*(n-1))]<-offspring[i,(11+2*(n-1))]
        offspring_map[i,(11+2*(n-1))]<-offspring[i,(12+2*(n-1))]
      }
    }
	
	for (j in 1:neut_loci)
    {
      names(offspring_map)[10+(2*FT_loci)+(2*(j-1))]<-paste("neut",n,"a",sep="")
      names(offspring_map)[11+(2*FT_loci)+(2*(j-1))]<-paste("neut",n,"b",sep="")
    }
	
    for (i in 1:pop_size){
      for (j in 1:neut_loci)
      {
        offspring_map[i,10+(2*FT_loci)+(2*(j-1))]<-offspring[i,11+(2*FT_loci)+(2*(j-1))]
        offspring_map[i,11+(2*FT_loci)+(2*(j-1))]<-offspring[i,12+(2*FT_loci)+(2*(j-1))]
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
    
    write.csv(offspring_map, paste("offspring_map", g, sep="_")
              , row.names=F)
    
  }
  
  if (g<no_gen){
    offspring->parents
    coancestry_o->coancestry_p}
  
}