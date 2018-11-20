library(dplyr)

##########Set-up directories##########
dir.create(path=paste(getwd(), paste("para_set", x, sep="_"), sep="/"), recursive=T)

setwd(paste(getwd(), paste("para_set", para_set, sep="_"), sep="/"))

r=1
dir.create(path=paste(getwd(), paste("model_run", r, sep="_"), sep="/"), recursive=T)

#dir.create(path=paste(getwd(), paste("para_set", para_set, sep="_"), paste("model_run", r, sep="_"), sep="/"), recursive=T)
#directory is created using the para_set number read in from the parameter script and parameters csv and the model_run number from the commander script

setwd(paste(getwd(), paste("model_run", r, sep="_"), sep="/"))
#setwd(paste(getwd(), paste("para_set", para_set, sep="_"), paste("model_run", r, sep="_"), sep="/"))
#directory is set relative to the current working directory, which is set in the commander script

output = as.data.frame(matrix(nrow=no_gen, ncol=11))
#no_gen is taken from the parameters script
names(output) = c("gen","FTmean","FTvar","FTh2","IBmean","pA","pB","pC","FisA","FisB","FisC")


##########Tables holding several values per generation##########

#TABLE TO HOLD FLOWERING TIMES BY GENERATION
onset_tab = (matrix(nrow=pop_size, ncol=no_gen))

#TABLE FOR A-LOCUS GENOTYPES BY GENERATION
Ageno_tab = matrix(nrow=pop_size, ncol=no_gen)

#TABLE FOR B-LOCUS GENOTYPES BY GENERATION
Bgeno_tab = matrix(nrow=pop_size, ncol=no_gen)

#TABLE FOR C-LOCUS GENOTYPES BY GENERATION
Cgeno_tab = matrix(nrow=pop_size, ncol=no_gen)

sum.N = as.data.frame(matrix(nrow=no_gen, ncol=3))
self.df = as.data.frame(matrix(nrow=no_gen, ncol=1))


##########Parental data-frame##########
parents = as.data.frame(matrix(nrow=pop_size, ncol=(10+2*FT_loci)+4+2*neut_loci))
#2 columns are needed for each locus so that each allele per locus can be stored

names(parents)[1:10] = c("Mother","Father","IB_coef","FLday","A1","A2","B1","B2","C1","C2")
#A, B and C are three additional neutral loci used primarily for heatmap purposes
#A can also be set to be under selection using the parameter .csv

for (n in 1:FT_loci)
{
  names(parents)[11+2*(n-1)] = paste("loc",n,"a",sep="")
  names(parents)[12+2*(n-1)] = paste("loc",n,"b",sep="")
}

names(parents)[(ncol(parents)-3): ncol(parents)] = c("FTnoise","position","X_pos","Y_pos")

for (j in 1:neut_loci)
{
  names(parents)[11+(2*FT_loci)+(2*(j-1))] = paste("neut",j,"a",sep="")
  names(parents)[12+(2*FT_loci)+(2*(j-1))] = paste("neut",j,"b",sep="")
}

#name assignments are set up so column positions are relative to the number of neutral loci specified

for (n in 1:FT_loci) {
  parents[,11+2*(n-1)] = sample(FT_alleles, size=pop_size, replace=T)
  parents[,12+2*(n-1)] = sample(FT_alleles, size=pop_size, replace=T)
}

#potential alleles are drawn from FT_alleles, which are specified in the parameters script
#to compare to neutral loci, best to make each FT locus biallelic

parents$FTnoise = rnorm(n=pop_size, mean=0, sd=fl_noise)

parents$FLday = floor(fl_mean + rowSums(parents[,11:(10+2*FT_loci)])+parents$FTnoise)

parents$A1 = sample(c("A","a"), size=pop_size, prob=c(freqA, 1-freqA), replace=T)
parents$A2 = sample(c("A","a"), size=pop_size, prob=c(freqA, 1-freqA), replace=T)

parents$B1 = sample(c("B","b"), size=pop_size, prob=c(freqB, 1-freqB), replace=T)
parents$B2 = sample(c("B","b"), size=pop_size, prob=c(freqB, 1-freqB), replace=T)

parents$C1 = sample(c("C","c"), size=pop_size, prob=c(freqC, 1-freqC), replace=T)
parents$C2 = sample(c("C","c"), size=pop_size, prob=c(freqC, 1-freqC), replace=T)

#starts population at HWE
#allele frequencies are typically starting at 0.5 (can be adjusted to alter starting conditions)

freqD = 0.5

for (j in 1:neut_loci)
{
  parents[,11+(2*FT_loci)+(2*(j-1))] = sample(c("D","d"), size=pop_size, prob=c(freqD, 1-freqD), replace=T)
  parents[,12+(2*FT_loci)+(2*(j-1))] = sample(c("D","d"), size=pop_size, prob=c(freqD, 1-freqD), replace=T)
}

#assign parents a space
parents$X_pos = xy.coor[,1]
parents$Y_pos = xy.coor[,2]

#coancestry
coancestry_p = matrix(nrow=pop_size, ncol=pop_size)

coancestry_p[,] = 0
for (i in 1:pop_size) {
  coancestry_p[i,i] = 0.5
}

parents$IB_coef = 0


##########Offspring generations##########

for (g in 1:no_gen) {
  
  ##########Record summary data from previous generation##########
  output$gen[g] = g
  output$FTmean[g] = mean(parents$FLday)
  output$FTvar[g] = var(parents$FLday)
  output$FTh2[g] = var(rowSums(parents[,11:(10+2*FT_loci)]))/var(parents$FLday)
  output$IBmean[g] = mean(parents$IB_coef)
  output$pA[g] = length(which(c(parents$A1, parents$A2)=="A"))/(2*pop_size)
  output$pB[g] = length(which(c(parents$B1, parents$B2)=="B"))/(2*pop_size)
  output$pC[g] = length(which(c(parents$C1, parents$C2)=="C"))/(2*pop_size)
  
  output$FisA[g] = ((2*output$pA[g]*(1-output$pA[g])) - (length(which(Ageno_tab[,g]=="aA"))/pop_size))/(2*output$pA[g]*(1-output$pA[g]))
  output$FisB[g] = ((2*output$pB[g]*(1-output$pB[g])) - (length(which(Bgeno_tab[,g]=="bB"))/pop_size))/(2*output$pB[g]*(1-output$pB[g]))
  output$FisC[g] = ((2*output$pC[g]*(1-output$pC[g])) - (length(which(Cgeno_tab[,g]=="cC"))/pop_size))/(2*output$pC[g]*(1-output$pC[g]))
  
  
  ##########Generate mating probability matrix##########
  #Calculating mating probabilities based on flowering time
  days = c(min(parents$FLday):(max(parents$FLday)+(duration-1)))
  
  flowers = as.data.frame(matrix(nrow=pop_size, ncol=length(days))) #flowers matrix is as wide as there are days in the flowering season, so an individual's flowers on each day are recorded
  names(flowers) = paste("d", days, sep="")
  
  for (i in 1:pop_size){
    flowers[i, which(days==parents$FLday[i]):(which(days==parents$FLday[i])+duration-1)] = fl_seq
  } #applies flowering schedule from parameters script to this range of days for an individual
  
  flowers[is.na(flowers)] = 0 #removes NA days

  if (phenology== "yes") {
    #Set maternal matrix as the proportion of flowers made by a single individual relative to those made by the total population (both over entire season)
    mat = as.matrix(flowers/sum(flowers))
    #dimensions are pop_size x length(days)
    
    #pat flower are adjusted for other flowers open on each day of an individual's flowering period
    pat = matrix(nrow=pop_size, ncol=length(days))
    for (t in 1:ncol(pat))	{	
      pat[,t] = flowers[,t]/sum(flowers[,t]) #t is a day in the flowering season as the columns are days in flowering season
    }
    
    #Fix pat matrix to deal with days where NO flowers are produced in the population (this removes NaN)
    if (any(colSums(flowers)==0)){
      pat[,which(colSums(flowers)==0)] = 0
    }
    
    #Place moms as rows and dads as columns
    mat_opp = mat %*% t(pat)
    
  } else if (phenology=="no") {
    mat_opp = matrix(nrow=pop_size, ncol=pop_size)
    mat_opp[,] = 1/(pop_size^2)
  } #if phenology is not affecting mating probabilities, then everyone has equal opportunity so prob for individuals i and j is just 1/n^2
  
  if (compatibility == "yes") {
    S_allele_matching = matrix(nrow=pop_size, ncol=pop_size)
    for (m in 1:pop_size) {
      for (p in 1:pop_size) {
        S_allele_matching[m, p] =  if (m == p) {0} else if (m != p) {1}
      } #assuming an infinite number of SI alleles, such that selfing cannot occur but individual i can mate with any other individual in the population (technically)
    }
  } else if (compatibility=="no") {
    S_allele_matching = matrix(nrow=pop_size, ncol=pop_size)
    S_allele_matching[,] = 1
  } #if SI is not included, then everyone just gets a 1
  
  mat_opp_adj = mat_opp * S_allele_matching
  
  #Re-adjust everything so that matrix sums to 1 (these need to be probability matrices)
  mat_opp_adj = mat_opp_adj/sum(mat_opp_adj)
  
  
  #6a1) Run through pollen dispersal function
  if (space_pollen=="yes"){ #does pollen have limited dispersal?
    
    distance1 = matrix(nrow=pop_size, ncol=pop_size)
    distance2 = matrix(nrow=pop_size, ncol=pop_size)
    
    for (m in 1:pop_size) {
      for (p in 1:pop_size) {
        distance1[m,p] = sqrt((parents$X_pos[m]-parents$X_pos[p])^2 + (parents$Y_pos[m] - parents$Y_pos[p])^2)	
      } #distance is just the Euclidean distance between two points
    }
    distance1 = distance1*m_per_unit
    distance2[,] = 0 #when not on a torus, distance is not considered in the second direction
    
    distance_for_mating1 = distance1
    distance_for_mating2 = distance2
    
    pollen_dispersal1 = (((0.375)^2)/(2*pi))*exp(-0.375*distance_for_mating1) #if space_pollen is no, everyone gets the same value
    pollen_dispersal2 = 0 #everyone gets the same value here as well
    
  } else if (space_pollen=="no") {
    pollen_dispersal1 = matrix(nrow=pop_size, ncol=pop_size)
    pollen_dispersal2 = matrix(nrow=pop_size, ncol=pop_size)
    pollen_dispersal1[,] = 1
    pollen_dispersal2[,] = 0
  }
  
  mat_opp_adj2 = mat_opp_adj * (pollen_dispersal1 + pollen_dispersal2) #because pollen can be dispersed in either direction 1 OR direction 2, the probabilities are additive
  mat_opp_adj2 = mat_opp_adj2/sum(mat_opp_adj2) #readjusts the matrix to a probability

  offspring = as.data.frame(matrix(nrow=pop_size, ncol=ncol(parents)))
  names(offspring) = names(parents)
  
  offspring$position = parents$position
  offspring$X_pos = parents$X_pos
  offspring$Y_pos = parents$Y_pos
  
  for (o in 1:pop_size) {
    
    #Calculate dispersal of seeds to location of offspring[o] as a weight on relative fitness and distance
    if (space_seed == "yes"){
      seed_dispersal1 = (exp(-(seed_lambda*distance1[o,])))
      seed_dispersal = seed_dispersal1
    } else if (space_seed == "no") {
      seed_dispersal1 = rep(1, times=pop_size)
      seed_dispersal = seed_dispersal1
    } 
    
    #Choose a mom based on 1) distance to seed location and 2) relative amount of pollen received
    dist.prob = (seed_dispersal)/sum(seed_dispersal)
    pollen.prob = rowSums(mat_opp_adj2)/sum(rowSums(mat_opp_adj2))
    joint.prob = (dist.prob*pollen.prob)/sum(dist.prob*pollen.prob)
    #Choose a mom based on relative fitness and chance of getting a seed to the specified location

    (mom = sample(c(1:pop_size), size=1, prob=joint.prob))
    
    #If there are no sires for chosen mom, choose another mom
    while (rowSums(mat_opp_adj2)[mom]==0) {
      mom = sample(c(1:pop_size), size=1, prob=joint.prob)
    }
    
    #Choose sire for selected mom
    pat_prob = mat_opp_adj2[mom,]/rowSums(mat_opp_adj2)[mom]
    dad = sample(x=c(1:pop_size), size=1, prob=pat_prob)
    
    #9b5) Record mom, dad, and inbreeding coefficient
    offspring$Mother[o] = mom
    offspring$Father[o] = dad
    offspring$IB_coef[o] = coancestry_p[mom,dad]
    
    #A-locus inhertiance (one from mom, one from dad)
    offspring$A1[o] = as.character(sample(x=parents[mom, 5:6], size=1))
    offspring$A2[o] = as.character(sample(x=parents[dad, 5:6], size=1))
    
    #B-locus inhertiance (one from mom, one from dad)
    offspring$B1[o] = as.character(sample(x=parents[mom, 7:8], size=1))
    offspring$B2[o] = as.character(sample(x=parents[dad, 7:8], size=1))
    
    #C-locus inhertiance (one from mom, one from dad)
    offspring$C1[o] = as.character(sample(x=parents[mom, 9:10], size=1))
    offspring$C2[o] = as.character(sample(x=parents[dad, 9:10], size=1))
    
    #Mendellian inheritance of flowering time loci 
    for (n in 1:(FT_loci)) {
      offspring[o,11 + 2*(n-1)] = sample(x=parents[mom,(11+2*(n-1)):(12 + 2*(n-1))], size=1)
      offspring[o,12 + 2*(n-1)] = sample(x=parents[dad,(11+2*(n-1)):(12 + 2*(n-1))], size=1)
    } #close loop over the flowering time loci
    
    #Mendellian inheritance of remaining neutral loci
    for (j in 1:(neut_loci)) {
      offspring[o,11+(2*FT_loci)+(2*(j-1))] = as.character(sample(x=parents[mom,(11+(2*FT_loci)+(2*(j-1))):(12+(2*FT_loci)+(2*(j-1)))], size=1))
      offspring[o,12+(2*FT_loci)+(2*(j-1))] = as.character(sample(x=parents[dad,(11+(2*FT_loci)+(2*(j-1))):(12+(2*FT_loci)+(2*(j-1)))], size=1))
    } #close loop over remaining neutral loci

  } #close loop over offspring
  
  offspring$FTnoise = rnorm(n=pop_size, mean=0, sd=fl_noise)
  offspring$FLday = floor(fl_mean + rowSums(offspring[,11:(10+2*FT_loci)]) + offspring$FTnoise) 
  
  coancestry_o = matrix(nrow=pop_size, ncol=pop_size)
  
  for (i in 2:pop_size){
    for (j in 1:(i-1)){	
      Mi = offspring$Mother[i]
      Mj = offspring$Mother[j]
      Fi = offspring$Father[i]
      Fj = offspring$Father[j]
      
      pathA = (1/4)*coancestry_p[Mi, Mj]
      pathB = (1/4)*coancestry_p[Fi, Mj]
      pathC = (1/4)*coancestry_p[Mi, Fj]
      pathD = (1/4)*coancestry_p[Fi, Fj]
      
      coancestry_o[i, j] = (pathA + pathB + pathC + pathD)
      coancestry_o[j, i] = (pathA + pathB + pathC + pathD)
    }
  }
  
  for (i in 1:pop_size){
    coancestry_o[i,i] = (1/2)*(1 + offspring$IB_coef[i])
  }
  
  if ((g == 1)||(g == 10)||(g == 20)||(g == 30)||(g == 40)||(g == 50)||(g == 60)||(g == 70)||(g == 80)||(g == 90)||(g == 100)||(g == 150)||(g == 200)||(g == 250)||(g == 300)||(g == 350)||(g == 400)||(g == 450)||(g == 500)){
    
    offspring_map = as.data.frame(matrix(nrow=pop_size, ncol=(11+(2*FT_loci)+(2*neut_loci))))
    names(offspring_map)[1:11] = c("FLday","X_pos","Y_pos","genotypeA","mapA","genotypeB","mapB","genotypeC","mapC","Mother","Father")
    
    offspring_map$FLday = offspring$FLday
    offspring_map$X_pos = offspring$X_pos
    offspring_map$Y_pos = offspring$Y_pos
    offspring_map$Mother = offspring$Mother
    offspring_map$Father = offspring$Father
    
    for (n in 1:FT_loci)
    {
      names(offspring_map)[12+2*(n-1)] = paste("loc",n,"a",sep="")
      names(offspring_map)[13+2*(n-1)] = paste("loc",n,"b",sep="")
    }
    
    for (i in 1:pop_size){
      for (n in 1:FT_loci)
      {
        offspring_map[i,(12+2*(n-1))] = offspring[i,(11+2*(n-1))]
        offspring_map[i,(13+2*(n-1))] = offspring[i,(12+2*(n-1))]
      }
    }
    
    for (j in 1:neut_loci)
    {
      names(offspring_map)[12+(2*FT_loci)+(2*(j-1))] = paste("neut",j,"a",sep="")
      names(offspring_map)[13+(2*FT_loci)+(2*(j-1))] = paste("neut",j,"b",sep="")
    }
    
    for (i in 1:pop_size){
      for (j in 1:neut_loci)
      {
        offspring_map[i,12+(2*FT_loci)+(2*(j-1))] = offspring[i,12+(2*FT_loci)+(2*(j-1))]
        offspring_map[i,13+(2*FT_loci)+(2*(j-1))] = offspring[i,13+(2*FT_loci)+(2*(j-1))]
      }
    }
    
    for (i in 1:pop_size){
      offspring_map$genotypeA[i] = paste(offspring$A1[i], offspring$A2[i], sep="")	
      offspring_map$mapA[i] = if (offspring_map$genotypeA[i]=="aa") {
        0} else if (offspring_map$genotypeA[i]=="AA"){2} else {1}
    }
    
    for (i in 1:pop_size){
      offspring_map$genotypeB[i] = paste(offspring$B1[i], offspring$B2[i], sep="")	
      offspring_map$mapB[i] = if (offspring_map$genotypeB[i]=="bb") {
        0} else if (offspring_map$genotypeB[i]=="BB"){2} else {1}
    }
    
    for (i in 1:pop_size){
      offspring_map$genotypeC[i] = paste(offspring$C1[i], offspring$C2[i], sep="")	
      offspring_map$mapC[i] = if (offspring_map$genotypeC[i]=="cc") {
        0} else if (offspring_map$genotypeC[i]=="CC"){2} else {1}
    }
    
    write.csv(offspring_map, paste("paraset_",x,"_offspring_map_", g, ".csv",sep="")
              , row.names=F)
    
  }
  #end loop over g if
  
  if (g<no_gen){
    offspring->parents
    coancestry_o->coancestry_p}
}
