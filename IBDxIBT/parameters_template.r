#Parameters 
para_set<-para$set[x]

#POPULATION PARAMETERS
no_runs<-para$no_runs[x]	#number of times to repeat this version of the model
no_gen<-para$no_gen[x]	#number of generations to run
pop_size<-para$pop_size[x]	#number of plants in the population (try 3 levels to affect ratio of allele number to pop size)
m_per_unit<-para$m_per_unit[x]	#scales the pairwise distance matrix... if m_per_unit =1, then plants are 1m apart from one another.

neut_loci<-para$neut_loci[x]
torus<-para$torus[x]
mut<-para$mut[x]
mut_rate<-para$u[x]
selfing<-para$selfing[x]
#FLOWERING TIME PARAMETERS
fl_mean<-100	#seed flowering time mean for parent generation

fl_var<-para$fl_var[x]	#seed flowering time variance for parents (higher variance means higher phenological assortative mating)

fl_noise<-1.5 #if (para$fl_noise[x]=="none") {0} else if (para$fl_noise[x]=="low") {1.5} else if (para$fl_noise[x]=="high") {5}

FT_loci<-para$FT_loci[x]

FT_alleles<-if (fl_var=="high") {c(-1, 1)} else if (fl_var=="low") {c(-0.5, 0.5)}	#quantitative effect of alleles at each of FT_loci (assumes same allele # and effect size at each locus)

fl_seq<-c(2, 3, 5, 8, 11, 13, 15, 17, 15, 13, 11, 8, 5, 3, 2)	# of flowers per plant for each day of flowering (here, 131 per plant, distributed over 15 days)

duration<-length(fl_seq)

phenology_matrix<-para$phenology_matrix[x]	#do we calculate effect of phenology on mating opportunities, yes or no?


#S-LOCUS PARAMETERS
no_S_alleles<-para$no_S_alleles[x]	#number of SI alleles

#probability of choosing each allele  (modified to accommodate 4 sets of allele probabilities, for scenarios when alleles do not start w/ random distribution over time)

allele_prob<-para$allele_prob[x]
	
stigma_expr<-0
pollen_expr<-0

#rec_rate<-para$rec_rate[x]	#recombination rate between S and the flowering time locus 1 (assume S1 is linked to loc1a, and S2 to loc1b).  ranges from 0.5 (no linkage, pure mendellian segregation) to 0 (100% linked, no recombination between S & loc1)


#THINGS AFFECTING MATING OPPORTUNITIES
compatibility<-para$compatibility[x]	#do we calculate effect of s-alleles on mating probabilities?  yes or no?

space_pollen<-para$space_pollen[x]  #do we calculate effect of space on mating probabilities? **SET TO NO FOR A-LOCUS VERSION OF THINGS!!!!

phenology<-para$phenology[x] #do we calculate effect of phenology on mating probabilities?

LA_inbreeding_flowers<-para$LA_inbreeding_flowers[x]	#Does inbreeding affect flower production? (affects M & F fitness) (yes or no)
fl_id_weight<-para$fl_id_weight[x]	#Penalty of inbreeding for flower production (fl = fl*fl_id_weight*(1-IBcoef)

LA_inbreeding_fecundity<-para$LA_inbreeding_fecundity[x]	#Does inbreeding affect relative female fitness? (yes or no)

EA_inbreeding<-para$EA_inbreeding[x]	#Does inbreeding affect probability of a seed making it to next generation? (yes or no)
EA_id_threshold<-para$EA_id_threshold[x]	#Maximum tolerated inbreeding coefficient before seed is rejected if EA_inbreeding == "yes"


#A LOCUS PARAMETERS (locus under selection)
freqA<-0.5	#frequency of big A allele.  (freq little a = 1-freqA)

selection_aa<-para$selection_aa[x]	#selection differential on FLday for aa genotype
selection_AA<-para$selection_AA[x]	#Should prob be exact opposite of selection_aa, though I could try making stronger or weaker.
selection_Aa<-0

#B LOCUS
freqB<-0.5

#C LOCUS
freqC<-0.5

freqD<-0.5


#SEED DISPERSAL
space_seed<-para$space_seed[x]
seed_lambda<-2
