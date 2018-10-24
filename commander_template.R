####### Phenological and spatial assortment, mutation, torus, neutral loci
#######M. Peters 26 April 2017
library(gtools)

args = commandArgs(trailingOnly = TRUE)
x = as.integer(args[1])

#COMMANDER 

base_wd<-getwd()

#Enter the proper names of the required files below, DO NOT INCLUDE PATH TO FILES
para_file<-"Test_Para"  #DO NOT end this one with .csv

para_script<-"parameters_template.r"  #DO end this one with .R

model_script<-"model_template.r"  #DO end this one with .R

#space_50by50_sq_grid<-read.csv(paste(getwd(),"space_50by50_sq_grid.csv",sep="/"))
#space_50by50_sq_grid<-space_50by50_sq_grid[,-1]

#Bring in table of parameter sets for model iterations
para<-read.csv(as.character(paste(base_wd, paste(para_file,".csv", sep=""), sep="/")))

#Start loop over rows of parameter table
#for (x in 1:nrow(para)){
  
#read in the parameters
source(file=paste(base_wd, para_script, sep="/"))
xy.coor = permutations(sqrt(pop_size),2,repeats.allowed=FALSE)
diagonal.coor = bind_cols(as.data.frame(1:30),as.data.frame(1:30))
names(diagonal.coor) = c('V1','V2')
xy.coor = bind_rows(diagonal.coor,as.data.frame(xy.coor))
##Run the model # of times specified

#Start by re-setting the wd
setwd(base_wd)
  
#run the model
source(file=paste(base_wd, model_script, sep="/"))
  
#Save files

write.csv(output, "data_by_generation.csv")

write.csv(onset_tab,
          "Flowering_time_onsets_over_generations.csv", row.names=F)

write.csv(Ageno_tab,
          "A_locus_genotypes_over_generations.csv", row.names=F)

write.csv(Bgeno_tab,
          "B_locus_genotypes_over_generations.csv", row.names=F)

write.csv(Cgeno_tab,
          "C_locus_genotypes_over_generations.csv", row.names=F)


#}#close loop over parameter sets
  
