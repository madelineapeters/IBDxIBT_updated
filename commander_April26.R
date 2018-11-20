####### Phenological and spatial assortment, mutation, torus, neutral loci
#######M. Peters 26 April 2017

#COMMANDER 

base_wd<-"C:/Users/Madeline/Desktop/Weis lab/EEB498"
###This needs to be the working directory where all files are going to be saved.
#Note- this does not mean that all files used for this version of commander have to be stored here.
#Shared subfiles folders are still ok.


#Enter the proper names of the required files below, DO NOT INCLUDE PATH TO FILES
para_file<-"April_26_para40"  #DO NOT end this one with .csv

para_script<-"parameters_April26.R"  #DO end this one with .R

model_script<-"model_April26.R"  #DO end this one with .R


space_50by50_sq_grid<-read.csv("C:/Users/Madeline/Desktop/Weis lab/EEB498/space_50by50_sq_grid.csv")
space_50by50_sq_grid<-space_50by50_sq_grid[,-1]


#Bring in table of parameter sets for model iterations
para<-read.csv(as.character(paste(base_wd, paste(para_file,".csv", sep=""), sep="/")))

#Start loop over rows of parameter table
for (x in 1:nrow(para)){
  
#read in the parameters
source(file=paste(base_wd, para_script, sep="/"))


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


}#close loop over runs
  
