####### IBD x IBT simulation - stats commander
#######M. Peters 5 June 2017

#Set working directories
base_wd<-"G:"
#Needs to be the working directory where all files are going to be saved and where scripts are stored

sub_wd<-"final sets"
#This is the sub-folder holding the parameter sets

#Set parameters
pop_size<-400
no_runs<-10
no_neutral<-10
FLmean<-100

#Set values for looping over parameter sets
  #Non-torus
  starting.set1<-42
  ending.set1<-53
  
  #Torus
  starting.set2<-55
  ending.set2<-64

#Start by re-setting the wd
setwd(base_wd)
  
#Run 2D Mantel statistics
source(file=paste(base_wd, "sapply_2D_Mantel_neutral.R", sep="/")) #updated
source(file=paste(base_wd, "sapply_2D_Mantel_neutral_torus.R", sep="/")) #updated
source(files=paste(base_wd, "sapply_2D_Mantel_neutral_FLday.R", sep="/")) #updated
source(files=paste(base_wd, "sapply_2D_Mantel_neutral_FLday_torus.R", sep="/")) #updated

#Run 3D Mantel statistics
source(file=paste(base_wd, "sapply_Mantel_neutral.R", sep="/")) #updated
source(file=paste(base_wd, "sapply_Mantel_neutral_torus.R", sep="/")) #updated

#Run neutral variance and frequency statistics
source(file=paste(base_wd, "Neutral p.R", sep="/")) #updated
source(file=paste(base_wd, "Neutral var.R", sep="/")) #updated

#Run pairwise Fst statistics (already averaged)
source(file=paste(base_wd, "sapply_Fst_neutral.R", sep="/")) #updated

#Run averaged statistics
source(file=paste(base_wd, "Average FL mantel statistics.R", sep="/")) #updated
source(file=paste(base_wd, "Average 2D mantel statistics.R", sep="/")) #updated
source(file=paste(base_wd, "Average 3D mantel statistics.R", sep="/")) #updated
source(file=paste(base_wd, "Average neutral dev p.R", sep="/")) #updated
source(file=paste(base_wd, "Average neutral variance.R", sep="/")) #updated

#Make individual plots
source(file=paste(base_wd, "FL mantel plots.R", sep="/")) #updated

#Make facetted plots (do not save, so run individually and save manually OR include save script here)
source(file=paste(base_wd, "Total dev p plot.R", sep="/")) #updated
source(file=paste(base_wd, "Total Fst plot.R", sep="/")) #updated
source(file=paste(base_wd, "Total 2D mantel plot.R", sep="/")) #updated
source(file=paste(base_wd, "Total 3D mantel plot.R", sep="/")) #updated
source(file=paste(base_wd, "Total FL mantel plot.R", sep="/"))


