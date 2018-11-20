#####Set-up#####
#load required packages
library(tidyr)
library(dplyr)

#set base working directory
base_wd<-"~/498"

#set parameters
para<-78
gen<-50
pop_size<-400
dispersal_distance<-1
d<-(1/(dispersal_distance)^2) #density for neighbourhood size equation

#make dataframe to store output
N.gen<-as.data.frame(matrix(nrow=gen, ncol=2))
names(N.gen)<-c("gen", "size")
N.gen$gen<-(1:gen)


#####Set up positions dataframe#####
#From the matrix of available positions in space, develop a list of positions from which to sample
pos<-as.data.frame(matrix(nrow=pop_size))
space_20by20_sq_grid<-read.csv("~/498/space_20by20_sq_grid.csv")
space_20by20_sq_grid<-space_20by20_sq_grid[,-1]

for (y in 1:ncol(space_20by20_sq_grid)) {
  if (y ==1 ) {
    pos_vector<-c(as.character(space_20by20_sq_grid[,y])) }
  else {
    pos_vector<-c(pos_vector, as.character(space_20by20_sq_grid[,y])) 
  }
}

#Give everyone a position from the set of possible positions, and then extract their x and y coordinates
pos$position<-(pos_vector)	
pos$X_pos<-sapply(X=pos$position, FUN=function(X) {
  as.numeric(strsplit(X, split=",")[[1]][1])
})
pos$Y_pos<-sapply(X=pos$position, FUN=function(X) {
  as.numeric(strsplit(X, split=",")[[1]][2])
})	


#####Start loop over generations#####
for (g in 1:gen) {
  #read in parentage dataframe
  NH<-read.csv(paste(base_wd, paste("para_set", para, sep="_"), "model_run_1", paste("NH_", g, ".csv", sep=""), sep="/"))
  
  #fill in offspring coordinates
  NH.2<-as.data.frame(matrix(nrow=pop_size, ncol=6))
  names(NH.2)<-c("O.x", "O.y", "M.x", "M.y", "P.x", "P.y")
  NH.2$O.x<-sapply(X=pos$position, FUN=function(X) {
    as.numeric(strsplit(X, split=",")[[1]][1])
  })
  NH.2$O.y<-sapply(X=pos$position, FUN=function(X) {
    as.numeric(strsplit(X, split=",")[[1]][2])
  })
  
  #fill in mom coordinates
  for (m in 1:pop_size) {
    mom.x<-NH$Mother[m]
    NH.2$M.x[m]<-sapply(X=pos$position[mom.x], FUN=function(X) {
      as.numeric(strsplit(X, split=",")[[1]][1])
    })
    NH.2$M.y[m]<-sapply(X=pos$position[mom.x], FUN=function(X) {
      as.numeric(strsplit(X, split=",")[[1]][2])
    })
    
  }
  
  #fill in dad coorindates
  for (p in 1:pop_size) {
    dad.x<-NH$Father[p]
    NH.2$P.x[p]<-sapply(X=pos$position[dad.x], FUN=function(X) {
      as.numeric(strsplit(X, split=",")[[1]][1])
    })
    NH.2$P.y[p]<-sapply(X=pos$position[dad.x], FUN=function(X) {
      as.numeric(strsplit(X, split=",")[[1]][2])
    })
    
  }
  
  #calculate Euclidean distances between parent and offspring
  NH.3<-bind_cols(NH.2, NH)
  for (m in 1:pop_size){
    NH.3$mom.dist[m]<-sqrt((NH.3$M.x[m]-NH.3$O.x[m])^2+(NH.3$M.y[m]-NH.3$O.y[m])^2)
  }
  for (m in 1:pop_size){
    NH.3$dad.dist[m]<-sqrt((NH.3$P.x[m]-NH.3$O.x[m])^2+(NH.3$P.y[m]-NH.3$O.y[m])^2)
  }
  
  #calculate the mean distance
  mean.dist<-sum(mean(NH.3$mom.dist), mean(NH.3$dad.dist))/2
  
  #calculate variance for each offspring
  NH.3<-mutate(NH.3, ind.var=(mom.dist-mean.dist)^2+(dad.dist-mean.dist)^2)
  
  #calculate population variance
  var<-sum(NH.3$ind.var)/(2*pop_size)
  
  #calculate neighbourhood size
  N.size<-4*3.14*var*d
  
  #store generation neighbourhood size
  N.gen$size[g]<-N.size

  }#next generation

library(ggplot2)
N.plot<-ggplot(data=N.gen, aes(x=gen, y=size))+geom_line()+ylim(0, 400)
plot(N.plot)
