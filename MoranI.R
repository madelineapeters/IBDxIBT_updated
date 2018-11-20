library(ape)

para_set = 8
run = 1
g = 100

g.list = c(1,10,20,30,49,50,60,70,80,90,100)

#Parameters for calculating temporal distance
fl_seq = c(2, 3, 5, 8, 11, 13, 15, 17, 15, 13, 11, 8, 5, 3, 2)	# of flowers per plant for each day of flowering (here, 131 per plant, distributed over 15 days)
duration = length(fl_seq)

for (g in g.list){
  
  df = read.csv(paste(getwd(),'/para_set_',para_set,'/model_run_',run,'/paraset_',para_set,'_offspring_map_',g,'.csv',sep=""))
  
  ####Spatial distance####
  spat.dists = as.matrix(dist(cbind(df$X_pos, df$Y_pos)))
  spat.dists = spat.dists/max(spat.dists)
  spat.dists.inv = 1/spat.dists
  diag(spat.dists.inv) = 0
  
  ####Spatial and temporal distance####
  temp.dists = matrix(nrow=nrow(df),ncol=nrow(df))
  for (f in 1:nrow(df)){
  temp.dists[f,] = abs(df$FLday[f] - df$FLday)
  }
  temp.dists = temp.dists/max(temp.dists)
  temp.dists.inv = 1/temp.dists
  temp.dists.inv[temp.dists.inv == Inf] = 0

  spat.temp.dists = sqrt((spat.dists^2)+(temp.dists^2))
  spat.temp.dists.inv = 1/spat.temp.dists
  diag(spat.temp.dists.inv) = 0
  
  spat.Moran.I = Moran.I(df$mapB,spat.dists.inv)
  spat.temp.Moran.I = Moran.I(df$mapB,temp.dists)
  
}
df = read.csv(paste(getwd(),'/para_set_',para_set,'/model_run_',run,'/paraset_',para_set,'_offspring_map_',g,'.csv',sep=""))

spat.dists = as.matrix(dist(cbind(df$X_pos, df$Y_pos)))
spat.dists.inv = 1/spat.dists
diag(spat.dists.inv) = 0

Moran.I(df$mapC,spat.dists.inv)
