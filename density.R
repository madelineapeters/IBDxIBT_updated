NH<-subset.data.frame(offspring, select=c(FLday, Mother, Father, X_pos, Y_pos))

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
ind.df<-as.data.frame(matrix(nrow=1, ncol=1))
ind.df.t<-as.data.frame(matrix(nrow=1, ncol=1))
#start loop over sampled offspring
for (w in off.samp){
  a<-NH$X_pos[w]
  b<-NH$Y_pos[w]
  FLday.o<-NH$FLday[w]
  p.x<-round(sqrt(var))+a #round so can translate to cartesian coordinates
  n.x<-a-round(sqrt(var))
  for (x in n.x:p.x){
    p.y<-round(sqrt(var-(x-a)^2))+b
    n.y<-b-round(sqrt(var-(x-a)^2))
    for (z in n.y:p.y) {
      NH.temp<-filter(NH, X_pos==x, Y_pos==z)
      FLday.t<-NH.temp$FLday[1]
      dif.t<-abs(FLday.o-FLday.t)
      dif.t<-(15-dif.t)/15
      ind.df.t[1,1]<-dif.t
      ind.df<-bind_rows(ind.df, ind.df.t)
    } #next y position
  } #next x position
  off.df[match(c(w),off.samp),]<-sum(ind.df[,1], na.rm=TRUE)/length(ind.df[,1])
} #next offspring