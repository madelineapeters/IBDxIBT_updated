library(dplyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)
library(factoextra)
library(cluster)
library(NbClust)
library(mclust)
library(rgl)

# Set parameters
para_set = 9
run = 8
g = 800

# Read in data and set up for k-means analysis
df = read.csv(paste(getwd(),'/IBDxIBT','/para_set_',para_set,'/model_run_',run,'/paraset_',para_set,'_offspring_map_',g,'.csv',sep=""))
neutral.df = df %>% select(.,FLday,X_pos,Y_pos,mapA,mapB,mapC,loc1a:loc5b,neut1a:neut24b)
neutral.df[] = lapply(neutral.df, as.character)
neutral.df[neutral.df == 'D'] = 1; 
neutral.df[neutral.df[,] == 'd'] = 0
neutral.df[] = lapply(neutral.df, as.numeric)
neutral.df = neutral.df %>% 
  mutate(.,F1 = loc1a+loc1b) %>% 
  mutate(.,F2 = loc2a+loc2b) %>% 
  mutate(.,F3 = loc3a+loc3b) %>% 
  mutate(.,F4 = loc4a+loc4b) %>% 
  mutate(.,F5 = loc5a+loc5b) %>% 
  mutate(.,map1 = neut1a+neut1b) %>% 
  mutate(.,map2 = neut2a+neut2b) %>% 
  mutate(.,map3 = neut3a+neut3b) %>% 
  mutate(.,map4 = neut4a+neut4b) %>% 
  mutate(.,map5 = neut5a+neut5b) %>% 
  mutate(.,map6 = neut6a+neut6b) %>% 
  mutate(.,map7 = neut7a+neut7b) %>% 
  mutate(.,map8 = neut8a+neut8b) %>% 
  mutate(.,map9 = neut9a+neut9b) %>% 
  mutate(.,map10 = neut10a+neut10b) %>% 
  mutate(.,map11 = neut11a+neut11b) %>% 
  mutate(.,map12 = neut12a+neut12b) %>% 
  mutate(.,map13 = neut13a+neut13b) %>% 
  mutate(.,map14 = neut14a+neut14b) %>% 
  mutate(.,map15 = neut15a+neut15b) %>% 
  mutate(.,map16 = neut16a+neut16b) %>% 
  mutate(.,map17 = neut17a+neut17b) %>% 
  mutate(.,map18 = neut18a+neut18b) %>% 
  mutate(.,map19 = neut19a+neut19b) %>% 
  mutate(.,map20 = neut20a+neut20b) %>% 
  mutate(.,map21 = neut21a+neut21b) %>% 
  mutate(.,map22 = neut22a+neut22b) %>% 
  mutate(.,map23 = neut23a+neut23b) %>% 
  mutate(.,map24 = neut24a+neut24b)

ind.neutral.df = neutral.df %>% select(.,FLday,X_pos,Y_pos,mapA:mapC,map1:map24)

neutral.df = neutral.df %>% select(.,mapA:mapC,map1:map24)
df.scaled = scale(neutral.df)
scaled.matrix = as.matrix(df.scaled)

################################
## How k-means clustering works
################################
#1 The process begins with k centroids initialised at random.
#2 These centroids are used to assign points to its nearest cluster.
#3 The mean of all points within the cluster is then used to update the position of the centroids.
#4 The above steps are repeated until the values of the centroids stabilise.

#################################
## How the elbow method works
#################################
# Looks at the percentage of variance explained as a function of the number of clusters
# Should choose a number of clusters so that adding another cluster doesn’t give much better modeling of the data
# If one plots the percentage of variance explained by the clusters against the number of clusters, the first clusters will add much information (explain a lot of variance)
# At some point the marginal gain will drop, giving an angle in the graph
# The number of clusters is chosen at this point, hence the “elbow criterion"
# This “elbow” cannot always be unambiguously identified.

set.seed(123)

comment.out = function(){
  # Compute and plot wss for k = 2 to k = 15
  k.max = 15 # Maximal number of clusters
  data = df.scaled
  wss = sapply(1:k.max, 
               function(k){kmeans(data, k, iter.max=20, nstart=25)$tot.withinss})
  
  plot(1:k.max, wss,
       type="b", pch = 19, frame = FALSE, 
       xlab="Number of clusters K",
       ylab="Total within-clusters sum of squares")
}

#################################
## How silhouette analysis works
#################################
# Used to study the separation distance between the resulting clusters
# The silhouette plot displays a measure of how close each point in one cluster is to points in the neighboring clusters
# Thus provides a way to assess parameters like number of clusters visually
# This measure has a range of [-1, 1]

# Silhouette coefficients near +1 indicate that the sample is far away from the neighboring clusters
# A value of 0 indicates that the sample is on or very close to the decision boundary between two neighboring clusters and negative values indicate that those samples might have been assigned to the wrong cluster

comment.out = function(){
  k.max = 15
  data = df.scaled
  sil = rep(0, k.max)
  
  # Compute the average silhouette width for 
  # k = 2 to k = 15
  for(i in 2:k.max){
    km.res = kmeans(data, centers = i, iter.max=30, nstart = 25)
    ss = silhouette(km.res$cluster, dist(data))
    sil[i] = mean(ss[, 3])
  }
  
  # Plot the  average silhouette width
  plot(1:k.max, sil, type = "b", pch = 19, 
       frame = FALSE, xlab = "Number of clusters k")
  abline(v = which.max(sil), lty = 2)
  k.best = which.max(sil)
  
  fviz_nbclust(df.scaled, hcut, method = "silhouette",
               hc_method = "complete")
}

#################################
## BIC for k-means
#################################
# Determine the optimal model and number of clusters according to the Bayesian Information Criterion for expectation-maximization, initialized by hierarchical clustering for parameterized Gaussian mixture models
# Set the modelNames parameter to mclust.options(“emModelNames”) so that it includes only those models for evaluation where the number of observations is greater than the dimensions of the dataset here 402>5. 

d_clust = Mclust(scaled.matrix, G=1:4, 
                  modelNames = mclust.options("emModelNames"))

BIC.plot = plot(d_clust, what = c("BIC"), 
               dimens = NULL, xlab = NULL, ylab = NULL, ylim = NULL, xlim = NULL,
               addEllipses = TRUE, main = FALSE)
title(main = "Spatial clustering of flowering time")
mtext(paste("BIC, generation", g, sep=" "))

############ K-means clustering ###############
k.means = d_clust$G

km.res = kmeans(df.scaled, k.means, iter.max = 20,nstart = 25)

# Visualize k-means clusters
fviz_cluster(km.res, data = df.scaled, geom = "point",
             stand = FALSE, ellipse.type = "norm")+theme_classic()+theme(legend.position='none')+
  ggtitle(paste('Neutral genetic clusters in two-dimensional space (k = ',d_clust$G,')',sep=""))+
  geom_point(aes(col=factor(ind.neutral.df$FLday)))

ggplot()+geom_tile(data=df,aes(x=X_pos,y=Y_pos,fill=km.res$cluster))+theme_classic()+theme(legend.position='none')+xlab('X position')+ylab('Y position')+ggtitle(paste('Neutral genetic clusters in two-dimensional space (k = ',d_clust$G,')',sep=""))+scale_fill_gradient(low='white',high='red')

ggplot()+geom_tile(data=df,aes(x=X_pos,y=Y_pos,fill=FLday))+theme_classic()+theme(legend.position='none')+xlab('X position')+ylab('Y position')+ggtitle(paste('Flowering day'))

# Spinning 3d Scatterplot
library(rgl)
c = km.res$cluster
cols = rainbow(n=7)[c]
plot3d(ind.neutral.df$X_pos, ind.neutral.df$Y_pos, ind.neutral.df$FLday, col=cols, size=3,
       xlab='',ylab='',zlab='')

