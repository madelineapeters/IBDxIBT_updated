library(dplyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)
library(factoextra)
library(cluster)
library(NbClust)


para_set = 16
run = 1
g = 500

df = read.csv(paste(getwd(),'/para_set_',para_set,'/model_run_',run,'/paraset_',para_set,'_offspring_map_',g,'.csv',sep=""))
neutral.df = df %>% select(.,FLday,X_pos,Y_pos,mapA,mapB,mapC,neut1a:neut24b)
neutral.df[] = lapply(neutral.df, as.character)
neutral.df[neutral.df == 'D'] = 1; 
neutral.df[neutral.df[,] == 'd'] = 0
neutral.df[] = lapply(neutral.df, as.numeric)
neutral.df = neutral.df %>% 
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

############ K-means clustering ###############
set.seed(123)
km.res = kmeans(df.scaled, 2, nstart = 25)
# k-means group number of each observation
km.res$cluster
# Visualize k-means clusters
fviz_cluster(km.res, data = df.scaled, geom = "point",
             stand = FALSE, ellipse.type = "norm")+theme_classic()+theme(legend.position='none')

############ Hierarchical clustering ###############
# Compute pairewise distance matrices
dist.res = dist(df.scaled, method = "euclidean")

# Hierarchical clustering results
hc = hclust(dist.res, method = "complete")
# Visualization of hclust
plot(hc, labels = FALSE, hang = -1)
# Add rectangle around 3 groups
rect.hclust(hc, k = 2, border = 2:4) 

############ Optimal number of clusters ###############
# Goal: define clusters such that the total intra-cluster variation (known as total within-cluster variation or total within-cluster sum of square) is minimized

# Elbow-method
    #1 Compute clustering algorithm (e.g., k-means clustering) for different values of k. For instance, by varying k from 1 to 10 clusters
    #2 For each k, calculate the total within-cluster sum of square (wss)
    #3 Plot the curve of wss according to the number of clusters k.
    #4 The location of a bend (knee) in the plot is generally considered as an indicator of the appropriate number of clusters.

set.seed(123)
# Compute and plot wss for k = 2 to k = 15
k.max = 15 # Maximal number of clusters
data = df.scaled
wss = sapply(1:k.max, 
              function(k){kmeans(data, k, iter.max=20, nstart=25)$tot.withinss})

plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
abline(v = 3, lty =2)


# Average silhouette method
    #1 Compute clustering algorithm (e.g., k-means clustering) for different values of k. For instance, by varying k from 1 to 10 clusters
    #2 For each k, calculate the average silhouette of observations (avg.sil)
    #3 Plot the curve of avg.sil according to the number of clusters k.
    #4 The location of the maximum is considered as the appropriate number of clusters.

k.max = 15
data = df.scaled
sil = rep(0, k.max)

# Compute the average silhouette width for 
# k = 2 to k = 15
for(i in 2:k.max){
  km.res = kmeans(data, centers = i, iter.max=20, nstart = 25)
  ss = silhouette(km.res$cluster, dist(data))
  sil[i] = mean(ss[, 3])
}

# Plot the  average silhouette width
plot(1:k.max, sil, type = "b", pch = 19, 
     frame = FALSE, xlab = "Number of clusters k")
abline(v = which.max(sil), lty = 2)

fviz_nbclust(df.scaled, hcut, method = "silhouette",
             hc_method = "complete")

