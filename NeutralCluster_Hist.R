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
run = 1
g = 800

for (s in 17:24){
  # Read in data and set up for k-means analysis
  df = read.csv(paste(getwd(),'/para_set_',s,'/model_run_',run,'/paraset_',s,'_offspring_map_',g,'.csv',sep=""))
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
  
  k.means = 2#d_clust$G
  
  km.res = kmeans(df.scaled, k.means, iter.max = 20,nstart = 25)
  
  hist.df = df %>% bind_cols(.,as.data.frame(km.res$cluster))
  names(hist.df)[ncol(hist.df)] = 'Cluster'
  hist.df$paraset = s
  
  if (s == 17){hist.joint = hist.df} else {hist.joint = bind_rows(hist.joint,hist.df)}
}

hist.joint$grouping = 'Selfing'
hist.joint$grouping[(1+nrow(hist.joint)/2):nrow(hist.joint)] = 'No selfing'
hist.joint$paraset[(hist.joint$paraset == 17)|(hist.joint$paraset == 21)] = 'Random'
hist.joint$paraset[(hist.joint$paraset == 18)|(hist.joint$paraset == 22)] = 'IBT'
hist.joint$paraset[(hist.joint$paraset == 19)|(hist.joint$paraset == 23)] = 'IBD'
hist.joint$paraset[(hist.joint$paraset == 20)|(hist.joint$paraset == 24)] = 'IBDxIBT'
names(hist.joint)[ncol(hist.joint)-1] = 'Isolation'

ggplot(data=hist.joint,aes(FLday)) + geom_histogram(aes(fill=factor(Cluster)),position='dodge')+guides(fill=guide_legend(title="Neutral cluster"))+theme_classic()+ylab('Count')+facet_grid(Isolation~grouping)
