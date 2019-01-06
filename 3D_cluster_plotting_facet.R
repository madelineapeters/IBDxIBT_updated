library(dplyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)
library(factoextra)
library(cluster)
library(NbClust)
library(mclust)
library(rgl)

for (s in 17:24){
      
    r = 1  
    gen = 800
    
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
    
    ################################
    ## How k-means clustering works
    ################################
    #1 The process begins with k centroids initialised at random.
    #2 These centroids are used to assign points to its nearest cluster.
    #3 The mean of all points within the cluster is then used to update the position of the centroids.
    #4 The above steps are repeated until the values of the centroids stabilise.
    
    ############ K-means clustering ###############
    k.means = 2#d_clust$G
    
    km.res = kmeans(df.scaled, k.means, iter.max = 30,nstart = 25)
    
    # Visualize k-means clusters
    #fviz_cluster(km.res, data = df.scaled, geom = "point",
                 #stand = FALSE, ellipse.type = "norm")+theme_classic()+theme(legend.position='none')+ggtitle("")+
      #ggtitle(paste('Neutral genetic clusters in two-dimensional space (k = ',k.means,')',sep=""))+
      #geom_point(aes(col=factor(ind.neutral.df$FLday)))
    
    hist.df = df %>% bind_cols(.,as.data.frame(km.res$cluster))
    names(hist.df)[ncol(hist.df)] = 'Cluster'
    
    clust.avg = hist.df %>% group_by(FLday) %>% summarise(.,avg.clust = mean(Cluster))
    clust.sd = hist.df %>% group_by(FLday) %>% summarise(.,sd.clust = sd(Cluster))
    clust.summary = as.data.frame(right_join(clust.avg,clust.sd,"FLday"))
    clust.summary$Set = s
    
    if (s == 17){clust.joint = clust.summary} else {clust.joint = bind_rows(clust.joint,clust.summary)}
    
}

clust.joint$grouping = 'Selfing'
clust.joint$grouping[(1+nrow(clust.joint)/2):nrow(clust.joint)] = 'No selfing'
clust.joint$Set[(clust.joint$Set == 17)|(clust.joint$Set == 21)] = 'Random'
clust.joint$Set[(clust.joint$Set == 18)|(clust.joint$Set == 22)] = 'IBT'
clust.joint$Set[(clust.joint$Set == 19)|(clust.joint$Set == 23)] = 'IBD'
clust.joint$Set[(clust.joint$Set == 20)|(clust.joint$Set == 24)] = 'IBDxIBT'
names(clust.joint)[4] = 'Isolation'


ggplot()+geom_point(data=filter(clust.joint,Isolation %in% c('IBD','IBDxIBT')),aes(x=FLday,y=avg.clust))+theme_test()+ylim(1,2)+ylab("Average cluster assigned based on neutral genetic variation")+xlab("Middle day of individual's flowering schedule")+facet_grid(Isolation~grouping)+scale_y_continuous(name="Average cluster assigned based on neutral genetic variation",breaks=c(1,1.5,2))

# Spinning 3d Scatterplot
library(rgl)
c = km.res$cluster
cols = rainbow(n=7)[c]
plot3d(ind.neutral.df$X_pos, ind.neutral.df$Y_pos, ind.neutral.df$FLday, col=cols, size=3,
       xlab='',ylab='',zlab='')

hist.df = df %>% bind_cols(.,as.data.frame(km.res$cluster))
names(hist.df)[ncol(hist.df)] = 'Cluster'

ggplot(data=hist.df,aes(FLday)) + geom_histogram(aes(fill=factor(km.res$cluster)),position='dodge')+guides(fill=guide_legend(title="Neutral cluster"))+theme_classic()+ylab('Count')
