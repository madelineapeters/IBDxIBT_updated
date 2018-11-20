library(dplyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)

para_set = 4 #parameter set
run = 1 #model run
g.list = c(1,10,20,30,49,50,60,70,80,90,100)

#Frequency plots
for (g in g.list){
  df = read.csv(paste(getwd(),'/para_set_',para_set,'/model_run_',run,'/paraset_',para_set,'_offspring_map_',g,'.csv',sep=""))
  flower.df = df %>% select(.,loc1a:loc5b)
  flower.df[flower.df == -1] = 0
  flower.df[flower.df == 1] = 0.5
  flower.df = flower.df %>% 
    mutate(.,map1 = loc1a+loc1b) %>% 
    mutate(.,map2 = loc2a+loc2b) %>% 
    mutate(.,map3 = loc3a+loc3b) %>% 
    mutate(.,map4 = loc4a+loc4b) %>% 
    mutate(.,map5 = loc5a+loc5b)
  flower.df = flower.df %>% select(.,map1:map5)
  flower.freq = colSums(flower.df)/(nrow(flower.df))
  
  if (g == 1){
    map.df = as.data.frame(flower.freq) %>% rownames_to_column()
    map.df$Generation = g
    names(map.df) = c('Locus','Frequency','Generation')
  } else {
    temp.df = as.data.frame(flower.freq) %>% rownames_to_column()
    temp.df$Generation = g
    names(temp.df) = c('Locus','Frequency','Generation')
    map.df = bind_rows(map.df,temp.df)
  }
}

ggplot()+geom_line(data=map.df,aes(x=Generation,y=Frequency,col=Locus))+
  theme_bw()+ylim(0,1)

greencols = brewer.pal(9, 'Greens')
newcol = colorRampPalette(greencols)
ncols = 100
greencols2 = newcol(ncols)#apply the function to get 100 colours


#Heatmaps
for (g in g.list){
  df = read.csv(paste(getwd(),'/para_set_',para_set,'/model_run_',run,'/paraset_',para_set,'_offspring_map_',g,'.csv',sep=""))
  
  heatmap.matFlower = matrix(nrow=sqrt(nrow(df)),ncol=sqrt(nrow(df)))
 
  
  for (x in 1:nrow(df)){
    heatmap.matFlower[df$X_pos[x],df$Y_pos[x]] = df$FLday[x]
  }
  
  heatmap(heatmap.matFlower,col=greencols2,Rowv=NA,Colv=NA,labRow=NA,labCol=NA,main=paste('Generation',g,sep=" "))
  
}

