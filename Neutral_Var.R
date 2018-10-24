library(dplyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)

run = 1 #model run
g.list = c(1,10,20,30,50,60,70,80,90,100,150,200,250,300,350,400,450,500)

Mantel.obs = as.data.frame(matrix(nrow=length(g.list),ncol=8))
names(Mantel.obs) = sapply(1:8, function(X) paste('paraset',X,sep="_"))

##############################################################
## Calculate spatial autocorrelation statistics
##############################################################
for (s in 1:8){
  for (g in 1:length(g.list)){
    
    gen = g.list[g]
    
  df = read.csv(paste(getwd(),'/para_set_',s,'/model_run_',run,'/paraset_',s,'_offspring_map_',gen,'.csv',sep=""))
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
  
  neutral.freq = colSums(neutral.df)/(2*nrow(neutral.df))
  
  if (g == 1){
    map.df = as.data.frame(neutral.freq) %>% rownames_to_column()
    map.df$Generation = g
    names(map.df) = c('Locus','Frequency','Generation')
  } else {
    temp.df = as.data.frame(neutral.freq) %>% rownames_to_column()
    temp.df$Generation = g
    names(temp.df) = c('Locus','Frequency','Generation')
    map.df = bind_rows(map.df,temp.df)
  }

  }
  
  if (s == 1){
    var.final = map.df
    var.final$paraset = s
  } else {
    var.temp = map.df
    var.temp$paraset = s
    
    var.final = bind_rows(var.final,var.temp)
  }
}
var.final$grouping[1:1944] = 'Selfing'
var.final$grouping[1945:3888] = 'No selfing'
var.final$paraset[(var.final$paraset == 1)|(var.final$paraset == 5)] = 'Null'
var.final$paraset[(var.final$paraset == 2)|(var.final$paraset == 6)] = 'IBT'
var.final$paraset[(var.final$paraset == 3)|(var.final$paraset == 7)] = 'IBD'
var.final$paraset[(var.final$paraset == 4)|(var.final$paraset == 8)] = 'IBDxIBT'
names(var.final)[4] = 'Isolation'
summary.df = var.final %>% group_by(.,Generation,grouping,Isolation) %>% summarise(var(Frequency))

ggplot()+geom_line(data=summary.df,aes(x=Generation,y=`var(Frequency)`,col=Isolation))+
  theme_classic()+ylim(0,0.05)+ylab('Variance')+facet_grid(.~grouping)+ggtitle('Variance in allele frequency at neutral loci')
