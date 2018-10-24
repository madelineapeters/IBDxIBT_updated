library(dplyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)

para_set = 4 #parameter set
run = 1 #model run
g.list = c(1,10,20,30,49,50,60,70,80,90,100,150,200,250,300,350,400,450,500)

#Frequency plots
for (g in g.list){
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

ggplot()+geom_line(data=map.df,aes(x=Generation,y=Frequency,col=Locus))+
  theme_bw()+ylim(0,1)+theme(legend.position="none")

summary.df = map.df %>% group_by(.,Generation) %>% summarise(var(Frequency))
ggplot()+geom_line(data=summary.df,aes(x=Generation,y=`var(Frequency)`))+
  theme_bw()+ylim(0,0.05)+ylab('Variance in allele frequency')

#Heatmaps
greencols = brewer.pal(9, 'Greens')
newcol = colorRampPalette(greencols)
ncols = 100
greencols2 = newcol(ncols)#apply the function to get 100 colours

for (g in g.list){
  df = read.csv(paste(getwd(),'/para_set_',para_set,'/model_run_',run,'/paraset_',para_set,'_offspring_map_',g,'.csv',sep=""))
  df$Generation = g
  
  if (g == 1){joint.dfs = df} else {joint.dfs = bind_rows(joint.dfs,df)}
  
}

ggplot()+geom_tile(data=joint.dfs,aes(x=X_pos,y=Y_pos,fill=factor(mapA)))+
  theme_classic()+theme(legend.position="none")+scale_fill_brewer(palette='Reds')+
  xlab('X position')+ylab('Y position')+ggtitle('Neutral locus 1')+
  facet_wrap(Generation~.)

ggplot()+geom_tile(data=joint.dfs,aes(x=X_pos,y=Y_pos,fill=factor(mapB)))+
  theme_classic()+theme(legend.position="none")+scale_fill_brewer(palette='Blues')+
  xlab('X position')+ylab('Y position')+ggtitle('Neutral locus 2')+
  facet_wrap(Generation~.)

ggplot()+geom_tile(data=joint.dfs,aes(x=X_pos,y=Y_pos,fill=factor(mapC)))+
  theme_classic()+theme(legend.position="none")+scale_fill_brewer(palette='Oranges')+
  xlab('X position')+ylab('Y position')+ggtitle('Neutral locus 3')+
  facet_wrap(Generation~.)

ggplot()+geom_tile(data=joint.dfs,aes(x=X_pos,y=Y_pos,fill=FLday))+
  theme_classic()+theme(legend.position="none")+scale_fill_gradient(low='lightgreen',high='darkgreen')+
  xlab('X position')+ylab('Y position')+ggtitle('Flowering day - by generation')+
  facet_wrap(Generation~.)

FLday.mapA = df %>% group_by(.,FLday) %>% summarise(.,mean(mapA))
FLday.mapB = df %>% group_by(.,FLday) %>% summarise(.,mean(mapB))
FLday.mapC = df %>% group_by(.,FLday) %>% summarise(.,mean(mapC))
ggplot()+geom_point(data=FLday.mapA,aes(x=FLday,y=`mean(mapA)`))+ylim(0,2)
ggplot()+geom_point(data=FLday.mapB,aes(x=FLday,y=`mean(mapB)`))+ylim(0,2)
ggplot()+geom_point(data=FLday.mapC,aes(x=FLday,y=`mean(mapC)`))+ylim(0,2)

hist(df$FLday)

#####Ordinal regression
require(foreign)
require(ggplot2)
require(MASS)
require(Hmisc)
require(reshape2)

df$mapA = factor(df$mapA, levels=c("0", "1", "2"), ordered=TRUE)
df$mapB = factor(df$mapB, levels=c("0", "1", "2"), ordered=TRUE)
df$mapC = factor(df$mapC, levels=c("0", "1", "2"), ordered=TRUE)

m = polr(mapA ~ FLday + X_pos + Y_pos + X_pos*Y_pos, data = df, Hess=TRUE)
summary(m)
ctable = coef(summary(m))
p = pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
ctable = cbind(ctable, "p value" = p)
ci = confint(m)
exp(cbind(OR = coef(m), ci))

mean(as.numeric(offspring_map$mapB))
