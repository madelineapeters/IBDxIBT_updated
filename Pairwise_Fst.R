library(dplyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)

r.list = c(1:2) #model run
g.list = c(1,10,20,30,40,50,60,70,80,90,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800)

##############################################################
## Fst function
##############################################################
fst.fun = function(i){
  
  p.1 = ind.1[1,i]; q.1 = 1 - p.1
  p.2 = ind.2[1,i]; q.2 = 1 - p.2
  
  H.exp1 = 1 - sum(c(p.1^2,q.1^2))
  H.exp2 = 1 - sum(c(p.2^2,q.2^2))
  HS = sum(c(H.exp1,H.exp2))/2
  
  p.T = (p.1+p.2)/2
  q.T = 1 - p.T
  
  HT = 1 - sum(c(p.T^2,q.T^2))

  fst = (HT - HS)/HT
  
  if(is.na(fst)){fst = 0}
  
  return(fst)
}

##############################################################
## Calculate pairwise fst and regress over distance
##############################################################
Set.Pairwise.Fst = as.data.frame(matrix(nrow=300,ncol=4))
names(Set.Pairwise.Fst) = c('Distance','Fst','Run','Set')

for (s in 17:24){
  for (r in r.list){

      gen = 800
      
      df = read.csv(paste(getwd(),'/para_set_',s,'/model_run_',r,'/paraset_',s,'_offspring_map_',gen,'.csv',sep=""))
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

      for (i in 1:300){
        
        ind.1 = ind.neutral.df[sample(1:nrow(ind.neutral.df),1),]
        ind.2 = ind.neutral.df[sample(1:nrow(ind.neutral.df),1),]
        
        ind.1[1,4:ncol(ind.1)] = ind.1[1,4:ncol(ind.1)]/2
        ind.2[1,4:ncol(ind.1)] = ind.2[1,4:ncol(ind.1)]/2
        
        dist = sqrt(((ind.1$X_pos-ind.2$X_pos)^2)+((ind.1$Y_pos-ind.2$Y_pos)^2))
        
        Fst = mean(sapply(4:ncol(ind.1),fst.fun))
        
        comment.out = function(){
          #Find the allele frequencies for each subpopulation.
          p.1 = sum(ind.1[1,4:ncol(ind.1)])/(ncol(ind.1)-3)
          p.2 = sum(ind.2[1,4:ncol(ind.2)])/(ncol(ind.2)-3)
          q.1 = 1-p.1
          q.2 = 1-p.2
          
          #Find the average allele frequencies for the total population.
          p.T  = sum(p.1,p.2)/2
          q.T = 1-p.T
          
          #Calculate the heterozygosity (2pq) for each subpopulation.
          H.1 = 2*p.1*q.1
          H.2 = 2*p.2*q.2
          
          #Calculate the average of these subpopulation heterozygosities. This is HS.
          HS = sum(H.1, H.2)/2
          
          #Calculate the heterozygosity based on the total population allele frequencies. This is HT.
          HT = 2*p.T*q.T
          
          #Finally, calculate FST=(HT-HS)/HT.
          Fst = (HT-HS)/HT
        }
        
        Set.Pairwise.Fst$Distance[i] = dist
        Set.Pairwise.Fst$Fst[i] = Fst
        
      }
      
      Set.Pairwise.Fst$Run = r
      Set.Pairwise.Fst$Set = s
      
      if ((r == 1)&(s == 17)){Pairwise.Fst = Set.Pairwise.Fst} else {Pairwise.Fst = bind_rows(Pairwise.Fst,Set.Pairwise.Fst)}
    
  }
}


Pairwise.Fst$grouping[1:(nrow(Pairwise.Fst)/2)] = 'Selfing'
Pairwise.Fst$grouping[(1+nrow(Pairwise.Fst)/2):nrow(Pairwise.Fst)] = 'No selfing'
Pairwise.Fst$Set[(Pairwise.Fst$Set == 17)|(Pairwise.Fst$Set == 21)] = 'Random'
Pairwise.Fst$Set[(Pairwise.Fst$Set == 18)|(Pairwise.Fst$Set == 22)] = 'IBT'
Pairwise.Fst$Set[(Pairwise.Fst$Set == 19)|(Pairwise.Fst$Set == 23)] = 'IBD'
Pairwise.Fst$Set[(Pairwise.Fst$Set == 20)|(Pairwise.Fst$Set == 24)] = 'IBDxIBT'
names(Pairwise.Fst)[4] = 'Isolation'

ggplot()+geom_line(data=Pairwise.Fst,aes(x=Distance,y=Fst,col=Isolation))+facet_grid(grouping~.)
ggplot()+geom_smooth(data=filter(Pairwise.Fst,Run==1),aes(x=Distance,y=Fst,color=Isolation))+facet_grid(grouping~.)
ggplot()+geom_smooth(data=filter(Pairwise.Fst,Run==2),aes(x=Distance,y=Fst,color=Isolation))+facet_grid(grouping~.)
