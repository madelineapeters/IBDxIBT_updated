detach("package:plyr")
detach("package:genetics")
detach("package:tibble")
detach("package:dplyr")
library(dplyr)
library(tibble)
library(genetics)
library(ggplot2)
library(RColorBrewer)

r.list = c(1:3) #model run
g.list = c(1,10,20,30,40,50,60,70,80,90,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800)


##############################################################
## Calculate spatial autocorrelation statistics
##############################################################
for (r in r.list){
  for (s in 17:24){
    
    for (g in 1:length(g.list)){
      
      gen = g.list[g]
      
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
        mutate(.,map24 = neut24a+neut24b) %>% 
        mutate(.,haplotype = paste(neut1a,neut1b,
                                   neut2a,neut2b,
                                   neut3a,neut3b,
                                   neut4a,neut4b,
                                   neut5a,neut5b,sep=""))
                                   neut6a,neut6b,
                                   neut7a,neut7b,
                                   neut8a,neut8b,
                                   neut9a,neut9b,
                                   neut10a,neut10b,
                                   neut11a,neut11b,
                                   neut12a,neut12b,
                                   neut13a,neut13b,
                                   neut14a,neut14b,
                                   neut15a,neut15b,
                                   neut16a,neut16b,
                                   neut17a,neut17b,
                                   neut18a,neut18b,
                                   neut19a,neut19b,
                                   neut20a,neut20b,
                                   neut21a,neut21b,
                                   neut22a,neut22b,
                                   neut23a,neut23b,
                                   neut24a,neut24b,
                                   sep=""))
      
      neut.loci = neutral.df %>% select(.,neut1a:neut24b)
      neut.haplo = neutral.df %>% select(.,haplotype)

      ind.neutral.df = neutral.df %>% select(.,FLday,X_pos,Y_pos,mapA:mapC,map1:map24)
      
      neutral.df = neutral.df %>% select(.,mapA:mapC,map1:map24)
      neutral.freq = colSums(neutral.df)/(2*nrow(neutral.df))
        
      obs.freq = matrix(nrow=24,ncol=24)
      exp.freq = neutral.freq%*%t(neutral.freq)    

            