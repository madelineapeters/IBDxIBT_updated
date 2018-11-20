library(dplyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)
library(psych)
library(tidyr)
library(plyr)

para_set = 9 #parameter set
r.list = 1:4 #model run
g.list = c(1,10,20,30,40,50,60,70,80,90,100,150,200,250,300,350,400,450,500)

Ne.df = as.data.frame(matrix(nrow=length(g.list),ncol=3))
names(Ne.df) = c('mean','var','Ne')

for (r in r.list){
  for (p in 9:16){
    for (g in 1:length(g.list)){
      
      gen = g.list[g]
      
      df = read.csv(paste(getwd(),'/para_set_',p,'/model_run_',r,'/paraset_',p,'_offspring_map_',gen,'.csv',sep=""))
      
      #Mothers
      Mother.df = bind_cols(as.data.frame(1:nrow(df)),as.data.frame(rep(0,nrow(df))))
      names(Mother.df) = c('Mother','Offspring')
      
      Mother.count = count(df,vars='Mother')
      names(Mother.count) = c('Mother','Offspring')
      
      anti.Mother = anti_join(Mother.df,Mother.count,by='Mother')
      
      Mother.count = bind_rows(Mother.count,anti.Mother)
      names(Mother.count) = c('Parent','Offspring1')
      
      #Fathers
      Father.df = bind_cols(as.data.frame(1:nrow(df)),as.data.frame(rep(0,nrow(df))))
      names(Father.df) = c('Father','Offspring')
      
      Father.count = count(df,vars='Father')
      names(Father.count) = c('Father','Offspring')
      Father.count$Father = as.numeric(Father.count$Father)
      
      anti.Father = anti_join(Father.df,Father.count,by='Father')
      
      Father.count = bind_rows(Father.count,anti.Father)
      names(Father.count) = c('Parent','Offspring2')
      
      #Parents
      Parent.count = left_join(Mother.count,Father.count,by='Parent') %>% mutate(.,Offspring = Offspring1+Offspring2)
      
      Parent.count = Parent.count %>% select(.,Parent,Offspring)
      Parent.count$Gen = gen
      Parent.count$paraset = p
      Parent.count$Run = r
      
      mean = mean(Parent.count$Offspring)
      var = var(Parent.count$Offspring)
      
      Ne.df$mean[g] = mean
      Ne.df$var[g] = var
      Ne.df$Ne[g] = 4*nrow(df)/(mean+var)
      Ne.df$Gen[g] = gen
      
      if ((p == 9)&(g==1)&(r==1)){Parent.joint = Parent.count} else {Parent.joint = bind_rows(Parent.joint,Parent.count)}
      
    }
    
    Ne.df$paraset = p
    
    if (p == 9){Ne.full = Ne.df} else {Ne.full = bind_rows(Ne.full,Ne.df)}
  }
  Ne.full$grouping[1:(nrow(Ne.full)/2)] = 'Selfing'
  Ne.full$grouping[(1+nrow(Ne.full)/2):nrow(Ne.full)] = 'No selfing'
  Ne.full$paraset[(Ne.full$paraset == 9)|(Ne.full$paraset == 13)] = 'Random'
  Ne.full$paraset[(Ne.full$paraset == 10)|(Ne.full$paraset == 14)] = 'IBT'
  Ne.full$paraset[(Ne.full$paraset == 11)|(Ne.full$paraset == 15)] = 'IBD'
  Ne.full$paraset[(Ne.full$paraset == 12)|(Ne.full$paraset == 16)] = 'IBDxIBT'
  names(Ne.full)[5] = 'Isolation'
  Ne.full$run = r
  if (r == 1){Ne.joint = Ne.full} else {Ne.joint = bind_rows(Ne.joint,Ne.full)}
}

detach(package:plyr)
detach(package:dplyr)
library(dplyr)

Parent.joint$grouping[1:(nrow(Parent.joint)/2)] = 'Selfing'
Parent.joint$grouping[(1+nrow(Parent.joint)/2):nrow(Parent.joint)] = 'No selfing'
Parent.joint$paraset[(Parent.joint$paraset == 9)|(Parent.joint$paraset == 13)] = 'Random'
Parent.joint$paraset[(Parent.joint$paraset == 10)|(Parent.joint$paraset == 14)] = 'IBT'
Parent.joint$paraset[(Parent.joint$paraset == 11)|(Parent.joint$paraset == 15)] = 'IBD'
Parent.joint$paraset[(Parent.joint$paraset == 12)|(Parent.joint$paraset == 16)] = 'IBDxIBT'
names(Parent.joint)[4] = 'Isolation'

Ne.avg = Ne.joint %>% group_by(grouping,Gen,Isolation) %>% summarize(.,avg = mean(Ne))
Var.avg = Ne.joint %>% group_by(grouping,Gen,Isolation) %>% summarize(.,avg = mean(var))
ggplot()+geom_line(data=filter(Ne.joint,run==1),aes(x=Gen,y=Ne,col=Isolation),alpha=0.25)+
  geom_line(data=filter(Ne.joint,run==2),aes(x=Gen,y=Ne,col=Isolation),alpha=0.25)+
  geom_line(data=filter(Ne.joint,run==3),aes(x=Gen,y=Ne,col=Isolation),alpha=0.25)+
  geom_line(data=filter(Ne.joint,run==4),aes(x=Gen,y=Ne,col=Isolation),alpha=0.25)+
  geom_line(data=Ne.avg,aes(x=Gen,y=avg,col=Isolation),size=1)+
  theme_classic()+ylab("Effective population size")+xlab("Generation")+facet_grid(.~grouping)+geom_hline(yintercept=2500,linetype='dashed')

ggplot()+geom_line(data=filter(Ne.joint,run==1),aes(x=Gen,y=var,col=Isolation),alpha=0.25)+
  geom_line(data=filter(Ne.joint,run==2),aes(x=Gen,y=var,col=Isolation),alpha=0.25)+
  geom_line(data=filter(Ne.joint,run==3),aes(x=Gen,y=var,col=Isolation),alpha=0.25)+
  geom_line(data=filter(Ne.joint,run==4),aes(x=Gen,y=var,col=Isolation),alpha=0.25)+
  geom_line(data=Var.avg,aes(x=Gen,y=avg,col=Isolation),size=1)+
  theme_classic()+ylab("Variation in offspring per parent")+xlab("Generation")+facet_grid(.~grouping)+geom_hline(yintercept=2,linetype='dashed')

Parent.filtered = filter(Parent.joint,Gen %in% c(1,50,150,500))#,grouping=="Selfing")
ggplot(data=filter(Parent.joint,Gen %in% c(1,50,150,500),grouping=="Selfing"),aes(Offspring,fill=Isolation))+geom_histogram(binwidth=0.5)+facet_grid(Gen~Isolation,labeller=labeller(Gen=c(`1`='Gen 1',`50`='Gen 50',`150`='Gen 150',`500`='Gen 500')))+theme_classic()+xlab("Number gametes to next generation per parent")+ylab('Count')

                                                                                                                                                     