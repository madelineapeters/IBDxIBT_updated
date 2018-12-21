library(dplyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)
library(psych)
library(tidyr)
library(plyr)

para_set = 9 #parameter set
r.list = c(1:3) #model run
g.list = c(1,10,20,30,40,50,60,70,80,90,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800)


for (r in r.list){
  for (p in 17:24){
    
    Ne.df = as.data.frame(matrix(nrow=length(g.list),ncol=3))
    names(Ne.df) = c('mean','var','Ne')
    
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
      
      if ((p == 17)&(g==1)&(r==1)){Parent.joint = Parent.count} else {Parent.joint = bind_rows(Parent.joint,Parent.count)}
      
    }
    
    Ne.df$paraset = p
    
    if (p == 17){Ne.full = Ne.df} else {Ne.full = bind_rows(Ne.full,Ne.df)}
  }
  Ne.full$grouping[1:(nrow(Ne.full)/2)] = 'Selfing'
  Ne.full$grouping[(1+nrow(Ne.full)/2):nrow(Ne.full)] = 'No selfing'
  Ne.full$paraset[(Ne.full$paraset == 17)|(Ne.full$paraset == 21)] = 'Random'
  Ne.full$paraset[(Ne.full$paraset == 18)|(Ne.full$paraset == 22)] = 'IBT'
  Ne.full$paraset[(Ne.full$paraset == 19)|(Ne.full$paraset == 23)] = 'IBD'
  Ne.full$paraset[(Ne.full$paraset == 20)|(Ne.full$paraset == 24)] = 'IBDxIBT'
  names(Ne.full)[5] = 'Isolation'
  Ne.full$run = r
  if (r == 1){Ne.joint = Ne.full} else {Ne.joint = bind_rows(Ne.joint,Ne.full)}
}

detach(package:plyr)
detach(package:dplyr)
library(dplyr)

Parent.joint$grouping[1:(nrow(Parent.joint)/2)] = 'Selfing'
Parent.joint$grouping[(1+nrow(Parent.joint)/2):nrow(Parent.joint)] = 'No selfing'
Parent.joint$paraset[(Parent.joint$paraset == 17)|(Parent.joint$paraset == 21)] = 'Random'
Parent.joint$paraset[(Parent.joint$paraset == 18)|(Parent.joint$paraset == 22)] = 'IBT'
Parent.joint$paraset[(Parent.joint$paraset == 19)|(Parent.joint$paraset == 23)] = 'IBD'
Parent.joint$paraset[(Parent.joint$paraset == 20)|(Parent.joint$paraset == 24)] = 'IBDxIBT'
names(Parent.joint)[4] = 'Isolation'

Ne.avg = Ne.joint %>% group_by(grouping,Gen,Isolation) %>% summarize(.,avg = mean(Ne)) 
Ne.sd = Ne.joint %>% group_by(grouping,Gen,Isolation) %>% summarize(.,SD = sd(Ne))
Ne.plot = bind_cols(Ne.avg,as.data.frame(Ne.sd$SD))
names(Ne.plot)[5] = 'SD'

Var.avg = Ne.joint %>% group_by(grouping,Gen,Isolation) %>% summarize(.,avg = mean(var))
Var.sd = Ne.joint %>% group_by(grouping,Gen,Isolation) %>% summarize(.,SD = sd(var))
Var.plot = bind_cols(Var.avg,as.data.frame(Var.sd$SD))
names(Var.plot)[5] = 'SD'

ggplot()+geom_line(data=filter(Ne.joint,run==1),aes(x=Gen,y=Ne,col=Isolation),alpha=0.25)+
  geom_line(data=filter(Ne.joint,run==2),aes(x=Gen,y=Ne,col=Isolation),alpha=0.25)+
  geom_line(data=filter(Ne.joint,run==3),aes(x=Gen,y=Ne,col=Isolation),alpha=0.25)+
  geom_line(data=filter(Ne.joint,run==4),aes(x=Gen,y=Ne,col=Isolation),alpha=0.25)+
  geom_line(data=filter(Ne.joint,run==5),aes(x=Gen,y=Ne,col=Isolation),alpha=0.25)+
  geom_line(data=filter(Ne.joint,run==6),aes(x=Gen,y=Ne,col=Isolation),alpha=0.25)+
  geom_line(data=filter(Ne.joint,run==8),aes(x=Gen,y=Ne,col=Isolation),alpha=0.25)+
  geom_line(data=Ne.avg,aes(x=Gen,y=avg,col=Isolation),size=1)+
  theme_classic()+ylab("Effective population size")+xlab("Generation")+facet_grid(.~grouping)+geom_hline(yintercept=2500,linetype='dashed')

ggplot(data=Ne.plot,aes(x=Gen))+
  geom_line(aes(y=avg,col=Isolation),size=1)+
  geom_ribbon(aes(ymin=avg-SD,ymax=avg+SD,fill=Isolation),alpha=0.2)+
  theme_classic()+ylab("Effective population size")+xlab("Generation")+facet_grid(.~grouping)+geom_hline(yintercept=2500,linetype='dashed')

ggplot()+geom_line(data=filter(Ne.joint,run==1),aes(x=Gen,y=var,col=Isolation),alpha=0.25)+
  geom_line(data=filter(Ne.joint,run==2),aes(x=Gen,y=var,col=Isolation),alpha=0.25)+
  geom_line(data=filter(Ne.joint,run==3),aes(x=Gen,y=var,col=Isolation),alpha=0.25)+
  geom_line(data=Var.avg,aes(x=Gen,y=avg,col=Isolation),size=1)+
  theme_classic()+ylab("Variation in offspring per parent")+xlab("Generation")+facet_grid(.~grouping)+geom_hline(yintercept=2,linetype='dashed')

ggplot(data=Var.plot,aes(x=Gen))+
  geom_line(aes(y=avg,col=Isolation),size=1)+
  geom_ribbon(aes(ymin=avg-SD,ymax=avg+SD,fill=Isolation),alpha=0.2)+
  theme_classic()+ylab("Variation in offspringer per parent")+xlab("Generation")+facet_grid(.~grouping)+geom_hline(yintercept=2,linetype='dashed')

Parent.filtered = filter(Parent.joint,Gen %in% c(1,50,150,500),grouping=="No selfing")
ggplot(data=filter(Parent.joint,Gen %in% c(1,50,150,500),grouping=="Selfing"),aes(Offspring,fill=Isolation))+geom_histogram(binwidth=0.5)+facet_grid(Gen~Isolation,labeller=labeller(Gen=c(`1`='Gen 1',`50`='Gen 50',`150`='Gen 150',`500`='Gen 500')))+theme_classic()+xlab("Number gametes to next generation per parent")+ylab('Count')

                                                                                                                                                     