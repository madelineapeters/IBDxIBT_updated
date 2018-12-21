library(dplyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)
library(psych)

para_set = 9 #parameter set
r.list = c(1:3) #model run
g.list = c(1,10,20,30,40,50,60,70,80,90,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800)


  for (p in 17:24){
    for (r in r.list){
      df = read.csv(paste(getwd(),'/para_set_',p,'/model_run_',r,'/data_by_generation.csv',sep=""))
      
      if ((r == 1)&(p==17)){
        joint.IB = select(df,gen,IBmean)
        joint.IB$run = r
        joint.IB$paraset = p
      } else{
        temp.IB = select(df,gen,IBmean)
        temp.IB$run = r
        temp.IB$paraset = p
        
        joint.IB = bind_rows(joint.IB,temp.IB)
      }
      
  }
  
}


  joint.IB$grouping[1:(nrow(joint.IB)/2)] = 'Selfing'
  joint.IB$grouping[(1+nrow(joint.IB)/2):nrow(joint.IB)] = 'No selfing'
  joint.IB$paraset[(joint.IB$paraset == 17)|(joint.IB$paraset == 21)] = 'Random'
  joint.IB$paraset[(joint.IB$paraset == 18)|(joint.IB$paraset == 22)] = 'IBT'
  joint.IB$paraset[(joint.IB$paraset == 19)|(joint.IB$paraset == 23)] = 'IBD'
  joint.IB$paraset[(joint.IB$paraset == 20)|(joint.IB$paraset == 24)] = 'IBDxIBT'
  names(joint.IB)[4] = 'Isolation'   
  names(joint.IB)[5] = 'Selfing status'
IB.avg = joint.IB %>% group_by(`Selfing status`,Isolation,gen) %>% summarise(.,avg=mean(IBmean))  
IB.var = joint.IB %>% group_by(`Selfing status`,Isolation,gen) %>% summarise(.,var=var(IBmean)) 
ggplot()+geom_line(data=filter(IB.avg,gen %in% g.list),aes(x=gen,y=avg,linetype=`Selfing status`,col=Isolation),size=1)+theme_classic()+xlab("Generation")+ylab("Average inbreeding coefficient across runs")
+geom_line(data=filter(IB.avg,grouping=="No selfing",gen %in% g.list),aes(x=gen,y=avg,col=Isolation))

ggplot()+geom_line(data=filter(joint.IB,run==1),aes(x=gen,y=IBmean,col=Isolation))+
  geom_line(data=filter(joint.IB,run==2),aes(x=gen,y=IBmean,col=Isolation))+
  theme_classic()+facet_grid(.~grouping)

