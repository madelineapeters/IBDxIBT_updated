library(dplyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)
library(tidyr)

r.list = c(1,2)
g.list = c(1,10,20,30,50,60,70,80,90,100,150,200,250,300,350,400,450,500)

selfing.df = as.data.frame(matrix(nrow=length(g.list),ncol=8))
names(selfing.df) = sapply(9:16, function(X) paste('paraset',X,sep="_"))

##############################################################
## Calculate spatial autocorrelation statistics
##############################################################
for (r in r.list){
  for (s in 9:16){
    for (g in 1:length(g.list)){
      
      gen = g.list[g]
      
      df = read.csv(paste(getwd(),'/para_set_',s,'/model_run_',r,'/paraset_',s,'_offspring_map_',gen,'.csv',sep=""))
      
      selfing.df[g,s-8] = nrow(filter(df,Mother==Father))/nrow(df)
      
    }
  }  
  
  selfing.plot = gather(selfing.df,'paraset','rate')
  selfing.plot$generation = rep(g.list,nrow(selfing.plot)/length(g.list))
  
  selfing.plot$grouping[1:(nrow(selfing.plot)/2)] = 'Selfing'
  selfing.plot$grouping[(1+nrow(selfing.plot)/2):nrow(selfing.plot)] = 'No selfing'
  selfing.plot[(selfing.plot == 'paraset_9')|(selfing.plot == 'paraset_13')] = 'Null'
  selfing.plot[(selfing.plot == 'paraset_10')|(selfing.plot == 'paraset_4')] = 'IBT'
  selfing.plot[(selfing.plot == 'paraset_11')|(selfing.plot == 'paraset_15')] = 'IBD'
  selfing.plot[(selfing.plot == 'paraset_12')|(selfing.plot == 'paraset_16')] = 'IBDxIBT'
  names(selfing.plot)[1] = 'Isolation'
  selfing.plot$Run = r
  
  if (r == 1){self.plot.full = selfing.plot} else {self.plot.full = bind_rows(self.plot.full,selfing.plot)}
}


Rate.avg = self.plot.full %>% group_by(grouping,generation,Isolation) %>% summarize(.,avg = mean(rate))
ggplot()+geom_line(data=filter(self.plot.full,grouping=="Selfing",Run==1),aes(x=generation,y=100*rate,col=Isolation),alpha=0.25)+
  geom_line(data=filter(self.plot.full,grouping=="Selfing",Run==2),aes(x=generation,y=100*rate,col=Isolation),alpha=0.25)+
  geom_line(data=filter(Rate.avg,grouping=="Selfing"),aes(x=generation,y=100*avg,col=Isolation),size=1)+
  xlab('Generation')+ylab('Selfing rate (% self-matings)')+
  theme_classic() 
