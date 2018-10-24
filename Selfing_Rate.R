library(dplyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)

g.list = c(1,10,20,30,50,60,70,80,90,100,150,200,250,300,350,400,450,500)

selfing.df = as.data.frame(matrix(nrow=length(g.list),ncol=8))
names(selfing.df) = sapply(1:8, function(X) paste('paraset',X,sep="_"))

##############################################################
## Calculate spatial autocorrelation statistics
##############################################################
for (s in 1:8){
  for (g in 1:length(g.list)){
    
    gen = g.list[g]

    df = read.csv(paste(getwd(),'/para_set_',s,'/model_run_',run,'/paraset_',s,'_offspring_map_',gen,'.csv',sep=""))
  
  selfing.df[g,s] = nrow(filter(df,Mother==Father))/nrow(df)
 
  }
}
 
selfing.plot = gather(selfing.df,'paraset','rate')
selfing.plot$generation = g.list

selfing.plot$grouping[1:72] = 'Selfing'
selfing.plot$grouping[73:144] = 'No selfing'
selfing.plot[(selfing.plot == 'paraset_1')|(selfing.plot == 'paraset_5')] = 'Null'
selfing.plot[(selfing.plot == 'paraset_2')|(selfing.plot == 'paraset_6')] = 'IBT'
selfing.plot[(selfing.plot == 'paraset_3')|(selfing.plot == 'paraset_7')] = 'IBD'
selfing.plot[(selfing.plot == 'paraset_4')|(selfing.plot == 'paraset_8')] = 'IBDxIBT'
names(selfing.plot)[1] = 'Isolation'

ggplot()+geom_line(data=filter(selfing.plot,grouping=="Selfing"),aes(x=generation,y=100*rate,col=Isolation))+
  xlab('Generation')+ylab('Selfing rate (% self-matings)')+
  theme_classic() 
