library(dplyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)

g.list = c(1,10,20,30,50,60,70,80,90,100,150,200,250,300,350,400,450,500)

FLvar.df = as.data.frame(matrix(nrow=length(g.list),ncol=8))
names(FLvar.df) = sapply(1:8, function(X) paste('paraset',X,sep="_"))

##############################################################
## Calculate spatial autocorrelation statistics
##############################################################
for (s in 1:8){
    
    df = read.csv(paste(getwd(),'/para_set_',s,'/model_run_',run,'/data_by_generation.csv',sep=""))
    
    if (s == 1){
      FLvar.df = df %>% select(.,gen,FTvar) %>% mutate(.,paraset = s)
    } else {
      FLvar.temp = df %>% select(.,gen,FTvar) %>% mutate(.,paraset = s)
      FLvar.df = FLvar.df %>% 
        bind_rows(.,FLvar.temp)
    }  
}

FLvar.df$grouping[1:2000] = 'Selfing'
FLvar.df$grouping[2001:4000] = 'No selfing'
FLvar.df$paraset[(FLvar.df$paraset == 1)|(FLvar.df$paraset == 5)] = 'Null'
FLvar.df$paraset[(FLvar.df$paraset == 2)|(FLvar.df$paraset == 6)] = 'IBT'
FLvar.df$paraset[(FLvar.df$paraset == 3)|(FLvar.df$paraset == 7)] = 'IBD'
FLvar.df$paraset[(FLvar.df$paraset == 4)|(FLvar.df$paraset == 8)] = 'IBDxIBT'
names(FLvar.df)[3] = 'Isolation'
ggplot()+geom_line(data=FLvar.df,aes(x=gen,y=FTvar,col=Isolation))+
  facet_wrap(grouping~.)+
  xlab('Generation')+ylab('Variance')+ggtitle('Variance in flowering time')+
  theme_classic()
