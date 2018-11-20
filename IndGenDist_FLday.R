library("tidyr")
library("dplyr")
library("poppr")
library("pegas")
library("ape")
library("adegenet")
library("ade4")

run = 1 #model run
g.list = c(1,10,20,30,50,60,70,80,90,100,150,200,250,300,350,400,450,500)

Mantel.obs = as.data.frame(matrix(nrow=length(g.list),ncol=8))
names(Mantel.obs) = sapply(1:8, function(X) paste('paraset',X,sep="_"))

Mantel.p = as.data.frame(matrix(nrow=length(g.list),ncol=8))
names(Mantel.p) = sapply(1:8, function(X) paste('paraset',X,sep="_"))

##############################################################
## Calculate spatial autocorrelation statistics
##############################################################
<<<<<<< HEAD
for (s in 3:8){
=======
for (s in 1:8){
>>>>>>> origin/master
  for (g in 1:length(g.list)){
    
    gen = g.list[g]
    print(paste('Working on generation',gen,'of parameter set',s,sep=" "))
    
    ##############################################################
    ## Convert the data
    ##############################################################
    Mydata = read.csv(paste(getwd(),'/para_set_',s,'/model_run_',run,'/paraset_',s,'_offspring_map_',gen,'.csv',sep=""))
    Mydata1 = Mydata %>% select(.,loc1a:loc5b)
    Mydata1[Mydata1 == 1] = 'D'
    Mydata1[Mydata1 == -1] = 'd'
<<<<<<< HEAD
    locus = Mydata1
=======
    locus = Mydata1 %>% select(.,)
>>>>>>> origin/master
    locus$genotype1 = paste(Mydata1$loc1a,Mydata1$loc1b,sep="")
    locus$genotype2 = paste(Mydata1$loc2a,Mydata1$loc2b,sep="")
    locus$genotype3 = paste(Mydata1$loc3a,Mydata1$loc3b,sep="")
    locus$genotype4 = paste(Mydata1$loc4a,Mydata1$loc4b,sep="")
    locus$genotype5 = paste(Mydata1$loc5a,Mydata1$loc5b,sep="")
    
    #Convert matrix to a  genind object (for the package adegenet). #The genind object can then easily be converted into loci objects (package pegas) (i.e. Mydata2)
    Mydata1 = df2genind(locus, ploidy = 2, sep="")
<<<<<<< HEAD
    #Mydata2 = genind2loci(Mydata1)
=======
    Mydata2 = genind2loci(Mydata1)
>>>>>>> origin/master
    
    ##############################################################
    ## Individual genetic distance: euclidean distance
    ##############################################################
    #Use the function dist() from adegenet 
    #Use euclidean distance among vector of allele frequencies
    #distgenEUCL = dist(Mydata1, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
    #hist(distgenEUCL)
    
    ##############################################################
    ## Individual genetic distance: num. loci at which ind. differ
    ##############################################################
    #distgenDIFF = dist.gene(Mydata2, method="pairwise", pairwise.deletion = TRUE, variance = FALSE)
    #hist(distgenDIFF)
    
    ##############################################################
    ## Number allelic differences between two individuals
    ##############################################################
    distgenDISS = diss.dist(Mydata1, percent = FALSE, mat = FALSE)
    #hist(distgenDISS)
    
    ##############################################################
    ## Comparison of different distance measures
    ##############################################################
    #boxplot(distgenEUCL, distgenDIFF, distgenDISS)
    
    ##############################################################
    ## Matrix of spatial distances
    ##############################################################
    dist.mat = matrix(nrow=sqrt(nrow(Mydata)), ncol=sqrt(nrow(Mydata)))
    
    X_pos = Mydata$X_pos
    Y_pos = Mydata$Y_pos
    
    for (m in 1:sqrt(nrow(Mydata))) {
      for (p in 1:sqrt(nrow(Mydata))) {
        dist.mat[m,p]<-sqrt((X_pos[m]-X_pos[p])^2 + (Y_pos[m] - Y_pos[p])^2)	
      } #distance is just the Euclidean distance between two points
    }
    dist.mat = dist(cbind(Mydata$X_pos,Mydata$Y_pos))
<<<<<<< HEAD
    mantel.out = mantel.rtest(distgenDISS, dist.mat, nrepet = 200)
=======
    mantel.out = mantel.rtest(distgenDISS, dist.mat, nrepet = 500)
>>>>>>> origin/master
    
    Mantel.obs[g,s] = mantel.out$obs
    Mantel.p[g,s] = mantel.out$pvalue
    
    write.csv(Mantel.obs,'FLday.Mantel.obs.csv',row.names=FALSE)
<<<<<<< HEAD
    write.csv(Mantel.p,'FLday.Mantel.p.csv',row.names=FALSE)
=======
    write.csv(Mantel.p,'FLday.Mantel.p.csv.',row.names=FALSE)
>>>>>>> origin/master
    
  }  
}

##############################################################
## Plot spatial autocorrelation statistics
##############################################################
comment.out = function(){
  library(ggplot2)

<<<<<<< HEAD
  Mantel.obs = read.csv(paste(getwd(),'FLday.Mantel.obs.csv',sep="/"))
  Mantel.obs.plot = gather(Mantel.obs,'paraset','correlation')
  Mantel.obs.plot$generation = g.list
  Mantel.p = read.csv(paste(getwd(),'FLday.Mantel.p.csv',sep="/"))
=======
  Mantel.obs.plot = gather(Mantel.obs,'paraset','correlation')
  Mantel.obs.plot$generation = g.list
  
>>>>>>> origin/master
  Mantel.p.plot = gather(Mantel.p,'paraset','significance')
  Mantel.p.plot$generation = g.list
  
  Mantel.plot = bind_cols(Mantel.obs.plot,Mantel.p.plot) %>% select(.,paraset,generation,correlation,significance)
  Mantel.plot$grouping[1:72] = 'Selfing'
  Mantel.plot$grouping[73:144] = 'No selfing'
  Mantel.plot[(Mantel.plot == 'paraset_1')|(Mantel.plot == 'paraset_5')] = 'Null'
  Mantel.plot[(Mantel.plot == 'paraset_2')|(Mantel.plot == 'paraset_6')] = 'IBT'
  Mantel.plot[(Mantel.plot == 'paraset_3')|(Mantel.plot == 'paraset_7')] = 'IBD'
  Mantel.plot[(Mantel.plot == 'paraset_4')|(Mantel.plot == 'paraset_8')] = 'IBDxIBT'
  names(Mantel.plot)[1] = 'Isolation'
  
  ggplot()+geom_line(data=Mantel.plot,aes(x=generation,y=correlation,col=Isolation))+
    geom_point(data=filter(Mantel.plot,significance<0.05),aes(x=generation,y=correlation,col=Isolation),shape=8)+
    facet_wrap(grouping~.)+
<<<<<<< HEAD
    xlab('Generation')+ylab('Correlation')+ggtitle('Correlation between spatial distance and mean flowering day')+
=======
    xlab('Generation')+ylab('Correlation')+ggtitle('Correlation between spatial distance and genetic distance (number allelic differences between individuals)')+
>>>>>>> origin/master
    theme_classic()
}
