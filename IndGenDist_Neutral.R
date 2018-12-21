library("tidyr")
library("dplyr")
library("poppr")
library("pegas")
library("ape")
library("adegenet")
library("ade4")

run = 1 #model run
g.short = c(1,seq(100,800,100))

Mantel.obs = as.data.frame(matrix(nrow=length(g.short),ncol=8))
names(Mantel.obs) = sapply(17:24, function(X) paste('paraset',X,sep="_"))

Mantel.p = as.data.frame(matrix(nrow=length(g.short),ncol=8))
names(Mantel.p) = sapply(17:24, function(X) paste('paraset',X,sep="_"))

##############################################################
## Calculate spatial autocorrelation statistics
##############################################################
for (s in 17:24){
  for (g in 1:length(g.short)){
    
    gen = g.short[g]
    print(paste('Working on generation',gen,'of parameter set',s,sep=" "))
    
    ##############################################################
    ## Convert the data
    ##############################################################
    Mydata = read.csv(paste(getwd(),'/para_set_',s,'/model_run_',run,'/paraset_',s,'_offspring_map_',gen,'.csv',sep=""))
    locus = Mydata %>% select(.,genotypeA,genotypeB,genotypeC)
    locus$genotype1 = paste(Mydata$neut1a,Mydata$neut1b,sep="")
    locus$genotype2 = paste(Mydata$neut2a,Mydata$neut2b,sep="")
    locus$genotype3 = paste(Mydata$neut3a,Mydata$neut3b,sep="")
    locus$genotype4 = paste(Mydata$neut4a,Mydata$neut4b,sep="")
    locus$genotype5 = paste(Mydata$neut5a,Mydata$neut5b,sep="")
    locus$genotype6 = paste(Mydata$neut6a,Mydata$neut6b,sep="")
    locus$genotype7 = paste(Mydata$neut7a,Mydata$neut7b,sep="")
    locus$genotype8 = paste(Mydata$neut8a,Mydata$neut8b,sep="")
    locus$genotype9 = paste(Mydata$neut9a,Mydata$neut9b,sep="")
    locus$genotype10 = paste(Mydata$neut10a,Mydata$neut10b,sep="")
    locus$genotype11 = paste(Mydata$neut11a,Mydata$neut11b,sep="")
    locus$genotype12 = paste(Mydata$neut12a,Mydata$neut12b,sep="")
    locus$genotype13 = paste(Mydata$neut13a,Mydata$neut13b,sep="")
    locus$genotype14 = paste(Mydata$neut14a,Mydata$neut14b,sep="")
    locus$genotype15 = paste(Mydata$neut15a,Mydata$neut15b,sep="")
    locus$genotype16 = paste(Mydata$neut16a,Mydata$neut16b,sep="")
    locus$genotype17 = paste(Mydata$neut17a,Mydata$neut17b,sep="")
    locus$genotype18 = paste(Mydata$neut18a,Mydata$neut18b,sep="")
    locus$genotype19 = paste(Mydata$neut19a,Mydata$neut19b,sep="")
    locus$genotype20 = paste(Mydata$neut20a,Mydata$neut20b,sep="")
    locus$genotype21 = paste(Mydata$neut21a,Mydata$neut21b,sep="")
    locus$genotype22 = paste(Mydata$neut22a,Mydata$neut22b,sep="")
    locus$genotype23 = paste(Mydata$neut23a,Mydata$neut23b,sep="")
    locus$genotype24 = paste(Mydata$neut24a,Mydata$neut24b,sep="")
    
    #Convert matrix to a  genind object (for the package adegenet). #The genind object can then easily be converted into loci objects (package pegas) (i.e. Mydata2)
    Mydata1 = df2genind(locus, ploidy = 2, sep="")
    Mydata2 = genind2loci(Mydata1)
    
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
    mantel.out = mantel.rtest(distgenDISS, dist.mat, nrepet = 500)
    
    Mantel.obs[g,s-16] = mantel.out$obs
    Mantel.p[g,s-16] = mantel.out$pvalue
    
    write.csv(Mantel.obs,'Neutral.2D.Mantel.obs.r.csv',row.names=FALSE)
    write.csv(Mantel.p,'Neutral.2D.Mantel.p.csv',row.names=FALSE)
  }  
}


##############################################################
## Plot spatial autocorrelation statistics
##############################################################
comment.out = function(){
  library(ggplot2)

  Mantel.obs = read.csv(paste(getwd(),'Neutral.2D.Mantel.obs.csv',sep="/"))
  Mantel.obs.plot = gather(Mantel.obs,'paraset','correlation')
  Mantel.obs.plot$generation = g.list
  
  Mantel.p = read.csv(paste(getwd(),'Neutral.2D.Mantel.p.csv',sep="/"))
  Mantel.p.plot = gather(Mantel.p,'paraset','significance')
  Mantel.p.plot$generation = g.list
  
  Mantel.plot = bind_cols(Mantel.obs.plot,Mantel.p.plot) %>% select(.,paraset,generation,correlation,significance)
  Mantel.plot$grouping[1:(nrow(Mantel.plot)/2)] = 'Selfing'
  Mantel.plot$grouping[(1+nrow(Mantel.plot)/2):nrow(Mantel.plot)] = 'No selfing'
  Mantel.plot[(Mantel.plot == 'paraset_1')|(Mantel.plot == 'paraset_5')] = 'Null'
  Mantel.plot[(Mantel.plot == 'paraset_2')|(Mantel.plot == 'paraset_6')] = 'IBT'
  Mantel.plot[(Mantel.plot == 'paraset_3')|(Mantel.plot == 'paraset_7')] = 'IBD'
  Mantel.plot[(Mantel.plot == 'paraset_4')|(Mantel.plot == 'paraset_8')] = 'IBDxIBT'
  names(Mantel.plot)[1] = 'Isolation'
  
  ggplot()+geom_line(data=Mantel.plot,aes(x=generation,y=correlation,col=Isolation))+
    geom_point(data=filter(Mantel.plot,significance<0.05),aes(x=generation,y=correlation,col=Isolation),shape=8)+
    facet_wrap(grouping~.)+
    xlab('Generation')+ylab('Correlation')+ggtitle('Correlation between spatial distance and genetic distance (number allelic differences between individuals)')+
    theme_classic()
}
