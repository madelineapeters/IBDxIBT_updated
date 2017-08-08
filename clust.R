library(mclust)

seq<-c(1,100, 200, 300, 400, 500)

for (p in 55){
  for (r in 1){
    start_wd<-(paste("~/498", paste("para_set", p, sep="_"), paste("model_run", r, sep="_"), sep="/"))
    setwd(start_wd)
    
    clust.df<-as.data.frame(matrix(nrow=20, ncol=length(seq)))
    names(clust.df)<-c("gen1", "gen250", "gen500")
    
    for (z in seq) {
      offspring_map <- read.csv(paste(getwd(), paste("offspring_map", z, sep="_"), sep="/"))
      
      fit <- Mclust(offspring_map[,1:3])
      
      BIC.plot<-plot(fit, what = c("BIC"), 
                     dimens = NULL, xlab = NULL, ylab = NULL, ylim = NULL, xlim = NULL,
                     addEllipses = TRUE, main = FALSE)
      title(main = "Spatial clustering of flowering time")
      mtext(paste("BIC, generation", z, sep=" "))
      dev.copy(png, paste(z, "BIC.png", sep="."))
      dev.off()
      
      class.plot<-plot(fit, what = c("classification"), 
                       dimens = c(2,3), xlab = TRUE, ylab = TRUE, ylim = NULL, xlim = NULL,
                       addEllipses = TRUE, main = FALSE)
      title(main = "Spatial flowering time clustering")
      mtext(paste("Clasffication, generation", z, sep=" "))
      dev.copy(png, paste(z, "class.png", sep="."))
      dev.off()
      
      uncert.plot<-plot(fit, what = c("uncertainty"), 
                        dimens = c(2,3) , xlab = NULL, ylab = NULL, ylim = NULL, xlim = NULL,
                        addEllipses = TRUE, main = FALSE)
      title(main = "Spatial flowering time clustering")
      mtext(paste("Uncertainty, generation", z, sep=" "))
      dev.copy(png, paste(z, "uncert.png", sep="."))
      dev.off()
      
      #fill in df
      clust.df[1,match(z,seq)]<-fit$modelName
      clust.df[2,match(z,seq)]<-fit$G
      G<-as.numeric(fit$G)
      clust.df[3,match(z,seq)]<-fit$bic
      clust.df[4,match(z,seq)]<-fit$loglik
      
      class.df<-as.data.frame(fit$classification)
      names(class.df)<-"class"
      for (g in 1:G) {
        clust.df[(4+g),match(z,seq)]<-dim(subset.data.frame(class.df, class == g))[1]
      }
      
      offspring_map$color <- factor(offspring_map$FLday,
                                    levels=unique(offspring_map$FLday),
                                    labels=pal(length(unique(offspring_map$FLday))))
      par(mfrow=c(1,3))
      plot(offspring_map$X_pos, offspring_map$Y_pos, type = 'p', col=as.character(offspring_map$color), pch=15, 
           xlab = "X position", ylab = "Y position", cex = 2.5)
      title(main = "Flowering time phenotype")
      plot(fit, what = c("classification"), 
           dimens = c(2,3) , xlab = "X position", ylab = "Y position", ylim = NULL, xlim = NULL,
           addEllipses = TRUE, main = FALSE)
      title(main = "Classification")
      plot(fit, what = c("uncertainty"), 
           dimens = c(2,3) , xlab = "X position", ylab = "Y position", ylim = NULL, xlim = NULL,
           addEllipses = TRUE, main = FALSE)
      title(main = "Uncertainty")
      dev.copy(png, paste(z, "comb.png", sep="."), width = 800, height = 300)
      dev.off()
      
    } #end z
    
    # now to write it out:
    write.csv(clust.df, # reorder columns to put CLUST first
              file="clust.df.csv",
              row.names=FALSE,                
              quote=FALSE)
    
  } #end r
} #end p
