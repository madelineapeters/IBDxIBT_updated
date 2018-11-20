library(readr)
library(ggplot2)

for (p in 82:85){
  for (r in 1){
    start_wd<-(paste(base_wd, paste("para_set", p, sep="_"), paste("model_run", r, sep="_"), sep="/"))
    setwd(start_wd)
    
    for (z in c(1, seq(from = 10, to = 100, by = 10))) {
    offspring_map <- read.csv(paste(getwd(), paste("NH", z, "csv", sep="."), sep="/"))
    
    
    FLday_hist<-hist(offspring_map$FLday,
                     main="",
                     xlab="Flowering day",
                     ylab="Count",
                     border="black", 
                     col="black",
                     xlim=c(75,125),
                     ylim=c(0,100),
                     las=1, 
                     breaks=40)
    
    heatmap_FLday<-ggplot(data = offspring_map, aes(x=X_pos, y=Y_pos, fill=FLday)) + 
      geom_tile()+scale_fill_gradient(name = "Flowering day", low= "darkolivegreen1" ,high="green4")+xlab("X position")+ylab("Y position")
    print(heatmap_FLday)
    
    png(filename=paste(z, "FLday.png", sep="_"))
    plot(heatmap_FLday)
    dev.off()
    
    png(filename=paste(z, "FLhist.png", sep="_"))
    plot(FLday_hist)
    dev.off()

} #end loop over generations
} #end loop over runs
} #end loop over parameter set
