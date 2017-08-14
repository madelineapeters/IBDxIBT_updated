library(readr)
library(ggplot2)

for (p in 109:111){
  for (r in 3){
    start_wd<-(paste("C:/Users/Madeline/Desktop/Weis lab/EEB498", paste("para_set", p, sep="_"), paste("model_run", r, sep="_"), sep="/"))
    setwd(start_wd)
    
    for (z in c(1, seq(from = 50, to = 500, by = 50))) {
    offspring_map <- read.csv(paste(getwd(), paste("offspring_map", z, sep="_"), sep="/"))
    
    
    FLday_hist<-hist(offspring_map$FLday,
                     main="",
                     xlab="Flowering day",
                     ylab="Count",
                     border="black", 
                     col="black",
                     xlim=c(75,125),
                     ylim=c(0,400),
                     las=1, 
                     breaks=40)
    FLday_hist<-ggplot(data=offspring_map, aes(x=FLday))+geom_histogram(binwidth=1, aes(fill=..x..))+
      xlab("Flowering day")+ylab("Count")+xlim(85,115)+ylim(0,200)+
      scale_fill_gradientn(
        colours = c("lightgreen", "olivedrab2", "green4", "darkgreen"),
        breaks=c(85, 95, 105, 115, Inf), limits = c(85,115),
        na.value = "red"
      ) + theme_bw()
    plot(FLday_hist)
    heatmap_FLday<-ggplot(data = offspring_map, aes(x=X_pos, y=Y_pos, fill=FLday)) + 
      geom_tile()+scale_fill_gradient(name = "Flowering day", low= "darkolivegreen1" ,high="green4")+xlab("X position")+ylab("Y position")
    print(heatmap_FLday)
    
    heatmap_Alocus<-ggplot(data = offspring_map, aes(x=X_pos, y=Y_pos, fill=mapA)) + 
      geom_tile() + scale_fill_gradient(name = "", low= "violet" ,high="violetred4") +
      xlab("X position") + ylab("Y position") + ggtitle(paste("A locus", " - Generation", z, sep=" ")) +
      theme(
        plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12))
    print(heatmap_Alocus)
    
    heatmap_Blocus<-ggplot(data = offspring_map, aes(x=X_pos, y=Y_pos, fill=mapB)) + 
      geom_tile() + scale_fill_gradient(name = "", low= "orange2" ,high="orangered4") +
      xlab("X position") + ylab("Y position") + ggtitle(paste("B locus", " - Generation", z, sep=" ")) +
      theme(
        plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12))
    print(heatmap_Blocus)
    
    heatmap_Clocus<-ggplot(data = offspring_map, aes(x=X_pos, y=Y_pos, fill=mapC)) + 
      geom_tile() + scale_fill_gradient(name = "", low= "steelblue1" ,high="steelblue4") +
      xlab("X position") + ylab("Y position") + ggtitle(paste("C locus", " - Generation", z, sep=" ")) +
      theme(
        plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12))
    print(heatmap_Clocus)
    
    png(filename=paste(z, "FLday.png", sep="_"))
    plot(heatmap_FLday)
    dev.off()
    
    png(filename=paste(z, "Alocus.png", sep="_"))
    plot(heatmap_Alocus)
    dev.off()
    
    png(filename=paste(z, "Blocus.png", sep="_"))
    plot(heatmap_Blocus)
    dev.off()
    
    png(filename=paste(z, "Clocus.png", sep="_"))
    plot(heatmap_Clocus)
    dev.off()
    
    png(filename=paste(z, "FLhist.png", sep="_"))
    plot(FLday_hist)
    dev.off()

} #end loop over generations
} #end loop over runs
} #end loop over parameter set
