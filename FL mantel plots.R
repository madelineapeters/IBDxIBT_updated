#load required packages
library(tidyr)
library(dplyr)
library(data.table)
library(doParallel)

cl<- makeCluster(2)
registerDoParallel(cl)

#make data frame for plotting
plot.df<-as.data.frame(c(1, seq(from=50, to=500, by=50)))
names(plot.df)<-("gen")

foreach (e = starting.set1:ending.set1) %dopar% {
  
  library(ggplot2)
  
  for (r in 1:no_runs) {
    #Set the working directory
    start_wd<-(paste(base_wd, sub_wd, paste("para_set", e, sep="_"), paste("model_run_", r, sep=""), sep="/"))
    setwd(start_wd)
    
    #ensure all dataframes are loaded into the environment
    assign(df.FL.r, read.csv(paste(start_wd, "/2D.mantel.FL.r", r, ".csv", sep=""), header = TRUE))
    df.FL.r<-subset(df.FL.r[,2])
    final.FL<-bind_cols(df.FL.r, plot.df)

    FL_plot<-ggplot()+geom_line(data=df.FL.r, aes(x=gen, y=V1))+
      xlab("Generation")+ylab("Spatial correlation coefficient r")+ylim(0,0.25)
    png(filename=paste(r, "2D mantel FL.png", sep="_"))
    plot(FL_plot)
    dev.off()
    
    if (r==1) {
      total.FL<-final.FL} else {
        total.FL<-bind_cols(total.FL, df.FL.r)
        names(total.FL[,r+1])<-paste("V", r, sep="")
        } #end loop over if/else
    
  } #end loop over runs

#reset working directory
start_wd<-paste(base_wd, sub_wd, paste("para_set", e, sep="_"), sep="/")
setwd(start_wd)

#load in average FL statistics
assign(df.FL.avg, read.csv(paste(start_wd, "FL.mantel.r.avg.csv", sep=""), header = TRUE))

df.FL.avg<-bind_cols(plot.df, df.FL.avg)
names(df.FL.avg)<-c("gen", "avg")

#make total average plot
avg.var_plot<-ggplot()+geom_line(data=total.FL, aes(x=gen, y=V1))+
  geom_line(data=total.FL, aes(x=gen, y=V2))+
  geom_line(data=total.FL, aes(x=gen, y=V3))+
  geom_line(data=total.FL, aes(x=gen, y=V4))+
  geom_line(data=total.FL, aes(x=gen, y=V5))+
  geom_line(data=total.FL, aes(x=gen, y=V6))+
  geom_line(data=total.FL, aes(x=gen, y=V7))+
  geom_line(data=total.FL, aes(x=gen, y=V8))+
  geom_line(data=total.FL, aes(x=gen, y=V9))+
  geom_line(data=total.FL, aes(x=gen, y=V10))+
  geom_line(data=df.FL.avg, aes(x=gen, y=avg), color="red", size=1.5)+
  xlab("Generation")+ylab("Spatial correlation coefficient r")+ylim(0,0.25)
png(filename=paste(e, "Avg 2D mantel FL.png", sep="_"))
plot(avg.var_plot)
dev.off()

} #end loop over parameter sets