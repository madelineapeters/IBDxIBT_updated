#load required packages
library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)


for (p in starting.set1:ending.set2) {
  #Set the working directory
  start_wd<-(paste(base_wd, sub_wd, paste("para_set", p, sep="_"), sep="/"))
  setwd(start_wd)
  
  #ensure all dataframes are loaded into the environment
  df.Fst<-paste("Fst.", p, sep="")
  assign(df.Fst, read.csv(paste("2.Fst.", p, ".csv", sep=""), header = TRUE))

} #end loop over parameter sets

rm(df.Fst)

#put dataframes into lists
dflist.Fst.avg<-lapply(ls(pattern = "Fst.*"), get)

#bind data frames and rename columns
total.Fst.avg<-rbind(dflist.Fst.avg[[1]], dflist.Fst.avg[[2]])
for (w in 3:length(dflist.Fst.avg)){
  total.Fst.avg<-rbind(total.Fst.avg, dflist.Fst.avg[[w]])
}

names(total.Fst.avg)<-c("index","Fst")

#create dataframe to be added with generation labels and parameter set labels
plot.df<-as.data.frame(matrix(nrow=(length(starting.set1:ending.set2)*11), ncol=4))
names(plot.df)<-c("set", "row", "gen", "col")
for (o in 1:length(starting.set1:ending.set2)){
  plot.df[(((o-1)*11+1):(o*11)),1]<-(starting.set1-1+o)
}
for (i in 1:(length(starting.set1:ending.set2)*11)) {
  plot.df[i,2]<- if ((plot.df[i,1] == 42)||(plot.df[i,1] == 48)||(plot.df[i,1] == 43)||(plot.df[i,1] == 49)||(plot.df[i,1] == 54)||(plot.df[i,1] == 60)||(plot.df[i,1] == 55)||(plot.df[i,1] == 61)){
    1} else if ((plot.df[i,1] == 44)||(plot.df[i,1] == 50)||(plot.df[i,1] == 45)||(plot.df[i,1] == 51)||(plot.df[i,1] == 56)||(plot.df[i,1] == 62)||(plot.df[i,1] == 57)||(plot.df[i,1] == 63)){
      2} else {3}
  plot.df[i,4]<- if ((plot.df[i,1] == 42)||(plot.df[i,1] == 44)||(plot.df[i,1] == 46)){
    1} else if ((plot.df[i,1] == 48)||(plot.df[i,1] == 50)||(plot.df[i,1] == 52)){
      2} else if ((plot.df[i,1] == 43)||(plot.df[i,1] == 45)||(plot.df[i,1] == 47)){
        3} else if ((plot.df[i,1] == 49)||(plot.df[i,1] == 51)||(plot.df[i,1] == 53)){
          4} else if ((plot.df[i,1] == 54)||(plot.df[i,1] == 56)||(plot.df[i,1] == 58)){
            5} else if ((plot.df[i,1] == 60)||(plot.df[i,1] == 62)||(plot.df[i,1] == 64)){
              6} else if ((plot.df[i,1] == 55)||(plot.df[i,1] == 57)||(plot.df[i,1] == 59)){
                7} else {8}
}
for (i in 1:length(starting.set1:ending.set2)) {
plot.df[((i-1)*11+1):(i*11),3]<-c(1, seq(from = 50, to = 500, by  = 50))
}

#attach data frames together
total.Fst.avg<-data.frame(total.Fst.avg, plot.df)

#make plot
Fst_plot<-ggplot()+geom_point(data=total.Fst.avg, aes(x=gen, y=Fst))+
  xlab("Generation")+ylab("Averaged pairwise Fst")+ylim(-0.01,0.1)+facet_grid(row~col)
print(Fst_plot)

