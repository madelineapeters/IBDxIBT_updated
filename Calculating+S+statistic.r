
#Set parameters
fl_mean<-100 #mean flowering day in generation 1
fl_seq<-c(2, 3, 5, 8, 11, 13, 15, 17, 15, 13, 11, 8, 5, 3, 2) #flowers produced per day
duration<-length(fl_seq) #length of flowering period
no_run<-5 #number of runs
seq<-c(1, seq(from=10, to=100, by=10)) #generations to test
seed_lambda<-1 #input into seed dispersal equation - larger numbers produce lower dispersal distances
pop_size<-400
space<-"yes" #are we calculating S and number mating pools based on space?
time<-"yes" #are we calculating S and number mating pools based on time?
torus<-"no" #were the parameter sets run on a torus?

#Start loop over parameter sets
for (e in 85) {
  
  #set metres per unit based on parameter set
  m_per_unit<-if ((e == 86) | (e == 87)) {0.5} else {1}
  
  #Start loop over runs
  for (r in 1:no_runs) {

S.N.df<-as.data.frame(matrix(nrow=length(seq), ncol=3)) #creates dataframe with one row per generation tested
names(S.N.df)<-c("gen", "S", "num.pool")
S.N.df[1]<-seq #fills in first column with generation labels from seq

#Start loop over generations
for (g in seq) {

  #Set working directory
  start_wd<-(paste("C:/Users/Madeline/Desktop/Weis lab/EEB498", paste("para_set", e, sep="_"), paste("model_run", r, sep="_"), sep="/"))
  setwd(start_wd)

  #read in generation dataframe
  offspring_map <- read.csv(paste(getwd(), paste("NH", g, "csv", sep="."), sep="/"))

if (time == "yes") {
    #set days of flowering schedule
    days<-c(min(offspring_map$FLday):(max(offspring_map$FLday)+(duration-1)))

    #fill in flowering schedule using parental generation abundances
    flowers<-as.data.frame(matrix(nrow=pop_size, ncol=length(days)))
    names(flowers)<-paste("d", days, sep="")
    for (i in 1:pop_size){
      flowers[i, which(days==offspring_map$FLday[i]):(which(days==offspring_map$FLday[i])+duration-1)]<-fl_seq
    }
    flowers[is.na(flowers)]<-0

    #set maternal matrix as the proportion of flowers made by a single individual relative to those made by the total population (both over entire season)
    mat<-as.matrix(flowers/sum(flowers))

    #set paternal matrix as the proportion of flowers made by a single individual relative to the entire population
    pat<-matrix(nrow=nrow(flowers), ncol=ncol(flowers))
    for (t in 1:ncol(pat)) {
      pat[,t]<-flowers[,t]/sum(flowers[,t])
    }

    #fix pat matrix to deal with days where NO flowers are produced in the population (this removes NaN)
    if (any(colSums(flowers)==0)){
        pat[,which(colSums(flowers)==0)]<-0
    }

    #place moms as rows and dads as columns
    mat_opp<-mat %*% t(pat)
} else {
    mat_opp<-matrix(nrow=pop_size, ncol=pop_size)
    mat_opp[,]<-1/(pop_size^2)
}

if (space == "yes") {

    #create dispersal mating matrices
    distance1<-matrix(nrow=pop_size, ncol=pop_size)
    distance2<-matrix(nrow=pop_size, ncol=pop_size)
    for (m in 1:pop_size) {
      for (p in 1:pop_size) {
        distance1[m,p]<-sqrt((offspring_map$X_pos[m]-offspring_map$X_pos[p])^2 + (offspring_map$Y_pos[m] - offspring_map$Y_pos[p])^2)
        distance2[m,p]<-sqrt(((sqrt(pop_size)+1)-abs(offspring_map$X_pos[m]-offspring_map$X_pos[p]))^2 + ((sqrt(pop_size)+1)-abs(offspring_map$Y_pos[m]-offspring_map$Y_pos[p]))^2)
      }
    }

    #adjust for metres per unit
    distance1<-distance1*m_per_unit
    distance2<-distance2*m_per_unit

    #run through pollen dispersal function
    distance_for_mating1<-distance1
    distance_for_mating2<-distance2
    
    pollen_dispersal1<-(((0.375)^2)/(2*pi))*exp(-0.375*distance_for_mating1)
    pollen_dispersal2<-(((0.375)^2)/(2*pi))*exp(-0.375*distance_for_mating2)


    #recalibrate mating opportunity matrix to account for space
    if (torus == "yes") {
        mat_opp_adj<-mat_opp * (pollen_dispersal1 + pollen_dispersal2)
        mat_opp_adj<-mat_opp_adj/sum(mat_opp_adj)
      } else {
        mat_opp_adj<-mat_opp * (pollen_dispersal1)
        mat_opp_adj<-mat_opp_adj/sum(mat_opp_adj)
      }
    
    
  } else {mat_opp_adj<-mat_opp} #end if statement

#eigenvalue extration and calculation of S statistic
EV<-eigen(mat_opp_adj, sym=FALSE) #calculate eigenvalues
ev.vector<-EV$values #get vector of eigenvalues
ev.1<-ev.vector[1] #store first eigenvalue
ev.sum<-sum(ev.vector) #store sum of all eigenvalues
S<-(ev.1/sum(ev.vector)) #calculate S
num_pools<-ev.sum/ev.1 #calculate number of mating pools

 #store output
  if (g == 1) {
    S.N.df[g,2]<-as.double(S)
    S.N.df[g,3]<-as.double(num_pools)
  } else {
    S.N.df[(g/10)+1,2]<-as.double(S)
    S.N.df[(g/10)+1,3]<-as.double(num_pools)
  }

} #end loop over g

write.csv(S.N.df, paste("S.space.df", r, "csv", sep="."))

} #end loop over runs
} #end loop over parameter sets
