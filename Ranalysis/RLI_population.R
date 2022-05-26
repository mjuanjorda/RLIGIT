#--------------------------------------------------------
#--------------------------------------------------------
# R script for calculation the Red List index of populations for oceanic tunas, billfishes and sharks
# Juan Jorda et al xxxxxxxx


RLI_population<- function (dat)
{
  
  #dat<- the data.frame with the yearly Red List status from 1950-2019 of all populations
  
  #create a column with the number of runs (number of iterations for each populations)
  number_runs<- 900 # number of model output runs per population and year
  number_stocks<-length(unique(dat$st))
  number_of_regressions <-dim(dat)[1]/number_runs # total number of regressions carried out for all the stocks and years included (which is the number of RL status calculations for all stocks and years included)
  dat$runs<-rep(seq(1:number_runs),number_of_regressions)
  # convert threat code to numeric
  dat$threat_code<-as.numeric(as.character(dat$threat_code))
  
  
  # calculate sum threat score, which is the weight of category for stocks at each year 
  rli_data<-aggregate(dat$threat_code,by=list(runs=dat$runs,years=dat$yearlast),sum)
  
  names(rli_data)[1]<-c("runs")
  names(rli_data)[2]<-c("years")
  names(rli_data)[3]<-c("sum_threat_score")
  
  # calculate max threat score
  rli_data$max_threat_score<-number_stocks*5 # number_stocks * 5 which is the maximum threat equivalent to extinct .
  
  #calculate RLI
  rli_data$RLI<-(rli_data$max_threat_score-rli_data$sum_threat_score)/rli_data$max_threat_score
  
  summary(rli_data)
  
  return(rli_data[,c(1,2,5)])
}

