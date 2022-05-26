#--------------------------------------------------------
#--------------------------------------------------------
# R script for calculation the Red List index of species for oceanic tunas, billfishes and sharks
# Juan Jorda et al xxxxxxxx
# 
# 

IGNORE for now

RLI_species<- function (dat)
{
  
  #dat<- the data.frame with the yearly RL threat status for all populations (52 populations)
  
  #create columns with the number of runs (number of iterations for each stock)
  number_runs<- 900 # number of model output runs per population and year
  
  #check total number of species
  number_species<-length(unique(dat$Species))
  
  #number of regresions (which is the number of Threat Status calculations for all stocks and years included)	
  number_of_regressions <-dim(dat)[1]/number_runs # total number of regressions carried out (for all the populations and years included)
  
  dat$runs<-rep(seq(1:number_runs),number_of_regressions)
  
  #calculate the WEIGHTED tpercent change for each SPECIES and year. Later I will have to re-asign a new Threat Category and threat code.
  library(plyr)
  global_data_speciesW<-ddply(dat,.(Taxa,Species,Fstatus_species_meanFgl10,Fstatus_species_meanFgl20,Fstatus_species_meanFgl30,Fstatus_species_meanFgl40,Fstatus_species_meanFgl50,yearlast,runs),summarise, wm = weighted.mean(tpercent_change3GL,FinalWeights))
  
  #assiged CriteriaA1 or CriteriaA2 based on F
  
  global_data_speciesW10<-global_data_speciesW %>%
    mutate(threat_categoryFgl10 = 
             ifelse(Fstatus_species_meanFgl10=="Not Overfishing" & (wm < -90), "Critically Endangered", 
                    ifelse(Fstatus_species_meanFgl10=="Not Overfishing" & (wm >= -90 & wm < -70 ), "Endangered", 
                           ifelse(Fstatus_species_meanFgl10=="Not Overfishing" & (wm >= -70 & wm < -50 ), "Vulnerable", 
                                  ifelse(Fstatus_species_meanFgl10=="Not Overfishing" & (wm >= -50 & wm < -40 ), "Near Threatened", 
                                         ifelse(Fstatus_species_meanFgl10=="Not Overfishing" & (wm >= -40), "Least Concern",
                                                ifelse(Fstatus_species_meanFgl10=="Overfishing" & (wm < -80), "Critically Endangered", 
                                                       ifelse(Fstatus_species_meanFgl10=="Overfishing" & (wm >= -80 & wm < -50 ), "Endangered", 
                                                              ifelse(Fstatus_species_meanFgl10=="Overfishing" & (wm >= -50 & wm < -30 ), "Vulnerable", 
                                                                     ifelse(Fstatus_species_meanFgl10=="Overfishing" & (wm >= -30 & wm < -20 ), "Near Threatened", 
                                                                            ifelse(Fstatus_species_meanFgl10=="Overfishing" & (wm >= -20), "Least Concern",NA)))))))))))
  
  
  global_data_speciesW20<-global_data_speciesW10 %>% mutate(threat_categoryFgl20 = 
                                                              ifelse(Fstatus_species_meanFgl20=="Not Overfishing" & (wm < -90), "Critically Endangered", 
                                                                     ifelse(Fstatus_species_meanFgl20=="Not Overfishing" & (wm >= -90 & wm < -70 ), "Endangered", 
                                                                            ifelse(Fstatus_species_meanFgl20=="Not Overfishing" & (wm >= -70 & wm < -50 ), "Vulnerable", 
                                                                                   ifelse(Fstatus_species_meanFgl20=="Not Overfishing" & (wm >= -50 & wm < -40 ), "Near Threatened", 
                                                                                          ifelse(Fstatus_species_meanFgl20=="Not Overfishing" & (wm >= -40), "Least Concern",
                                                                                                 ifelse(Fstatus_species_meanFgl20=="Overfishing" & (wm < -80), "Critically Endangered", 
                                                                                                        ifelse(Fstatus_species_meanFgl20=="Overfishing" & (wm >= -80 & wm < -50 ), "Endangered", 
                                                                                                               ifelse(Fstatus_species_meanFgl20=="Overfishing" & (wm >= -50 & wm < -30 ), "Vulnerable", 
                                                                                                                      ifelse(Fstatus_species_meanFgl20=="Overfishing" & (wm >= -30 & wm < -20 ), "Near Threatened", 
                                                                                                                             ifelse(Fstatus_species_meanFgl20=="Overfishing" & (wm >= -20), "Least Concern",NA)))))))))))
  
  global_data_speciesW30<-global_data_speciesW20 %>%
    mutate(threat_categoryFgl30 = 
             ifelse(Fstatus_species_meanFgl30=="Not Overfishing" & (wm < -90), "Critically Endangered", 
                    ifelse(Fstatus_species_meanFgl30=="Not Overfishing" & (wm >= -90 & wm < -70 ), "Endangered", 
                           ifelse(Fstatus_species_meanFgl30=="Not Overfishing" & (wm >= -70 & wm < -50 ), "Vulnerable", 
                                  ifelse(Fstatus_species_meanFgl30=="Not Overfishing" & (wm >= -50 & wm < -40 ), "Near Threatened", 
                                         ifelse(Fstatus_species_meanFgl30=="Not Overfishing" & (wm >= -40), "Least Concern",
                                                ifelse(Fstatus_species_meanFgl30=="Overfishing" & (wm < -80), "Critically Endangered", 
                                                       ifelse(Fstatus_species_meanFgl30=="Overfishing" & (wm >= -80 & wm < -50 ), "Endangered", 
                                                              ifelse(Fstatus_species_meanFgl30=="Overfishing" & (wm >= -50 & wm < -30 ), "Vulnerable", 
                                                                     ifelse(Fstatus_species_meanFgl30=="Overfishing" & (wm >= -30 & wm < -20 ), "Near Threatened", 
                                                                            ifelse(Fstatus_species_meanFgl30=="Overfishing" & (wm >= -20), "Least Concern",NA)))))))))))
  
  
  global_data_speciesW40<-global_data_speciesW30 %>%
    mutate(threat_categoryFgl40 = 
             ifelse(Fstatus_species_meanFgl40=="Not Overfishing" & (wm < -90), "Critically Endangered", 
                    ifelse(Fstatus_species_meanFgl40=="Not Overfishing" & (wm >= -90 & wm < -70 ), "Endangered", 
                           ifelse(Fstatus_species_meanFgl40=="Not Overfishing" & (wm >= -70 & wm < -50 ), "Vulnerable", 
                                  ifelse(Fstatus_species_meanFgl40=="Not Overfishing" & (wm >= -50 & wm < -40 ), "Near Threatened", 
                                         ifelse(Fstatus_species_meanFgl40=="Not Overfishing" & (wm >= -40), "Least Concern",
                                                ifelse(Fstatus_species_meanFgl40=="Overfishing" & (wm < -80), "Critically Endangered", 
                                                       ifelse(Fstatus_species_meanFgl40=="Overfishing" & (wm >= -80 & wm < -50 ), "Endangered", 
                                                              ifelse(Fstatus_species_meanFgl40=="Overfishing" & (wm >= -50 & wm < -30 ), "Vulnerable", 
                                                                     ifelse(Fstatus_species_meanFgl40=="Overfishing" & (wm >= -30 & wm < -20 ), "Near Threatened", 
                                                                            ifelse(Fstatus_species_meanFgl40=="Overfishing" & (wm >= -20), "Least Concern",NA)))))))))))
  
  
  global_data_speciesW50<-global_data_speciesW40 %>%
    mutate(threat_categoryFgl50 = 
             ifelse(Fstatus_species_meanFgl50=="Not Overfishing" & (wm < -90), "Critically Endangered", 
                    ifelse(Fstatus_species_meanFgl50=="Not Overfishing" & (wm >= -90 & wm < -70 ), "Endangered", 
                           ifelse(Fstatus_species_meanFgl50=="Not Overfishing" & (wm >= -70 & wm < -50 ), "Vulnerable", 
                                  ifelse(Fstatus_species_meanFgl50=="Not Overfishing" & (wm >= -50 & wm < -40 ), "Near Threatened", 
                                         ifelse(Fstatus_species_meanFgl50=="Not Overfishing" & (wm >= -40), "Least Concern",
                                                ifelse(Fstatus_species_meanFgl50=="Overfishing" & (wm < -80), "Critically Endangered", 
                                                       ifelse(Fstatus_species_meanFgl50=="Overfishing" & (wm >= -80 & wm < -50 ), "Endangered", 
                                                              ifelse(Fstatus_species_meanFgl50=="Overfishing" & (wm >= -50 & wm < -30 ), "Vulnerable", 
                                                                     ifelse(Fstatus_species_meanFgl50=="Overfishing" & (wm >= -30 & wm < -20 ), "Near Threatened", 
                                                                            ifelse(Fstatus_species_meanFgl50=="Overfishing" & (wm >= -20), "Least Concern",NA)))))))))))
  
  
  # convert threat category to factor
  
  summary(global_data_speciesW50)
  
  
  global_data_speciesW50$threat_categoryFgl10<-factor(global_data_speciesW50$threat_categoryFgl10,levels=c("Critically Endangered", "Endangered","Vulnerable","Near Threatened","Least Concern"))
  
  global_data_speciesW50$threat_categoryFgl20<-factor(global_data_speciesW50$threat_categoryFgl20,levels=c("Critically Endangered", "Endangered","Vulnerable","Near Threatened","Least Concern"))
  
  global_data_speciesW50$threat_categoryFgl30<-factor(global_data_speciesW50$threat_categoryFgl30,levels=c("Critically Endangered", "Endangered","Vulnerable","Near Threatened","Least Concern"))
  
  global_data_speciesW50$threat_categoryFgl40<-factor(global_data_speciesW50$threat_categoryFgl40,levels=c("Critically Endangered", "Endangered","Vulnerable","Near Threatened","Least Concern"))
  
  global_data_speciesW50$threat_categoryFgl50<-factor(global_data_speciesW50$threat_categoryFgl50,levels=c("Critically Endangered", "Endangered","Vulnerable","Near Threatened","Least Concern"))
  
  summary(global_data_speciesW50)
  
  
  
  
  
  ### threatcode
  
  
  global_data_speciesWcode10<-global_data_speciesW50 %>% mutate(threat_codeFgl10=case_when(threat_categoryFgl10=="Critically Endangered" ~ 4,threat_categoryFgl10=="Endangered" ~ 3,threat_categoryFgl10=="Vulnerable" ~ 2,threat_categoryFgl10=="Near Threatened" ~ 1,threat_categoryFgl10=="Least Concern" ~ 0))
  
  
  global_data_speciesWcode20<-global_data_speciesWcode10 %>% mutate(threat_codeFgl20=case_when(threat_categoryFgl20=="Critically Endangered" ~ 4,threat_categoryFgl20=="Endangered" ~ 3,threat_categoryFgl20=="Vulnerable" ~ 2,threat_categoryFgl20=="Near Threatened" ~ 1,threat_categoryFgl20=="Least Concern" ~ 0))
  
  global_data_speciesWcode30<-global_data_speciesWcode20 %>% mutate(threat_codeFgl30=case_when(threat_categoryFgl30=="Critically Endangered" ~ 4,threat_categoryFgl30=="Endangered" ~ 3,threat_categoryFgl30=="Vulnerable" ~ 2,threat_categoryFgl30=="Near Threatened" ~ 1,threat_categoryFgl30=="Least Concern" ~ 0))
  
  global_data_speciesWcode40<-global_data_speciesWcode30 %>% mutate(threat_codeFgl40=case_when(threat_categoryFgl40=="Critically Endangered" ~ 4,threat_categoryFgl40=="Endangered" ~ 3,threat_categoryFgl40=="Vulnerable" ~ 2,threat_categoryFgl40=="Near Threatened" ~ 1,threat_categoryFgl40=="Least Concern" ~ 0))
  
  
  global_data_speciesWcode50<-global_data_speciesWcode40 %>% mutate(threat_codeFgl50=case_when(threat_categoryFgl50=="Critically Endangered" ~ 4,threat_categoryFgl50=="Endangered" ~ 3,threat_categoryFgl50=="Vulnerable" ~ 2,threat_categoryFgl50=="Near Threatened" ~ 1,threat_categoryFgl50=="Least Concern" ~ 0))
  
  
  # Fgl10 - calculate sum threat score, which is the weight of category for stocks at each year 
  
  rli_dataFgl10<-aggregate(global_data_speciesWcode50$threat_codeFgl10,by=list(runs=global_data_speciesWcode50$runs,years=global_data_speciesWcode50$yearlast),sum)
  
  head(rli_dataFgl10)
  names(rli_dataFgl10)[1:3]<-c("runs","years","sum_threat_score")
  # calculate max threat score
  rli_dataFgl10$max_threat_score<-number_species*5 # number_species * 5 which is the maximum threat equivalent to extinct .
  
  #calculate RLI ---Fgl10
  rli_dataFgl10$RLIFgl10<-(rli_dataFgl10$max_threat_score-rli_dataFgl10$sum_threat_score)/rli_dataFgl10$max_threat_score
  
  
  # Fgl20 - calculate sum threat score, which is the weight of category for stocks at each year 
  
  rli_dataFgl20<-aggregate(global_data_speciesWcode50$threat_codeFgl20,by=list(runs=global_data_speciesWcode50$runs,years=global_data_speciesWcode50$yearlast),sum)
  
  head(rli_dataFgl20)
  names(rli_dataFgl20)[1:3]<-c("runs","years","sum_threat_score")
  # calculate max threat score
  rli_dataFgl20$max_threat_score<-number_species*5 # number_species * 5 which is the maximum threat equivalent to extinct .
  
  #calculate RLI ---Fgl20
  rli_dataFgl20$RLIFgl20<-(rli_dataFgl20$max_threat_score-rli_dataFgl20$sum_threat_score)/rli_dataFgl20$max_threat_score
  
  
  
  # Fgl30 - calculate sum threat score, which is the weight of category for stocks at each year 
  
  rli_dataFgl30<-aggregate(global_data_speciesWcode50$threat_codeFgl30,by=list(runs=global_data_speciesWcode50$runs,years=global_data_speciesWcode50$yearlast),sum)
  
  head(rli_dataFgl30)
  names(rli_dataFgl30)[1:3]<-c("runs","years","sum_threat_score")
  # calculate max threat score
  rli_dataFgl30$max_threat_score<-number_species*5 # number_species * 5 which is the maximum threat equivalent to extinct .
  
  #calculate RLI ---Fgl30
  rli_dataFgl30$RLIFgl30<-(rli_dataFgl30$max_threat_score-rli_dataFgl30$sum_threat_score)/rli_dataFgl30$max_threat_score
  
  # Fgl40 - calculate sum threat score, which is the weight of category for stocks at each year 
  
  rli_dataFgl40<-aggregate(global_data_speciesWcode50$threat_codeFgl40,by=list(runs=global_data_speciesWcode50$runs,years=global_data_speciesWcode50$yearlast),sum)
  
  head(rli_dataFgl40)
  names(rli_dataFgl40)[1:3]<-c("runs","years","sum_threat_score")
  # calculate max threat score
  rli_dataFgl40$max_threat_score<-number_species*5 # number_species * 5 which is the maximum threat equivalent to extinct .
  
  #calculate RLI ---Fgl40
  rli_dataFgl40$RLIFgl40<-(rli_dataFgl40$max_threat_score-rli_dataFgl40$sum_threat_score)/rli_dataFgl40$max_threat_score
  
  # Fgl50 - calculate sum threat score, which is the weight of category for stocks at each year 
  
  rli_dataFgl50<-aggregate(global_data_speciesWcode50$threat_codeFgl50,by=list(runs=global_data_speciesWcode50$runs,years=global_data_speciesWcode50$yearlast),sum)
  
  head(rli_dataFgl50)
  names(rli_dataFgl50)[1:3]<-c("runs","years","sum_threat_score")
  # calculate max threat score
  rli_dataFgl50$max_threat_score<-number_species*5 # number_species * 5 which is the maximum threat equivalent to extinct .
  
  #calculate RLI ---Fgl50
  rli_dataFgl50$RLIFgl50<-(rli_dataFgl50$max_threat_score-rli_dataFgl50$sum_threat_score)/rli_dataFgl50$max_threat_score
  
  
  #Save outputs
  
  global_data_speciesW2<-global_data_speciesWcode50
  
  return(list(global_data_speciesW2,rli_dataFgl10[,c(1,2,5)],rli_dataFgl20[,c(1,2,5)],rli_dataFgl30[,c(1,2,5)],rli_dataFgl40[,c(1,2,5)],rli_dataFgl50[,c(1,2,5)]))
  
}


