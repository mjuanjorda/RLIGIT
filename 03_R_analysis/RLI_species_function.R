#--------------------------------------------------------
#--------------------------------------------------------
# R scripts for calculation the Red List Index (population and species level) for oceanic tunas, billfishes and sharks
# Juan-Jordá MJ, Murua H, Arrizabalaga H, Merino G, Pacoureau N, Dulvy NK. Science 378, eabj0211 (2022).DOI: 10.1126/science.abj0211
#READ THE FULL ARTICLE AT https://doi.org/10.1126/science.abj0211
#Corresponding author: Maria José Juan Jordá, mjuan.jorda@ieo.csic.es
#--------------------------------------------------------

#--------------------------------------------------------
# This script assissts in the calculation of the Red List Index of "populations" for oceanic tunas, billfishes and sharks
#--------------------------------------------------------



RLI_species_function<- function (dat){ #dat is the data.frame with the yearly RL threat status for all populations (52 populations)
  
  #create columns with the number of runs (number of iterations for each stock)
  number_runs<- 900 # number of model output runs per population and year
  
  #check total number of species
  number_species<-length(unique(dat$Species))
  
  #number of regresions (which is the number of Threat Status calculations for all stocks and years included)	
  number_of_regressions <-dim(dat)[1]/number_runs # total number of regressions carried out (for all the populations and years included)
  
  dat$runs<-rep(seq(1:number_runs),number_of_regressions)
  
  #calculate the WEIGHTED percent change for each SPECIES and year.
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
  
  
    #convert threat category to factor
  summary(global_data_speciesW50)
  global_data_speciesW50$threat_categoryFgl10<-factor(global_data_speciesW50$threat_categoryFgl10,levels=c("Critically Endangered", "Endangered","Vulnerable","Near Threatened","Least Concern"))
  

  ### threat code
  global_data_speciesWcode10<-global_data_speciesW50 %>% mutate(threat_codeFgl10=case_when(threat_categoryFgl10=="Critically Endangered" ~ 4,threat_categoryFgl10=="Endangered" ~ 3,threat_categoryFgl10=="Vulnerable" ~ 2,threat_categoryFgl10=="Near Threatened" ~ 1,threat_categoryFgl10=="Least Concern" ~ 0))
  
  # Fgl10 - calculate sum threat score, which is the weight of category for stocks at each year 
  
  rli_dataFgl10<-aggregate(global_data_speciesWcode50$threat_codeFgl10,by=list(runs=global_data_speciesWcode50$runs,years=global_data_speciesWcode50$yearlast),sum)
  
  head(rli_dataFgl10)
  names(rli_dataFgl10)[1:3]<-c("runs","years","sum_threat_score")
  
  # calculate max threat score
  rli_dataFgl10$max_threat_score<-number_species*5 # number_species * 5 which is the maximum threat score, equivalent to Extinct
  
  #calculate Red List Index --- Fgl10
  rli_dataFgl10$RLIFgl10<-(rli_dataFgl10$max_threat_score-rli_dataFgl10$sum_threat_score)/rli_dataFgl10$max_threat_score
  
  
  #Save outputs
  
  global_data_speciesW2<-global_data_speciesWcode50
  
  return(list(global_data_speciesW2,rli_dataFgl10[,c(1,2,5)]))
  
}


