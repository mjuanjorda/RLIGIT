
# STEP7 estimates the total percent population reduction over three generation lengths (or 10 years, whichever is longer) for each species by combining the information for all available subpopulations. For those species with multiple populations, we weighted the estimated total percent change in biomass of each population by their maximum sustainable yield to account for the contribution of different population sizes to the global species


rm(list=ls())

setwd("/Users/mjuanjorda/Dropbox/POSTDOC/MarieCurieProject/Redlistindex/RLI_january2015/RLI_2020/Threat_status_stocks_A1_A2/model_output_A1A2_varying_over_time_meanFgl")

load(file="Model_output_A1_A2_varying_over_time_based_meanFgl.RData")



head(table_save_all_stocks)

length(unique(table_save_all_stocks$st)) # Now 52 stocks

length(unique(table_save_all_stocks$sp)) #Now 18 sp

unique(table_save_all_stocks$st)

unique(table_save_all_stocks$threat_category)

levels(table_save_all_stocks$threat_category)


table_save_all_stocks$threat_category<-factor(table_save_all_stocks$threat_category,levels=c("Critically Endangered", "Endangered","Vulnerable","Near Threatened","Least Concern"))


#### Need to ADD columns with taxonomic information
#### ----------------------------------------------------
table_save_all_stocks$sp<-as.factor(table_save_all_stocks$sp)

library("dplyr")

table_save_all_stocks<-arrange(table_save_all_stocks, sp)

levels(table_save_all_stocks$sp)

#"ALB" "BET" "BFT" "BLM" "BSH" "BUM" "FAL" "MLS" "OCS" "PBF" "POR" "SAI" "SBT" "SKJ" "SMA" "SWO" "WHM" "YFT"

table_save_all_stocks$Taxa<-table_save_all_stocks$sp

levels(table_save_all_stocks$Taxa)<-c("Tunas","Tunas","Tunas","Billfishes","Sharks","Billfishes","Sharks","Billfishes","Sharks","Tunas","Sharks","Billfishes","Tunas","Tunas","Sharks","Billfishes","Billfishes","Tunas")

summary(table_save_all_stocks)


#-----------------------------------------------
# Combine data sets with stock weights
#-----------------------------------------------

setwd("/Users/mjuanjorda/Dropbox/POSTDOC/MarieCurieProject/Redlistindex/RLI_january2015/RLI_2020/RLI_GLOBAL")

datosWeights<-read.csv(file="Species_stocks_weight_to_calculate_SpeciesRLI_Global.csv",sep=",")
head(datosWeights)
str(datosWeights)
unique(datosWeights$SpStock)

unique(table_save_all_stocks$st)

#merge data sets

Model_output_A1_A2_varying_over_time_with_weights<-merge(table_save_all_stocks,datosWeights,by.x="st",by.y="SpStock",all=TRUE)


head(Model_output_A1_A2_varying_over_time_with_weights)


summary(Model_output_A1_A2_varying_over_time_with_weights)

names(Model_output_A1_A2_varying_over_time_with_weights)

Model_output_A1_A2_varying_over_time_with_weights<-subset(Model_output_A1_A2_varying_over_time_with_weights, select=-c(Group,Ocean,taxo))

names(Model_output_A1_A2_varying_over_time_with_weights)

setwd("/Users/mjuanjorda/Dropbox/POSTDOC/MarieCurieProject/Redlistindex/RLI_january2015/RLI_2020/Threat_status_stocks_A1_A2/model_output_A1A2_varying_over_time_meanFgl")


save(Model_output_A1_A2_varying_over_time_with_weights,file="Model_output_A1_A2_varying_over_time_meanFgl_with_weights.RData")



#-----------------------------------------------
# Estimate the weighted mean of the posteriors
#-----------------------------------------------

setwd("/Users/mjuanjorda/Dropbox/POSTDOC/MarieCurieProject/Redlistindex/RLI_january2015/RLI_2020/Threat_status_stocks_A1_A2/model_output_A1A2_varying_over_time_meanFgl")

load(file="Model_output_A1_A2_varying_over_time_meanFgl_with_weights.RData")

#create columns with the number of runs (number of iterations for each stock)
number_runs<- 900 # number of runs per stocks 

#number of species
number_species<-length(unique(Model_output_A1_A2_varying_over_time_with_weights$Species))

#number of regresions (which is the number of Threat Status calculations for all stocks and years included) 
number_of_regressions <-dim(Model_output_A1_A2_varying_over_time_with_weights)[1]/number_runs # total number of regressions carried out (for all the stocks and years included)

Model_output_A1_A2_varying_over_time_with_weights$runs<-rep(seq(1:number_runs),number_of_regressions)

#calculate the WEIGHTED tpercent change for each SPECIES and year. Later I will have to re-asign a new Threat Category and threat code.
library(plyr)
library(dplyr)

head(Model_output_A1_A2_varying_over_time_with_weights)

global_data_speciesW<-ddply(Model_output_A1_A2_varying_over_time_with_weights,.(Taxa,sp,yearlast,Fstatus_species_meanFgl10,Fstatus_species_meanFgl20,Fstatus_species_meanFgl30,Fstatus_species_meanFgl40,Fstatus_species_meanFgl50,runs),summarise, wm = weighted.mean(tpercent_change3GL,FinalWeights))

head(global_data_speciesW)
summary(global_data_speciesW)


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

global_data_speciesW20<-global_data_speciesW10 %>%
  mutate(threat_categoryFgl20 = 
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


summary(global_data_speciesW50)


global_data_speciesW50$threat_categoryFgl10<-factor(global_data_speciesW50$threat_categoryFgl10,levels=c("Critically Endangered", "Endangered","Vulnerable","Near Threatened","Least Concern"))

global_data_speciesW50$threat_categoryFgl20<-factor(global_data_speciesW50$threat_categoryFgl20,levels=c("Critically Endangered", "Endangered","Vulnerable","Near Threatened","Least Concern"))

global_data_speciesW50$threat_categoryFgl30<-factor(global_data_speciesW50$threat_categoryFgl30,levels=c("Critically Endangered", "Endangered","Vulnerable","Near Threatened","Least Concern"))

global_data_speciesW50$threat_categoryFgl40<-factor(global_data_speciesW50$threat_categoryFgl40,levels=c("Critically Endangered", "Endangered","Vulnerable","Near Threatened","Least Concern"))

global_data_speciesW50$threat_categoryFgl50<-factor(global_data_speciesW50$threat_categoryFgl50,levels=c("Critically Endangered", "Endangered","Vulnerable","Near Threatened","Least Concern"))

summary(global_data_speciesW50)

