#--------------------------------------------------------
#--------------------------------------------------------
# R scripts for calculation the Red List Index (population and species level) for oceanic tunas, billfishes and sharks
# Juan-Jordá MJ, Murua H, Arrizabalaga H, Merino G, Pacoureau N, Dulvy NK. Science 378, eabj0211 (2022).DOI: 10.1126/science.abj0211
#READ THE FULL ARTICLE AT https://doi.org/10.1126/science.abj0211
#Corresponding author: Maria José Juan Jordá, mjuan.jorda@ieo.csic.es
#--------------------------------------------------------


#--------------------------------------------------------
# This script calculates the continuous Red List Index of oceanic predatory fishes between 1950 and 2019 using the estimated Red List extinction risk of the 18 species of tunas, billfishes, and sharks. We applied the IUCN Criteria A "population reduction" to calculate the Red List extinction risk status for all species. Criterion A measures the extinction risk based on exceeding thresholds of population decline over a time frame of three generation lengths (or a minimum of 10 years, whichever is longer).  


# These steps are followed below:
# STEP1 loads the population trends of each species (and subpopulations) 
# STEP2 loads the generation lengths for each species (and subpopulations)
# STEP3 loads the fishing mortality trends of the species (and subpopulations) for understanding whether a population is being sustainably fished (or not) and also whether a species, as a whole, is being sustainably fished (or not) according to IUCN Red List criteria
# STEP4 estimates the total percent population reduction over three generation lengths (or 10 years, whichever is longer) for each subpopulation
# STEP5 applies Criterion A population reduction A1 or A2 thresholds for the determination of the IUCN Red List categories for each subpopulation
# STEP6 calculates the Red List Index of populations (including 52 populations of tunas, billfishes and sharks)
# STEP7 calculates the Red List Index of species (including 18 species of tunas, billfishes and sharks)


#---------------


# LOAD REQUIRED PACKAGES ####
library(rstudioapi)
library(plyr)
library(dplyr)
library(ggplot2)
library(rstan)
library(readbulk)
library(here)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#-------------------------------------------------------------------------
# STEP1 loads the population trends of each species (and subpopulations) 
#-------------------------------------------------------------------------
# Read file with time-series of B/Bmsy and then get them ready for analysis

data<-read.csv(file="../02_data_input/time_series_biomass_ratio_extrapolated.csv",sep=",",header=TRUE)

unique(data$SpStock) #check total number of population (52 pop)
data$SpStock_ModelRun<-as.factor(paste(data$SpStock,data$ModelRun))
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
# STEP2 loads the generation lengths for each species (and subpopulations)
#-------------------------------------------------------------------------
data$three_gl<-round(3*data$GL) #Calculate 3GL for each population
data$three_gl<-replace(data$three_gl, data$three_gl<10,10) #if 3GL is less than 10 years, then replace by 10 years following IUCN guidelines
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
# STEP3 loads the fishing mortality trends of the species (and subpopulations) for understanding whether a population is being sustainably fished (or not) and also whether a species, as a whole, is being sustainably fished (or not) according to IUCN Red List criteria
#-------------------------------------------------------------------------
# Read file with time-series of F/Fmsy and then get them ready for analysis

FF<-read.csv("../02_data_input/time_series_fishing_mortality_ratio_Fstatus_stock_and_species_level.csv")
FF$Ocean<-as.factor(FF$Ocean)
FF$Species<-as.factor(FF$Species)
FF$SpStock<-as.factor(FF$SpStock)
FF$Fstatus_meanFgl<-as.factor(FF$Fstatus_meanFgl) #yearly exploitation status of each population based on the average ratio of the fishing mortality relative to FMSY for the greater of one generation length or five years prior to the focal assessment year. This is used to determine if a population is sustainably managed in each focal assessment year in order to apply eaither A1 or A2 thresholds.
FF$Fstatus_species_meanFgl10<-as.factor(FF$Fstatus_species_meanFgl10) #yearly exploitation status of each species used to determined if a species is sustainably managed. Based on the IUCN Red List Criteria, a species needs to be sustainably managed in at least 90% of its range to justify the use of A1 thresholds, otherwise the A2 thresholds apply. 
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
# STEP4 estimates the total percent population reduction over three generation lengths (or 10 years, whichever is longer) for each subpopulation
# AND
# STEP5 applies Criterion A population reduction A1 or A2 thresholds for the determination of the IUCN Red List categories for each subpopulation
#-------------------------------------------------------------------------

# We use the following three-step estimation procedure to calculate the total percent change in biomass over three generation lengths for each focal assessment year. 
# 1. Conversion of the raw time-series of biomass of each population (and their multiple model runs) to successive annual rates of change
# 2. Estimation, for each population and focal assessment year, of the average of the annual rates of change across all the years over the prior three generation lengths using an intercept-only hierarchical Bayesian model. For  populations having multiple population biomass trends from several population assessment models (and model runs), we used the following two-level intercept-only hierarchical model (stan model file:modelRATE_multipleruns_hiearchical.stan) to estimate the average annual rate of change across all years i (at level 1) over three generation lengths window and nested within different model runs (at level 2). For  populations having only one or two time-series of biomass from one/two population assessment model and run, we used a one level intercept-only model to estimate the average annual rate of change across all years i (at level 1) over the three generations length window (stan model file: model_rate_onetuna_pooled.stan). 
# 3. Calculation of the total percent change in biomass over three generation lengths 

# Prepare the code to calculate annual rates of change
data<-data[!is.na(data$Value),] #Remove NaNs of Value Columns
data<-arrange(data,SpStock,SpStock_ModelRun,Year) #Rearrange time-series in order by Year 
data$Value.log<-log(data$Value) #log transform the biomass ratio values
data_diff<-ddply(data,.(RFMO, Ocean,Species,Stock,SpStock,ModelRun,SpStock_ModelRun,GL,three_gl),summarize,Year=head(Year,-1)+1,d1=diff(Value.log,1)) # Use the difference to get the annual rate of change


# For each population, estimates the average annual rate of change over the prior 3GL yearly, then calculates the total percent change over the 3GL yearly, and then A1 or A2 thresholds are applied for assigning a Red List Category 
data_diff$stockid<-as.numeric(factor(data_diff$SpStock)) #prepare population ID for following loop
data_diff<-arrange(data_diff, stockid,SpStock_ModelRun,Year)


for(spp in 1:length(unique(data_diff$stockid))){  
  table_save2<-data.frame(ocean=c(),sp=c(), st=c(), Fstatus_species_meanFgl10=c(), yearlast=c(), threat_category=c(), threat_code=c()) #data_frame where all the outputs for each population individually will be saved
  mystock<-subset(data_diff,data_diff$stockid==spp)
  n_regressions<-length(unique(mystock$Year))-unique( mystock$three_gl)+1 #For each stock, calculate the number of times the extent of change will be calculated across the the whole time-series based on the length of the time-series and the GL of the stock
  GLX3<-unique(mystock$three_gl)
  
  ## Calculates (1) average annual rate of change over 3GL (2) total percent change over 3GL, and (3) assigns Red List Category for each time step (yearly)
  for(i in 1:n_regressions){
    # identify the first year and the last year of the time window to estimate extent of change
    firstyear<-sort(unique(mystock$Year))[i]
    indice<-i+GLX3-1
    lastyear<-sort(unique(mystock$Year))[indice]
    
    # data within time window
    data_subset_by_y<-subset(mystock,mystock$Year>=firstyear & mystock$Year<=lastyear) 
    
    #prepare data to run the Stan model ( intercept-only hierarchical Bayesian model)
    data_subset_by_y$ModelRunId<-as.numeric(factor(data_subset_by_y$ModelRun))
    numberruns<-length(unique(as.numeric(factor(data_subset_by_y$ModelRun))))
    data.list<-list(N=length(data_subset_by_y$Year), dbio=data_subset_by_y$d1,N_group=length(unique(as.numeric(factor(data_subset_by_y$ModelRun)))),modelruns=data_subset_by_y$ModelRunId) #data input for the STAN model
    
    #Calculate number of model runs/scenarios within each stock
    number_modelruns<-length(unique(as.numeric(factor(data_subset_by_y$ModelRun))))
    
    #if number of model runs/scenarios is only 1 or 2, use a one level intercept-only model "model_rate_one_tuna_pooled.stan"
    if(number_modelruns<=3){
      model1 <- stan_model("model_rate_onetuna_pooled.stan") #single run model 
      model_output <- sampling(model1, data = data.list,chains=3, iter=6000, warmup=3000, thin=10,control = list(adapt_delta = 0.8),cores=6)#Runs smaller number of iterations to decrease the calculation time
      posterior <- as.array(model_output)
      list_of_draws <- extract(model_output)
      beta<-list_of_draws$beta  #model-average annual rates of change across all years in the time window nested in all model runs, which is used in step 3 to estimate the total extent of change within 3GL
      
    }else{
      #if number of model runs/scenarios is > 3, use a two-level intercept-only hierarchical model "modelRATE_multipleruns_hiearchical.stan"
      model1 <- stan_model("modelRATE_multipleruns_hiearchical.stan") #model with multiple runs, so it needs hiearchical model
      model_output <- sampling(model1, data = data.list,chains=3, iter=6000, warmup=3000, thin=10,control = list(adapt_delta = 0.8),cores=6)#Run smaller number of iterations to decrease the calculation time
      posterior <- as.array(model_output)
      list_of_draws <- extract(model_output)
      beta<-list_of_draws$mu_a1 #model-average annual rates of change across all years in the time window nested in all model runs, which is used in step 3 to estimate the total extent of change within 3GL
    }
    
    #calculating the total percent population change over three generation lengths 
    AverageAnnualRC<-round((exp(beta)-1)*100,3)
    tpercent_change3GL<-round((exp(beta*GLX3)-1)*100,3)
    tpercent_change3GL<-as.numeric(tpercent_change3GL)
    
    #Prepare data frame to store all the outputs
    ocean<-unique(data_subset_by_y$Ocean)
    sp<-unique(data_subset_by_y$Species)
    st<-unique(data_subset_by_y$SpStock)
    firstyear<-unique(data_subset_by_y$Year)[1]
    lastyear<-unique(data_subset_by_y$Year)[GLX3]
    
    #Choose A1 or A2 thresholds based whether a population is being sustainably managed or not (based on the yearly exploitation status of each population)
    Fextract<-FF %>% filter (SpStock==st & Year==lastyear)
    Fstatus_meanFgl<-Fextract$Fstatus_meanFgl #population level Fstatus_meanFgl
    Fstatus_species_meanFgl10<-Fextract$Fstatus_species_meanFgl10 #species level Fstatus_meanFgl10
    
    if(length(Fstatus_meanFgl)==0){
      Fstatus_meanFgl<-c("Not Overfishing")
      Fstatus_species_meanFgl10<-c("Not Overfishing")
    }
    
    if(Fstatus_meanFgl=="Overfishing"){ #apply A2 thesholds and assing Red List status category
      threat_category<-cut(tpercent_change3GL, c(-110,-80,-50,-30,-20,Inf),labels = c("Critically Endangered", "Endangered", "Vulnerable","Near Threatened","Least Concern")) #allocates the category that corresponds to the percent change estimated
      threat_code<-cut(tpercent_change3GL, c(-110,-80,-50,-30,-20,Inf),labels = c(4,3,2,1,0))
    }else if(Fstatus_meanFgl=="Not Overfishing"){ 
      #apply A1 thesholds and assing Red List status category
      threat_category<-cut(tpercent_change3GL, c(-110,-90,-70,-50,-40,Inf),labels = c("Critically Endangered", "Endangered", "Vulnerable","Near Threatened","Least Concern")) #allocates the category that corresponds to the percent change estimated
      threat_code<-cut(tpercent_change3GL, c(-110,-90,-70,-50,-40,Inf),labels = c(4,3,2,1,0))
    }else{
      threat_category<-cut(tpercent_change3GL, c(-110,-90,-70,-50,-40,Inf),labels = c("Critically Endangered", "Endangered", "Vulnerable","Near Threatened","Least Concern")) #allocates the category that corresponds to the percent change estimated
      threat_code<-cut(tpercent_change3GL, c(-110,-90,-70,-50,-40,Inf),labels = c(4,3,2,1,0))
    }
    
    #Prepare data to save in data frame
    ocean<-rep(ocean,length(tpercent_change3GL))
    sp=rep(sp,length(tpercent_change3GL))
    st=rep(st,length(tpercent_change3GL))
    Fstatus_meanFgl=rep(Fstatus_meanFgl,length(tpercent_change3GL))
    Fstatus_species_meanFgl10=rep(Fstatus_species_meanFgl10,length(tpercent_change3GL))
    yearlast=rep(lastyear,length(tpercent_change3GL))
    threat_category=threat_category
    threat_code=threat_code
    table_results<-data.frame(ocean,sp,st,Fstatus_meanFgl,Fstatus_species_meanFgl10,yearlast,threat_category,threat_code)
    table_save<-table_results
    table_save2<-rbind(table_save2,table_save) #save individual file ouputs for each stock
  }
  
  #For the calculation of the Red List Index, all the populations need to have a yearly Red List status category between 1950 and 2019. Therefore, we assigned a Least Concern category for these early years without having a Red List status Category assigned when needed.
  minyear<-min(table_save2$yearlast)
  if(minyear>1950){
    minyear<-min(table_save2$yearlast)-1
    yearlast_series<-c(1950:minyear)
    n_runs<-dim(table_save)[1]
    ocean<-rep(unique(table_save$ocean),length(yearlast_series)*n_runs)
    sp<-rep(unique(table_save$sp),length(yearlast_series)*n_runs)
    st<-rep(unique(table_save$st),length(yearlast_series)*n_runs)
    Fstatus_meanFgl<-gl(1,length(yearlast_series)*n_runs,label="Not Overfishing")
    Fstatus_species_meanFgl10<-gl(1,length(yearlast_series)*n_runs,label=unique(Fstatus_species_meanFgl10))
    yearfirst<-rep(NA,length(yearlast_series)*n_runs) 
    yearlast<-rep(yearlast_series,times=n_runs)
    AverageAnnualRC<-rep(0,length(yearlast_series)*n_runs) 
    tpercent_change3GL<-rep(0,length(yearlast_series)*n_runs) 
    threat_category<-gl(1,length(yearlast_series)*n_runs,label="Least Concern")
    threat_code<-gl(1,length(yearlast_series)*n_runs,label="0")
    newyears<-data.frame(ocean,sp,st,Fstatus_meanFgl, Fstatus_species_meanFgl10,yearlast,threat_category,threat_code)
    table_save2<-rbind(newyears,table_save2)
  }else{
    table_save2<-table_save2
  }
  
  #Outputs for each stock in individual files
  filename<-paste("Model_output_",unique(data_subset_by_y$SpStock),".csv", sep="") 
  path_output_stocks <- "04_Output" 
  
  write.csv(table_save2,paste0(path_output_stocks, filename)) 
} #end main loop



#-----------------------------------------------------
# STEP6 calculates the Red List Index (RLI) of populations (including 52 populations of tunas, billfishes and sharks)
#-----------------------------------------------------
#The file "Yearly_RL_status_populations.RData" contains all the model outputs with the Red List status estimates for all 52 populations from 1950 to 2019.

 load(file="../04_Output/Yearly_RL_status_populations.RData")
#Prepare the data for the Red List calculation
names(threat_data)
length(unique(threat_data$st)) #Red List status for 52 populations
length(unique(threat_data$sp)) #Red List status for 18 species
threat_data$st<-as.factor(threat_data$st)
threat_data$threat_category<-as.factor(threat_data$threat_category)
unique(threat_data$threat_category)
levels(threat_data$threat_category)
threat_data$threat_category<-factor(threat_data$threat_category,levels=c("Critically Endangered", "Endangered","Vulnerable","Near Threatened","Least Concern"))

#Use function RLI_population.R to calculate the population Red List Index
source("RLI_population_function.R")  
RLI_global<-RLI_population(threat_data)
summary(RLI_global)

#Calculate mean and percentiles and plot
require(plyr)

RLIq<-ddply(RLI_global,.(years),summarise, X0 = quantile(RLI, probs = 0),X2.5 = quantile(RLI, probs = 0.025),X10 = quantile(RLI, probs = 0.10),X20 = quantile(RLI, probs = 0.20), X25 = quantile(RLI, probs = 0.25), X50 = quantile(RLI, probs = 0.50), Median=median(RLI),X75 = quantile(RLI, probs = 0.75),X80 = quantile(RLI, probs = 0.80), X90 = quantile(RLI, probs = 0.90), X97.5 = quantile(RLI, probs = 0.975) ,X100 = quantile(RLI, probs = 1))

require(ggplot2)

ggplot(data=RLIq, aes(x=years, y=Median)) + geom_line(aes(x=years, y=Median)) +geom_ribbon(data=RLIq, aes(x=years, ymin=X2.5, ymax=X97.5),fill="gray30", alpha=0.2) + ggtitle("Population-level Red List Index") +  xlab("Years") + ylab("Median and 50th and 95th percentiles") 

#-----------------------------------------------------
#STEP7 calculates the Red List Index of species (including 18 species of tunas, billfishes and sharks)
#-----------------------------------------------------

#Load the data
#The file "Yearly_RL_status_populations_with_weights.RData" contains all the model outputs with the red list status estimates for all 52 populations between 1950 to 2019 with information about the weights of each population to inform the species-level calculation of the red list index.
load(file="../04_Output/Yearly_RL_status_populations_with_weights.RData")

summary(threat_data_with_weights)

#Use function RLI_species_function.R to calculate the species level Red List Index
#First, it estimates the total percent population reduction over three generation lengths (or 10 years, whichever is longer) for each species by combining the information for all available subpopulations. For those species with multiple populations, we weighted the estimated total percent change in biomass of each population by their maximum sustainable yield to account for the contribution of different population sizes to the global species
#Second, it applies Criterion A population reductions A1 or A2 thresholds for the determination of the IUCN Red List categories for each species
#Third, it calculates the species-level Red List Index

source("RLI_species_function.R")
RLI_global<-RLI_species(threat_data_with_weights) #This output contains two lists
summary(RLI_global)

# LIST (1) Extract threat status at the species level
global_data_speciesW<-RLI_global[[1]]
head(global_data_speciesW)


# LIST (2) Extract Red List Index estimation
RLI_species_global_Fgl10<-RLI_global[[2]]
str(RLI_species_global_Fgl10)

#Calculate mean and percentiles
head(RLI_species_global_Fgl10)
summary(RLI_species_global_Fgl10)

library(plyr)
RLI_species_Fgl10<-ddply(RLI_species_global_Fgl10,.(years),summarise, X0 = quantile(RLIFgl10, probs = 0),X2.5 = quantile(RLIFgl10, probs = 0.025),X10 = quantile(RLIFgl10, probs = 0.10),X20 = quantile(RLIFgl10, probs = 0.20), X25 = quantile(RLIFgl10, probs = 0.25), X50 = quantile(RLIFgl10, probs = 0.50), Median=median(RLIFgl10),X75 = quantile(RLIFgl10, probs = 0.75),X80 = quantile(RLIFgl10, probs = 0.80), X90 = quantile(RLIFgl10, probs = 0.90), X97.5 = quantile(RLIFgl10, probs = 0.975) ,X100 = quantile(RLIFgl10, probs = 1))


head(RLI_species_Fgl10)
summary(RLI_species_Fgl10)

RLI_species_Fgl10<- RLI_species_Fgl10 %>% filter(years>=1950)


ggplot(data=RLI_species_Fgl10, aes(x=years, y=Median)) + geom_line(aes(x=years, y=Median)) +geom_ribbon(data=RLI_species_Fgl10, aes(x=years, ymin=X2.5, ymax=X97.5),fill="gray30", alpha=0.2) + ggtitle("RLI SPECIES GLOBAL") +  xlab("Years") + ylab("Median and 50th and 100th percentiles")