#--------------------------------------------------------
#--------------------------------------------------------
# R scripts for calculation the Red List index of stock and species for oceanic tunas, billfishes and sharks
# Juan Jorda et al xxxxxxxx ##NP: + add email of corresponding author. -> instead of all author, for script, you usually just write the developper(s) of the code
#--------------------------------------------------------
#--------------------------------------------------------
# This script calculates a continuous Red List Index for oceanic predatory fishes between 1950 and 2019 using the estimated extinction risk of the 18 species of tunas, billfishes, and sharks. We applied the IUCN Criteria A "population reduction" to calculate the extinction risk (Red List Threat Status) for 18 species of tunas, billfishes, and sharks. Criterion A measures the extinction risk based on exceeding thresholds of population decline over a time frame of three generation lengths (or a minimum of 10 years, whichever is longer).  ##NP: no need here but keep it for readme file

# The R code below calculates: 
# 1. The Red List threat status the 52 populations of tunas, billfishes and sharks which are then used to calculate the population-level Red List Index .
# 2. The Red list threat status for the 18 species of tunas, billfishes and sharks which then are used to calculate the species-level Red List Index. 

# These steps are followed below:
# STEP1 reads the population trends of the species (and subpopulations) ##NP: replace 'reads'by 'load' everywhere in the doc
# 
# STEP2 reads the generation length for the species (and subpopulations)
# 
# STEP3 reads the fishing mortality trends of the species (and subpopulations) for understanding of whether a population is sustainably fished and also whether a species, as a whole, is sustainably fished
# 
# STEP4 estimates the total percent population reduction over three generation lengths (or 10 years, whichever is longer) for each subpopulation
# 
# STEP5  applies Criterion A population reductions A1 or A2 thresholds for the determination of the IUCN Red List categories for each subpopulation
# 
# STEP6 calculates the Red List Index of populations (including 52 populations of tunas, billfishes and sharks)
# 
# STEP7 estimates the total percent population reduction over three generation lengths (or 10 years, whichever is longer) for each species by combining the information for all available subpopulations. For those species with multiple populations, we weighted the estimated total percent change in biomass of each population by their maximum sustainable yield to account for the contribution of different population sizes to the global species
# 
# STEP8 applies Criterion A population reductions A1 or A2 thresholds for the determination of the IUCN Red List categories for each species
#
# STEP9 calculates the Red List Index of species (including 18 species of tunas, billfishes and sharks)
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

# set the working directory to the folder containing this script: ##NP: delete this line and the two below, people will do this by themself
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # it sucks this does not work. So I need to set my WD below: ##NP: delete, people will do this by themself
setwd("/Users/mjuanjorda/Dropbox/POSTDOC/MarieCurieProject/Redlistindex/RLI_january2015/RLI_2020/Data_code_sharing/Ranalysis") ##NP: delete, people will do this by themself

## replace with this
## http://jenrichmond.rbind.io/post/how-to-use-the-here-package/
## https://github.com/jennybc/here_here

here::here()
# in my case this looks like this : [1] "/Users/nickdulvy/Dropbox/Data_code_sharing"

#-------------------------------------------------------------------------
# STEP1 reads the population trends of the species (and subpopulations)
#-------------------------------------------------------------------------
# Read file with time series of B/Bmsy and then get them ready for analysis

###-------------------------------------------------------
# generic data ingest using here()
datos <- read.csv(
  here("data_input", "time_series_biomass_ratio_extrapolated.csv"))
#-------------------------------------------------------------------------
datos<-read.csv(file="../data_input/time_series_biomass_ratio_extrapolated.csv",sep=",",header=TRUE)

unique(datos$SpStock) #check total number of population (52 pop)
datos$SpStock_ModelRun<-as.factor(paste(datos$SpStock,datos$ModelRun))
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
# STEP2 reads the generation length for the species (and subpopulations)
#-------------------------------------------------------------------------
datos$three_gl<-round(3*datos$GL) #Calculate 3xGL for each population
datos$three_gl<-replace(datos$three_gl, datos$three_gl<10,10) #if 3GL is less than 10 years, then replace by 10 years following IUCN guidelines
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
# STEP3 reads the fishing mortality trends of the species (and subpopulations) for understanding of whether a population is sustainably fished and also whether a species, as a whole, is sustainably fished
#-------------------------------------------------------------------------
# Read file with time series of F/Fmsy and then get them ready for analysis
FF <- read.csv(
  here("data_input", "time_series_fishing_mortality_ratio_Fstatus_stock_and_species_level.csv"))

## FF<-read.csv("../data_input/time_series_fishing_mortality_ratio_Fstatus_stock_and_species_level.csv")
FF$Ocean<-as.factor(FF$Ocean)
FF$Species<-as.factor(FF$Species)
FF$SpStock<-as.factor(FF$SpStock)
FF$Fstatus_meanFgl<-as.factor(FF$Fstatus_meanFgl) #This is the yearly exploitation status of each population based on the average ratio of the fishing mortality relative to FMSY for the greater of one generation length or five years prior to the focal assessment year. This is used to determine if a population is sustainably managed n each focal assessment year in order to apply eaither A1 or A2 thresholds.
FF$Fstatus_species_meanFgl10<-as.factor(FF$Fstatus_species_meanFgl10) #This is the yearly exploitation status of each species used to determined if a species is sustainably managed. Based on the IUCN Red List Categories and Criteria state, a species needs to be sustainably managed in at least 90% of its range to justify the use of A1 thresholds, otherwise the A2 thresholds apply. 
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
# STEP4 estimates the total percent population reduction over three generation lengths (or 10 years, whichever is longer) for each subpopulation
# AND
# STEP5 applies Criterion A population reductions A1 or A2 thresholds for the determination of the IUCN Red List categories for each subpopulation
#-------------------------------------------------------------------------

# We use the following three-step estimation procedure to calculate the total percent change in biomass over three generation lengths for each focal assessment year. 
# 1. we converted the raw time-series of biomass of each population (and their multiple model runs) to successive annual rates of change
# 2. we estimated for each population and focal assessment year the average of the annual rates of change across all the years over the prior three generation lengths using an intercept-only hierarchical Bayesian model. ##NP: I would delete everything after this ## The regression intercept-only model implicitly assumes that the true population follows a constant deterministic trend, and it assumes that all the deviations from this trend are attributable to statistically-independent errors. The hierarchical model allows to account for the nested structure of the data, as the majority of the populations have multiple time-series of population biomass derived from multiple assessment models (and multiple model runs). For those populations having multiple trends in population biomass from several population assessment models (and model runs), we used the following two-level intercept-only hierarchical model (file: xxxxxxx stan model) to estimate the average annual rate of change across all years i (at level 1) over three generation lengths window and nested within different model runs (at level 2). For those populations having only one or two time-series of biomass from one/two population assessment model and run, we used a one level intercept-only model (Equation 4) to estimate the average annual rate of change across all years i (at level 1) over the three generations length window (file: xxxxxxx stan model). 
# 3. we calculated the total % change in biomass over three generation lengths 
# comment : A1 thresholds are applied when a population is determined to be sustainably managed, and A2 thresholds are applied when a population is determined to be unsustainably managed  

# Set Working directory file, where outputs for each stock in individual files will be stored ##NP: to add
path_output_stocks <- "..." ##NP: to add

## Suggest using here()


# Prepare the code to calculate annual rates of change
datos<-datos[!is.na(datos$Value),] #Remove NaNs of Value Columns
datos<-arrange(datos,SpStock,SpStock_ModelRun,Year) #make sure the time series are in order by Year 
datos$Value.log<-log(datos$Value) #log transform the biomass ratio values
datos_diff<-ddply(datos,.(RFMO, Ocean,Species,Stock,SpStock,ModelRun,SpStock_ModelRun,GL,three_gl),summarize,Year=head(Year,-1)+1,d1=diff(Value.log,1)) # take the difference to get the annual rate of change


# For each population, estimates the average annual rate of change over the prior 3GL yearly, then calculates the total percent change over the 3GL yearly, and then A1 or A2 thresholds are applied for assigning a Red List Category 

datos_diff$stockid<-as.numeric(factor(datos_diff$SpStock)) ## prepare population ID for the loop
datos_diff<-arrange(datos_diff, stockid,SpStock_ModelRun,Year)


for(spp in 1:length(unique(datos_diff$stockid))){  
  table_save2<-data.frame(ocean=c(),sp=c(), st=c(), Fstatus_species_meanFgl10=c(), yearlast=c(), threat_category=c(), threat_code=c()) #data_frame where all the outputs for each population individually will be saved
  mystock<-subset(datos_diff,datos_diff$stockid==spp)
  n_regressions<-length(unique(mystock$Year))-unique( mystock$three_gl)+1     #For each stock, calculate the number of times the extent of change will be calculated across the the whole time series based on the length of the time series and the GL of the stock
  GLX3<-unique(mystock$three_gl)
  
  ## Calculates (1) average annual rate of change over 3GL (2) total percent change over 3GL and (3) assigns Red List Category for each time step (yearly)
  for(i in 1:n_regressions){
    # identify the first year and the last year of the time window to estimate extent of change
    firstyear<-sort(unique(mystock$Year))[i]
    indice<-i+GLX3-1
    lastyear<-sort(unique(mystock$Year))[indice]
    
    # data within time window
    data_subset_by_y<-subset(mystock,mystock$Year>=firstyear & mystock$Year<=lastyear) 
    
    #prepare data to run the Stan model (the intercept-only hierarchical Bayesian model)
    data_subset_by_y$ModelRunId<-as.numeric(factor(data_subset_by_y$ModelRun))
    numberruns<-length(unique(as.numeric(factor(data_subset_by_y$ModelRun))))
    data.list<-list(N=length(data_subset_by_y$Year), dbio=data_subset_by_y$d1,N_group=length(unique(as.numeric(factor(data_subset_by_y$ModelRun)))),modelruns=data_subset_by_y$ModelRunId) # this is the data input for the STAN model.
    
    #Calculate number of model runs/scenarios within each stock
    number_modelruns<-length(unique(as.numeric(factor(data_subset_by_y$ModelRun))))
    
    #if number of model runs/scenarios is only 1 or 2, use a one level intercept-only model "model_rate_one_tuna_pooled.stan"
    if(number_modelruns<=3){
      model1 <- stan_model("model_rate_onetuna_pooled.stan") # single run model 
      #model_output <- sampling(model1, data = data.list,chains=3, iter=12000, warmup=3000, thin=10,control = list(adapt_delta = 0.8),cores=6)  ##NP: to delete? keep this one or the one below, preferably the one you use for the real analysis
      # model_output <- sampling(model1, data = data.list,chains=3, iter=6000, warmup=3000, thin=10,control = list(adapt_delta = 0.8),cores=6)#Run smaller number of iterations to decrease the calculation time
      model_output <- sampling(model1, data = data.list,chains=1, iter=600, warmup=300, thin=1)#Run smaller number of iterations to decrease the calculation time
      
      posterior <- as.array(model_output)
      list_of_draws <- extract(model_output)
      beta<-list_of_draws$beta  #this parameter is the model-average annual rates of change across all years in the time window nested in all model runs, which we use in step 3 to estimate the total extent of change within 3GL
      
    }else{
      #if number of model runs/scenarios is > 3, use a two-level intercept-only hierarchical model "modelRATE_multipleruns_hiearchical.stan"
      model1 <- stan_model("modelRATE_multipleruns_hiearchical.stan") #model with multiple runs, so it needs hiearchical model
      # model_output <- sampling(model1, data = data.list,chains=3, iter=6000, warmup=3000, thin=10,control = list(adapt_delta = 0.8),cores=6) #Run smaller number of iterations to decrease the calculation time
      model_output <- sampling(model1, data = data.list,chains=1, iter=600, warmup=300, thin=1) #Run smaller number of iterations to decrease the calculation time
      posterior <- as.array(model_output)
      list_of_draws <- extract(model_output)
      beta<-list_of_draws$mu_a1 #this parameter is the model-average annual rates of change across all years in the time window nested in all model runs, which we use in step 3 to estimate the total extent of change within 3GL
    }
    
    #calculating the total % population change over three generation lengths 
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
    Fstatus_meanFgl<-Fextract$Fstatus_meanFgl #this is population level Fstatus_meanFgl
    Fstatus_species_meanFgl10<-Fextract$Fstatus_species_meanFgl10 #this is the species level Fstatus_meanFgl10
    
    if(length(Fstatus_meanFgl)==0){
      Fstatus_meanFgl<-c("Not Overfishing")
      Fstatus_species_meanFgl10<-c("Not Overfishing")
    }
    
    if(Fstatus_meanFgl=="Overfishing"){ #apply A2 thesholds and assing Red List status category
      threat_category<-cut(tpercent_change3GL, c(-110,-80,-50,-30,-20,Inf),labels = c("Critically Endangered", "Endangered", "Vulnerable","Near Threatened","Least Concern")) ### allocates the category that correspond to the percent change estimated
      threat_code<-cut(tpercent_change3GL, c(-110,-80,-50,-30,-20,Inf),labels = c(4,3,2,1,0))
    }else if(Fstatus_meanFgl=="Not Overfishing"){ 
      #apply A1 thesholds and assing Red List status category
      threat_category<-cut(tpercent_change3GL, c(-110,-90,-70,-50,-40,Inf),labels = c("Critically Endangered", "Endangered", "Vulnerable","Near Threatened","Least Concern")) ### allocates the category that correspond to the percent change estimated
      threat_code<-cut(tpercent_change3GL, c(-110,-90,-70,-50,-40,Inf),labels = c(4,3,2,1,0))
    }else{
      threat_category<-cut(tpercent_change3GL, c(-110,-90,-70,-50,-40,Inf),labels = c("Critically Endangered", "Endangered", "Vulnerable","Near Threatened","Least Concern")) ### allocates the category that correspond to the percent change estimated
      threat_code<-cut(tpercent_change3GL, c(-110,-90,-70,-50,-40,Inf),labels = c(4,3,2,1,0))
    }
    
    #Prepare data to save in data.frame
    ocean<-rep(ocean,length(tpercent_change3GL))
    sp=rep(sp,length(tpercent_change3GL))
    st=rep(st,length(tpercent_change3GL))
    Fstatus_meanFgl=rep(Fstatus_meanFgl,length(tpercent_change3GL))
    Fstatus_species_meanFgl10=rep(Fstatus_species_meanFgl10,length(tpercent_change3GL))
    #yearfirst=rep(firstyear,length(tpercent_change3GL)) ##NP: to delete?
    yearlast=rep(lastyear,length(tpercent_change3GL))
    #AverageAnnualRC=AverageAnnualRC ##NP: to delete?
    #tpercent_change3GL=tpercent_change3GL ##NP: to delete?
    threat_category=threat_category
    threat_code=threat_code
    table_results<-data.frame(ocean,sp,st,Fstatus_meanFgl,Fstatus_species_meanFgl10,yearlast,threat_category,threat_code)
    table_save<-table_results
    table_save2<-rbind(table_save2,table_save) # save individual file ouputs for each stock
  }
  
  #for some populations the first year of the starts after 1950s with a threat status of Least concern. For the calculation of the Red List Index, all the populations need to have a yearly Red List status category between 1950 and 2019. Therefore, we assigned a Least Concern category for these early year without having a Red List status Category assigned when needed.
  minanno<-min(table_save2$yearlast)
  if(minanno>1950){
    minyear<-min(table_save2$yearlast)-1
    yearlast_series<-c(1950:minyear)
    n_runs<-dim(table_save)[1]
    ocean<-rep(unique(table_save$ocean),length(yearlast_series)*n_runs)
    sp<-rep(unique(table_save$sp),length(yearlast_series)*n_runs)
    st<-rep(unique(table_save$st),length(yearlast_series)*n_runs)
    Fstatus_meanFgl<-gl(1,length(yearlast_series)*n_runs,label="Not Overfishing")
    Fstatus_species_meanFgl10<-gl(1,length(yearlast_series)*n_runs,label=unique(Fstatus_species_meanFgl10))
    #yearfirst<-rep(NA,length(yearlast_series)*n_runs) ##NP: to delete?
    yearlast<-rep(yearlast_series,times=n_runs)
    #AverageAnnualRC<-rep(0,length(yearlast_series)*n_runs) ##NP: to delete?
    #tpercent_change3GL<-rep(0,length(yearlast_series)*n_runs) ##NP: to delete?
    threat_category<-gl(1,length(yearlast_series)*n_runs,label="Least Concern")
    threat_code<-gl(1,length(yearlast_series)*n_runs,label="0")
    newyears<-data.frame(ocean,sp,st,Fstatus_meanFgl, Fstatus_species_meanFgl10,yearlast,threat_category,threat_code)
    table_save2<-rbind(newyears,table_save2)
  }else{
    table_save2<-table_save2
  }
  
  ### outputs for each stock in individual files
  filename<-paste("Model_output_",unique(data_subset_by_y$SpStock),".csv", sep="") 
  path <- "/Users/mjuanjorda/Dropbox/POSTDOC/MarieCurieProject/Redlistindex/RLI_january2015/RLI_2020/Data_code_sharing/output"  ##NP: to delete
  write.csv(table_save2,file.path(path, filename))  ##NP: to delete
  write.csv(table_save2,paste0(path_output_stocks, filename))  ##NP: to add
} #end main loop


##-------------------------------------------------------------------  ##NP: to delete
## Example plot: plot the  Red List status over time for one population

setwd("/Users/mjuanjorda/Dropbox/POSTDOC/MarieCurieProject/Redlistindex/RLI_january2015/RLI_2020/Data_code_sharing/Ranalysis") #BORRAR AL FINAL ##NP: to delete

example<-read.csv(file="../output/Model_output_FAL Western Pacific.csv",sep=",",header=TRUE)
summary(example)

detach("package:plyr", unload=TRUE)
library(dplyr)
example2 <- example %>% 
  group_by(sp,st,yearlast,threat_category) %>% 
  summarise(count=n()) %>% 
  mutate(perc=count/sum(count))

example2$threat_category<-factor(example2$threat_category,levels=c("Critically Endangered", "Endangered","Vulnerable","Near Threatened","Least Concern"))

example2<-as.data.frame(example2)
example2$perc<-as.numeric(example2) ##NP: not working


#dev.new()
require(ggplot2)
ggplot(example2, aes(x = factor(yearlast), y = perc*100, fill = factor(threat_category))) + geom_bar(stat="identity", width = 0.7)+ labs(x = "Year", y = "Probability of being\nclassified in Red\nList Categories", fill = "threat_category") + scale_fill_manual(name="Threat category",values = c("Critically Endangered"="#D81E05","Endangered"="#FC7F3F","Vulnerable"="#F9E814","Near Threatened"="#CCE226","Least Concern"="#60C659","No evaluation"="#FFFFFF"))+ ggtitle("Threat status over time") + scale_x_discrete(breaks = seq(1930,2020,10)) +  guides(fill = guide_legend(override.aes = list(colour = "black")))+theme(legend.position="bottom")+ facet_wrap(~ st,ncol=1)

##-------------------------------------------------------------  ##NP: to delete
## The output folder contains yearly RL status for each population. Read and combine all the population into one data frame in preparation of the RLI calculation  ##NP: move lines 260 to 270 in the next step

require(readbulk)
threat_data <- read_bulk(directory = "../output", extension = ".csv") #read and combines all the individual files into a dataframe

threat_data$sp<-as.factor(threat_data$sp)
levels(threat_data$sp)
threat_data$taxo<-threat_data$sp
levels(threat_data$taxo)<-c("Tunas","Tunas","Tunas","Billfishes","Sharks","Billfishes","Sharks","Billfishes","Sharks","Tunas","Sharks","Billfishes","Tunas","Tunas","Sharks","Billfishes","Billfishes","Tunas")

save(threat_data,file="../output/Yearly_RL_status_populations.RData")
##-------------------------------------------------------------  ##NP: to add


#-----------------------------------------------------
# STEP6 calculates the Red List Index of populations (including 52 populations of tunas, billfishes and sharks)
#-----------------------------------------------------
# load(file="../output/Yearly_RL_status_populations.RData")  ##NP: to delete?
# Prepare the data for the RL calculation
length(unique(threat_data$st)) # RL status for 52 populations
length(unique(threat_data$sp)) # RL status for 18 species
threat_data$st<-as.factor(threat_data$st)
threat_data$threat_category<-as.factor(threat_data$threat_category)
unique(threat_data$threat_category)
levels(threat_data$threat_category)
threat_data$threat_category<-factor(threat_data$threat_category,levels=c("Critically Endangered", "Endangered","Vulnerable","Near Threatened","Least Concern"))

# This function (RLI_population.R) calculates de RLI
summary(threat_data)
path_script_RLI_population <- '...'   ##NP: to add
source("/Users/mjuanjorda/Dropbox/POSTDOC/MarieCurieProject/Redlistindex/RLI_january2015/RLI_2020/Data_code_sharing/Ranalysis/RLI_population.R")  ##NP: to delete?
source(paste0(path_script_RLI_population, "RLI_population.R"))  ##NP: to add
RLI_global<-RLI_population(threat_data)
summary(RLI_global)
dim(RLI_global)  ##NP: to delete?

#Calculate mean and percentiles and plot
require(plyr)

RLIq<-ddply(RLI_global,.(years),summarise, X0 = quantile(RLI, probs = 0),X2.5 = quantile(RLI, probs = 0.025),X10 = quantile(RLI, probs = 0.10),X20 = quantile(RLI, probs = 0.20), X25 = quantile(RLI, probs = 0.25), X50 = quantile(RLI, probs = 0.50), Median=median(RLI),X75 = quantile(RLI, probs = 0.75),X80 = quantile(RLI, probs = 0.80), X90 = quantile(RLI, probs = 0.90), X97.5 = quantile(RLI, probs = 0.975) ,X100 = quantile(RLI, probs = 1))

require(ggplot2)

ggplot(data=RLIq, aes(x=years, y=Median)) + geom_line(aes(x=years, y=Median)) +geom_ribbon(data=RLIq, aes(x=years, ymin=X2.5, ymax=X97.5),fill="gray30", alpha=0.2) + ggtitle("Population-level Red List Index") +  xlab("Years") + ylab("Median and 50th and 95th percentiles") 








