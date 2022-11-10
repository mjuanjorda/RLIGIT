# Seventy years of tuna, billfish and sharks as sentinels of global ocean health


This Red List Index analysis published in Science (Juan-Jordá et al. 2022) develops a robust global indicator (the continuous Red List Index of oceanic predatory fishes) for monitoring ocean health in oceanic marine ecosystems & tracking global sustainability & biodiversity targets.

These files contain the main data sources and RCode for conducting the global Red List Index calculation.

Data_raw: 
This folder contains the raw time series of biomass and fishing mortalty rates for 52 populations of tunas, billfishes and sharls. We compiled the most recent (as of June 2020) fish stock assessments for 52 populations (18 species) of tunas, billfishes, and sharks from the five tuna-Regional Fisheries Management Organizations (RFMOs) (Table 1). For each fish stock assessment, we extracted the following: (a) the estimated time-series of biomass or time-series of biomass relative to the biomass that produces the Maximum Sustainable Yield (MSY) [B/BMSY], (b) the estimated time-series of fishing mortality relative to the fishing mortality that produces the MSY [F/FMSY], and (c) the standard biological reference points used to determine population status, generally the current adult biomass relative to the adult biomass producing MSY (Bcurrent/BMSY) and current fishing mortality rate relative to the fishing mortality that maintains MSY (Fcurrent/FMSY). This data was extracted from the assessment models (and model runs) used to determine population status and provide management advice by the Scientific Committees of each of the tuna-RFMOs.

Data_input: This folder contains the time series of biomass for the 52 populations (18 species) of tunas, billfishes, and sharks extrapolated backward and forward to ensure that Criterion A can be applied to assign a Red List category to each popualtion and year between 1950 and 2009. In order to assess a population against Criterion A for a given focal assessment year, a time-series of biomass with a length corresponding to three times the generation length (or 10 years, whichever is longer) is needed prior to the focal assessment year. We found that in some cases the time-series of biomass were not long enough to calculate the total percentage change in biomass over three generations lengths in yearly time steps from 1950 to 2019. Therefore, we had to extrapolate backwards and forward some time-series of biomass as needed. We extrapolated the time-series of biomass backwards in time until sufficient number of years were available to apply Criteria A starting in 1950 for every population as need it. The length of the backward extrapolation depended on the generation length of the population, with populations having larger generation length usually requiring a larger number of years of extrapolation.

Ranalysis: 
This folder contains the Rcode for calculating the extinction risk for 18 species of tunas, billfishes and sharks based on the IUCN Red List Categories and Criteria.  All species of oceanic tunas, billfishes, and sharks were assessed under IUCN Red List Criterion A “population reduction”. Criterion A was applied to both, the taxonomic unit of population and the taxonomic unit of species, to assign a Red List Category to each population and species of tunas, billfishes, and sharks between 1950 and 2019. 

The folder also contains the R code for calculating the global Red List Index of oceanic predatory fishes. We calculated a continuous Red List Index for oceanic predatory fishes between 1950 and 2019 using the estimated extinction risk of the 18 species of tunas, billfishes, and sharks. 

Output:This folder contains the estimated Red List categories between 1950 and 2009 for each population of tuna, billfish and shark, which are used to estimate calculate the global Red List Index of oceanic predatory species.


MAIN SOURCE
Cite this article as M.J. Juan-Jordá et al., Science 378, eabj0211 (2022).
DOI: 10.1126/science.abj0211
READ THE FULL ARTICLE AT
https://doi.org/10.1126/science.abj0211




