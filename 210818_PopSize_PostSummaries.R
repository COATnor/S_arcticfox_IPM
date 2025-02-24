#############################################################
#### INTEGRATED POPULATION MODEL FOR SVALBARD ARCTIC FOX ####
#############################################################

library(coda)
library(reshape)
library(ggplot2)
library(viridis)
library(plyr)
library(gridExtra)

## Load workspace containing all data
load('200228_AF_IPM_Data.RData')

## Load posterior samples
load('200429_AF_IPM_VersionB4.RData')

## Re-arrange data
out.mat <- as.matrix(AF.IPM.varB4)

## Make an empty data frame
N.data <- data.frame(
  Year = c(1997:2019),
  localNtot_mean = NA,
  localNtot_median = NA,
  localNtot_lCI = NA,
  localNtot_uCI = NA,
  Ntot_mean = NA,
  Ntot_median = NA,
  Ntot_lCI = NA,
  Ntot_uCI = NA
)

## Make posterior summaries for local and total population size in each year
for(t in 1:23){
  localNtot_t <- out.mat[,paste0('N.tot[', t, ']')] - out.mat[,paste0('Imm[', t, ']')]
  N.data$localNtot_mean[t] <- mean(localNtot_t) 
  N.data[t,c('localNtot_lCI', 'localNtot_median', 'localNtot_uCI')] <- quantile(localNtot_t, probs = c(0.025, 0.5, 0.975))
  
  Ntot_t <- out.mat[,paste0('N.tot[', t, ']')]
  N.data$Ntot_mean[t] <- mean(Ntot_t) 
  N.data[t,c('Ntot_lCI', 'Ntot_median', 'Ntot_uCI')] <- quantile(Ntot_t, probs = c(0.025, 0.5, 0.975))
}

## Write csv
write.csv(N.data, '210818_SvalbardArcticFox_PopEstimates.csv', row.names = F)
