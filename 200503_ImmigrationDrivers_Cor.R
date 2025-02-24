#############################################################
#### INTEGRATED POPULATION MODEL FOR SVALBARD ARCTIC FOX ####
#############################################################

library(coda)
library(reshape)

## Load workspace containing all data
load('200228_AF_IPM_Data.RData')


## Load posterior samples
load('200429_AF_IPM_VersionB4.RData')
ls()

## Re-arrange data
out.mat <- as.matrix(AF.IPM.varB4)
data <- melt(out.mat)
colnames(data) <- c('index', 'parameter', 'value')

#-------------------------------------------------------------#
# FUNCTION FOR POST-HOC CALCULATION OF IMMIGRATION CORRELATES #
#-------------------------------------------------------------#

Imm.posthocTest <- function(MCMC.mat, Tmax, RepAgeClass, i){ # i = sample number
  
  ## Extract time-series of immigrant numbers
  Imm <- MCMC.mat[i, paste("Imm[", 1:Tmax, "]", sep = "")]
  Imm[1] <- NA
  
  ## Extract time-series of population sizes
  Ntot <- MCMC.mat[i, paste("N.tot[", 1:Tmax, "]", sep = "")] # all individuals
  Nad <- Ntot - MCMC.mat[i, paste("N[1, ", 1:Tmax, "]", sep = "")] # adults only
  Nloc <- Ntot - MCMC.mat[i, paste("Imm[", 1:Tmax, "]", sep = "")] # all local individuals
  
  ## Extract time-series of local recruits
  Rtot <- MCMC.mat[i, paste("R.tot[", 1:Tmax, "]", sep = "")]
    
  ## Extract time-series of relevant vital rates
  mHj <- c(MCMC.mat[i, paste("mH[2, ", 1:(Tmax-1), "]", sep = "")], NA)
  mHa <- c(MCMC.mat[i, paste("mH[1, ", 1:(Tmax-1), "]", sep = "")], NA)	
  mOj <- c(MCMC.mat[i, paste("mO[2, ", 1:(Tmax-1), "]", sep = "")], NA)
  mOa <- c(MCMC.mat[i, paste("mO[1, ", 1:(Tmax-1), "]", sep = "")], NA)

  
  ## Arrange vital rates and population sizes into a data frame
  data <- data.frame(Imm = Imm, Ntot = Ntot, Ntot_prev = c(NA, Ntot[-Tmax]), Nad = Nad, Nloc = Nloc, Rtot = Rtot, mHj = mHj, mHa = mHa, mOj = mOj, mOa = mOa)
  
  
  ## Use linear models to calculate coefficients
  results <- data.frame(sample = i, 
  beta_Ntot.Imm = cor(data$Imm, data$Ntot_prev, use = "complete.obs"),
  beta_Nad.Imm = cor(data$Imm, data$Nad, use = "complete.obs"),
  beta_Nloc.Imm = cor(data$Imm, data$Nloc, use = "complete.obs"),
  
  beta_Rtot.Imm = cor(data$Imm, data$Rtot, use = "complete.obs"),
  
  beta_mHj.Imm = cor(data$Imm, data$mHj, use = "complete.obs"),
  beta_mHa.Imm = cor(data$Imm, data$mHa, use = "complete.obs"),
  beta_mOj.Imm = cor(data$Imm, data$mOj, use = "complete.obs"),
  beta_mOa.Imm = cor(data$Imm, data$mOa, use = "complete.obs")
  )
  
  rownames(results) <- NULL
  
  ## Return results
  return(results)
}


#----------------------------------------------------#
# CALCULATE COEFFICIENTS OF DD AND EXTRACT QUANTILES #
#----------------------------------------------------#

## Running the function
Immtest <- do.call("rbind", sapply(1:nrow(out.mat), FUN = function(i) Imm.posthocTest(out.mat, 23, 3, i), simplify = FALSE))

## Save results
save(Immtest, file = '200503_Immtest_Results_Cor.RData')

## Extract quantiles
quantile(Immtest$beta_Ntot.Imm, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
quantile(Immtest$beta_Nad.Imm, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
quantile(Immtest$beta_Nloc.Imm, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
# --> All 95% CI's overlap 0
# --> The 50% CI's for effects of adult/local populationsdo not include 0, and the effects are negative (= less immigrants when local population is larger)

quantile(Immtest$beta_Rtot.Imm, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
# --> 95% CI overlaps 0
# --> 50% CI does not overlap 0, and effect is negative (= less immigrants when more recruits were produced locally)

quantile(Immtest$beta_mHj.Imm, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
quantile(Immtest$beta_mHa.Imm, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
quantile(Immtest$beta_mOj.Imm, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
quantile(Immtest$beta_mOa.Imm, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
# --> All 95% CI's overlap 0
# --> 50% CI's for harvest do not include 0 and indicate a negative effect (less immigrants when more harvest is going on) -> counter to what we would expect under the population sink hypothesis


