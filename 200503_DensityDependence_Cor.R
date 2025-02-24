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

#---------------------------------------------------------#
# FUNCTION FOR POST-HOC CALCULATION OF COEFFICIENTS OF DD #
#---------------------------------------------------------#

DD.posthocTest <- function(MCMC.mat, Tmax, RepAgeClass, i){ # i = sample number
  
  ## Extract time-series of population sizes
  Ntot <- MCMC.mat[i, paste("N.tot[", 1:Tmax, "]", sep = "")] # all individuals
  Nad <- Ntot - MCMC.mat[i, paste("N[1, ", 1:Tmax, "]", sep = "")] # adults only
  Nloc <- Ntot - MCMC.mat[i, paste("Imm[", 1:Tmax, "]", sep = "")] # all local individuals
  
  ## Extract time-series of vital rates
  mHj <- c(MCMC.mat[i, paste("mH[2, ", 1:(Tmax-1), "]", sep = "")], NA)
  mHa <- c(MCMC.mat[i, paste("mH[1, ", 1:(Tmax-1), "]", sep = "")], NA)	
  mOj <- c(MCMC.mat[i, paste("mO[2, ", 1:(Tmax-1), "]", sep = "")], NA)
  mOa <- c(MCMC.mat[i, paste("mO[1, ", 1:(Tmax-1), "]", sep = "")], NA)
  
  Psi <- MCMC.mat[i, paste("Psi[", RepAgeClass, ", ", 1:Tmax, "]", sep = "")]
  rho <- MCMC.mat[i, paste("rho[", RepAgeClass, ", ", 1:Tmax, "]", sep = "")]
  m0 <- MCMC.mat[i, paste("m0t[", 1:Tmax, "]", sep = "")]
  
  ## Arrange vital rates and population sizes into a data frame
  data <- data.frame(Ntot = Ntot, Ntot_prev = c(NA, Ntot[-Tmax]), Ntot_next = c(Ntot[-1], NA), Nad = Nad, Nloc = Nloc, mHj = mHj, mHa = mHa, mOj = mOj, mOa = mOa, Psi = Psi, rho = rho, m0 = m0)
  
  ## Calculate population growth rate
  data$lambda <- data$Ntot_next/data$Ntot

  ## Use linear models to calculate coefficient of DD
  results <- data.frame(sample = i, 
  beta_Ntot.lambda = cor(data$lambda, data$Ntot, use = "complete.obs"),
  
  beta_Ntot.mHj = cor(data$mHj, data$Ntot, use = "complete.obs"),
  beta_Ntot.mHa = cor(data$mHa, data$Ntot, use = "complete.obs"),
  beta_Ntot.mOj = cor(data$mOj, data$Ntot, use = "complete.obs"),
  beta_Ntot.mOa = cor(data$mOa, data$Ntot, use = "complete.obs"),
  
  beta_Nloc.mHj = cor(data$mHj, data$Nloc, use = "complete.obs"),
  beta_Nloc.mHa = cor(data$mHa, data$Nloc, use = "complete.obs"),
  beta_Nloc.mOj = cor(data$mOj, data$Nloc, use = "complete.obs"),
  beta_Nloc.mOa = cor(data$mOa, data$Nloc, use = "complete.obs"),
  
  beta_Ntot.Psi = cor(data$Psi, data$Ntot_prev, use = "complete.obs"),
  beta_Ntot.rho = cor(data$rho, data$Ntot_prev, use = "complete.obs"),
  
  beta_Nad.Psi = cor(data$Psi, data$Nad, use = "complete.obs"),
  beta_Nad.rho = cor(data$rho, data$Nad, use = "complete.obs"),
  beta_Nad.m0 = cor(data$m0, data$Nad, use = "complete.obs")
  )
  
  rownames(results) <- NULL
  
  ## Return results
  return(results)
}


#----------------------------------------------------#
# CALCULATE COEFFICIENTS OF DD AND EXTRACT QUANTILES #
#----------------------------------------------------#

## Running the function
DDtest <- do.call("rbind", sapply(1:nrow(out.mat), FUN = function(i) DD.posthocTest(out.mat, 23, 3, i), simplify = FALSE))

## Save results
save(DDtest, file = '200503_DDtest_Results_Cor.RData')

## Extract quantiles
quantile(DDtest$beta_Ntot.lambda, probs = c(0.05, 0.5, 0.95))
# --> Negative, 90% CI does not overlap 0 (lower growth rate when population is larger)
# --> Median = -0.63, indicating fairly strong correlation

quantile(DDtest$beta_Ntot.mHj, probs = c(0.05, 0.5, 0.95))
quantile(DDtest$beta_Ntot.mHa, probs = c(0.05, 0.5, 0.95))
quantile(DDtest$beta_Nloc.mHj, probs = c(0.05, 0.5, 0.95))
quantile(DDtest$beta_Nloc.mHa, probs = c(0.05, 0.5, 0.95))
# --> All CI's overlap 0

quantile(DDtest$beta_Ntot.mOj, probs = c(0.05, 0.5, 0.95))
quantile(DDtest$beta_Ntot.mOa, probs = c(0.05, 0.5, 0.95))
quantile(DDtest$beta_Nloc.mOj, probs = c(0.05, 0.5, 0.95))
quantile(DDtest$beta_Nloc.mOa, probs = c(0.05, 0.5, 0.95))
# --> Positive, 90% CI's do not overlap 0 (higher natural mortality when population size is larger)
# --> Higher median for juveniles than adults, and for total population than local population only 

quantile(DDtest$beta_Ntot.Psi, probs = c(0.05, 0.5, 0.95))
quantile(DDtest$beta_Nad.Psi, probs = c(0.05, 0.5, 0.95))
# --> All CI's overlap 0
# --> Effect size larger for total population the previous year than adult population in current year

quantile(DDtest$beta_Ntot.rho, probs = c(0.05, 0.5, 0.95))
quantile(DDtest$beta_Nad.rho, probs = c(0.05, 0.5, 0.95))
# --> All CI's overlap 0
# --> Effect size larger for adult population in current year than total population the previous year

quantile(DDtest$beta_Nad.m0, probs = c(0.05, 0.5, 0.95))
# --> CI overlaps 0


# --> Overall, this seems to indicate that there may be negative DD in population growth rate which may be mediated primarily via natural mortality (but immigration remains to be checked)





quantile(DDtest$beta_Ntot.mHj, probs = c(0.05, 0.25, 0.5, 0.75, 0.97))
quantile(DDtest$beta_Ntot.mHa, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
quantile(DDtest$beta_Nloc.mHj, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
quantile(DDtest$beta_Nloc.mHa, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
# --> All CI's overlap 0, also for 50% CI's


quantile(DDtest$beta_Ntot.Psi, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
quantile(DDtest$beta_Nad.Psi, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
# --> 50% CI for effect of Ntot does not include 0 and is negative (more individuals present before winter leads to less breeding / lower proportion breeding in the next spring)

quantile(DDtest$beta_Ntot.rho, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
quantile(DDtest$beta_Nad.rho, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
# --> All CI's overlap 0, also for 50% CI's

quantile(DDtest$beta_Nad.m0, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
# --> All CI's overlap 0
