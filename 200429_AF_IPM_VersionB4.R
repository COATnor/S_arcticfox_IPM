#############################################################
#### INTEGRATED POPULATION MODEL FOR SVALBARD ARCTIC FOX ####
#############################################################

# NOTE: All of the models in this code will use immigration model 2 (constrained by a time-varying mean)


#library(MASS)
#library(plyr)
#library(dclone)
#library(reshape2)
#library(ggplot2)

library(coda)
library(nimble)
nimbleOptions(disallow_multivariate_argument_expressions = FALSE)

mySeed <- 0

##########################
#### DATA PREPARATION ####
##########################

## Load workspace containing all data
load('Data/200228_AF_IPM_Data.RData')
IPM.data <- readRDS('AF_IPM_Data.rds')

## Data Object
str(IPM.data)

## Data for Age-at-Harvest module (based on object 'carcass_core')
IPM.data$C # Age-at-Harvest matrix (5 age classes, 22 seasons)
IPM.data$A # Oldest age class (5+)
IPM.data$TmaxC # Number of years in harvest carcass data
IPM.data$pAgeSex # Annual proportion core area harvests with sex and age information
IPM.data$pLoc # Annual proportion of harvests with location information

## Data for Age-at-Death module
IPM.data$M # Age-at-Death matrix (5 age classes, 20 seasons)
IPM.data$TmaxM # Number of years in other carcass data
IPM.data$pAgeSex.M # Proportion core area recoveries with sex and age information
IPM.data$pLoc.M # Proportion of recoveries with location information

## Data for Mark-Recovery module (based on object 'fox')
IPM.data$y.CMR # Brownie-type capture histories (142 individuals, 9 years)
IPM.data$initsCH.CMR # Initial values for latent states in Brownie model
IPM.data$n.occasions # Number of years in MR data
IPM.data$nind # Number of individuals in MR data
IPM.data$firstyear # matrix indicating life stage of individuals in all years (1 = juvenile, 0 = adult)
IPM.data$first # First year for each individual

## Additional data for Mark-Recapture-Recovery module extension
IPM.data$y.CMRR # Multistate capture histories incl. recaptures and harvest recoveries 
IPM.data$initsCH.CMRR # Initial values for latent states in multistate models (recaptures + harvest recoveries)
IPM.data$y.CMRR2 # Multistate capture histories incl. recaptures and all recoveries 
IPM.data$y.CMRR3 # Multistate capture histories incl. recaptures and all recoveries (except other cause recoveries NOT included in carcass data)
IPM.data$initsCH.CMRR # Initial values for latent states in multistate models (recaptures + harvest recoveries)

## Data for Placental Scar module (based on object 'plac_scar2')
IPM.data$P1 # Vector of observed numbers of placental scars (across individuals and years)
IPM.data$P1_age # Ages associated with elements of P1
IPM.data$P1_year # Years associated with elements of P1
IPM.data$X1 # Length of P1
IPM.data$P2 # Vector of observed presence/absence of placental scars (across individuals and years)
IPM.data$P2_age # Ages associated with elements of P2
IPM.data$P2_year # Years associated with elements of P2
IPM.data$X2 # Length of P2

## Data for Den Survey module (based on objects 'denS_sum' and 'denS_pups')
IPM.data$NoOcc # Vector of the number of occupied dens in each year
IPM.data$NoMon # Vector of the number of monitored dens in each year
IPM.data$NoPups # Vector of the minimum number of pups observed on each den in each year
IPM.data$DS_year # Years associated with NoPups
IPM.data$X3 # Length of NoPups
IPM.data$k.Dens # Total number of known dens in the study area


## Re-arrange into data and constants for nimble
fox.data <- list(
  C=IPM.data$C,
  y=IPM.data$y.CMR, 
  P1=IPM.data$P1,
  P2=IPM.data$P2,
  NoOcc=IPM.data$NoOcc, NoMon=IPM.data$NoMon, 
  NoPups=IPM.data$NoPups, 
  RdCarcass = IPM.data$RdCarcass, 
  WinterTemp = IPM.data$MeanWinterTemp, 
  GooseCount = IPM.data$GooseCount, 
  GooseRep = IPM.data$GooseJuvProp,
  #GooseRep0 = c(IPM.data$GooseJuvProp[-23],0), # Missing value set to 0 for testing
  #PtarmDens = IPM.data$PtarmiganDens,
  #PtarmDens0 = c(0, 0, 0, IPM.data$PtarmiganDens[4:22], 0), # Missing values set to 0 for testing
  #PtarmRep = IPM.data$PtarmiganJuvProp,
  #PtarmRep0 = c(0, 0, IPM.data$PtarmiganJuvProp[3:23]), # Missing values set o 0 for testing
  #PbBCI = IPM.data$PolarBearBCI, 
  #AprSeaIceSv = IPM.data$AprSeaIceSvalbard,
  #AprSeaIceSv0 = c(IPM.data$AprSeaIceSvalbard[-23], 0), 
  #AprMaySeaIceIsfj = IPM.data$AprMaySeaIceIsfj,
  JanJunSeaIceIsfj = IPM.data$JanJunSeaIceIsfj,
  WinterAO = IPM.data$WinterAO,
  HPeriod = c(rep(0, 12), rep(1, 11)))


fox.constants <- list(
  A=IPM.data$A, Tmax=IPM.data$Tmax, TmaxC = IPM.data$TmaxC, TmaxD = IPM.data$TmaxD,
  n.occasions=IPM.data$n.occasions, nind=IPM.data$nind, first=IPM.data$first, firstyear=IPM.data$firstyear, 
  X1=IPM.data$X1,  P1_age=IPM.data$P1_age, P1_year=IPM.data$P1_year,
  X2=IPM.data$X2, P2_age=IPM.data$P2_age, P2_year=IPM.data$P2_year,
  X3=IPM.data$X3,  DS_year=IPM.data$DS_year, k.Dens=IPM.data$k.Dens,
  pAgeSex=IPM.data$pAgeSex, pLoc=IPM.data$pLoc)

############################
#### PARAMETER OVERVIEW ####
############################

## Population Model
# N[a,t] = Number of individuals in age class a at time t
# R[a,t] = Number of recruits produced by age class a at time t
# L[a,t] = Number of offspring conceive by age class a at time t
# Surv[a,t] = Survival of age class a individuals at time t
# S0 = Denning survival of offspring (conception to emergence from the den)
# rho.age[a,t] = number of placental scars of age class a at time t
# Psi.age[a,t] = pregnancy rate of age class a at time t

## AaH and MR modules
# mH[a] = age-specific harvest mortality hazard rate
# mO[a] = age-specific background mortality hazard rate
# S[a] = age-specific survival probability
# alpha[a] = age-specific proportion of deaths due to harvest
# h[a] = age-specific probability of dying from harvest

# (index: 1 = adult, 2 = juvenile)

## PC module
# mean.rho = baseline number of placental scars
# alpha1 = linear age effect on the number of placental scars
# alpha2 = quadratic age effect on the number of placental scars
# par.a = pregancy rate for old females
# par.b = slope for age effect on logit(pregnany rate)
# par.c = age when pregancy rate = par.a/2

## DS module
# tot.B[t] = total size of the breeding population at time t
# OR[t] = den occupancy rate at time t
# meanLS[t] = mean number of pups present on dens at time t
# u.Dens = number of unknown dens in the area


#######################################
#### SeaIce + RdCarcass + GooseRep ####
#######################################
# Includes random year variation in mH, Psi, and rho
# Sea ice covariate: January-June Isfjorden sea ice
# Goose covariate: count -> denning mortality, reproduction -> juv. mortality
# Trends centered around t=12 (year 2008)

fox.code.varB4 <- nimbleCode({
  
  ##########################  
  #### POPULATION MODEL ####
  ##########################
  
  ### Likelihood (age classes: 1, 2, 3, 4, 5+)
  
  ## Survival
  
  # Age class 1
  #N[1,2:Tmax] <- sum(R[2:A,2:Tmax]) + Imm[2:Tmax]     
  # NOTE: Should be possible to vectorize the above calculation, but a function for making matrix rowsums is needed
  
  for(t in 1:(Tmax-1)){ 
    
    # Age class 1
    N[1,t+1] <- sum(R[2:A,t+1]) + Imm[t+1] # Assuming all immigrants are juveniles
    
    
    # Age classes 2 to 4          				
    for(a in 1:(A-2)){ 
      
      N[a+1,t+1] ~ dbin(Surv[a,t], N[a,t])
    } 
    
    # Age class 5+
    N[A,t+1] ~ dbin(Surv[A,t], N[A-1,t] + N[A,t])
  }
  
  
  ## Reproduction
  
  # Age class 1 (young of the year --> do not reproduce)
  B[1,1:Tmax] <- 0
  L[1,1:Tmax] <- 0
  R[1,1:Tmax] <- 0
  Imm[1] <- 0
  
  # Age classes 2 - 5+    	    
  for(t in 1:Tmax){        				
    for(a in 2:A){
      
      # Breeding Population Size
      #B[a,t] ~ dbin(Psi.age[a], N[a,t])
      B[a,t] ~ dbin(Psi[a,t], N[a,t])
      
      # Litter Size
      #L[a,t] ~ dpois(B[a,t]*rho.age[a]*0.5)
      L[a,t] ~ dpois(B[a,t]*rho[a,t]*0.5)
      
      # Number Recruits
      R[a,t] ~ dbin(S0t[t], L[a,t])
    } 
  }
  
  
  ### Priors and constraints

  for(t in 1:(Tmax-1)){
    Surv[1,t] <-  S[2,t]
    Surv[2:A,t] <-  S[1,t]
  }
    
  for(t in 1:Tmax){
    S0t[t] <- exp(-m0t[t])
    log(m0t[t]) <- log(-log(S0)) + betaY.m0*(t-12) + betaRC.m0*RdCarcass[t] + betaSI.m0*JanJunSeaIceIsfj[t] + epsilon.m0[t]
    epsilon.m0[t] ~ dnorm(0, sd = sigma.m0)
  }
  
  S0 ~ dunif(0, 1)
  
  betaRC.m0 ~ dunif(-5, 5)
  betaSI.m0 ~ dunif(-5, 5)
  betaY.m0 ~ dunif(-5, 5)
  
  sigma.m0 ~ dunif(0, 5)
  
  # Discrete uniform distribution
  
  for(t in 2:Tmax){
    Imm[t] ~ dpois(ImmT[t])
    ImmT[t] ~ dnorm(avgImm, sd = sigma.Imm)
    # In theory, ImmT could be negative (non-truncated normal distribution)
    # --> this may end up causing issues at some point (does not now)
  }
  
  avgImm ~ dunif(0, 400)
  sigma.Imm ~ dunif(0, 200)
  
  for(a in 1:A){
    N[a,1] ~ dcat(DU.prior[1:400])
  }
  
  DU.prior[1:400] <- 1/400
  
  
  ### Additional calculations
  
  for(t in 1:Tmax){
    N.tot[t] <- sum(N[1:A,t])
    R.tot[t] <- sum(R[1:A,t])		
    B.tot[t] <- sum(B[1:A,t])
  }
  # NOTE: Should be possible to vectorize the above calculation, but a function for making matrix rowsums is needed
  
  
  ###############################
  #### AGE-AT-HARVEST MODULE ####
  ###############################
  
  ### Parameters:
  # N = number of individiuals in a given age class at a given time
  # h = time-dependent probability of dying from hunting ([1] = adults, [2] = juveniles)
  
  ### Data:
  # C = age-at-harvest matrix
  
  ### Likelihood
  
  for(t in 1:TmaxC){
    
    # Age class 1 (juveniles)
    #C[1,t] ~ dbin(h[2,t], N[1,t])
    C[1,t] ~ dbin(h[2,t]*pAgeSex[t]*pLoc[t], N[1,t]) # Added year-spec. proportions
    
    # Age classes 2 to 5+ (adults)
    for(a in 2:A){
      
      #C[a,t] ~ dbin(h[1,t], N[a,t])
      C[a,t] ~ dbin(h[1,t]*pAgeSex[t]*pLoc[t], N[a,t]) # Added year-spec. proportions
    }
  }
  
  ##############################
  #### MARK-RECOVERY MODULE ####
  ##############################
  
  ### Parameters
  # mH = time-dependent hunting mortality hazard rate ([1,] = adults, [2,] = juveniles)
  # mean.mO = time-dependent background mortality hazard rate ([1,] = adults, [2,] = juveniles)
  # S = annual survival
  # alpha = annual proportion dying from hunting (relative to all other causes)
  # h = annual probability of dying from hunting = (S-1)*alpha
  
  ### Data
  # y = matrix of capture histories
  # z = matrix of known states
  
  ### Likelihood
  
  for(i in 1:nind){
    
    # Latent state at first capture
    z[i,first[i]] <- 1
    
    for(t in (first[i]+1):n.occasions){
      
      # State process
      z[i,t] ~ dbern(S[firstyear[i,t-1]+1,t-1]*z[i,t-1])
      
      # Observation process
      y[i,t] ~ dbern(alpha[firstyear[i,t-1]+1,t-1]*(z[i,t-1]-z[i,t])) # alpha is used here because we assume that the reporting rate is 1
    }
  }
  
  ### Priors and constraints
  
  for(a in 1:2){ # Index 1 = adults, index 2 = juveniles    
    Mu.mH[a] ~ dunif(0.01, 3)
    Mu.mO[a] ~ dunif(0.01, 3)
  }
  
  for(t in 1:(Tmax-1)){
    
    #mH[1:2,t] <- Mu.mH[1:2]
    #mO[1:2,t] <- Mu.mO[1:2]
    
    mH[1:2,t] <- exp(log(Mu.mH[1:2]) + betaHP.mH*HPeriod[t] + epsilon.mH[t])
    
    mO[1:2,t] <- exp(log(Mu.mO[1:2]) + betaY.mO*(t-12) + betaT.mO*WinterTemp[t+1] + betaRC.mO*RdCarcass[t+1] + betaG.mO*Juv[1:2]*GooseRep[t] + betaSI.mO*JanJunSeaIceIsfj[t+1] + epsilon.mO[t])
    
    #mO[1:2,t] <- exp(log(Mu.mO[1:2]) + betaY.mO*(t-12) + betaT.mO*WinterTemp[t+1] + betaRC.mO*RdCarcass[t+1] + betaG.mO*GooseRep[t] + betaSI.mO*JanJunSeaIceIsfj[t+1] + epsilon.mO[t])
    
    
    S[1:2,t] <- exp(-(mH[1:2,t]+mO[1:2,t]))
    alpha[1:2,t] <- mH[1:2,t]/(mH[1:2,t]+mO[1:2,t])
    h[1:2,t] <- (1 - S[1:2,t])*alpha[1:2,t]
    
    epsilon.mH[t] ~ dnorm(0, sd = sigma.mH)
    epsilon.mO[t] ~ dnorm(0, sd = sigma.mO)
  }
  
  Juv[1] <- 0
  Juv[2] <- 1
  
  betaHP.mH ~ dunif(-5, 5)
  betaY.mO ~ dunif(-5, 5)
  #betaT.mO ~ dunif(-5, 5)
  betaRC.mO ~ dunif(-5, 5)
  betaG.mO ~ dunif(-5, 5)  
  betaSI.mO ~ dunif(-5, 5)
  
  sigma.mH ~ dunif(0, 5)
  sigma.mO ~ dunif(0, 5)
  
  #betaHP.mH <- 0
  #betaY.mO <- 0
  betaT.mO <- 0
  #betaRC.mO <- 0
  #sigma.mH <- 0    
  
  ###############################
  #### PLACENTAL SCAR MODULE ####
  ###############################
  
  ### Parameters:
  # rho = expected number of placental scars (fetuses)
  # a.eff1 = log-linear age effect on rho
  # a.eff2 = log-quadratic age effect on rho
  
  # Psi = pregnancy rate
  # par.a = pregnancy rate of old females
  # par.b = slope of age effect on logit(Psi)
  # par.c = age when Psi = a/2
  # betaRC.Psi = effect of reindeer carcass availability on pregnancy rate
  
  ## Data:
  # P1 = individual placental scar counts
  # P1_age = individual ages associated with P1
  
  # P2 = individual presence/absence of placental scars
  # P2_age = individual ages associated with P2
  
  
  ### Likelihood (litter size)
  
  for(x in 1:X1){
    P1[x] ~ dpois(rho[P1_age[x], P1_year[x]])
  }
  
  ### Likelihood (pregnancy rate)
  
  for(x in 1:X2){
    P2[x] ~ dbern(Psi[P2_age[x],P2_year[x]])
  }
  
  ### Priors and constraints (litter size)
  
  
  ## Log-linear model for litter size
  #log(rho[1:A]) <- log(mean.rho) + a.eff1*(1:A)
  
  ## Logit-linear model for pregnancy rate
  for(t in 1:Tmax){

  	## Log-linear model for litter size
  	log(rho[1:A,t]) <- log(mean.rho) + a.eff1*(1:A) + betaY.rho*(t-12) + betaSI.rho*JanJunSeaIceIsfj[t] + betaRC.rho*RdCarcass[t] + betaT.rho*WinterTemp[t] + epsilon.rho[t]
      
    #logit(eta[1:A,t]) <- par.b*(par.c - (1:A))
    logit(eta[1:A,t]) <- par.b*(par.c - (1:A)) + betaY.Psi*(t-12) + betaRC.Psi*RdCarcass[t] +  betaT.Psi*WinterTemp[t] + betaSI.Psi*JanJunSeaIceIsfj[t] + epsilon.Psi[t]
    
    Psi[1:A,t] <- par.a*eta[1:A,t] 
    
    epsilon.Psi[t] ~ dnorm(0, sd = sigma.Psi)
    epsilon.rho[t] ~ dnorm(0, sd = sigma.rho) 	
  }
  
  
  mean.rho ~ dunif(0, 16) # 16 is the maximum number ever observed (the mean is 6.4 in the data)
  a.eff1 ~ dnorm(0, 1)
   
  betaSI.rho ~ dunif(-5, 5)
  betaRC.rho ~ dunif(-5, 5)
  #betaT.rho ~ dunif(-5, 5)
  betaY.rho ~ dunif(-5, 5)
  
  #betaSI.rho <- 0 
  #betaRC.rho <- 0 
  betaT.rho <- 0 
  
  ## Priors and constraints (pregnancy rate)
  
  par.a ~ dunif(0.5, 1)
  par.b ~ dunif(-5, 5)
  par.c ~ dunif(1, 5)
  
  betaRC.Psi ~ dunif(-5, 5) 
  betaY.Psi ~ dunif(-5, 5)
  
  #betaT.Psi ~ dunif(-5, 5) # Winter temperature
  betaSI.Psi ~ dunif(-5, 5) # Sea ice
  
  betaT.Psi <- 0
  
  sigma.Psi ~ dunif(0, 5)
  sigma.rho ~ dunif(0, 5)
  
  
  ###########################
  #### DEN SURVEY MODULE ####
  ###########################
  
  ### Parameters:
  # B: Number of breeding individuals (per age class)
  # Pmon: Annual proportion of dens that are monitored (relative to all dens in the area)
  # u.Dens: Number of unknown breeding dens
  # meanLS: Average annual litter size (males and females)
  
  ### Data:
  # NoMon: Number of monitored dens in each year
  # NoOcc: Number of occupied dens in each year
  # NoPups: Number of pups (males and females) observed per den
  # k.Dens: Number of known breeding dens
  
  ### Likelihood (Breeding population size)
  
  Pmon[1:TmaxD] <- NoMon[1:TmaxD]/(k.Dens + u.Dens)
  
  for(t in 1:TmaxD){ 
    NoOcc[t] ~ dbin(Pmon[t], sum(B[2:A,t]))
  }
  
  ### Likelihood (Number of pups)
  
  for(x in 1:X3){
    NoPups[x] ~ dpois(meanLS[DS_year[x]])
  }
  
  ### Priors and constraints 
  
  for(t in 1:TmaxD){
    meanLS[t] <- (sum(R[2:A,t])*2)/sum(B[2:A,t])
  }
  # NOTE: Should be possible to vectorize the above calculation, but a function for making matrix rowsums is needed
  
  u.Dens ~ dpois(4)
  
  #########################
  #### COVARIATE MODEL ####
  #########################
  
  # Currently, the last value for the GooseRep covariate is NA
  # Here, we estimate it assuming it follows the same distribution as the remaining values
  
  GooseRep[23] ~ dnorm(0, sd = 1) # Mean = 0 and sd = 1 because the covariate is standardized
  
})

## Setting initial values
IPM.inits <- function() {list(
  z = IPM.data$initsCH.CMR, 
  Mu.mH = runif(2,0.1, 0.2), Mu.mO = c(runif(1, 0.4, 0.6),runif(1, 1.1, 1.3)),
  betaHP.mH = runif(1,-0.2,0.2),
  #betaT.mO = runif(1,-0.2,0.2),
  betaRC.mO = runif(1,-0.2,0),
  betaY.mO = 0,
  betaG.mO = 0,
  betaSI.mO = 0,
  mean.rho = runif(1, 3, 6),
  a.eff1 = runif(1,-0.2,0.2),
  #betaT.rho = 0,
  betaRC.rho = 0, 
  betaSI.rho = 0, 
  betaY.rho = 0,
  #a.eff2 = runif(1,-0.02,0.02),
  par.a = runif(1, 0.8, 1),
  par.b = runif(1,-2, 2),
  par.c = runif(1, 2, 4),
  betaRC.Psi = runif(1, 0, 0.2),
  betaY.Psi = runif(1, 0, 0.2),
  #betaT.Psi = 0,
  betaSI.Psi = 0,
  S0 = runif(1, 0.8, 0.9),
  betaRC.m0 = 0,
  betaSI.m0 = 0,
  betaY.m0 = 0,
  sigma.mH = runif(1,0.05,0.5),
  epsilon.mH = rep(0, fox.constants$Tmax-1),
  sigma.mO = runif(1,0.05,0.5),
  epsilon.mH = rep(0, fox.constants$Tmax-1),
  sigma.Psi = runif(1,0.05,0.5),
  epsilon.Psi = rep(0, fox.constants$Tmax),
  sigma.rho = runif(1,0.05,0.5),
  epsilon.rho = rep(0, fox.constants$Tmax),
  sigma.m0 = runif(1,0.05,0.5),
  epsilon.m0 = rep(0, fox.constants$Tmax),
  avgImm = runif(1, 50, 100),
  sigma.Imm = (runif(1, 10, 40)),
  GooseRep = c(rep(NA, 22), 0),
  
  u.Dens = rpois(1, 4)
  #Imm = ,
  # B = ,
  # L = ,  
  # R = , 
  # NOTE: NIMBLE runs without these initial values, and seems to do quite well at filling in the blanks
  
  #N =  N.inits # --> N.inits from bottom of code
  
)} 


initValsFull <- list(IPM.inits(), IPM.inits(), IPM.inits())


## Setting parameters to monitor
IPM.params <- c(
"Mu.mH", "betaHP.mH", "sigma.mH", 
"Mu.mO", "betaY.mO", "betaRC.mO", "betaG.mO", "betaSI.mO", "sigma.mO",  
"S0", "betaRC.m0", "betaSI.m0", "betaY.m0", "sigma.m0", 
"mean.rho", "a.eff1", "betaRC.rho", "betaSI.rho", "betaY.rho", "sigma.rho", 
"par.a", "par.b", "par.c", "betaY.Psi", "betaRC.Psi", "betaSI.Psi", "sigma.Psi", 
"N.tot", "N", "R.tot", "B.tot", "B", "meanLS", 
"avgImm", "sigma.Imm", "Imm", 
"epsilon.mH", "epsilon.mO", "epsilon.m0", "epsilon.Psi", "epsilon.rho",
"m0t", "mH", "mO", "Psi", "rho",
"u.Dens")


## Test run
niter <- 20000
nburnin <- 5000
nchains <- 3
nthin <- 1


#AF.IPM.varB4 <- nimbleMCMC(fox.code.varB4, fox.constants, fox.data, initValsFull, monitors = IPM.params, niter = 2, nburnin = 0, nchains = nchains, thin = 1, setSeed = mySeed, samplesAsCodaMCMC = TRUE)

AF.IPM.varB4 <- nimbleMCMC(fox.code.varB4, fox.constants, fox.data, initValsFull, monitors = IPM.params, niter = niter, nburnin = nburnin, nchains = nchains, thin = nthin, setSeed = mySeed, samplesAsCodaMCMC = TRUE)


save(AF.IPM.varB4, file = '200429_AF_IPM_VersionB4.RData')



