#############################################################
#### INTEGRATED POPULATION MODEL FOR SVALBARD ARCTIC FOX ####
#############################################################


#***************#
# DATA OVERVIEW #
#***************#

## General
# A = Number of age classes
# Tmax = Number of years in the study period
# HPeriod[t] = indicator of harvest regulation period for year t (0 = first period, 1 = second period)
# RdCarcass[t] = availability of reindeer carcasses in the winter of year t-1 to t (standardized)
# JanJunSeaIceIsfj[t] = sea ice extent in the winter of year t-1 to t (standardized)
# GooseRep[t] = proportion of juveniles in the goose population during the count in fall of year t


## Age-at-harvest data
# C[a,t] = Number of individuals in age class a harvested during the interval t to t+1
# TmaxC = Number of years in the age-at-harvest data
# pAgeSex[t] = Proportion of individuals harvested in the study area during the interval t to t+1 which could be aged
# pLoc[t] = Proportion of individuals harvested in all of Svalbard during the interval t to t+1 which could be assigned a location (inside vs. outside study area) 


## Mark-recovery data
# y[i,t] = Capture history of individual i spanning all t from marking to the end of the study period
# n.occasions = Number of capture occasions in the mark-recovery data
# nind = Number of individuals in the mark-recovery data
# first[i] = Time index of the marking event of individual i
# firstyear[i,t] = Life-stage of individual i in year t (1 = juvenile/age class 0, 0 = adult)


## Placental scar data
# P1[x] = Number of placental scars counted for carcass of individual x (given any placental scars were present)
# X1 = Length of P1
# P1_age[x] = Age of individual x at harvest
# P1_year[x] = Year in which individual x was harvested
# P2[x] = Presence (1) or absence (0) of placental scars for carcass of individual x 
# X2 = Length of P2
# P2_age[x] = Age of individual x at harvest
# P2_year[x] = Year in which individual x was harvested


## Den survey data
# NoOcc[t] = Number of dens that were observed to be occupied in year t
# k.Dens = Total number of individual dens monitored throughout the study period
# TmaxD = Number of years in the den survey data
# NoPups[x] = Number of pups observed for a given den-year combination x
# X3 = Length of NoPups
# DS_year[x] = Year in which observation NoPups[x] was made



#********************#
# PARAMETER OVERVIEW #
#********************#

## Population Model
# N[a,t] = Number of individuals in age class a at time t
# R[a,t] = Number of recruits produced by age class a at time t
# L[a,t] = Number of offspring conceive by age class a at time t
# Surv[a,t] = Survival of age class a individuals at time t
# S0 = Denning survival of offspring (conception to emergence from the den)
# rho.age[a,t] = Number of placental scars of age class a at time t
# Psi.age[a,t] = Pregnancy rate of age class a at time t


## Age-at-harvest and mark-recovery modules
# mH[a] = Age-specific harvest mortality hazard rate
# mO[a] = Age-specific background mortality hazard rate
# S[a] = Age-specific survival probability
# alpha[a] = Age-specific proportion of deaths due to harvest
# h[a] = Age-specific probability of dying from harvest

# (index a: 1 = adult, 2 = juvenile)


## Placental scar module
# mean.rho = Baseline number of placental scars
# alpha1 = Linear age effect on the number of placental scars
# alpha2 = Quadratic age effect on the number of placental scars
# par.a = Pregancy rate for old females
# par.b = Slope for age effect on logit(pregnany rate)
# par.c = Age when pregancy rate = par.a/2


## Den survey module
# tot.B[t] = Total size of the breeding population at time t
# OR[t] = Den occupancy rate at time t
# meanLS[t] = Mean number of pups present on dens at time t
# u.Dens = Number of unknown dens in the area


## Time-variation
# betaRdCarcass.X = Effect of reindeer carcass availability on vital rate X (on the relevant link scale)
# betaSI.X = Effect of sea ice extent on vital rate X (on the relevant link scale)
# betaG.X = Effect of goose reproduction on vital rate X (on the relevant link scale)
# betaHP.X = Effect of shifted harvest regulation in the second time period on vital rate X (on the relevant link scale)
# betaY.X = Time trend in vital rate X (on the relevant link scale)
# epsilon.X[t] = random year effect on vital rate X in year t (on the relevant link scale)
# sigma.X = standard deviation of random year effects on vital rate X


#************#
# MODEL CODE #
#************#

## NOTE:
# --> Age classes are defined as 0, 1, 2, 3, and 4+ in the paper, but indices in the model code are age class + 1 (e.g. age class 0 has index 1)


fox.code <- nimbleCode({
  
  ##########################  
  #### POPULATION MODEL ####
  ##########################
  
  ### Likelihood (age classes: 0, 1, 2, 3, 4+)
  
  ## Survival
    
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
  
  # Age class 0 (young of the year --> do not reproduce)
  B[1,1:Tmax] <- 0
  L[1,1:Tmax] <- 0
  R[1,1:Tmax] <- 0
  Imm[1] <- 0
  
  # Age classes 1 - 4+    	    
  for(t in 1:Tmax){        				
    for(a in 2:A){
      
      # Breeding Population Size
      B[a,t] ~ dbin(Psi[a,t], N[a,t])
      
      # Litter Size
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
  
  
  for(t in 2:Tmax){
    Imm[t] ~ dpois(ImmT[t])
    ImmT[t] ~ dnorm(avgImm, sd = sigma.Imm)
  }
  
  avgImm ~ dunif(0, 400)
  sigma.Imm ~ dunif(0, 200)
  
  for(a in 1:A){
    N[a,1] ~ dcat(DU.prior[1:400]) # = discrete uniform prior
  }
  
  DU.prior[1:400] <- 1/400
  
  
  ### Additional calculations
  
  for(t in 1:Tmax){
    N.tot[t] <- sum(N[1:A,t])
    R.tot[t] <- sum(R[1:A,t])		
    B.tot[t] <- sum(B[1:A,t])
  }
  
  
  ###############################
  #### AGE-AT-HARVEST MODULE ####
  ###############################
  
  ### Likelihood
  
  for(t in 1:TmaxC){
    
    # Age class 0 (juveniles)
    C[1,t] ~ dbin(h[2,t]*pAgeSex[t]*pLoc[t], N[1,t]) 
    
    # Age classes 1 to 4+ (adults)
    for(a in 2:A){
      
      C[a,t] ~ dbin(h[1,t]*pAgeSex[t]*pLoc[t], N[a,t])
    }
  }
  
  
  ##############################
  #### MARK-RECOVERY MODULE ####
  ##############################

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

    
    mH[1:2,t] <- exp(log(Mu.mH[1:2]) + betaHP.mH*HPeriod[t] + epsilon.mH[t])
    
    mO[1:2,t] <- exp(log(Mu.mO[1:2]) + betaY.mO*(t-12) + betaRC.mO*RdCarcass[t+1] + betaG.mO*Juv[1:2]*GooseRep[t] + betaSI.mO*JanJunSeaIceIsfj[t+1] + epsilon.mO[t])
    
       
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
  betaRC.mO ~ dunif(-5, 5)
  betaG.mO ~ dunif(-5, 5)  
  betaSI.mO ~ dunif(-5, 5)
  
  sigma.mH ~ dunif(0, 5)
  sigma.mO ~ dunif(0, 5)
   
  
  ###############################
  #### PLACENTAL SCAR MODULE ####
  ###############################
  
  ### Likelihood (litter size)
  
  for(x in 1:X1){
    P1[x] ~ dpois(rho[P1_age[x], P1_year[x]])
  }
  
  ### Likelihood (pregnancy rate)
  
  for(x in 1:X2){
    P2[x] ~ dbern(Psi[P2_age[x],P2_year[x]])
  }
  
  ### Priors and constraints (litter size)
 
  for(t in 1:Tmax){

  	## Log-linear model for litter size
  	log(rho[1:A,t]) <- log(mean.rho) + a.eff1*(1:A) + betaY.rho*(t-12) + betaSI.rho*JanJunSeaIceIsfj[t] + betaRC.rho*RdCarcass[t] + epsilon.rho[t]
      
    ## Logit-linear model for pregnancy rate
    logit(eta[1:A,t]) <- par.b*(par.c - (1:A)) + betaY.Psi*(t-12) + betaRC.Psi*RdCarcass[t] + betaSI.Psi*JanJunSeaIceIsfj[t] + epsilon.Psi[t]
    
    Psi[1:A,t] <- par.a*eta[1:A,t] 
    
    epsilon.Psi[t] ~ dnorm(0, sd = sigma.Psi)
    epsilon.rho[t] ~ dnorm(0, sd = sigma.rho) 	
  }
  
  
  mean.rho ~ dunif(0, 16) # 16 is the maximum number ever observed (the mean is 6.4 in the raw data)
  a.eff1 ~ dnorm(0, 1)
   
  betaSI.rho ~ dunif(-5, 5)
  betaRC.rho ~ dunif(-5, 5)
  betaY.rho ~ dunif(-5, 5)

  
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
  
  u.Dens ~ dpois(4)

  
})
