#' Write NIMBLE code for arctic fox IPM
#'
#' @return an R call object specifying the model structure for the red fox IPM. 
#' @export
#'
#' @examples

writeCode_arcticfoxIPM <- function(){
  
  arcticfox.code <- nimbleCode({
    
    
    ##########################  
    #### POPULATION MODEL ####
    ##########################
    
    ### Likelihood (age classes: 1, 2, 3, 4, 5+)
    
    ## Survival
    
    for(t in 1:(Tmax-1)){
      # Age class 1 (young-of-the-year): local reproduction & immigration
      N[1, t+1] <- sum(R[2:Amax, t+1]) + Imm[t+1]
      
      # Age classes 2 to 4 (1, 2, and 3-year olds): age classes 1, 2, and 3 survivors    
      for(a in 1:(Amax-2)){
        N[a+1, t+1] ~ dbin(S[a, t], N[a, t])
      }			
      
      # Age class 5+ (index = Amax = 5): age class 4 and 5+ survivors
      N[Amax, t+1] ~ dbin(S[Amax, t], N[Amax-1, t] + N[Amax, t])
    }
    
    
    ## Reproduction

    # Age class 1 (young of the year --> do not reproduce in year of birth)
    B[1, 1:Tmax] <- 0
    L[1, 1:Tmax] <- 0
    R[1, 1:Tmax] <- 0
    
    # Age classes 1 to 3+    	    
    for(t in 1:Tmax){        				
      for(a in 2:Amax){
        
        # Breeding Population Size: Number of females that reproduce
        B[a, t] ~ dbin(Psi[a, t], N[a, t])
        
        # Litter Size (in utero): Number of pups produced by females of age class a
        L[a, t] ~ dpois(B[a, t]*rho[a, t]*0.5)
        
        # Number Recruits: Number of pups surviving to emerge from the den
        R[a, t] ~ dbin(S0[t], L[a, t])
      } 
    }
    
    #===============================================================================================
    
    
    
    ############################
    #### DERIVED QUANTITIES ####
    ############################
    
    for(t in 1:Tmax){
      N.tot[t] <- sum(N[1:Amax, t])
      R.tot[t] <- sum(R[1:Amax, t])		
      B.tot[t] <- sum(B[1:Amax, t])
    }
    
    #===============================================================================================
    
    
    
    ###############################
    #### AGE-AT-HARVEST MODULE ####
    ###############################
    
    ### Parameters:
    # N = number of individuals in a given age class at a given time (start of June)
    # h = age- and time-dependent probability of dying from hunting
    
    ### Data:
    # C_w = age-at-harvest matrix
    # pAgeSex = annual proportion of carcasses from core area with complete age/sex info
    # pLoc = annual proportion of carcasses with information on location
    
    ### Likelihood
    
    for(t in 1:TmaxC){
      
      # Age class 1 (juveniles)
      C[1, t] ~ dbin(h[2, t]*pAgeSex[t]*pLoc[t], N[1, t])
      
      # Age classes 2 to 5+ (adults)
      for(a in 2:Amax){
        
        C[a, t] ~ dbin(h[1, t]*pAgeSex[t]*pLoc[t], N[a, t])
      }
    }
    
    #===============================================================================================
    
    
    
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
    
    #===============================================================================================
    
    
    ##############################
    #### MARK-RECOVERY MODULE ####
    ##############################
    
    ### Parameters
    # S = annual survival
    # alpha = annual proportion dying from hunting (relative to all other causes)
    # h = annual probability of dying from hunting = (S-1)*alpha
    
    ### Data
    # y = matrix of capture histories
    # z = matrix of known states
    
    ### Likelihood
    
    for(i in 1:nind){
      
      # Latent state at first capture
      z[i, first[i]] <- 1
      
      for(t in (first[i]+1):n.occasions){
        
        # State process
        z[i, t] ~ dbern(S[lifestage[i, t-1], t-1]*z[i, t-1])
        
        # Observation process
        y[i, t] ~ dbern(alpha[lifestage[i, t-1], t-1]*(z[i, t-1]-z[i, t])) # alpha is used here because we assume that the reporting rate is 1
      }
    }
    
    #===============================================================================================
    
    
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
      NoOcc[t] ~ dbin(Pmon[t], sum(B[2:Amax, t]))
    }
    
    ### Likelihood (Number of pups)
    
    for(x in 1:X3){
      NoPups[x] ~ dpois(meanLS[DS_year[x]])
    }
    
    
    
    #===============================================================================================
    
    
    ################################
    #### PRIORS AND CONSTRAINTS ####
    ################################
    
    ## Survival and mortality
    
    for(t in 1:(Tmax-1)){ 
      
      # Harvest mortality hazard rate
      log(mH[1:Amax, t]) <- log(Mu.mH[1:Amax]) + betaHP.mH*HPeriod[t] + epsilon.mH[t]
      
      # Other (natural) mortality hazard rate
      log(mO[1:Amax, t]) <- log(Mu.mO[1:Amax]) + betaY.mO*(t-12) + betaRC.mO*RdCarcass[t+1] + betaG.mO*Juv[1:Amax]*GooseRep[t] + betaSI.mO*JanJunSeaIceIsfj[t+1] + epsilon.mO[t]
      
      # Survival probability
      S[1:Amax, t] <- exp(-(mH[1:Amax, t] + mO[1:Amax,t]))
      
      # Proportion winter harvest mortality
      alpha[1:Amax, t] <- mH[1:Amax, t]/(mH[1:Amax, t] + mO[1:Amax, t])
      
      # Winter harvest rate
      h[1:Amax, t] <- (1-S[1:Amax, t])*alpha[1:Amax, t]
    }
    
    # Median mortality hazard rates
    
    for(a in 1:2){
      Mu.mH[a] ~ dunif(0, 5)
      Mu.mO[a] ~ dunif(0, 5)
    }
    
    Mu.mH[3:Amax] <- Mu.mH[2]
    Mu.mO[3:Amax] <- Mu.mO[2]
    
    
    # Covariate effects
    betaHP.mH ~ dunif(-5, 5)
    betaY.mO ~ dunif(-5, 5)
    betaRC.mO ~ dunif(-5, 5)
    betaG.mO ~ dunif(-5, 5)  
    betaSI.mO ~ dunif(-5, 5)
    
    
    #---------------------------------------------------------------------------------------------
    
    
    ## Pregnancy rate
    
    for(t in 1:Tmax){

      logit(eta[1:Amax, t]) <- par.b*(par.c - (1:Amax)) + betaY.Psi*(t-12) + betaRC.Psi*RdCarcass[t] +  betaSI.Psi*JanJunSeaIceIsfj[t] + epsilon.Psi[t]
      Psi[1:Amax,t] <- par.a*eta[1:Amax,t] 
    }
    
    par.a ~ dunif(0.5, 1)
    par.b ~ dunif(-5, 5)
    par.c ~ dunif(1, 5)
    
    betaRC.Psi ~ dunif(-5, 5) 
    betaY.Psi ~ dunif(-5, 5)
    betaSI.Psi ~ dunif(-5, 5)
    
    
    #---------------------------------------------------------------------------------------------
    
    
    ## Litter size
    
    for(t in 1:(Tmax+1)){
      
      log(rho[1:Amax, t]) <- log(mean.rho) + a.eff1*(1:Amax) + betaY.rho*(t-12) + betaSI.rho*JanJunSeaIceIsfj[t] + betaRC.rho*RdCarcass[t] + epsilon.rho[t]
    }
    
    mean.rho ~ dunif(0, maxPups)
    a.eff1 ~ dnorm(0, 1)
    
    betaSI.rho ~ dunif(-5, 5)
    betaRC.rho ~ dunif(-5, 5)
    betaY.rho ~ dunif(-5, 5)
    
    #---------------------------------------------------------------------------------------------  
    
    
    ## Denning survival
    
    for(t in 1:Tmax){
      S0[t] <- exp(-m0t[t])
      log(m0t[t]) <- log(-log(Mu.S0)) + betaY.m0*(t-12) + betaRC.m0*RdCarcass[t] + betaSI.m0*JanJunSeaIceIsfj[t] + epsilon.m0[t]
    }
    
    Mu.S0 ~ dunif(0, 1)
    
    betaRC.m0 ~ dunif(-5, 5)
    betaSI.m0 ~ dunif(-5, 5)
    betaY.m0 ~ dunif(-5, 5)
    
    #---------------------------------------------------------------------------------------------
    
    
    ## Immigration
    
    # Lognormal prior for immigrant numbers
    Imm[1] <- 0
    
    for(t in 2:Tmax){
      Imm[t] <- round(ImmT[t])
      ImmT[t] ~ dlnorm(meanlog = log(Mu.Imm), sdlog = logsigma.Imm) 
    }
    
    Mu.Imm ~ dunif(1, uLim.Imm)
    logsigma.Imm ~ dunif(0, 10)
    
    # Derivation of immigration rates
    immR[1:Tmax] <- Imm[1:Tmax] / R.tot[1:Tmax]
      
  
    
    #---------------------------------------------------------------------------------------------
    
    
    ## Initial population size (discrete uniform prior) 
    for(a in 1:Amax){
      N[a, 1] ~ dcat(DU.prior.N[1:uLim.N]) 
    }
    
    DU.prior.N[1:uLim.N] <- 1/uLim.N
    
    ## Average litter size
    for(t in 1:TmaxD){
      meanLS[t] <- (sum(R[2:Amax, t])*2)/sum(B[2:Amax, t])
    }
    
    
    #---------------------------------------------------------------------------------------------
    
    
    ## Random year variation
    for(t in 1:(Tmax-1)){  
      epsilon.mH[t] ~ dnorm(0, sd = sigma.mH)
      epsilon.mO[t] ~ dnorm(0, sd = sigma.mO)
    }
    
    for(t in 1:Tmax){
      epsilon.Psi[t] ~ dnorm(0, sd = sigma.Psi)
      epsilon.rho[t] ~ dnorm(0, sd = sigma.rho) 
      epsilon.m0[t] ~ dnorm(0, sd = sigma.m0)
    }
    
    sigma.mH ~ dunif(0, 5)
    sigma.mO ~ dunif(0, 5)
    sigma.Psi ~ dunif(0, 5)
    sigma.rho ~ dunif(0, 5)
    sigma.m0 ~ dunif(0, 5)

    
    #---------------------------------------------------------------------------------------------
    
    ## Number of unknown dens
    u.Dens ~ dpois(nrDens_unknown)
    
    #===============================================================================================
    
    
    
    ##############################
    #### COVARIATE IMPUTATION ####
    ##############################

    for(t in 1:Tmax){
      RdCarcass[t] ~ dnorm(0, sd = 1)
      GooseRep[t] ~ dnorm(0, sd = 1)
      JanJunSeaIceIsfj[t] ~ dnorm(0, sd = 1)
      
    }
    
  })
  
  
  return(arcticfox.code)
}

