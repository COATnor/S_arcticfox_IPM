#' Simulate a complete set of initial values
#'
#' @param nim.data list of input objects representing data
#' @param nim.constants list of input objects representing constants
#' @param z_inits matrix of initial values for mark-recapture data latent states.
#' @param minN1 integer vector. Lower bound for initial population size (per age class).
#' @param maxN1 integer vector. Upper bound for initial population size (per age class).
#' @param minImm integer. Lower bound for the annual number of immigrants. 
#' @param maxImm integer. Upper bound for the annual number of immigrants.
#' 
#' @return a list containing a complete set of initial values for all parameters
#' in the IPM. 
#' @export
#'
#' @examples

simulateInitVals <- function(nim.data, nim.constants, z_inits, 
                             minN1, maxN1, minImm, maxImm){
  
  Amax <- nim.constants$Amax
  Tmax <- nim.constants$Tmax
  
  Juv <- nim.constants$Juv
  
  #-------------------------------------------------#
  # Set initial values for missing covariate values #
  #-------------------------------------------------#
  
  ## Reindeer carcasses
  RdCarcass <- nim.data$RdCarcass
  if(NA %in% RdCarcass){
    RdCarcass[which(is.na(RdCarcass))] <- mean(RdCarcass, na.rm = TRUE)
  }
  
  ## Goose reproduction
  GooseRep <- nim.data$GooseRep
  if(NA %in% GooseRep){
    GooseRep[which(is.na(GooseRep))] <- mean(GooseRep, na.rm = TRUE)
  }
  
  ## Sea ice availability
  SeaIceIsfj <- nim.data$SeaIceIsfj
  if(NA %in% SeaIceIsfj){
    SeaIceIsfj[which(is.na(SeaIceIsfj))] <- mean(SeaIceIsfj, na.rm = TRUE)
  }

  ## Harvest period
  HPeriod <- nim.data$HPeriod
  if(NA %in% HPeriod){
    stop("Covariate HPeriod contains NA values. HPeriod has to be defined a priori for all years.")
  }
  
  #---------------------------------------------------#
  # Set initial values for vital rate base parameters #
  #---------------------------------------------------#
  
  ## Harvest mortality
  Mu.mH.juv <- runif(1, 0.1, 0.2)
  Mu.mH.ad <- runif(1, 0.1, 0.2)
  Mu.mH <- c(Mu.mH.juv, rep(Mu.mH.ad, Amax - 1))
  
  ## Natural mortality
  Mu.mO.juv <- runif(1, 1.1, 1.3)
  Mu.mO.ad <- runif(1, 0.4, 0.6)
  Mu.mO <- c(Mu.mO.juv, rep(Mu.mO.ad, Amax - 1))
  
  ## Pregnancy rate
  par.a <- runif(1, 0.8, 1)
  par.b <- runif(1, -2, 2)
  par.c <- runif(1, 2, 4)
  
  ## Placental scars
  mean.rho <- runif(1, 3, 6)
  a.eff1 <- runif(1, -0.2, 0.2)
  
  ## Early survival
  Mu.S0 <- runif(1, 0.8, 0.9)
  
  ## Immigration
  avgImm <- runif(1, minImm, maxImm)
  sigma.Imm <- runif(1, 10, 40)
  
  # Mu.Imm <- runif(1, minImm, maxImm)
  # logsigma.Imm <- sd(log(truncnorm::rtruncnorm(10000, a = 0, b = maxImm, mean = Mu.Imm, sd = sigma.Imm)))
  
  
  ## Random effect standard deviations
  sigma.mH <- runif(1, 0.05, 0.5)
  sigma.mO <- runif(1, 0.05, 0.5)
  sigma.Psi <- runif(1, 0.05, 0.5)
  sigma.rho <- runif(1, 0.05, 0.5)
  sigma.m0 <- runif(1, 0.05, 0.5)
  
  ## Random effects (initialize at to 0)
  epsilon.mH <- rep(0, Tmax - 1)
  epsilon.mO <- rep(0, Tmax - 1)
  epsilon.Psi <- rep(0, Tmax)
  epsilon.rho <- rep(0, Tmax)
  epsilon.m0 <- rep(0, Tmax)
  
  ## Covariate effects
  
  betaHP.mH <- runif(1, -0.2, 0.2)
  
  betaRC.mO <- runif(1, -0.2, 0)
  betaSI.mO <- runif(1, -0.2, 0.2)
  betaG.mO <- runif(1, -0.2, 0)
  betaY.mO <- runif(1, -0.02, 0.02)
  
  betaRC.rho <- runif(1, 0, 0.2)
  betaSI.rho <- runif(1, -0.2, 0.2)
  betaY.rho <- runif(1, -0.02, 0.02)
  
  betaRC.Psi <- runif(1, 0, 0.2)
  betaSI.Psi <- runif(1, -0.2, 0.2)
  betaY.Psi <- runif(1, -0.02, 0.02)
  
  betaRC.m0 <- runif(1, -0.2, 0)
  betaSI.m0 <- runif(1, -0.2, 0.2)
  betaY.m0 <- runif(1, -0.02, 0.02)
  
  
  #-------------------------------------#
  # Calculate year-specific vital rates #
  #-------------------------------------#
  
  mH <- mO <- matrix(NA, nrow = Amax, ncol = Tmax - 1)
  eta <- rho <- matrix(NA, nrow = Amax, ncol = Tmax)
  m0t <- rep(NA, Tmax)
  
  for(t in 1:Tmax){
    
    if(t < Tmax){
      ## Harvest mortality hazard rate
      mH[1:Amax, t] <- exp(log(Mu.mH[1:Amax]) + betaHP.mH*HPeriod[t] + epsilon.mH[t])
      
      ## Other (natural) mortality hazard rate
      mO[1:Amax, t] <- exp(log(Mu.mO[1:Amax]) + betaY.mO*(t-12) + betaRC.mO*RdCarcass[t+1] + betaG.mO*Juv[1:Amax]*GooseRep[t] + betaSI.mO*SeaIceIsfj[t+1] + epsilon.mO[t])
    }
    
    ## Pregnancy rate
    eta[1:Amax, t] <- plogis(par.b*(par.c - (1:Amax)) + betaY.Psi*(t-12) + betaRC.Psi*RdCarcass[t] +  betaSI.Psi*SeaIceIsfj[t] + epsilon.Psi[t])
    
    ## Placental scars
    rho[1:Amax, t] <- exp(log(mean.rho) + a.eff1*(1:Amax) + betaY.rho*(t-12) + betaSI.rho*SeaIceIsfj[t] + betaRC.rho*RdCarcass[t] + epsilon.rho[t])
    
    ## Early mortality
    m0t[t] <- exp(log(-log(Mu.S0)) + betaY.m0*(t-12) + betaRC.m0*RdCarcass[t] + betaSI.m0*SeaIceIsfj[t] + epsilon.m0[t])
  }

  ## Survival probability
  S <- exp(-(mH + mO))
  
  ## Proportion harvest mortality
  alpha <- mH/(mH + mO)

  ## Harvest rate
  h <- (1 - S)*alpha

  ## Pregnancy rate
  Psi <- par.a*eta
  
  ## Early survival
  S0 <- exp(-m0t)
  
  ## Immigrant numbers
  Imm <- round(truncnorm::rtruncnorm(Tmax, a = 0, b = maxImm, 
                                     mean = avgImm, sd = sigma.Imm))
  Imm[1] <- 0
  

  #---------------------------------------------------------#
  # Simulate initial values for population-level quantities #
  #---------------------------------------------------------#
  
  ## Prepare empty vectors and matrices
  N <- B <- L <- R <- matrix(NA, nrow = Amax, ncol = Tmax)

  ## Set initial population sizes
  for(a in 1:Amax){
    N[a, 1] <- round(runif(1, minN1[a], maxN1[a]))
  }
  
  ## Set age class 0 reproductive contributions to 0
  B[1, 1:Tmax] <- L[1, 1:Tmax] <- R[1, 1:Tmax] <- 0
  
  
  for (t in 1:Tmax){

    # a) Project local survivors to the next year
    #---------------------------------------------
    
    if(t > 1){
      
      ## Age classes 1 to 3 (indices = 2, 3, 4): age classes 0, 1, and 2 survivors    
      for(a in 2:(Amax-1)){
        N[a, t] <- rbinom(1, size = N[a-1, t-1], prob = S[a-1, t-1])
      }			
      
      ## Age class 4+ (index = 5): age class 3 and 4+ survivors
      N[Amax, t] <- rbinom(1, size = N[Amax-1, t-1] + N[Amax, t-1], prob = S[Amax, t-1])
    }

      
    # b) Sample through reproductive season 
    #--------------------------------------
      
    for (a in 2:Amax){
        
      ## Breeding Population Size: Number of females that reproduce
      B[a, t] <- rbinom(1, size = N[a, t], prob = Psi[a, t])
        
      ## Litter Size (in utero): Number of pups produced by females of age class a
      L[a, t] <- rpois(1, lambda = B[a, t] * rho[a, t] * 0.5)
        
      ## Number Recruits: Number of pups surviving to emerge from the den
      R[a, t] <- rbinom(1, size = L[a, t], prob = S0[t])
    }
      
      
    # c) Add new recruits and immigrants
    #------------------------------------
    if(t > 1){
      N[1, t] <- sum(R[2:Amax, t]) + Imm[t]
    }  
 }
  
  
  # d) Check for years with more harvests than alive individuals
  #-------------------------------------------------------------
  
  if(any(nim.data$C > N[, 1:nim.constants$TmaxC], na.rm = TRUE)){
    stop('Simulation resulted in less alive than harvested (winter). Retry.')
  }
  
  # e) Derive average litter size
  meanLS <- colSums(R[2:Amax, 1:nim.constants$TmaxD, drop = FALSE]) * 2 / colSums(B[2:Amax, 1:nim.constants$TmaxD, drop = FALSE])
  
  
  #----------------------------------------------------#
  # Simulate initial values for observation parameters #
  #----------------------------------------------------#
  
  ## Number of unknown and monitored dens
  u.Dens <- rpois(1, nim.constants$nrDens_unknown)
  Pmon <- nim.data$NoMon/(nim.constants$k.Dens + u.Dens)
  
  
  #-----------------------------------------------------#
  # Assemble and return complete list of initial values #
  #-----------------------------------------------------#
  
  ## List all initial values
  InitVals <- list(
    N = N,
    B = B,
    L = L,
    R = R,
    Imm = Imm,
    ImmT = Imm,
    meanLS = meanLS,
    
    N.tot = colSums(N),
    B.tot = colSums(B),
    R.tot = colSums(R),
    
    Mu.mH = Mu.mH,
    Mu.mO = Mu.mO,
    par.a = par.a,
    par.b = par.b,
    par.c = par.c,
    mean.rho = mean.rho,
    a.eff1 = a.eff1,
    Mu.S0 = Mu.S0,
    avgImm = avgImm,
    
    sigma.mH = sigma.mH,
    sigma.mO = sigma.mO,
    sigma.Psi = sigma.Psi,
    sigma.rho = sigma.rho,
    sigma.m0 = sigma.m0,
    sigma.Imm = sigma.Imm,
    
    epsilon.mH = epsilon.mH,
    epsilon.mO = epsilon.mO,
    epsilon.Psi = epsilon.Psi,
    epsilon.rho = epsilon.rho,
    epsilon.m0 = epsilon.m0,
    
    betaHP.mH = betaHP.mH,
    betaRC.mO = betaRC.mO,
    betaSI.mO = betaSI.mO,
    betaG.mO = betaG.mO,
    betaY.mO = betaY.mO,
    betaRC.rho = betaRC.rho,
    betaSI.rho = betaSI.rho,
    betaY.rho = betaY.rho,
    betaRC.Psi = betaRC.Psi,
    betaSI.Psi = betaSI.Psi,
    betaY.Psi = betaY.Psi,
    betaRC.m0 = betaRC.m0,
    betaSI.m0 = betaSI.m0,
    betaY.m0 = betaY.m0,
    
    mH = mH,
    mO = mO,
    eta = eta,
    rho = rho,
    m0t = m0t,
    S = S,
    alpha = alpha,
    h = h,
    Psi = Psi,
    S0 = S0,
    
    RdCarcass = RdCarcass,
    GooseRep = GooseRep,
    SeaIceIsfj = SeaIceIsfj,
    
    u.Dens = u.Dens,
    Pmon = Pmon,
    
    z = z_inits,
    
    DU.prior.N = rep(1/nim.constants$uLim.N, nim.constants$uLim.N)
  )
  
  
  ## Return initial values
  return(InitVals)
}
