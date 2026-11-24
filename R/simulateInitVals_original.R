simulateInitVals_original <- function(z_inits, Tmax, Amax){
  
  initVals <- list(
    #z = IPM.data$initsCH.CMR, 
    z = z_inits,
    Mu.mH = c(runif(1, 0.1, 0.2), rep(runif(1, 0.1, 0.2), Amax - 1)),
    Mu.mO = c(runif(1, 1.1, 1.3), rep(runif(1, 0.4, 0.6), Amax - 1)),
    betaHP.mH = runif(1, -0.2, 0.2),
    betaRC.mO = runif(1, -0.2, 0),
    betaY.mO = 0,
    betaG.mO = 0,
    betaSI.mO = 0,
    mean.rho = runif(1, 3, 6),
    a.eff1 = runif(1, -0.2, 0.2),
    betaRC.rho = 0, 
    betaSI.rho = 0, 
    betaY.rho = 0,
    par.a = runif(1, 0.8, 1),
    par.b = runif(1, -2, 2),
    par.c = runif(1, 2, 4),
    betaRC.Psi = runif(1, 0, 0.2),
    betaY.Psi = runif(1, 0, 0.2),
    betaSI.Psi = 0,
    Mu.S0 = runif(1, 0.8, 0.9),
    betaRC.m0 = 0,
    betaSI.m0 = 0,
    betaY.m0 = 0,
    sigma.mH = runif(1, 0.05, 0.5),
    epsilon.mH = rep(0, Tmax - 1),
    sigma.mO = runif(1, 0.05, 0.5),
    epsilon.mH = rep(0, Tmax - 1),
    sigma.Psi = runif(1, 0.05, 0.5),
    epsilon.Psi = rep(0, Tmax),
    sigma.rho = runif(1, 0.05, 0.5),
    epsilon.rho = rep(0, Tmax),
    sigma.m0 = runif(1, 0.05, 0.5),
    epsilon.m0 = rep(0, Tmax),
    avgImm = runif(1, 50, 100),
    sigma.Imm = runif(1, 10, 40),
    #sigma.Imm = runif(1, 0.5, 2),
    GooseRep = c(rep(NA, 22), 0),
    
    u.Dens = rpois(1, 4)
    #Imm = ,
    # B = ,
    # L = ,  
    # R = , 
    # NOTE: NIMBLE runs without these initial values, and seems to do quite well at filling in the blanks
    
    #N =  N.inits # --> N.inits from bottom of code
    
  )
  
  return(initVals)
  
}