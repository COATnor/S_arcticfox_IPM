#' Set up model, data, and initial values for running MCMC
#'
#' @param modelCode string. Relative path to the model file to be used
#' @param nim.data list of input objects representing data
#' @param nim.constants list of input objects representing constants.
#' @param minN1 integer vector. Lower bound for initial population size 
#' (per age class) to use in initial value simulation.
#' @param maxN1 integer vector. Upper bound for initial population size 
#' (per age class) to use in initial value simulation. 
#' @param minImm integer. Lower bound for the annual number of immigrants to 
#' use in initial value simulation.  
#' @param maxImm integer. Upper bound for the annual number of immigrants to 
#' use in initial value simulation. 
#' @param z_inits matrix of initial values for mark-recapture data latent states.
#' @param niter integer. Number of MCMC iterations (default = 20000)
#' @param nthin integer. Thinning factor (default = 1)
#' @param nburn integer. Number of iterations to discard as burn-in (default = 5000)
#' @param nchains integer. Number of chains to run (default = 3).
#' @param testRun logical. If TRUE, sets up for a test run with 10 iterations,
#' no thinning, and no burn-in (default = FALSE)
#' @param initVals.seed integer. Seed to use for inital value simulation.
#'
#' @return list of list containing all components necessary for running model 
#' with `nimble::nimbleMCMC()`
#' @export
#'
#' @examples

setupModel <- function(modelCode,
                       nim.data, nim.constants,
                       minN1, maxN1, minImm, maxImm,
                       z_inits,
                       niter = 20000, nthin = 1, nburn = 5000, nchains = 3,
                       testRun = FALSE, initVals.seed){
  
  
  ## Set parameters to monitor in all model versions
  params <- c(
    "Mu.mH", "betaHP.mH", "sigma.mH", 
    "Mu.mO", "betaY.mO", "betaRC.mO", "betaG.mO", "betaSI.mO", "sigma.mO",  
    "Mu.S0", "betaRC.m0", "betaSI.m0", "betaY.m0", "sigma.m0", 
    "mean.rho", "a.eff1", "betaRC.rho", "betaSI.rho", "betaY.rho", "sigma.rho", 
    "par.a", "par.b", "par.c", "betaY.Psi", "betaRC.Psi", "betaSI.Psi", "sigma.Psi", 
    "N.tot", "N", "R.tot", "B.tot", "B", "meanLS", 
    "avgImm", "sigma.Imm", "Imm", "ImmT", 
    "epsilon.mH", "epsilon.mO", "epsilon.m0", "epsilon.Psi", "epsilon.rho",
    "m0t", "mH", "mO", "Psi", "rho",
    "u.Dens")
  
  ## Add additional parameters to monitor depending on model version
  # if(Condition){
  #   params <- c(params, "ExtraParam")
  # }
  
  ## Simulate initial values
  set.seed(initVals.seed)
  initVals <- list()
  for(c in 1:nchains){
    initVals[[c]] <- simulateInitVals(nim.data = nim.data,
                                      nim.constants = nim.constants,
                                      z_inits = z_inits,
                                      minN1 = minN1, maxN1 = maxN1,
                                      minImm = minImm, maxImm = maxImm)
    # initVals[[c]] <- simulateInitVals_original(z_inits = z_inits,
    #                                            Tmax = nim.constants$Tmax,
    #                                            Amax = nim.constants$Amax)
  }
  
  ## Adjust MCMC parameters if doing a test run
  if(testRun){
    niter <- 50
    nthin <- 1
    nburn <- 0
  }
  
  ## Collate model setup variables in a list
  setup <- list(
    modelCode = modelCode,
    modelParams = params,
    initVals = initVals,
    mcmcParams = list(niter = niter, nthin = nthin, 
                      nburn = nburn, nchains = nchains)
  )
  
  return(setup) 
}
