library(ggplot2)
library(nimble)
library(sf)
library(reshape2)
library(remotes)
library(ckanr)
library(purrr)
library(dplyr)
library(metafor)
library(patchwork)
library(coda)

options(tibble.width = Inf)

#**********#
# 0) SETUP #
#**********#

## Set seed
mySeed <- 10

## Set general parameters
Amax <- 5 # Number of age classes
Tmax <- 23  # Number of years
minYear <- 1997 # First year to consider
area_selection <- "Advent-/Sassendalen" # Areas to consider

## Set COAT dataset names, versions, and directories, and access
reindeer.dataset.name <-"S_ungulates_carcasses_adventdalen_summer_v3"
rodent.dataset.version <- 3

COAT_key <- Sys.getenv("COAT_API")

## Set filepaths for local datasheets
carcass.data.path <- c("Data/200108_carcass_trap.csv")
denSurvey.data.path <- c("Data/200214_DenSurvey.csv")
cmrr.data.path <- c("Data/SvalbardArcticFox_CaptureRecapture_1984-2016.csv")
gps.data.path <- c("Data/Survival_GPS.csv")
gps.metadata.path <- c("Data/Metadata_Survival_GPS.csv")


seaIce.data.path <- c("Data/isfjorden_1966-2019.csv")
goose.data.path <- c("Data/Data/SvalbardTerr_Covariates.csv")

## Source all functions in "R" folder
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir('R')

## Generate plotting directory if not available
if(!dir.exists("Plots")){
  dir.create("Plots")
}

## Set "switches" for running different model versions
# --> TBA as need arises


#*********************#
# 1) DATA PREPARATION #
#*********************#

# 1a) Load and reformat carcass data
#----------------------------------#

## Load carcass data
carcass.data.raw <- readr::read_csv(carcass.data.path)
# --> This gives a warning about some unexpected values that get set to NA. 
# --> It should not be an issue as it is not data that actually gets used in original analyses.
# --> We assume that this will be resolved once data is on COAT data portal. 

## Reformat carcass data
carcass.data <- reformatData_carcass(Amax = Amax, 
                                     minYear = minYear, maxYear = minYear + Tmax - 1,
                                     area_selection = area_selection,
                                     carcass.data.raw = carcass.data.raw)


# 1b) Age-at-Harvest data #
#--------------------------------#

## Extract AaH data from processed carcass data
AaH.data <- list(C = carcass.data$AaH.mat,
                 pLoc = carcass.data$pLoc,
                 pAgeSex = carcass.data$pAgeSex)



# 1c) Reproduction data #
#-----------------------#

## Extract reproduction data from processed carcass data
rep.data <- list(P1 = carcass.data$P1.data,
                 P2 = carcass.data$P2.data)


# 1d) Mark-recapture recovery data #
#----------------------------------#

## Load CMRR data
cmrr.data.raw <- readr::read_csv(cmrr.data.path)

## Prepare CMRR data
cmrr.data <- reformatData_cmrr(cmrr.data.raw = cmrr.data.raw,
                               minYear = minYear, maxYear = minYear + Tmax - 1,
                               area_selection = area_selection)


# 1e) Den survey data #
#---------------------#

## Prepareden survey data
denS.data <- wrangleData_denS(datapath = denSurvey.data.path,
                              minYear = minYear)


# 1f) GPS collar survival data #
#------------------------------#

## Load GPS collar data
gps.data <- readr::read_csv(gps.data.path)
gps.metadata <- readr::read_csv(gps.metadata.path)

# reformatData_gps(gps.data, gps.metadata, 
#                  minYear, nMonths_aggregate = 1,
#                  area_selection)

# 1g) Environmental data #
#------------------------#

## Download rodent data
rodent.data.raw <- downloadData_COAT(COAT_key = COAT_key, 
                                     COATdataset.name = rodent.dataset.name,
                                     COATdataset.version = rodent.dataset.version)

## Reformat rodent data
rodent.data <- reformatData_rodent(rodent.dataset = rodent.data.raw,
                                   minYear = minYear)

## Reformat reindeer data
reindeer.data <- reformatData_reindeer(minYear = minYear,
                                       Tmax = Tmax,
                                       reinCov.VarTana = reinCov.VarTana)


# 1h) Conceptual year information #
#---------------------------------#

YearInfo <- collate_yearInfo(minYear = minYear,
                             Tmax = Tmax)


#**********************#
# 2) PRIOR INFORMATION #
#**********************#

## Parameters/paths for making informative priors for survival based on meta-analysis of literature data
meta.datafile <- "Data/RedFox_LiteratureData.csv"
simulateSD <- TRUE

## Parameters/paths for making informative priors for natural mortality using Tom Porteus' Hoenig model approach
mu.t.max <- 22.61062
hoenig.datafile <- "Data/HoenigMod_Posteriors_fromTomPorteus.txt"
nsim <- 30

## Collate all prior information
surv.priors <- collate_priorInfo(meta.datafile = meta.datafile,
                                 simulateSD = simulateSD,
                                 hoenig.datafile = hoenig.datafile, 
                                 nsim = nsim, 
                                 mu.t.max = mu.t.max, 
                                 maxAge = maxAge_yrs,
                                 S0.mean.offset = S0.mean.offset,
                                 S0.sd.factor = S0.sd.factor)

## Define type of prior to use for annual survival
survPriorType <- definePriorType_AnnSurv(HoenigPrior = HoenigPrior, 
                                         sPriorSource = sPriorSource)

#****************#
# 3) MODEL SETUP #
#****************#

# 3a) Write model code #
#----------------------#

redfox.code <- writeCode_redfoxIPM(indLikelihood.genData = indLikelihood.genData)
#redfox.code <- writeCode_redfoxIndepMod(indLikelihood.genData = indLikelihood.genData)


# 3b) Assemble IPM input data #
#-----------------------------#

input.data <- assemble_inputData(Amax = Amax, 
                                 Tmax = Tmax, 
                                 minYear = minYear,
                                 maxPups = 14,
                                 uLim.N = 800,
                                 uLim.Imm = 3000,
                                 nLevels.rCov = nLevels.rCov,
                                 standSpec.rCov = standSpec.rCov,
                                 poolYrs.genData = poolYrs.genData,
                                 pImm.type = pImm.type,
                                 wAaH.data = wAaH.data, 
                                 sAaH.data = sAaH.data,
                                 rep.data = rep.data, 
                                 gen.data = gen.data,
                                 pup.data = pup.data,
                                 rodent.data = rodent.data, 
                                 reindeer.data = reindeer.data,
                                 hunter.data = hunter.data, 
                                 surv.priors = surv.priors,
                                 survPriorType = survPriorType)


# 3c) Set up for model run (incl. simulating initial values) #
#------------------------------------------------------------#

model.setup <- setupModel(modelCode = redfox.code, 
                          nim.data = input.data$nim.data, 
                          nim.constants = input.data$nim.constants, 
                          minN1 = c(600, 50, 50, 50, 50), 
                          maxN1 = c(800, 400, 400, 400, 400), 
                          minImm = 50, 
                          maxImm = 600,
                          fitCov.mH = fitCov.mH, 
                          fitCov.mO = fitCov.mO,
                          fitCov.Psi = fitCov.Psi, 
                          fitCov.rho = fitCov.rho,
                          fitCov.immR = fitCov.immR,
                          rCov.idx = rCov.idx,
                          mO.varT = mO.varT,
                          HoenigPrior = HoenigPrior,
                          imm.asRate = imm.asRate,
                          testRun = FALSE,
                          initVals.seed = mySeed
)


####################
# 4) MODEL FITTING #
####################

t1 <- Sys.time()
IPM.out <- nimbleMCMC(code = model.setup$modelCode,
                      data = input.data$nim.data, 
                      constants = input.data$nim.constants,
                      inits = model.setup$initVals, 
                      monitors = model.setup$modelParams,
                      nchains = model.setup$mcmcParams$nchains, 
                      niter = model.setup$mcmcParams$niter, 
                      nburnin = model.setup$mcmcParams$nburn, 
                      thin = model.setup$mcmcParams$nthin, 
                      samplesAsCodaMCMC = TRUE, 
                      setSeed = 0)
Sys.time() - t1


saveRDS(IPM.out, file = "RedFoxIPM_main_singleCensus_combHarvest2.rds") 

#MCMCvis::MCMCtrace(IPM.out)


########################
# 5) MODEL COMPARISONS #
########################

## Simplified models
compareModels(Amax = Amax, 
              Tmax = Tmax, 
              minYear = minYear, 
              post.filepaths = c("AF_IPM.rds", 
                                 "AF_IPM_logImm.rds"), 
              model.names = c("Original", 
                              "logNormal Immigration dist."), 
              plotFolder = "Plots/Comp_ImmDistributions")


###########################################
# 6) IPM RESULTS - STUDY PERIOD ESTIMATES #
###########################################

IPM.out <- readRDS("RedfoxIPM_main.rds")
#IPM.out <- readRDS("RedfoxIPM_ModelRun.rds")


## Plot basic IPM outputs (vital rate & population size estimates)
plotIPM_basicOutputs(MCMC.samples = IPM.out,
                     nim.data = input.data$nim.data,
                     Amax = Amax, Tmax = Tmax, minYear = minYear)

## Plot covariate relationships
plotIPM_covariateEffects(MCMC.samples = IPM.out,
                         rCov.idx = rCov.idx,
                         rodentMIN = -1.75, rodentMAX = 4,
                         AgeClass = 1) 

#########################
# 7) FOLLOW-UP ANALYSES #
#########################

## Extract parameter samples
paramSamples <- extractParamSamples(MCMC.samples = IPM.out,
                                    Amax = Amax, Tmax = Tmax)

## Calculate sensitivities and elasticities
sensitivities <- calculateSensitivities(paramSamples = paramSamples,
                                        Amax = Amax)

## Plot sensitivities
plotSensitivities(sensitivities = sensitivities,
                  Amax = Amax)


## Set LTRE options
HazardRates <- TRUE
PopStructure <- TRUE

## Run random design LTRE
randomLTRE <- runLTRE_randomDesign(paramSamples = paramSamples, 
                                   sensitivities = sensitivities, 
                                   Amax = Amax, Tmax = Tmax, 
                                   HazardRates = HazardRates, 
                                   PopStructure = PopStructure)

## Plot results from random design LTRE
plotLTRE_randomDesign(LTRE_results = randomLTRE,
                      Amax = Amax,
                      HazardRates = HazardRates,
                      PopStructure = PopStructure)

## Run fixed design LTRE
fixedLTRE <- runLTRE_fixedDesign_allYears(paramSamples = paramSamples, 
                                          Amax = Amax, Tmax = Tmax, 
                                          HazardRates = HazardRates, 
                                          PopStructure = PopStructure)

## Plot results from fixed design LTRE
plotLTRE_fixedDesign(LTRE_results = fixedLTRE, 
                     Amax = Amax, Tmax = Tmax, minYear = minYear, 
                     HazardRates = HazardRates, 
                     PopStructure = PopStructure)

## Plot decomposition of mO into covariates and random effect
plotVariance_comp_mO(MCMC.samples = IPM.out, 
                     Tmax = Tmax,
                     minYear = minYear)

## Calculate post-hoc parameter correlations to check for signs of density dependence
calculate_p.hoc_param.corr(MCMC.samples = IPM.out, 
                           Tmax = Tmax)
