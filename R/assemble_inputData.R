#' Assemble demographic and environmental data for running the IPM
#'
#' @param Amax integer. Number of age classes to consider in analyses.
#' @param Tmax integer. The number of years to consider in analyses.
#' @param minYear integer. First year to consider in analyses.
#' @param maxPups integer. Upper prior bound for average litter size.
#' @param uLim.N integer. Upper prior bound for initial number of individuals per age class.
#' @param nLevels.rCov integer. Number of levels of categorical rodent abundance to use.
#' @param standSpec.rCov logical. If TRUE, standardises rodent numbers per species before summing 
#' to offset catchability, If FALSE simple sums alls rodent numbers. 
#' @param poolYrs.genData integer. Whether or not genetic immigration data is pooled across years.
#' @param pImm.type character. Which type of individual-level data to use for immigration. 
#' "original" = p values as output by Geneclass 2. "rescaled" = p values as output
#' by Geneclass 2 and standardized so that the minimum immigrant probability = 0.
#' "LL-based" = log likelihood other / log likelihood other + log likelihood Varanger. 
#' @param uLim.Imm integer. Upper prior bound for annual number of immigrants. 
#' @param wAaH.data a list containing a winter Age-at-Harvest matrix (C) and a vector of
#' yearly proportions of individuals aged/included in Age-at-Harvest data (pData).
#' @param sAaH.data a list containing a summer Age-at-Harvest matrix (C) and a vector of
#' yearly proportions of individuals aged/included in Age-at-Harvest data (pData).
#' @param rep.data a list containing formatted reproduction data in two data 
#' frames: 'P1' (counts) and 'P2' (presences/absences).
#' @param gen.data a list containing relevant data on genetically determined 
#' probabilities of individuals being immigrants.
#' @param pup.data a list containing data on numbers of pups observed on dens. 
#' @param rodent.data a list containing rodent abundance data as a continuous variable (cont),
#' and categorical variable with two (cat2) and three (cat3) levels.
#' @param reindeer.data a list containing reindeer carcass abundance and proportion
#' of foxes with reindeer in stomachs.
#' @param hunter.data a dataframe containing original and scaled counts of successful 
#' hunters per year.
#' @param surv.priors a list of lists containing parameters to define informative priors
#' for early survival, age-specific annual survival, and juvenile/adult natural
#' mortality hazard rate.
#' @param survPriorType a list containing information on prior for annual survival.
#' @param save logical. If TRUE, saves assembled data as an .rds file in the 
#' working directory. Default = FALSE. 
#'
#' @return a list containing all data necessary for running the IPM. 
#' @export
#'
#' @examples

assemble_inputData <- function(Amax, Tmax, minYear,
                               AaH.data = AaH.data,
                               rep.data = rep.data,
                               cmrr.data = cmrr.data,
                               denS.data = denS.data,
                               reindeer.data = reindeer.data,
                               goose.data = goose.data,
                               seaIce.data = seaIce.data,
                               info.priors = info.priors,
                               YearInfo = YearInfo){
  
  ## Select relevant years from observational data
  
  # Set max year
  maxYear <- minYear + Tmax
  
  # Age-at-Harvest data
  C.start <- min(which(grepl(minYear, colnames(AaH.data$C))))
  C.end <- max(which(grepl((maxYear-1), colnames(AaH.data$C))))
  
  C <- AaH.data$C[,C.start:C.end]
  pLoc <- AaH.data$pLoc[C.start:C.end]
  pAgeSex <- AaH.data$pAgeSex[C.start:C.end]
  
  # Reproduction data
  P1 <- subset(rep.data$P1, P1_year %in% 1:Tmax)
  P2 <- subset(rep.data$P2, P2_year %in% 1:Tmax)
  
  # Mark-recovery data
  # --> Already contains only relevant years. 
  # TODO: May want to include an explicit filter in the reformatting function to enforce this also in cases where unorthodox time periods are used. 
  
  # Den survey data
  # --> Already contains only relevant years (filtered in reformatting function).
  
  
  ## Select relevant environmental covariates (incl. time periods)
  RdCarcass_orig <- subset(reindeer.data, year %in% minYear:maxYear)$reindeer_carcass
  GooseRep_orig <- subset(goose.data, year %in% minYear:maxYear)$prop_goose_juv
  SeaIceIsfj_orig <- subset(seaIce.data, year %in% minYear:maxYear)$IFSI.meanJanJun
  
  
  ## Detrend and scale environmental covariates
  RdCarcass_prep <- SpecsVerification::Detrend(RdCarcass_orig)
  GooseRep_prep <- SpecsVerification::Detrend(GooseRep_orig)
  SeaIceIsfj_prep <- SpecsVerification::Detrend(SeaIceIsfj_orig)
  
  RdCarcass <- as.numeric(scale(RdCarcass_prep))
  GooseRep <- as.numeric(scale(GooseRep_prep))
  SeaIceIsfj <- as.numeric(scale(SeaIceIsfj_prep))
  
  envCov.metadata <- data.frame(
    covariate = c("RdCarcass", "GooseRep", "SeaIceIsfj"),
    original_mean = c(mean(RdCarcass_orig, na.rm = TRUE), mean(GooseRep_orig, na.rm = TRUE), mean(SeaIceIsfj_orig, na.rm = TRUE)),
    original_sd = c(sd(RdCarcass_orig, na.rm = TRUE), sd(GooseRep_orig, na.rm = TRUE), sd(SeaIceIsfj_orig, na.rm = TRUE)),
    detrended_mean = c(mean(RdCarcass_prep, na.rm = TRUE), mean(GooseRep_prep, na.rm = TRUE), mean(SeaIceIsfj_prep, na.rm = TRUE)),
    detrended_sd = c(sd(RdCarcass_prep, na.rm = TRUE), sd(GooseRep_prep, na.rm = TRUE), sd(SeaIceIsfj_prep, na.rm = TRUE))
  )
  
  
  ## Make harvest period covariate
  
  
  ## List all relevant data (split into data and constants as used by NIMBLE)
  # Data
  nim.data <- list(
    C_w = C_w,
    pData_w = pData_w,
    
    C_s = C_s,
    pData_s = pData_s,
    
    P1 = P1$P1,
    
    P2 = P2$P2,
    
    NoPups = pup.data$NoPups,
    
    HarvestEffort = hunter.data$NHunters_std,
    RodentAbundance = RodentAbundance,
    RodentAbundance2 = RodentAbundance2,
    RodentIndex = RodentIndex,
    RodentIndex2 = RodentIndex2,
    Reindeer = Reindeer
  )
  
  # Constants
  nim.constants <- list(
    Amax = Amax,
    Tmax = Tmax,
    minYear = minYear,
    
    maxPups = maxPups,
    uLim.N = uLim.N,
    uLim.Imm = uLim.Imm,
    
    XsH = XsH,
    sH_year = sH_year,
    
    P1_age = P1$age_adj,
    P1_year = P1$RepYearIndex,
    X1 = length(P1$P1),
    
    P2_age = P2$age_adj,
    P2_year = P2$RepYearIndex,
    X2 = length(P2$P2),
    
    NoPups_year = pup.data$NoPups_year,
    X3 = pup.data$X3,
    
    nLevels.rCov = nLevels.rCov
  )
  
  
  ## Add relevant prior information
  nim.constants <- c(nim.constants, surv.priors$earlySurv)
  
  
  ## Combine data and constants in a list
  inputData <- list(nim.data = nim.data, 
                    nim.constants = nim.constants)
  
  ## Save (optional) and return data
  if(save){
    saveRDS("inputData_formatted.rds")
  }
  
  return(inputData)
  
}