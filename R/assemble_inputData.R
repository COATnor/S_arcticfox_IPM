#' Assemble demographic and environmental data for running the IPM
#'
#' @param Amax integer. Number of age classes to consider in analyses.
#' @param Tmax integer. The number of years to consider in analyses.
#' @param minYear integer. First year to consider in analyses.
#' @param AaH.data list containing formatted age-at-harvest data.
#' @param rep.data list containing formatted reproduction data from obductions.
#' @param cmrr.data list containing formatted capture-recapture-recovery data.
#' @param denS.data list containing formatted den survey data.
#' @param reindeer.data data frame containing reindeer carcass covariate data per year.
#' @param goose.data data frame containing goose covariate data per year.
#' @param seaIce.data data frame containing sea ice covariate data per year.
#' @param info.priors list containing constants needed for setting informative priors.
#' @param logN_Imm logical. Model toggle for immigration model.
#' @param YearInfo data frame summarizing years in analysis and corresponding indices. 
#' @param save logical. If TRUE, saves assembled data as an .rds file in the 
#' working directory. Default = FALSE. 
#'
#' @return a list containing all data necessary for running the IPM. 
#' @export
#'
#' @examples

assemble_inputData <- function(Amax, Tmax, minYear,
                               AaH.data, rep.data, cmrr.data, denS.data,
                               reindeer.data, goose.data, seaIce.data,
                               info.priors, logN_Imm, YearInfo,
                               save){
  
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
  # --> Already contains only relevant years (filtered in reformatting function). 

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
  IPM.data <- list(
    C = C,
    y = cmrr.data$rCH, 
    P1 = P1$P1,
    P2 = P2$P2,
    NoOcc = denS.data$NoOcc, 
    NoMon = denS.data$NoMon, 
    NoPups = denS.data$NoPups, 
    RdCarcass = RdCarcass, 
    GooseRep = GooseRep,
    #JanJunSeaIceIsfj = input.data.old$JanJunSeaIceIsfj,
    SeaIceIsfj = SeaIceIsfj,
    HPeriod = c(rep(0, 12), rep(1, 11)))
  
  
  IPM.constants <- list(
    Amax = Amax, 
    Juv = c(1, rep(0, Amax - 1)),
    Tmax = Tmax, 
    TmaxC = C.end, 
    TmaxD = denS.data$TmaxD,
    n.occasions = ncol(cmrr.data$rCH), 
    nind = nrow(cmrr.data$rCH), 
    first = cmrr.data$f, 
    #firstyear = cmrr.data$firstyear, 
    lifestage = cmrr.data$lifestage,
    X1 = nrow(P1),  
    P1_age = P1$P1_age, 
    P1_year = P1$P1_year,
    X2 = nrow(P2), 
    P2_age = P2$P2_age, 
    P2_year = P2$P2_year,
    X3 = denS.data$X3,
    DS_year = denS.data$DS_year, 
    k.Dens = denS.data$k.Dens,
    pAgeSex = pAgeSex, 
    pLoc = pLoc,
    logN_Imm = logN_Imm)
  
  
  ## Add relevant prior information
  IPM.constants <- c(IPM.constants, info.priors)
  
  ## Combine data and constants in a list
  inputData <- list(IPM.data = IPM.data, 
                    IPM.constants = IPM.constants)
  
  ## Save (optional) and return data
  if(save){
    saveRDS(inputData, "inputData_formatted.rds")
  }
  
  return(inputData)
  
}