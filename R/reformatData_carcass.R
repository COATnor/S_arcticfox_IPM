#' Reformat the carcass data to make Age at harvest matrix & P1var (nr of embryos/scars) & P2var (breeding or not breeding)
#'
#' @param Amax  integer. Number of age classes to consider in analyses.
#' @param minYear  integer. First year to consider in analyses.
#' @param maxYear  integer. Last year to consider in analyses.
#' @param area_selection a vector of study-area sub-area names to consider in the analyses: c("Inner", "BB",  "Tana")
#' @param carcass.dataset dataframe containing the dataset downloaded from the COAT dataportal.
#'
#' @return a list containing the age-at-harvest matrix and dataframes with embryo count (P1var) and placental scar presence-absence (P2var) data.
#' @export
#'
#' @examples


reformatData_carcass <- function (Amax, minYear, maxYear,
                                  area_selection, carcass.data) {
  
  #---------------------------#
  # CONSOLIDATION & FILTERING #
  #---------------------------#
  
  ## Check for and drop any entries with unknown trapping season
  if(any(is.na(carcass.data$Trapseason))){
    numberNA <- length(which(is.na(carcass.data$Trapseason)))
    carcass.data <- subset(carcass.data, !is.na(Trapseason))
    message(paste0("Dropped ", numberNA, " entries with missing trapping season."))
  }
  
  ## Subset data to time period of interest
  carcass.data <- carcass.data %>%
    dplyr::mutate(seasonStart_yr = as.numeric(stringr::str_split_i(Trapseason, pattern = "-", i = 1)),
                  seasonEnd_yr = as.numeric(stringr::str_split_i(Trapseason, pattern = "-", i = 2))) %>%
    dplyr::filter(seasonStart_yr >= minYear & seasonEnd_yr <= maxYear)
  
  ## Add integer trapping session
  carcass.data <- carcass.data %>%
    dplyr::mutate(Session = seasonStart_yr - minYear + 1)
  
  ## Consolidate age (class) information
  carcass.data <- carcass.data %>%
    dplyr::mutate(
      
      # Assign age = 1 to individuals noted as juveniles or with baby teeth
      Age = dplyr::case_when(Life_stage == "juvenile" & is.na(Tooth) ~ 1,
                             Tooth == "1 (0) åpen pulpa" ~ 1, # NOTE: With new versions of data, we need to double-check that there are no "new" spellings of this information
                             TRUE ~ as.numeric(Tooth)),
      
      ## Collapse ages 5 and up into one age class
      AgeClass = ifelse(Age >= Amax, Amax, Age)
    )
  
  ## Drop males
  carcass.dataF <- carcass.data %>%
    dplyr::filter(Sex == "F")
  
  ## Subset to relevant area(s) for AaH data
  if(area_selection == "Advent-/Sassendalen"){
    carcass.dataF_sub <- subset(carcass.dataF, Core_area == 1)
  }else{
    stop("Invalid area_selection provided. 
         For now, the code only handles 'Advent-/Sassendalen' alone.
         Once we move to new data, we should have an option to also do Ny-Ålesund/Kongsfjorden")
  }
  
  #---------------------#
  # AGE-AT-HARVEST DATA #
  #---------------------#

  ## Build age-at-harvest matrix
  AaH.mat <- carcass.dataF_sub %>%
    dplyr::filter(!is.na(AgeClass)) %>%
    dplyr::group_by(AgeClass, seasonStart_yr) %>%
    dplyr::count() %>%
    tidyr::pivot_wider(names_from = seasonStart_yr, values_from = n, values_fill = 0) %>%
    dplyr::ungroup() %>%
    dplyr::select(-AgeClass) %>%
    as.matrix()
  
  dimnames(AaH.mat)[[1]] <- 1:Amax
  dimnames(AaH.mat)[[2]] <- unique(carcass.dataF$Trapseason)

  ## Calculate annual proportions of harvested females that were aged
  aged <- carcass.dataF %>%
    dplyr::mutate(Aged = ifelse(is.na(AgeClass), 0, 1)) %>%
    dplyr::group_by(seasonStart_yr) %>%
    dplyr::summarise(prop = mean(Aged))
  
  pData <- aged$prop
  
  
  #---------------------#
  # PLACENTAL SCAR DATA #
  #---------------------#
  
  ## Filter data
  placScar.data <- carcass.dataF %>%
    
    # Subset to drop entries with missing placental scar or age info
    dplyr::filter(!is.na(Placental_scar) & !is.na(Age)) %>%
    
    # Drop all individuals in age class 1 (juveniles that should not have placental scars)
    dplyr::filter(Age > 1) %>%
    
    # Make binary reproduction variable
    dplyr::mutate(Rep = ifelse(Placental_scar > 0, 1, 0))
  
  ## Select relevant columns for #placental scars dataset (P1)
  P1.data <- placScar.data %>%
    dplyr::filter(Placental_scar > 0) %>%
    dplyr::select(Placental_scar, Age, Session, Core_area) %>%
    dplyr::rename(P1 = Placental_scar, 
                  P1_age = Age,
                  P1_year = Session,
                  P1_core = Core_area)
  
  ## Select relevant columns for reproduction yes/no dataset (P2)
  P2.data <- placScar.data %>%
    dplyr::select(Rep, Age, Session, Core_area) %>%
    dplyr::rename(P2 = Rep, 
                  P2_age = Age,
                  P2_year = Session,
                  P2_core = Core_area)
  
  
  #-----------#
  # COLLATION #
  #-----------#
  
  ## Combine data
  carcassData <- list(AaH.mat = AaH.mat,
                      pData = pData,
                      P1.data = P1.data,
                      P2.data = P2.data)
  
  ## Return data
  return(carcassData)
}
