#' Format and collate reindeer carcass data
#'
#' @param reindeer.data.raw data frame containing reindeer carcass data 
#' downloaded from COAT data portal. The earliest data are available for 2012.
#' @param reindeer.data.old data frame containing the reindeer carcass data 
#' time series used in the original analysis (Nater et al. 2021). Includes data
#' from 1997 - 2019.
#' @param prioritiseCOAT logical. If TRUE, the function will prioritise the COAT 
#' data over the old data. If FALSE, available values in the original time series
#' will be given priority. 
#'
#' @returns a dataframe containing covariate data on availability of reindeer
#' carcasses. 
#' @export
#'
#' @examples
#' 
reformatData_reindeerCov <- function(reindeer.data.raw, reindeer.data.old, prioritiseCOAT = TRUE){
  
  ## Summarize COAT data -> calculate counts of carcasses per year
  reindeer.counts <- reindeer.data.raw %>%
    dplyr::mutate(t_date = as.Date(t_date, format = "%Y-%m-%d"),  # convert to Date
                  year = lubridate::year(t_date)) %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(n = n(), .groups = "drop") %>%
    dplyr::filter(!is.na(year))
  
  ## Combine old and new data
  reindeer.data <- reindeer.data.old %>%
    dplyr::select(year, reindeer_carcass) %>%
    dplyr::full_join(reindeer.counts, by = "year") %>%
    dplyr::rename(reindeer_carcass_old = "reindeer_carcass",
                  reindeer_carcass_COAT = "n")
  
  ## Define variable for use in analyses
  if(prioritiseCOAT){
    
    reindeer.data <- reindeer.data %>%
      dplyr::mutate(reindeer_carcass = dplyr::case_when(
        !is.na(reindeer_carcass_COAT) ~ reindeer_carcass_COAT,
        is.na(reindeer_carcass_COAT) & !is.na(reindeer_carcass_old) ~ reindeer_carcass_old,
        TRUE ~ NA_real_))
    
  }else{
    
    reindeer.data <- reindeer.data %>%
      dplyr::mutate(reindeer_carcass = dplyr::case_when(
        !is.na(reindeer_carcass_old) ~ reindeer_carcass_old,
        is.na(reindeer_carcass_old) & !is.na(reindeer_carcass_COAT) ~ reindeer_carcass_COAT,
        TRUE ~ NA_real_))
  }
  
  ## Return data
  return(reindeer.data)
}
