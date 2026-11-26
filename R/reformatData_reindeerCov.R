#' Format and collate reindeer carcass data
#'
#' @param reindeer.data.raw data frame containing reindeer carcass data 
#' downloaded from COAT data portal.
#'
#' @returns a dataframe containing covariate data on availability of reindeer
#' carcasses. 
#' @export
#'
#' @examples
#' 
reformatData_reindeerCov <- function(reindeer.data.raw){
  
    ## Extract the old format and new format datasets
    reindeer.newFormat <- reindeer.data.raw[[1]]
    reindeer.oldFormat <- reindeer.data.raw[[2]] #Data values missing from data portal
    
    ## Merge old and new formats together
    reindeer.data <- reindeer.oldFormat %>%
      dplyr::bind_rows(reindeer.newFormat)
    
    ## Summarize COAT data -> calculate counts of carcasses per year
    reindeer.data <- reindeer.data %>%
      dplyr::rename(year = t_year) %>%
      dplyr::group_by(year) %>%
      dplyr::summarise(reindeer_carcass = n(), .groups = "drop") %>%
      dplyr::filter(!is.na(year))
    
  ## Return data
  return(reindeer.data)
}
