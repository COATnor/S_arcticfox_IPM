#' Format and collate goose data
#'
#' @param goose.data.raw data frame containing goose nestfate data downloaded
#' from the COAT data protal. 
#'
#' @returns a dataframe containing covariate data on goose nest success (proportion of nests hatching successfully). 
#' @export
#'
#' @examples
#' 
reformatData_gooseCov <- function(goose.data.raw){
  
  ## Select relevant columns
  goose.data <- goose.data.raw[[1]] %>%
    dplyr::rename(year = v_year) %>%
    dplyr::filter(sn_locality %in% c("Sassendalen", "Adventdalen")) %>%
    dplyr::mutate(nestfate = dplyr::case_when(v_nestfate == "hatched" ~ 1,
                                              v_nestfate == "unsuccessful" ~ 0, 
                                              TRUE ~ NA)) %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(nest_success = mean(nestfate, na.rm = TRUE))
  
  ## Return data
  return(goose.data)
  
}