#' Format and collate sea ice data
#'
#' @param seaIce.data.raw data frame containing summarised sea ice date derived
#' from MET's daily sea ice maps. The original data are polygons and have been 
#' processed with python (script TBA) by personell at NP.
#'
#' @returns a dataframe containing annual averages of very close drift ice and 
#' fast ice for the periods April-May and January-June.
#' @export
#'
#' @examples
 
reformatData_seaIceCov <- function(seaIce.data.raw){
  
  ## Summarise data into annual covariates
  
  # 1) April-May: Period where foxes are expected to prey on seal pups
  AprMay.sum <- seaIce.data.raw %>%
    dplyr::filter(month %in% c(4, 5)) %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(
      IFSI.meanAprMay = mean(average),
      .groups = "drop"
    )
  
  # 2) January-June: Entire season
  JanJun.sum <- seaIce.data.raw %>%
    dplyr::filter(month %in% c(1:6)) %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(
      IFSI.meanJanJun = mean(average),
      .groups = "drop"
    )
  
  ## Combine measures into one data frame
  IFSI.data <- merge(AprMay.sum, JanJun.sum, by = 'year', all = T)
  
  ## Return data
  return(IFSI.data)
}