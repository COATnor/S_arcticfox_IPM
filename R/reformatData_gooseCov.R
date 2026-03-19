#' Format and collate goose data
#'
#' @param goose.data.ipm data frame containing estimates from Jesper Madsen's
#' IPM for pink-footed geese.
#'
#' @returns a dataframe containing covariate data on goose abundance, juvenile
#' proportion, and brood size. 
#' @export
#'
#' @examples
#' 
reformatData_gooseCov <- function(goose.data.ipm){
  
  ## Select relevant columns
  goose.data <- goose.data.ipm %>%
    dplyr::select(year, goose_count, prop_goose_juv, goose_brood_size)
  
  ## Return data
  return(goose.data)
  
}