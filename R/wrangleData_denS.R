#' Loading and reformatting den survey data
#'
#' @param datapath character string with the path to the file containing the
#' data. 
#' @param minYear integer. The first year to consider in the analyses. 
#' @param maxYear integer. The last year to consider in the analyses.'
#' 
#' @return a list containing vector-format data on numbers of monitored and 
#' occupied dens, as well as observed numbers of pups across years. 
#' @export
#'
#' @examples
#' 
wrangleData_denS <- function(datapath, minYear, maxYear){
  
  ## Load data file
  denData <- readr::read_csv(datapath)
  
  ## Remove non-core area and secondary dens
  denS <- denData %>%
    dplyr::filter(CoreArea == 1 & MainDen == 1)
  
  ## Transform data into longitudinal format
  year_cols_idx <- grep("^(19|20)[0-9]{2}$", names(denS))
  year_cols_names <- colnames(denS)[year_cols_idx]
  
  denS2 <- denS %>%
    dplyr::select(ID, dplyr::all_of(year_cols_names)) %>%
    reshape2::melt(id.vars = c("ID"), measure.vars = year_cols_names) %>%
    dplyr::rename(Year = "variable", NoPups = "value")

  ## Add info on breeding activity and monitoring
  denS2 <- denS2 %>%
    dplyr::mutate(Year = as.numeric(as.character(Year)),
                  Breeding = ifelse(is.na(NoPups), NA, ifelse(NoPups > 0, 1, 0)),
                  Monitoring = ifelse(is.na(NoPups), 0, 1))

  ## Remove the year with missing data (1996) and subset to relevant time period
  denS2 <- denS2 %>%
    dplyr::filter(Year != 1996) %>%
    dplyr::filter(between(Year, minYear, maxYear))
  
  ## Calculate yearly summaries (no. of dens monitored, no. of dens occupied, mean no. of pups observed)
  denS_sum <- denS2 %>%
    dplyr::group_by(Year) %>%
    dplyr::summarise(TotMon = sum(Monitoring, na.rm = TRUE),
                     NoOcc = sum(Breeding, na.rm = TRUE),
                     MeanNoPups = mean(NoPups, na.rm = TRUE))
  
  ## Collate data on number of pups observed
  denS_pups <- denS2 %>%
    dplyr::filter(NoPups > 0) %>%
    dplyr::mutate(time = Year - minYear + 1)
  
  
  ## Return data
  return(list(TmaxD = nrow(denS_sum),
              NoOcc = denS_sum$NoOcc, 
              NoMon = denS_sum$TotMon, 
              NoPups = denS_pups$NoPups, 
              DS_year = denS_pups$time, 
              X3 = nrow(denS_pups),
              k.Dens = nrow(denS)))
  
}