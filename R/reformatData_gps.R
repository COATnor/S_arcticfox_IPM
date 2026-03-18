reformatData_gps <- function(gps.data, gps.metadata, 
                             minYear, nMonths_aggregate = 1,
                             area_selection){
  
  
  ## Assign "mark-recapture" observation state
  # 1 = confirmed alive (based on criteria)
  # 2 = confirmed dead (based on criteria)
  # 3 = no signal due to collar malfunction
  # 4 = signal but unknown alive/dead state (failure to assign based on criteria)
  
  gps.data <- gps.data %>%
    dplyr::mutate(obs_state = dplyr::case_when(
      alive == "Maybe" ~ 4,
      is.na(num_data_points) ~ 3,
      alive == "N" ~ 2,
      alive == "Y" ~ 1
    ))
  
  ## Optional: subset based on area 
  
  ## Assign sessions based on specified month aggregation
  
  ## Write multi-state capture histories
  
  ## Assemble auxiliary data (age class, sex, season, area)
  
  ## Generate latent state initial values
  
  ## Collate data and return as list
}
