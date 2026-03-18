makeCHs_fromLongitudinal <- function(data, session_name, id_name, state_name, age_name){
  
  data <- as.data.frame(data)
  
  Ncapture <- nrow(data)  
  Nanimal <- length(unique(data[, id_name]))  
  occasion <- data[, session_name] - min(data[, session_name]) + 1
  
  msCH <- matrix(0, Nanimal, max(occasion))
  ageCH <- matrix(0, Nanimal, max(occasion))
  
  for (i in 1:Ncapture) {
    msCH[data[i, id_name], occasion[i]] <- data[i, state_name]
    ageCH[data[i, id_name], occasion[i]] <- data[i, age_name]
  }
  
  return(list(ch = msCH, firstyear = ageCH))
}