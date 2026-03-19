
reformatData_cmrr <- function(cmrr.data.raw, minYear, maxYear, area_selection){
 
  ## Check minYear and stop if it's earlier than the defined CMRR sessions
  if(minYear < 1990){
    stop("minYear is < 1990. 
         CMRR data is available for 1984-1989, but capture sessions for that time period are not defined in the data.
         The current version of this function relies on defined capture session and therefore does not work for data before 1990.
         Support for earlier data will be added later.")
  }
  
  ## Remove NA events and recode events as states
  # alive (marking/recapture) = 1
  # dead from harvest = 2
  # dead from other causes = 3
  # alive (resighting) = 4
  # alive (telemetry) = 5
  fox <- cmrr.data.raw %>%
    dplyr::filter(!is.na(Event)) %>%
    dplyr::mutate(state = dplyr::case_when(
      Event %in% c("marking", "recapture") ~ 1,
      Event == "recovery" & Death_cause == "trapping" ~ 2,
      Event == "recovery" ~ 3,
      Event == "resighting" ~ 4,
      Event %in% c("telemetry", "last telemetry") ~ 5
    ))

  ## Remove resighting and telemetry events
  fox <- fox %>%
    dplyr::filter(state <= 3)
  
  ## Assign and filter on individual's marking location
  mark.outside <- fox %>%
    dplyr::filter(Event == "marking") %>%
    dplyr::mutate(mark_out = ifelse(Region %in% area_selection, FALSE, TRUE)) %>%
    dplyr::filter(mark_out)
  
  fox <- fox %>%
    dplyr::filter(!(Tag_ID %in% mark.outside$Tag_ID))
  
  ## Remove within-session recaptures (keep only one recapture per year)
  fox <- fox %>% 
    dplyr::filter(!is.na(Session_CMRR_adj))
  
  ## Assign an integer ID for each individual
  indIDs <- data.frame(Tag_ID = unique(fox$Tag_ID), no = 1:length(unique(fox$Tag_ID)))
  fox <- fox %>%
    dplyr::left_join(indIDs, by = "Tag_ID") 
  
  ## Add a column for life stage (juvenile = 1, adult = 2)
  fox <- fox %>%
    dplyr::mutate(lifestage = ifelse(Year_first == 1, 1, 2))
  
  ## Transform longitudinal data to multistate capture histories
  CH.data <- makeCHs_fromLongitudinal(data = fox, 
                                      session_name = "Session_CMRR_adj", 
                                      id_name = "no", 
                                      state_name = "state", 
                                      age_name = "lifestage")
  
  CH <- CH.data$ch
  lifestage <- CH.data$age.ch
  
  ## Compute vector with occasion of first capture
  get.first <- function(x) min(which(x!=0))
  f <- apply(CH, 1, get.first)
  
  ## Complete "lifestage" matrix
  # (all 0's after first capture are set to 2. Entries before first capture remain 0)
  for(i in 1:dim(lifestage)[1]){
    
    for(t in (f[i]+1):dim(lifestage)[2]){
      if(lifestage[i,t] == 0){
        lifestage[i,t] <- 2
      }
    }
  }
  
  ## Creating different capture histories, auxiliary data (known states), and initial values for MCMC
  AF_CMRR_data <- function(CH, f){
    
    # CAPTURE HISTORIES (MAIN DATA) - Recaptures and harvest recoveries
    rrCH <- CH
    rrCH[rrCH%in%c(0,4,5)] <- 3 # --> removing other cause deaths/resightings/telemetry
    
    # CAPTURE HISTORIES (MAIN DATA) - Recaptures and all recoveries
    rrCH2 <- CH
    rrCH2[rrCH2%in%c(0,5)] <- 4 # --> removing resightings/telemetry
    
    # CAPTURE HISTORIES FOR ANALYSIS (MAIN DATA) - Harvest recoveries only
    rCH <- CH
    for(i in 1:dim(rCH)[1]){
      for(t in (f[i]+1):dim(rCH)[2]){
        if(rCH[i,t] == 1){rCH[i,t] <- 3}
      }
    } # --> removing recaptures
    
    rCH[which(rCH==2)] <- 1
    rCH[rCH%in%c(3,4,5)] <- 0 # --> removing other cause deaths/resightings/telemetry
    
    # AUXILIARY DATA ON INDIVIDUAL STATE
    dataCH <- CH 
    
    # Resightings and telemetry = alive (1)
    dataCH[dataCH%in%c(4,5)] <- 1
    
    # Between 1 and 1, 2, or 3 = alive (1)
    mcap <- which(rowSums(dataCH) > 1)
    for(i in 1:length(mcap)){
      last <- max(which(dataCH[mcap[i],]!=0))
      if(last > (f[mcap[i]]+1)){dataCH[mcap[i], (f[mcap[i]]+1):(last-1)] <- 1}
    }
    
    # After a recovery (harvest or other) = dead (4)
    recov <- which(dataCH > 1, arr.ind = T)
    for(i in 1:nrow(recov)){
      if((recov[i,2]+1) <= dim(dataCH)[2]){dataCH[recov[i,1], (recov[i,2]+1):dim(dataCH)[2]] <- 4}
    }
    
    # All remaining 0's are unknown = NA
    dataCH[which(dataCH==0)] <- NA
    
    # Remove the information on the first capture (marking)
    for(i in 1:dim(dataCH)[1]){
      dataCH[i, f[i]] <- NA
    }
    
    # Split matrices for CMRR model with all recoveries, or harvest only
    dataCH2 <- dataCH
    dataCH[which(dataCH==4)] <- 3
    
    # Make a separate matrix for the mark-recovery only model
    BdataCH <- dataCH
    BdataCH[which(BdataCH==2)] <- 0
    BdataCH[which(BdataCH==3)] <- 0
    
    # INITIAL VALUES - based on auxiliary data
    initsCH <- dataCH2
    
    # Write initial values for remaining entries (assumption: individual alive, init = 1)
    for(i in 1:dim(initsCH)[1]){
      for(t in (f[i]+1):dim(initsCH)[2]){
        if(is.na(initsCH[i,t])){initsCH[i,t] <- 1}
      }
    }
    
    # Split intial value matrices for CMRR model with all recoveries, or harvest only
    initsCH2 <- initsCH
    initsCH[which(initsCH==4)] <- 3
    
    # Make a separate initial values matrix for the mark-recovery only model
    BinitsCH <-initsCH
    BinitsCH[which(BinitsCH%in%c(2,3))] <- 0
    
    # RETURN OUTPUT
    output <- list(rrCH = rrCH, 
                   rrCH2 = rrCH2, 
                   rCH = rCH, 
                   dataCH = dataCH, 
                   dataCH2 = dataCH2, 
                   initsCH = initsCH, 
                   initsCH2 = initsCH2, 
                   BdataCH = BdataCH, 
                   BinitsCH = BinitsCH)
    return(output)
  }
  
  
  CMRR.data <- AF_CMRR_data(CH, f)
  
  
  ## Making a second version of the capture histories excluding other-cause recoveries made by Nina (and not included in carcass data)
  # NOTE: With this formulation, the recovery rate r will be the same for the CMRR and age-at-death likelihoods
  
  # --> The foxes also in carcass data are:
  #		- RW137-OG
  #		- WB009-RR
  #		- YY016-YY
  #		- G-Y001 (not included because marked outside core area)
  
  target.ids <- unique(subset(fox, Tag_ID %in% c("RW137-OG", "WB009-RR", "YY016-YY", "G-Y001"))$no)
  
  rrCH3 <- CMRR.data$rrCH2
  rrCH3[which(rrCH3==3)] <- 4
  rrCH3[target.ids,] <- CMRR.data$rrCH2[target.ids,]
  
  CMRR.data$rrCH3 <- rrCH3
  
  
  ## Collate all outputs
  CMRR.data.out <- list(rrCH = CMRR.data$rrCH, 
                        rrCH2 = CMRR.data$rrCH2, 
                        rrCH3 = CMRR.data$rrCH3, 
                        rCH = CMRR.data$rCH, 
                        dataCH = CMRR.data$dataCH, 
                        dataCH2 = CMRR.data$dataCH2, 
                        initsCH = CMRR.data$initsCH, 
                        initsCH2 = CMRR.data$initsCH2, 
                        BdataCH = CMRR.data$BdataCH, 
                        BinitsCH = CMRR.data$BinitsCH,
                        lifestage = lifestage,
                        f = f)
  
  ## Return output
  return(CMRR.data.out)
  
  }
  