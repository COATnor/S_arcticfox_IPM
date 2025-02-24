library(plyr)
library(reshape)

###############################
#### 1) MARK-RECOVERY DATA #### 
###############################

#**************************************
# Jan 20 update: 
# Few errors corrected
# Code now removes recoveries/resigthings outside core area
#
# Feb 20 update:
# Code now also makes capture histories (and initial values / known state data) for other cause recoveries
#
#**************************************

## Read in data file
fox <- read.csv('/Users/chloern/Dropbox/Arctic_Fox_SUSTAIN/Data/Tagging/200110_AF_tagging_red.csv', sep = ',')
fox <- subset(fox, !is.na(Event))

## Coding events as numbers
# alive (marking/recapture) = 1
# dead from harvest = 2
# dead from other causes = 3
# alive (resighting) = 4
# alive (telemetry) = 5
unique(fox$Event)

fox$state <- NA

for(i in 1:nrow(fox)){
  
  if(fox$Event[i] %in% c('marking', 'recapture')){fox$state[i] <- 1}
  
  if(fox$Event[i] == 'recovery'){
    if(!is.na(fox$Death_cause[i]) & fox$Death_cause[i] == 'trapping'){fox$state[i] <- 2} else {fox$state[i] <- 3}
  }
  
  if(fox$Event[i] == 'resighting'){fox$state[i] <- 4}
  
  if(fox$Event[i]%in%c('telemetry', 'last telemetry')){fox$state[i] <- 5}
}

table(fox$Event)

## Remove individuals that were not marked in the core area
out.ids <- unique(subset(fox, Event=='marking' & Core_area!=1)$Tag_ID)
fox$mark.out <- ifelse(fox$Tag_ID%in%out.ids, 1, 0)
fox <- subset(fox, mark.out == 0)
table(fox$Event)
table(fox$Event, fox$Core_area)

## Remove recoveries and resightings outside the core area
fox <- subset(fox, Core_area == 1)
table(fox$Event, fox$Core_area)


## Remove non-informative sessions (marked as X)
fox <- subset(fox, Session_adj != 'X')
fox$Session <- as.numeric(as.character(fox$Session_adj))
table(fox$Event)

# NOTE: Here I loose most of the recaptures (as the majority happens in the marking year)
table(fox$Event, fox$Core_area)


## Assigning an id no. that starts at 1
id.rank = function(x){
  
  n = length(x)
  pp=c(1,rep(0,n-1))
  for(i in 2:n){
    if(x[i]==x[i-1]) pp[i]=pp[i-1]
    else pp[i]=pp[i-1]+1}
  pp}

fox = fox[order(fox$Tag_ID, fox$Session),] 
fox$no = id.rank(fox$Tag_ID)

## Function for transforming longitudinal data into multistate capture histories
fox.data=function(){
  
  # Function to create the capture histories and age at occasion (first year or not) for each individual
  cmr.sort = function(data,t,no,state,firstyear){   # t=time, no=animal number
    
    Ncapture = nrow(data)  #counts the number of capture entries
    Nanimal = data[Ncapture,no]  
    occasion = data[,t]-min(data[,t])+1
    
    cmr = matrix(0,Nanimal,max(occasion))
    first = matrix(0,Nanimal,max(occasion))
    
    for (i in 1:Ncapture) {
    	cmr[data[i,no],occasion[i]] <- data[i,state]
    	first[data[i,no],occasion[i]] <- data[i,firstyear]
    }
    return(list(cmr=cmr,first=first))
  }
  
  fox.data = cmr.sort(fox,"Session","no","state","firstyear")
  fox.cmr = fox.data$cmr
  fox.age = fox.data$first
  return(list(ch = fox.cmr, firstyear = fox.age))
}

## Write multistate capture histories
CH <- fox.data()$ch
firstyear <- fox.data()$firstyear

## Remove NA's from firstyear
firstyear[is.na(firstyear)] <- 0

## Compute vector with occasion of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

## Creating capture histories, auxiliary data (known states), and initial values for MCMC
AF_CMRR_data <- function(CH, f){
  
  # CAPTURE HISTORIES (MAIN DATA) - Recaptures and harvest recoveries
  rrCH <- CH
   rrCH[rrCH%in%c(0,4,5)] <- 3 # --> removing other cause deaths/resightings/telemetry
   
  # CAPTURE HISTORIES (MAIN DATA) - Recaptures and all recoveries
  rrCH2 <- CH
  rrCH2[rrCH2%in%c(0,5)] <- 4 # --> removing resightings/telemetry
  
  # CAPTURE HISTORIES FOR JAGS (MAIN DATA) - Harvest recoveries only
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
 output <- list(rrCH = rrCH, rrCH2 = rrCH2, rCH = rCH, dataCH = dataCH, dataCH2 = dataCH2, initsCH = initsCH, initsCH2 = initsCH2, BdataCH = BdataCH, BinitsCH = BinitsCH)
  return(output)
}

  
CMRR.data <- AF_CMRR_data(CH, f)


## Making a second version of the capture histories excluding other-cause recoveries made by Nina (and not included in carcass data)
# NOTE: With this formulation, the recovery rate r will be the same for the CMRR and age-at-death likelihoods

# --> The foxes also in carcass data are:
#		- No 69: RW137-OG
#		- No 76:  WB009-RR
#		- No 142: YY016-YY
#		- G-Y001 (not included because marked outside core area)
		
rrCH3 <- CMRR.data$rrCH2
rrCH3[which(rrCH3==3)] <- 4
rrCH3[69,] <- CMRR.data$rrCH2[69,]
rrCH3[76,] <- CMRR.data$rrCH2[76,]		
rrCH3[142,] <- CMRR.data$rrCH2[142,]

CMRR.data$rrCH3 <- rrCH3


#############################
#### AGE-AT-HARVEST DATA ####
#############################

#**************************************
# Nov 19 update: 
# New carcass data added 
#
# Jan 20 update:
# 1 error corrected
# Age estimated for some individuals based on characteristics (new column tooth_from_estimate)
# Missing core_area info assigned based on trapper name where possible 
# Code now also reads in proportions of carcasses with unknown age/sex or location
#**************************************

## Read in data file
carcass <- read.csv('/Users/chloern/Dropbox/Arctic_Fox_SUSTAIN/Data/Carcasses_Trap/Raw_data/200108_carcass_trap.csv')
# --> Required: Age-at-harvest matrix

## Assign age = 1 to individuals noted as juveniles or with baby teeth
carcass$Tooth[which(carcass$Life_stage == 'juvenile' & is.na(carcass$Tooth))] <- 1
carcass$Tooth[which(carcass$Tooth == '1/0 apen pulpa')] <- 1

## Remove all entries that were trapped in the 1996-97 period
carcass2 <- subset(carcass, Trapseason!='1996-1997')
carcass2$Tooth <- as.numeric(as.character(carcass2$Tooth))
unique(sort(carcass2$Tooth, na.last = TRUE))

# --> I drop the 1996-97 season because it's the first season, representing the population one year before the current analysis starts. For now, it's fine to drop this, but it may become useful to include if I decide to go for the hindcast. 

## Collapse ages 5 and up into one age classe
carcass2$AgeClass <- ifelse(carcass2$Tooth >= 5, 5, carcass2$Tooth)

## Assign the year
# --> trapping season t1-t2 gives information about survival from t1 to t2, and will therefore be assigned to t1 (= 1997)
str(carcass2$Trapseason)
carcass2$Session <- as.numeric(carcass2$Trapseason) - 1

## Data availability
table(carcass2$Core_area, carcass2$Sex, useNA = 'always')
# --> There are 17 individuals trapped in the core area and 1 in an unknown area (total = 18) that are of unknown sex 
# --> There are also 15+8 = 23 individuals for which we don't know if they were trapped in the core area or not
table(carcass2$Core_area, carcass2$Sex, useNA = 'always')

## Remove all entries with ambiguous age, sex, location and extract age-at-harvest matrix
carcass_coreAS <- subset(carcass2, Core_area == 1 & Tooth %in% c(1:16) & Sex == 'F')
C <- as.matrix(table(carcass_coreAS[,c('AgeClass','Session')]))
dimnames(C) <- NULL
C

## Load data on proportions of carcasses missing age/sex and location information
load('200131_MissingDataProps.RData')


#############################
#### PLACENTAL SCAR DATA ####
#############################

plac_scar <- subset(carcass2, !is.na(Placental_scar) & !is.na(AgeClass))
# --> Required: A) presence/absence of placental scars (individual-, age- and year-specific)
#               B) number of placental scars given presence (individual-, age- and year-specific)

## Assign the year
# --> placental scars from trapping season t1-t2 give information about reproduction in t1 (= 1997)
plac_scar$Session <- as.numeric(plac_scar$Trapseason) - 1

## Assign an individual id starting at 1
plac_scar$id <- c(1:nrow(plac_scar)) # this works, as every individual only appears once

## Remove unneccessary information
plac_scar2 <- plac_scar[,c('id','Session','Placental_scar', 'AgeClass', 'Core_area')]

## Add presence/absence data for placental scars
plac_scar2$PS_presabs <- ifelse(plac_scar2$Placental_scar == 0, 0, 1)

## Split data into vectors for passing them to JAGS
# Number of placental scars
P1 <- subset(plac_scar2, Placental_scar > 0)$Placental_scar

# Presence/absence of placental scars
P2 <- plac_scar2$PS_presabs

# Years and ages associated with placental scars
P1_age <- subset(plac_scar2, Placental_scar > 0)$AgeClass
P1_year <- subset(plac_scar2, Placental_scar > 0)$Session

P2_age <- plac_scar2$AgeClass
P2_year <- plac_scar2$Session

# Areas associated with placental scars
P1_core <- subset(plac_scar2, Placental_scar > 0)$Core_area
P2_core <- plac_scar2$Core_area

CarcassLoc <- list(P1_core = P1_core, P2_core = P2_core)

save(CarcassLoc, file = '200131_CarcassLocations.RData')


###########################
#### AGE-AT-DEATH DATA ####
###########################

#**************************************
# Feb 20 update:
# Added this section to also make an age-at-death matrix for the non-harvest carcasses
#**************************************

## Read in data file
carcassN <- read.csv('/Users/chloern/Dropbox/Arctic_Fox_SUSTAIN/Data/Carcasses_Natural/Raw_data/200228_Carcass_natural_red.csv')
# --> Required: Age-at-death matrix

## Assign age = 1 to individuals noted as pups, juveniles or with baby teeth
carcassN$Tooth[which(carcassN$Life_stage%in% c('juvenile', 'pup') & is.na(carcassN$Tooth))] <- 1
carcassN$Tooth[which(carcassN$Tooth%in% c('1 (0) åpen pulpa', '<1'))] <- 1
carcassN$Tooth <- as.numeric(as.character(carcassN$Tooth))
unique(sort(carcassN$Tooth, na.last = TRUE))

## Data availability
table(carcassN$Core_area, carcassN$Sex, useNA = 'always')
table(carcassN$Core_area, carcassN$Tooth, useNA = 'always')
table(carcassN$Sex, carcassN$Tooth, useNA = 'always')
# --> 2 individuals in the core area, 8 outside, and 1 in an unknown location are of unknown sex
# --> There are also 3 individuals for which we don't know if they were trapped in the core area or not
# --> Half of the missing sex individuals are also missing age (missing info not independent, as for harvest data)

## Calculate proportions missing sex&age of location data for other cause carcasses
pAgeSex.M <- nrow(subset(carcassN, !is.na(Sex) & !is.na(Tooth)))/nrow(carcassN)
pLoc.M <-  length(which(!is.na(carcassN$Core_area)))/nrow(carcassN)

## Remove all entries for which death season is unavailable, or prior to 1997-1998
# NOTE: So far, we only defined death season for carcasses from the core area (or with unknown location)
carcassN2 <- subset(carcassN, !(DeathSeason%in%c('1995-1996', '1996-1997', 'TBA', NA)))

# --> I drop the seasons prior to 1997 because that is the first season, representing the population one year before the current analysis starts. For now, it's fine to drop this, but it may become useful to include if I decide to go for the hindcast. 

## Collapse ages 5 and up into one age classe
carcassN2$AgeClass <- ifelse(carcassN2$Tooth >= 5, 5, carcassN2$Tooth)

## Assign the year
# --> death season t1-t2 gives information about survival from t1 to t2, and will therefore be assigned to t1 (= 1997)
str(carcassN2$DeathSeason)
seasonsM <- data.frame(DeathSeason = paste(c(1997:2018),'-',c(1998:2019), sep = ''), Session = c(1:22))
carcassN2 <- merge(carcassN2, seasonsM, by = 'DeathSeason', all.x = T)


## Reformat data to extract a complete table
carcassN2$Session <- as.factor(carcassN2$Session)
carcassN2$AgeClass <- as.factor(carcassN2$AgeClass)

## Remove all entries with ambiguous age, sex, location and extract age-at-death matrix
carcassN_coreAS <- subset(carcassN2, Core_area == 1 & !is.na(AgeClass) & Sex == 'F')
M <- as.matrix(table(carcassN_coreAS[,c('AgeClass','Session')]))
M

## Add "missing" years (when no carcasses were delivered)
m1 <- m5 <- m8 <- m15 <- rep(0, 5)
M <- cbind(m1, M[,(2:4)-1], m5, M[,(6:7)-2], m8, M[,(9:14)-3], m15, M[,(16:20)-4])
dim(M)

dimnames(M) <- NULL
M
# --> There is very, very little information in this...


#########################
#### DEN SURVEY DATA ####
#########################

#**************************************
# Feb 20 update:
# Added new data (up to and incl. 2019)
#**************************************

## Read in data file
denS <- read.csv('/Users/chloern/Dropbox/Arctic_Fox_SUSTAIN/Data/Den_Survey/200214_DenSurvey.csv', check.names = FALSE)

## Remove non-core area and secondary dens
denS <- subset(denS, CoreArea == 1 & MainDen == 1)

## Transform data into longitudinal format
denS2 <- melt(denS[,c(1,5:28)],id.vars = c("ID"), measure.vars = as.character(c(1996:2019)))
names(denS2) <- c('ID','Year','NoPups')

## Add info on breeding activity and monitoring
denS2$Breeding <- ifelse(is.na(denS2$NoPups), NA, ifelse(denS2$NoPups>0, 1, 0))
denS2$Monitoring <- ifelse(is.na(denS2$NoPups), 0, 1)

## Remove the year with missing data (1996)
denS2 <- subset(denS2, Year!=1996)

## Calculate yearly summaries (no. of dens monitored, no. of dens occupied, mean no. of pups observed)
denS_sum <- ddply(denS2, .(Year), summarise,
                  TotMon = sum(Monitoring, na.rm = T),
                  NoOcc = sum(Breeding, na.rm = T),
                  MeanNoPups = mean(NoPups, na.rm = T))
denS_sum

## Collate data on number of pups observed
denS_pups <- subset(denS2, NoPups > 0)
denS_pups$time <- as.numeric(as.character(denS_pups$Year)) - 1996


####################################################################
#### ORGANIZE DATA INTO A LIST TOGETHER WITH ENVIRONMENTAL DATA ####
####################################################################

## Write down the sessions and corresponding years
SessionYears <- cbind(Index = c(1:22), Interval = paste(c(1997:2018),'(Jun)-',c(1998:2019), '(May)', sep = ''), Census = paste(c(1997:2018),'(Jun)', sep = ''))

## Load environmental data
load('200207_EnvCov_forIPM.RData')


## Organize all data into a single list
IPM.data <- list(

  # General information
  SessionYears = SessionYears,
  Tmax = max(dim(C)[2], nrow(denS_sum)),
  	
  # Mark-recovery module
  y.CMR = CMRR.data$rCH, initsCH.CMR = CMRR.data$BinitsCH,
  y.CMRR = CMRR.data$rrCH, initsCH.CMRR = CMRR.data$initsCH,
  y.CMRR2 = CMRR.data$rrCH2, y.CMRR3 = CMRR.data$rrCH3, initsCH.CMRR2 = CMRR.data$initsCH2,
  first = f, firstyear = firstyear, nind = dim(CMRR.data$rCH)[1], n.occasions = dim(CMRR.data$rCH)[2],
  
  # Age-at-harvest module
  C = C, TmaxC = dim(C)[2], A = dim(C)[1],
  pAgeSex = MissingDataProps$pAgeSex,
  pLoc = MissingDataProps$pLoc,
  
  # Age-at-death module
  M = M, TmaxM = dim(M)[2],
  pAgeSex.M = pAgeSex.M, 
  pLoc.M = pLoc.M,
  
  # Placental scar module
  P1 = P1, P1_age = P1_age, P1_year = P1_year, X1 = length(P1), 
  P2 = P2, P2_age = P2_age, P2_year = P2_year, X2 = length(P2),
  
  # Den survey module
  TmaxD = nrow(denS_sum),
  NoOcc = denS_sum$NoOcc, NoMon = denS_sum$TotMon, 
  NoPups = denS_pups$NoPups, DS_year = denS_pups$time, X3 = nrow(denS_pups),
  k.Dens = nrow(denS),
  
  # Environmental data
  MeanWinterTemp = envCov$MeanWinterTemp.sc, 
  WholeYearPrecip = envCov$WholeYearPrecip.sc,
  AmountROS = envCov$AmountROS.sc,
  NoDaysROS = envCov$NoDaysROS.sc,
  PtarmiganBag = envCov$PtarmiganBag.sc,
  PtarmiganDens = envCov$PtarmiganDens.sc,
  PtarmiganJuvProp = envCov$PtarmiganJuvProp.sc,
  RdCarcass = envCov$RdCarcass.sc,
  GooseCount = envCov$GooseCount.sc,
  GooseJuvProp = envCov$GooseJuvProp.sc,
  GooseBroodSize = envCov$GooseBroodSize.sc,
  PolarBearBCI = envCov$PolarBearBCI.sc,
  AprSeaIceSvalbard = envCov$AprSeaIceSvalbard.sc,
  AprMaySeaIceIsfj = envCov$AprMaySeaIceIsfj.sc,
  JanJunSeaIceIsfj = envCov$JanJunSeaIceIsfj.sc,
  WinterAO = envCov$WinterAO.sc,
  SpringAO = envCov$SpringAO.sc)


save(IPM.data, file = '200228_AF_IPM_Data.RData')


