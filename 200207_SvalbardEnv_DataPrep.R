library(plyr)
library(reshape)
library(ggplot2)
library(viridis)
library(lubridate)
library(SpecsVerification)

###############################
#### BASIC ABIOTIC FACTORS ####
###############################

## Load data
envA <- read.csv('SvalbardAirport_Environment.csv', sep = ';')

## Rename month colums to numbers
names(envA)[4:15] <- c(1:12)

## Transform into longitudinal format
envA2 <- melt(envA, id.vars = c('Description','Variable','Year'), measure.vars = names(envA)[4:15])

names(envA2)[4:5] <- c('Month','Value')

## Add a continuous date to allow plotting
envA2$Date <- parse_date_time(paste(envA2$Year,'-',envA2$Month, sep = ''), "ym")

## Plot all data (high resolution)
ggplot(envA2, aes(x = Date, y = Value, group = Year)) + geom_line(aes(color = Year)) + facet_wrap(~Description, scales = 'free_y') + theme_bw() + scale_color_viridis()
# --> Trends are visible in all temperature measures (and esp. in deviation from norm)
# --> All measures of precipitation show similar patterns
# --> Snow data has a large gap early in the time series

varSelect <- c('TAM', 'TAMA', 'RA', 'RR')
ggplot(subset(envA2, Variable%in%varSelect), aes(x = Date, y = Value, group = Year)) + geom_line(aes(color = Year)) + facet_wrap(~Description, scales = 'free_y') + theme_bw() + scale_color_viridis()

## Plot temp/precip. by month and year
ggplot(subset(envA2, Variable=='TAM'), aes(x = Month, y = Value, group = Year)) + geom_line(aes(color = Year)) + theme_bw() + scale_color_viridis()

ggplot(subset(envA2, Variable=='TAM'), aes(x = Year, y = Value, group = Month)) + geom_line(aes(color = Month)) + theme_bw() + scale_color_viridis(discrete = T)
# --> Temperature in the winter months has been much more variable (and trending) than in the summer months

ggplot(subset(envA2, Variable=='RR'), aes(x = Month, y = Value, group = Year)) + geom_line(aes(color = Year)) + theme_bw() + scale_color_viridis()

ggplot(subset(envA2, Variable=='RR'), aes(x = Year, y = Value, group = Month)) + geom_line(aes(color = Month)) + theme_bw() + scale_color_viridis(discrete = T)


## Add the 'fox year' to the data
envA2$Month <- as.numeric(as.character(envA2$Month)) 
envA2$FoxYear <- ifelse(envA2$Month>5, envA2$Year+1, envA2$Year)

## Summarize mean winter temperature and total precipitation into yearly values

TempData.sum <- ddply(subset(envA2, Month%in%c(10:12,1:4) & Variable == 'TAM'), .(FoxYear), summarise, MeanWinterTemp = mean(Value))

PrecipData.sum <- ddply(subset(envA2, Variable == 'RR'), .(FoxYear), summarise, TotalPrecip = sum(Value))

## Organize into vectors that can be used as data for the IPM
MeanWinterTemp <- TempData.sum$MeanWinterTemp[1:23]
WholeYearPrecip <- PrecipData.sum$TotalPrecip[1:23]

#############################
#### RAIN-ON-SNOW EVENTS ####
#############################

# NOTE: Data and code obtained from Filippo (Oct 19)

klima <- read.table("ROS_science_Brage.txt",header=T)

klima$Date <- as.POSIXct(klima$Date, tz = "GMT", format = "%d.%m.%Y")
klima[["dy"]] <- as.numeric(format(klima$Date, "%d"))
klima[["mo"]] <- as.numeric(format(klima$Date, "%m"))
klima[["yr"]] <- as.numeric(format(klima$Date, "%Y"))
# date, daily prec, daily prec adj, daily T, daily T adj, day, month, year
df <- klima
years <- unique(klima$yr)
CL <- data.frame(YEAR = years)
CL$Temp <- NA 
CL$Prec <- NA
CL$Snow <- NA
CL$ROS <- NA
CL$ROS_days

# Hansen et al. 2013 : dec-mar
start.month = 12
end.month = 3
#
#start.month = 11
#end.month = 4

for(t in 2:length(years)) { 
 
  win <- with(df, df[(mo >= start.month & yr == years[t-1]) | (mo <= end.month & yr == years[t]),]) # seleziona range
 
  CL$Prec[t] <- with(win, sum(RRTA, na.rm=T)) # total Prec
  CL$ROS[t] <- sum(win[which(win$TAM >= 1),]$RRTA, na.rm=T) # ROS mm
  CL$Snow[t] <-  sum(win[which(win$TAM < 1),]$RRTA, na.rm=T)
  CL$ROS_days[t] <- nrow(win[which(win$TAM >= 1 & win$RRTA >= 1),])
  CL$Temp[t] <- round(with(win, mean(TAM, na.rm=TRUE)),2) # mean  Temp
 
}

## Organize into vectors that can be used as data for the IPM
AmountROS <- c(CL$ROS[41:62], NA) 
NoDaysROS <- c(CL$ROS_days[41:62], NA)


########################
#### BIOTIC FACTORS ####
########################

## Read in data
envB <- read.csv('SvalbardTerr_Covariates.csv')


## Organize into vectors that can be used as data for the IPM

# Harvest data (collected primarily in winter)
PtarmiganBag <- c(NA, as.vector(envB$ptarm_hun)[1:22])
PtarmiganJuvProp <- c(NA, as.vector(envB$ptarm_prop_juv)[1:22])


# Monitoring data (collected primarily in summer)
PtarmiganDens <- as.vector(envB$ptarm_dens)
RdCarcass <- as.vector(envB$reindeer_carcass) 

GooseCount <- as.vector(envB$goose_count)
GooseJuvProp <- as.vector(envB$prop_goose_juv)
GooseBroodSize <- as.vector(envB$goose_brood_size)


###################################
#### POLAR BEAR & SEA ICE DATA ####
###################################

## Read in data
envC <- read.csv('SeaIce_PolarBearBCI/191130_SeaIceCovariates.csv')
envD <- read.csv('SeaIce_PolarBearBCI/200207_IsfjordenSeaIceCovariates.csv')

## Organise into vectors that can be used as data for the IPM
PolarBearBCI <- envC$bciPB[6:28]
AprSeaIceSvalbard <- envC$SIApr[6:28]

AprMaySeaIceIsfj <- envD$IFSI.meanAprMay[32:54]
JanJunSeaIceIsfj <- envD$IFSI.meanJanJun[32:54]


#################################
#### ARCTIC OSCILLATION DATA ####
#################################

## Load arctic oscillation (AO) data
AOdata <- read.csv('/Users/chloern/Dropbox/Arctic_Fox_SUSTAIN/Data/Covariates/SeaIce_PolarBearBCI/200128_AOdata_raw.csv')

AOdata <- subset(AOdata, year%in%c(1997:2019))

## Calculate mean AO for the winter (Jan - Apr) and spring (Apr-Jun)
AOdata$wAO <- rowSums(AOdata[,2:5])/4
AOdata$sAO <- rowSums(AOdata[,5:7])/3

WinterAO <- AOdata$wAO
SpringAO <- AOdata$sAO

##################
#### OVERVIEW ####
##################

## Plot data series

pdf('200207_TimeTrends_Covariates.pdf', width = 15, height = 9)
par(mfrow = c(3,6))
plot(c(1997:2019), MeanWinterTemp, type = 'l', xlab = 'Year', ylab = 'Mean winter temperature', col = 'cornflowerblue')
plot(c(1997:2019), WholeYearPrecip, type = 'l', xlab = 'Year', ylab = 'Total annual precipitation', col = 'cornflowerblue')
plot(c(1997:2019), AmountROS, type = 'l', xlab = 'Year', ylab = 'Total amount ROS', col = 'cornflowerblue')
plot(c(1997:2019), NoDaysROS, type = 'l', xlab = 'Year', ylab = 'Number of days with ROS', col = 'cornflowerblue')
plot(c(1997:2019), PtarmiganBag, type = 'l', xlab = 'Year', ylab = 'Ptarmigans hunted', col = 'forestgreen')
plot(c(1997:2019), PtarmiganDens, type = 'l', xlab = 'Year', ylab = 'Ptarmigan density', col = 'forestgreen')
plot(c(1997:2019), PtarmiganJuvProp, type = 'l', xlab = 'Year', ylab = 'Proportion juvenile ptarmigan', col = 'forestgreen')
plot(c(1997:2019), RdCarcass, type = 'l', xlab = 'Year', ylab = 'Reindeer carcasses', col = 'forestgreen')
plot(c(1997:2019), GooseCount, type = 'l', xlab = 'Year', ylab = 'Number of geese', col = 'forestgreen')
plot(c(1997:2019), GooseJuvProp, type = 'l', xlab = 'Year', ylab = 'Proportion juvenile geese', col = 'forestgreen')
plot(c(1997:2019), GooseBroodSize, type = 'l', xlab = 'Year', ylab = 'Goose brood size', col = 'forestgreen')
plot(c(1997:2019), PolarBearBCI, type = 'l', xlab = 'Year', ylab = 'Polar bear BCI', col = 'purple')
plot(c(1997:2019), AprSeaIceSvalbard, type = 'l', xlab = 'Year', ylab = 'Whole Svalbard sea ice (Apr)', col = 'purple')
plot(c(1997:2019), AprMaySeaIceIsfj, type = 'l', xlab = 'Year', ylab = 'Isfjorden sea ice (Apr-May)', col = 'purple')
plot(c(1997:2019), JanJunSeaIceIsfj, type = 'l', xlab = 'Year', ylab = 'Isfjorden sea ice (Jan-Jun)', col = 'purple')
plot(c(1997:2019), WinterAO, type = 'l', xlab = 'Year', ylab = 'Arctic osciallation (winter)', col = 'purple')
plot(c(1997:2019), SpringAO, type = 'l', xlab = 'Year', ylab = 'Arctic osciallation (spring)', col = 'purple')
dev.off()


## Detrend all covariates  (residuals of linear regression)
MeanWinterTemp <- as.vector(Detrend(MeanWinterTemp))
WholeYearPrecip <- as.vector(Detrend(WholeYearPrecip))
AmountROS <- as.vector(Detrend(AmountROS))
NoDaysROS <- as.vector(Detrend(NoDaysROS))
PtarmiganBag <- as.vector(Detrend(PtarmiganBag))
PtarmiganDens <- as.vector(Detrend(PtarmiganDens))
PtarmiganJuvProp <- as.vector(Detrend(PtarmiganJuvProp))
RdCarcass <- as.vector(Detrend(RdCarcass))
GooseCount <- as.vector(Detrend(GooseCount))
GooseJuvProp <- as.vector(Detrend(GooseJuvProp))
GooseBroodSize <- as.vector(Detrend(GooseBroodSize))
PolarBearBCI <- as.vector(Detrend(PolarBearBCI))
AprSeaIceSvalbard <- as.vector(Detrend(AprSeaIceSvalbard))
AprMaySeaIceIsfj <- as.vector(Detrend(AprMaySeaIceIsfj))
JanJunSeaIceIsfj <- as.vector(Detrend(JanJunSeaIceIsfj))
WinterAO <- as.vector(Detrend(WinterAO))
SpringAO <- as.vector(Detrend(SpringAO))

pdf('200207_TimeTrends_CovariatesDetrend.pdf', width = 15, height = 9)
par(mfrow = c(3,6))
plot(c(1997:2019), MeanWinterTemp, type = 'l', xlab = 'Year', ylab = 'Mean winter temperature', col = 'cornflowerblue')
plot(c(1997:2019), WholeYearPrecip, type = 'l', xlab = 'Year', ylab = 'Total annual precipitation', col = 'cornflowerblue')
plot(c(1997:2019), AmountROS, type = 'l', xlab = 'Year', ylab = 'Total amount ROS', col = 'cornflowerblue')
plot(c(1997:2019), NoDaysROS, type = 'l', xlab = 'Year', ylab = 'Number of days with ROS', col = 'cornflowerblue')
plot(c(1997:2019), PtarmiganBag, type = 'l', xlab = 'Year', ylab = 'Ptarmigans hunted', col = 'forestgreen')
plot(c(1997:2019), PtarmiganDens, type = 'l', xlab = 'Year', ylab = 'Ptarmigan density', col = 'forestgreen')
plot(c(1997:2019), PtarmiganJuvProp, type = 'l', xlab = 'Year', ylab = 'Proportion juvenile ptarmigan', col = 'forestgreen')
plot(c(1997:2019), RdCarcass, type = 'l', xlab = 'Year', ylab = 'Reindeer carcasses', col = 'forestgreen')
plot(c(1997:2019), GooseCount, type = 'l', xlab = 'Year', ylab = 'Number of geese', col = 'forestgreen')
plot(c(1997:2019), GooseJuvProp, type = 'l', xlab = 'Year', ylab = 'Proportion juvenile geese', col = 'forestgreen')
plot(c(1997:2019), GooseBroodSize, type = 'l', xlab = 'Year', ylab = 'Goose brood size', col = 'forestgreen')
plot(c(1997:2019), PolarBearBCI, type = 'l', xlab = 'Year', ylab = 'Polar bear BCI', col = 'purple')
plot(c(1997:2019), AprSeaIceSvalbard, type = 'l', xlab = 'Year', ylab = 'Whole Svalbard sea ice (Apr)', col = 'purple')
plot(c(1997:2019), AprMaySeaIceIsfj, type = 'l', xlab = 'Year', ylab = 'Isfjorden sea ice (Apr-May)', col = 'purple')
plot(c(1997:2019), JanJunSeaIceIsfj, type = 'l', xlab = 'Year', ylab = 'Isfjorden sea ice (Jan-Jun)', col = 'purple')
plot(c(1997:2019), WinterAO, type = 'l', xlab = 'Year', ylab = 'Arctic osciallation (winter)', col = 'purple')
plot(c(1997:2019), SpringAO, type = 'l', xlab = 'Year', ylab = 'Arctic osciallation (spring)', col = 'purple')
dev.off()



## Write a generic scaling function that does not attach attributes
scale.simple = function(v){

	v.sc <- (v - mean(v, na.rm = T)) / sd (v, na.rm = T)	
	return(v.sc)
}

## Save all covariate vectors into a list

envCov <- list(	
	MeanWinterTemp = MeanWinterTemp,
	MeanWinterTemp.sc = scale.simple(MeanWinterTemp), 
	WholeYearPrecip = WholeYearPrecip,
	WholeYearPrecip.sc = scale.simple(WholeYearPrecip),
	AmountROS = AmountROS, 
	AmountROS.sc = scale.simple(AmountROS),
	NoDaysROS = NoDaysROS, 
	NoDaysROS.sc = scale.simple(NoDaysROS),
	PtarmiganBag = PtarmiganBag,
	PtarmiganBag.sc = scale.simple(PtarmiganBag),
	PtarmiganDens = PtarmiganDens,
	PtarmiganDens.sc = scale.simple(PtarmiganDens),
	PtarmiganJuvProp = PtarmiganJuvProp,
	PtarmiganJuvProp.sc = scale.simple(PtarmiganJuvProp),
	RdCarcass = RdCarcass, 
	RdCarcass.sc = scale.simple(RdCarcass),
	GooseCount = GooseCount, 
	GooseCount.sc = scale.simple(GooseCount),
	GooseJuvProp = GooseJuvProp, 
	GooseJuvProp.sc = scale.simple(GooseJuvProp),
	GooseBroodSize = GooseBroodSize, 
	GooseBroodSize.sc = scale.simple(GooseBroodSize), 
	PolarBearBCI = PolarBearBCI, 
	PolarBearBCI.sc = scale.simple(PolarBearBCI),
	AprSeaIceSvalbard = AprSeaIceSvalbard, 
	AprSeaIceSvalbard.sc = scale.simple(AprSeaIceSvalbard),
	AprMaySeaIceIsfj = AprMaySeaIceIsfj, 
	AprMaySeaIceIsfj.sc = scale.simple(AprMaySeaIceIsfj), 
	JanJunSeaIceIsfj = JanJunSeaIceIsfj, 
	JanJunSeaIceIsfj.sc = scale.simple(JanJunSeaIceIsfj), 
	WinterAO = WinterAO, 
	WinterAO.sc = scale.simple(WinterAO),
	SpringAO = SpringAO, 
	SpringAO.sc = scale.simple(SpringAO))

save(envCov, file = '200207_EnvCov_forIPM.RData')



