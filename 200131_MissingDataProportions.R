library(plyr)
library(reshape)


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

## Convert sex = 'F?' to 'F'
carcass2$Sex[which(carcass2$Sex=='F?')] <- 'F'
carcass2$Sex <- as.factor(as.character(carcass2$Sex))


## 1) Missing data and Sex
##************************

table(carcass2$Sex, useNA = 'always')
# --> 5% (122/2459) of all carcasses have no sex

table(carcass2$AgeClass, carcass2$Sex, useNA = 'always')
table(carcass2$AgeClass, carcass2$Sex, useNA = 'always')[,3]/rowSums(table(carcass2$AgeClass, carcass2$Sex, useNA = 'always'))
# --> Across age classes, the proportion of missing sex individual ranges from 0.8% to 4.1% (no substantial difference)
# --> Sex is clearly more likely to be missing for individuals that also do not have an age

# --> Calculate sexing probability across age classes (but exclude non-aged individuals)


table(carcass2$Core_area, carcass2$Sex, useNA = 'always')
table(carcass2$Core_area, carcass2$Sex, useNA = 'always')[,3]/rowSums(table(carcass2$Core_area, carcass2$Sex, useNA = 'always'))
# --> Outside the core area, 5.7% of carcasses have no sex
# --> Inside the core area, 2% of carcasses have no sex
# --> The difference of 3.7% is not substantial I would say 
# --> Clearly more unknown area carcasses also have no sex (26%)

# --> Calculate sexing probability either from core-area only individuals, or from core and not-core area individuals (but exclude unknown area individuals)



## 2) Missing data and Age
##************************

table(carcass2$AgeClass, useNA = 'always')
table(carcass2$AgeClass, useNA = 'always')[6]/sum(table(carcass2$AgeClass, useNA = 'always'))
# --> 8.7% (214/2459) of all carcasses have no age

table(carcass2$Sex, carcass2$AgeClass, useNA = 'always')
table(carcass2$Sex, carcass2$AgeClass, useNA = 'always')[,6]/rowSums(table(carcass2$Sex, carcass2$AgeClass, useNA = 'always'))
# --> 6.1 % of females have no age
# --> 5.4 % of males have no age
# --> The difference of 0.7% is negligible
# --> Clearly more unknown sex carcasses also have no age

# --> Calculate aging probability either from females only, or from females and males (but exclude unknown-age individuals)

table(carcass2$Core_area, carcass2$AgeClass, useNA = 'always')
table(carcass2$Core_area, carcass2$AgeClass, useNA = 'always')[,6]/rowSums(table(carcass2$Core_area, carcass2$AgeClass, useNA = 'always'))
# --> Outside the core area, 9.1% of carcasses have no age
# --> Inside the core area, 6.5% of carcasses have no age
# --> The difference of 2.6% is not substantial I would say 
# --> Clearly more unknown area carcasses also have no age (29%)

# --> Calculate aging probability either from core-area only individuals, or from core and not-core area individuals (but exclude unknown area individuals)


## 3) Missing data and Location
##*****************************

table(carcass2$Core_area, useNA = 'always')
table(carcass2$Core_area, useNA = 'always')[3]/sum(table(carcass2$AgeClass, useNA = 'always'))
# --> 1.0% (24/2459) of all carcasses have an unknown area

table(carcass2$Sex, carcass2$Core_area, useNA = 'always')
table(carcass2$Sex, carcass2$Core_area, useNA = 'always')[,3]/rowSums(table(carcass2$Sex, carcass2$Core_area, useNA = 'always'))
# --> 1.4 % of females have no age
# --> 0.6 % of males have no age
# --> 0.8 % of unknown sex individuals have no location
# --> Probablility of missing from the data due to no location information seems independent of sex and can thus be calculated for all sexes

table(carcass2$AgeClass, carcass2$Core_area, useNA = 'always')
table(carcass2$AgeClass, carcass2$Core_area, useNA = 'always')[,3]/rowSums(table(carcass2$AgeClass, carcass2$Core_area, useNA = 'always'))
# --> Across age classes, the proportion of unknown location individual ranges from 0.0% to 1.3% (no substantial difference)
# --> For unknown age individuals, the proportion is 0.9%

# --> Presence of location information does not seem to depend on age and probablility of missing from the data due to no location can thus be calculated for all age classes


## CONCLUSIONS:
#  --> The probabilities of being sexed and being aged are definitely not independent
#  --> Carcasses without sex and/or age info are more likely to lack information on location, but are not more likely to come from inside vs. outside the core area

## 	AGE-AT-HARVEST MATRIX OBSERVATION MODEL:
##  --> conditional on: 
#       (1) being harvested (prob = h[a,t])
#       (2) being sexed AND aged, given you from core area (prob = pAgeSex[t])
#
##  C[a,t] ~ dbin(h[a,t]*pAgeSex[t], N[a,t])
#
#       (3) being located (prob = pLoc[t])  
#
##  C[a,t] ~ dbin(h[a,t]*pAgeSex[t]*pLoc[t], N[a,t])


## ANNUAL PROPORTION OF CORE AREA CARCASSES WITH AGE AND SEX
carcass.core <- subset(carcass2, Core_area == 1)
carcass.core$MissingData <- ifelse(is.na(carcass.core$Sex) | is.na(carcass.core$AgeClass), 1, 0)
table(carcass.core$Session, carcass.core$MissingData, useNA = 'ifany')

pAgeSex <- unname(table(carcass.core$Session, carcass.core$MissingData, useNA = 'ifany')[,1]/rowSums(table(carcass.core$Session, carcass.core$MissingData, useNA = 'ifany')))
plot(c(1997:2018), pAgeSex, type = 'l')


## ANNUAL PROPORTION OF KNOWN LOCATION CARCASSES 
carcass2$MissingData <- ifelse(is.na(carcass2$Core_area), 1, 0)
table(carcass2$Session, carcass2$MissingData, useNA = 'ifany')

pLoc <- unname(table(carcass2$Session, carcass2$MissingData, useNA = 'ifany')[,1]/rowSums(table(carcass2$Session, carcass2$MissingData, useNA = 'ifany')))
plot(c(1997:2018), pLoc, type = 'l')

pdf('200131_MissingDataProportions.pdf', width = 8, height = 5)
plot(c(1997:2018), 1-pAgeSex, type = 'l', col = 'red', ylab = 'Proportion missing data', xlab = 'Trapping season (start year)', lwd = 2)
lines(c(1997:2018), 1-pLoc, col = 'blue', lwd = 2)
abline(v = c(1997:2018), lty = 3, col = 'grey80')
legend(2010, 0.5, legend = c('Sex and/or age', 'Location'), col = c('red', 'blue'), lwd = 2, bty = 'n')
dev.off()

MissingDataProps <- list(pAgeSex = pAgeSex, pLoc = pLoc)
save(MissingDataProps, file = '200131_MissingDataProps.RData')

