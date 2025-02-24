library(ggplot2)
library(viridis)
library(gridExtra)
library(plyr)


## Load updated Isfjorden sea ice data (from Mikhail Itkin via Eva)
SI <- read.csv('isfjorden_1966-2019.csv', sep = ';')


##########################
#### DATA PREPARATION ####
##########################

# Isfjorden sea ice data
# ----------------------

## Check data availability
table(SI$year, SI$month, useNA = 'ifany')
# --> Data is now complete for all years 1966-2019


## Plot within-year change in sea ice (all time)
ggplot(SI, aes(x = as.factor(month), y = average)) + geom_boxplot(aes(group = as.factor(month)), col = 'grey70') + geom_point(aes(col = year)) + theme_bw() + scale_color_viridis(discrete = F) + theme(panel.grid = element_blank())

ggplot(SI, aes(x = as.factor(month), y = average)) + geom_line(aes(group = as.factor(year), col = year)) +  theme_bw() + scale_color_viridis(discrete = F) + theme(panel.grid = element_blank())

# --> Sea ice starts to appear in November, peaks in March, and is usually gone by July


## Plot within-year change in sea ice (study period)
ggplot(subset(SI, year%in%c(1996:2019)), aes(x = as.factor(month), y = average)) + geom_boxplot(aes(group = as.factor(month)), col = 'grey70') + geom_point(aes(col = year)) + theme_bw() + scale_color_viridis(discrete = F) + theme(panel.grid = element_blank())

ggplot(subset(SI, year%in%c(1996:2019)), aes(x = as.factor(month), y = average)) + geom_line(aes(group = as.factor(year), col = year)) +  theme_bw() + scale_color_viridis(discrete = F) + theme(panel.grid = element_blank())

# --> Consistent with above, but it's quite clear that it's still very low in November and December



## Summarise data into annual covariates

# 1) April-May: Period where foxes are expected to prey on seal pups
AprMay.sum <- ddply(subset(SI, month%in%c(4,5)), .(year), summarise, IFSI.meanAprMay = mean(average))

# 2) January-June: Entire season
JanJun.sum <- ddply(subset(SI, month%in%c(1:6)), .(year), summarise, IFSI.meanJanJun = mean(average))


## Combine measures into one data frame
IFSI.data <- merge(AprMay.sum, JanJun.sum, by = 'year', all = T)


## Check correlation
cor.test(IFSI.data$IFSI.meanAprMay, IFSI.data$IFSI.meanJanJun)
# Correlation coefficient = 0.90

ggplot(IFSI.data, aes(x = IFSI.meanAprMay, y = IFSI.meanJanJun)) + geom_point(aes(col = year)) + geom_smooth(method = 'lm', col = 'grey55') + theme_bw() + scale_color_viridis(discrete = F)


## Saving data
write.csv(IFSI.data, '200207_IsfjordenSeaIceCovariates.csv')

