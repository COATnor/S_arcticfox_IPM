library(ggplot2)
library(viridis)
library(gridExtra)
library(lme4)
library(splines)
library(MuMIn)
library(plyr)

## Load data
PBdata <- read.table("TabBCI92_19.txt", sep="\t", header=TRUE)
head(PBdata)
str(PBdata)

## Take a subset including only male bears
mPBdata <- subset(PBdata, sex == 'm')


################################
#### BASIC DATA EXPLORATION ####
################################

## Dataset in numbers
nrow(mPBdata) # 578 observations
length(unique(mPBdata$id)) # 374 individuals
length(unique(mPBdata$age)) # 23 ages ranging...
min(mPBdata$age) # ... from 5...
max(mPBdata$age) # ... to 27
hist(mPBdata$age) # Quite evenly distributed, except for old individuals
length(unique(mPBdata$yr)) # 28 years ranging...
min(mPBdata$yr) # ... from 1992...
max(mPBdata$yr) # ... to 2019
hist(mPBdata$yr, breaks = 28) # Quite evenly distributed, except for the year 2000 (which had many more measurements)


## Distribution of BCI
hist(mPBdata$bci, breaks = 100)
# --> Looks fairly normal

## BCI ~ age
ggplot(mPBdata, aes(x = age, y = bci)) + geom_point(alpha = 0.5) + geom_smooth(col = 'blue', fill = 'blue') + theme_bw()
# --> There is a clear increase with age until ~11-12 years old and no relationship afterwards


## BCI ~ sampling date
ggplot(mPBdata, aes(x = DayM, y = bci)) + geom_point(alpha = 0.5) + geom_smooth(col = 'blue', fill = 'blue') + theme_bw()
# --> There may or may not be a relationship
# --> Seems likely that any pattern is driven by the very early and very late sampling occasions
table(subset(mPBdata, DayM < 30)$yr)
table(subset(mPBdata, DayM > 55)$yr)
ggplot(mPBdata, aes(x = as.factor(yr), y = DayM)) + geom_boxplot() + theme_bw()
# --> Sampling seasons are confounded with year, so the pattern here may be explained by year variation (not change in condition within year)

## BCI ~ year
ggplot(mPBdata, aes(x = yr, y = bci)) + geom_point(alpha = 0.5) + geom_smooth(col = 'blue', fill = 'blue') + theme_bw()
# --> There is clearly something going on here, particularly for the first half of the time-series

## BCI ~ sea ice
ggplot(mPBdata, aes(x = seaiceApr, y = bci)) + geom_point(alpha = 0.5) + geom_smooth(col = 'blue', fill = 'blue') + theme_bw()
# --> Not much going on here

## BCI ~ location
p1 <- ggplot(subset(mPBdata, long > 10), aes(x = long, y = bci)) + geom_point(aes(color = lat)) + geom_smooth(col = 'grey55', fill = 'grey55') + theme_bw() + scale_color_viridis()
p2 <- ggplot(mPBdata, aes(x = lat, y = bci)) + geom_point(aes(color = long)) + geom_smooth(col = 'grey55', fill = 'grey55') + theme_bw() + scale_color_viridis()
grid.arrange(p1, p2, ncol = 2)
# --> There is some patterning, but it's hard to interpret
# --> A group of high longitude low latitude measurements (quite low BCI) really stand out
subset(mPBdata, long > 35 & lat < 76)
# --> with the exception of one, these are all from the same year (1998)
# --> according to google maps, the locations are about halfway between Svalbard and Russia. Should ask Jon about these. 


############################
#### MIXED SPLINE MODEL ####
############################

## Removing outliers (Barents sea pelagic bears)
mPBdata <- subset(mPBdata, long < 35)


## Candidate models

# Possible combinations 
x1 <- c('none','lin','poly2','bs3','bs4')
x2 <- c('lin','poly2','bs3','bs4')
expand.grid(DayM=x1,age=x2)

# Linear age effect
mod1 <- lmer(bci ~ 0 + age + (1|id) + (1|yr), data = mPBdata, na.action = na.fail)
mod2 <- lmer(bci ~ 0 + DayM + age + (1|id) + (1|yr), data = mPBdata, na.action = na.fail)
mod3 <- lmer(bci ~ 0 + poly(DayM, 2) + age + (1|id) + (1|yr), data = mPBdata, na.action = na.fail)
mod4 <- lmer(bci ~ 0 + bs(DayM, df=3) + age + (1|id) + (1|yr), data = mPBdata, na.action = na.fail)
mod5 <- lmer(bci ~ 0 + bs(DayM, df=4) + age + (1|id) + (1|yr), data = mPBdata, na.action = na.fail)

# Quadratic age effect
mod6 <- lmer(bci ~ 0 + poly(age,2) + (1|id) + (1|yr), data = mPBdata, na.action = na.fail)
mod7 <- lmer(bci ~ 0 + DayM + poly(age,2) + (1|id) + (1|yr), data = mPBdata, na.action = na.fail)
mod8 <- lmer(bci ~ 0 + poly(DayM, 2) + poly(age,2) + (1|id) + (1|yr), data = mPBdata, na.action = na.fail)
mod9 <- lmer(bci ~ 0 + bs(DayM, df=3) + poly(age,2) + (1|id) + (1|yr), data = mPBdata, na.action = na.fail)
mod10 <- lmer(bci ~ 0 + bs(DayM, df=4) + poly(age,2) + (1|id) + (1|yr), data = mPBdata, na.action = na.fail)

# 3df spline age effect
mod11 <- lmer(bci ~ 0 + bs(age, df=3) + (1|id) + (1|yr), data = mPBdata, na.action = na.fail)
mod12 <- lmer(bci ~ 0 + DayM + bs(age, df=3) + (1|id) + (1|yr), data = mPBdata, na.action = na.fail)
mod13 <- lmer(bci ~ 0 + poly(DayM, 2) + bs(age, df=3) + (1|id) + (1|yr), data = mPBdata, na.action = na.fail)
mod14 <- lmer(bci ~ 0 + bs(DayM, df=3) + bs(age, df=3) + (1|id) + (1|yr), data = mPBdata, na.action = na.fail)
mod15 <- lmer(bci ~ 0 + bs(DayM, df=4) + bs(age, df=3) + (1|id) + (1|yr), data = mPBdata, na.action = na.fail)

# 4df spline age effect
mod16 <- lmer(bci ~ 0 + bs(age, df=4) + (1|id) + (1|yr), data = mPBdata, na.action = na.fail)
mod17 <- lmer(bci ~ 0 + DayM + bs(age, df=4) + (1|id) + (1|yr), data = mPBdata, na.action = na.fail)
mod18 <- lmer(bci ~ 0 + poly(DayM, 2) + bs(age, df=4) + (1|id) + (1|yr), data = mPBdata, na.action = na.fail)
mod19 <- lmer(bci ~ 0 + bs(DayM, df=3) + bs(age, df=4) + (1|id) + (1|yr), data = mPBdata, na.action = na.fail)
mod20 <- lmer(bci ~ 0 + bs(DayM, df=4) + bs(age, df=4) + (1|id) + (1|yr), data = mPBdata, na.action = na.fail)


## Model selection

# Comparison of all models
anova(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20)

# Comparison of all models significantly better than mod1
anova(mod2,mod4,mod5,mod6,mod7,mod8,mod9,mod14,mod15,mod17)

# Comparison of models significantly better than mod2
anova(mod6,mod8,mod9,mod14)

anova(mod8,mod14)


###########################
#### MODEL FIT TO DATA ####
###########################

## Make predictions under the top 3 models (average conditions)
dataAge <- data.frame(age = c(min(mPBdata$age):max(mPBdata$age)), DayM = mean(mPBdata$DayM), bci = 0, id = NA, yr = NA)
#dataAge <- data.frame(age = c(min(mPBdata$age):max(mPBdata$age)), DayM = median(mPBdata$DayM), bci = 0, id = NA, yr = NA)

dataDayM <- data.frame(age = c(mean(mPBdata$age)), DayM = c(min(mPBdata$DayM):max(mPBdata$DayM)), bci = 0, id = NA, yr = NA)
#dataDayM <- data.frame(age = c(median(mPBdata$age)), DayM = c(min(mPBdata$DayM):max(mPBdata$DayM)), bci = 0, id = NA, yr = NA)

predict.PBbci = function(model,data){
	
	mm <- model.matrix(terms(model),data)
	data$bci <- predict(model,data,re.form=NA)

	pvar1 <- diag(mm %*% tcrossprod(vcov(model),mm))
	tvar1 <- pvar1+VarCorr(model)$Subject[1] 
	cmult <- 1.96 
	data <- data.frame(data, plo = data$bci-cmult*sqrt(pvar1), phi = data$bci+cmult*sqrt(pvar1))

	return(data)
}

rAge.mod8 <- predict.PBbci(mod8, dataAge)
rAge.mod14 <- predict.PBbci(mod14, dataAge)

rDayM.mod8 <- predict.PBbci(mod8, dataDayM)
rDayM.mod14 <- predict.PBbci(mod14, dataDayM)


## Plot top 2 models over data

# BCI ~ age (scatterplot)
p1 <- ggplot(mPBdata, aes(x = age, y = bci)) + geom_point(alpha = 0.5) + theme_bw() + theme(panel.grid = element_blank())
p1 <- p1 + 
geom_line(data = rAge.mod8, aes(x = age, y = bci), col = viridis(3)[1]) + geom_ribbon(data = rAge.mod8, aes(ymin = plo, ymax = phi), fill = viridis(3)[1], alpha = 0.3) + 
geom_line(data = rAge.mod14, aes(x = age, y = bci), col = viridis(3)[2]) + geom_ribbon(data = rAge.mod14, aes(ymin = plo, ymax = phi), fill = viridis(3)[2], alpha = 0.3)

# BCI ~ day (scatterplot)
p2 <- ggplot(mPBdata, aes(x = DayM, y = bci)) + geom_point(alpha = 0.5) + theme_bw() + theme(panel.grid = element_blank())
p2 <- p2 + 
geom_line(data = rDayM.mod8, aes(x = DayM, y = bci), col = viridis(3)[1]) + geom_ribbon(data = rDayM.mod8, aes(ymin = plo, ymax = phi), fill = viridis(3)[1], alpha = 0.3) + 
geom_line(data = rDayM.mod14, aes(x = DayM, y = bci), col = viridis(3)[2]) + geom_ribbon(data = rDayM.mod14, aes(ymin = plo, ymax = phi), fill = viridis(3)[2], alpha = 0.3)

grid.arrange(p1, p2, nrow = 1)

# BCI ~ age (boxplot)
p3 <- ggplot(mPBdata, aes(x = age, y = bci, group = age)) + geom_boxplot(color = 'grey40') + theme_bw() + theme(panel.grid = element_blank())
p3 <- p3 + 
geom_line(data = rAge.mod8, aes(x = age, y = bci, group = 1), col = viridis(3)[1]) + geom_ribbon(data = rAge.mod8, aes(ymin = plo, ymax = phi, group = 1), fill = viridis(3)[1], alpha = 0.3) + 
geom_line(data = rAge.mod14, aes(x = age, y = bci, group = 1), col = viridis(3)[2]) + geom_ribbon(data = rAge.mod14, aes(ymin = plo, ymax = phi, group = 1), fill = viridis(3)[2], alpha = 0.3)

# BCI ~ day (boxplot)
p4 <- ggplot(mPBdata, aes(x = DayM, y = bci, group = DayM)) + geom_boxplot(color = 'grey40') + theme_bw() + theme(panel.grid = element_blank())
p4 <- p4 + 
geom_line(data = rDayM.mod8, aes(x = DayM, y = bci, group = 1), col = viridis(3)[1]) + geom_ribbon(data = rDayM.mod8, aes(ymin = plo, ymax = phi, group = 1), fill = viridis(3)[1], alpha = 0.3) + 
geom_line(data = rDayM.mod14, aes(x = DayM, y = bci, group = 1), col = viridis(3)[2]) + geom_ribbon(data = rDayM.mod14, aes(ymin = plo, ymax = phi, group = 1), fill = viridis(3)[2], alpha = 0.3)

pdf('ModelFit_toData.pdf', width = 8, height = 6)
grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()
# --> The fit with age looks quite good for all models
# --> For day on the other hand, it looks like models largely predict too high BCI on average...


################################
#### OBSERVED VS. PREDICTED ####
################################

## Make predictions for specific data points
mPBdata$bciPred.mod8 <- predict(mod8)
mPBdata$bciPred.mod14 <- predict(mod14)

## Plot observed vs. predicted
pp1 <- ggplot(mPBdata, aes(x = bci, y = bciPred.mod8)) + geom_point(alpha = 0.3) + theme_bw() + theme(panel.grid = element_blank()) + geom_abline(slope = 1, intercept = 0, linetype = 'dashed') + ylab('Predicted BCI') + xlab('Observed BCI') + ggtitle('Model 8') + geom_smooth(method = 'lm', fill = viridis(3)[1], color = viridis(3)[1])

pp2 <- ggplot(mPBdata, aes(x = bci, y = bciPred.mod14)) + geom_point(alpha = 0.3) + theme_bw() + theme(panel.grid = element_blank()) + geom_abline(slope = 1, intercept = 0, linetype = 'dashed') + ylab('Predicted BCI') + xlab('Observed BCI') + ggtitle('Model 14') + geom_smooth(method = 'lm', fill = viridis(3)[2], color = viridis(3)[2])


## Plot oserved vs. residuals
pp4 <- ggplot(mPBdata, aes(x = bciPred.mod8-bci, y = bciPred.mod8)) + geom_point(alpha = 0.3) + theme_bw() + theme(panel.grid = element_blank()) + geom_hline(yintercept = 0, linetype = 'dashed') + ylab('Residuals') + xlab('Observed BCI') + ggtitle('Model 8') + geom_smooth(method = 'lm', fill = viridis(3)[1], color = viridis(3)[1])

pp5 <- ggplot(mPBdata, aes(x = bciPred.mod14-bci, y = bciPred.mod14)) + geom_point(alpha = 0.3) + theme_bw() + theme(panel.grid = element_blank()) + geom_hline(yintercept = 0, linetype = 'dashed') + ylab('Residuals') + xlab('Observed BCI') + ggtitle('Model 14') + geom_smooth(method = 'lm', fill = viridis(3)[2], color = viridis(3)[2])


pdf('ModelDiagnostics.pdf', width = 10, height = 6)
grid.arrange(pp1,pp2,pp4,pp5, nrow = 2)
dev.off()
# --> This is really not good. All models overestimate low BCIs and underestimate high BCIs (to a similar degree)
# --> I checked (graphically) whether the lack of fit is correlated with variables in the dataframe (year, age, DayM, latitude, longitude, etc.) but did not find any obvious patterns
# --> So this means I have literally no idea what is causing the bad fit...

par(mfrow = c(1,2))
qqnorm(residuals(mod8))
qqnorm(residuals(mod14))
# --> No major problem with normality of residuals


#######################################################
#### EXTRACTION OF POLAR BEAR BCI ANNUAL COVARIATE ####
#######################################################

## Top model comparison with fixed year effects
mod8b <- lmer(bci ~ 0 + as.factor(yr) + poly(DayM, 2) + poly(age,2) + (1|id), data = mPBdata, na.action = na.fail)
mod14b <- lmer(bci ~ 0 + as.factor(yr) + bs(DayM, df=3) + bs(age, df=3) + (1|id), data = mPBdata, na.action = na.fail)
anova(mod8b, mod14b)
# --> This does not favour any model over the other, and since the estimates from all models are highly similar, I may go with either. I'll go with model 14 for now because it seems to better account for high ages and very early/late sampling days

## Extract year-specific parameters
bciPB <- data.frame(year = c(1992:2019), bciPB = unname(fixef(mod14b)[1:28]))


## Compare year-specific parameters to raw data
p7 <- ggplot(mPBdata, aes(x = yr, y = bci, group = yr)) + geom_boxplot(color = 'grey50') + theme_bw() + theme(panel.grid = element_blank()) + geom_line(data = bciPB, aes(x = year, y = bciPB, group = 1), col = viridis(3)[2], size = 1) + xlab('Year') + ylab('Polar bear BCI')

bciPB.raw <- ddply(mPBdata, .(yr), summarise, mean.bci = mean(bci))
bciPB$bci.raw <- bciPB.raw$mean.bci

p8 <- ggplot(bciPB, aes(x = bci.raw, y = bciPB)) + geom_smooth(method = 'lm', col = 'grey70', fill = 'grey70', linetype = 'dashed') + geom_point(col = viridis(3)[2], size = 2) + theme_bw() + theme(panel.grid = element_blank()) + ylab('Polar bear BCI proxy') + xlab('Population-average BCI (raw data)')

pdf('Proxy_vs_RawData.pdf', width = 6, height = 7)
grid.arrange(p7, p8, nrow = 2)
dev.off()

summary(lm(bciPB ~ bci.raw, data = bciPB))
# 81% of variance in our proxy is explained by variance in the raw data

cor.test(bciPB$bciPB, bciPB$bci.raw)
# the correlation coefficient it 0.90

# --> Overall, this seems quite good!

## Save bci data
write.csv(bciPB, '191130_PolarBearBCI_Covariate.csv')

