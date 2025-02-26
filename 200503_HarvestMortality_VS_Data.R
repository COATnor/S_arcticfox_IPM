#############################################################
#### INTEGRATED POPULATION MODEL FOR SVALBARD ARCTIC FOX ####
#############################################################

library(coda)
library(reshape2)
library(ggplot2)
library(viridis)
library(plyr)
library(gridExtra)
library(GGally)
library(corrplot)

## Load workspace containing data used for model fitting
#load('200228_AF_IPM_Data.RData')
IPM.data <- readRDS("AF_IPM_Data.rds") 

## Load additional preliminary data on harvest effort
effData <- read.csv('Data/200324_TrappingEffort.csv')

## Load posterior samples
#load('200429_AF_IPM_VersionB4.RData')
AF.IPM <- readRDS("AF_IPM.rds")

## Re-arrange data
out.mat <- as.matrix(AF.IPM)
data <- melt(out.mat)
colnames(data) <- c('index', 'parameter', 'value')

## Summarise posterior into median and 50% and 95% CI
data.sum <- ddply(data, .(parameter), summarise, median = median(value, na.rm = T), lCI_95 = quantile(value, probs = 0.05, na.rm = T), uCI_95 = quantile(value, probs = 0.95, na.rm = T), lCI_50 = quantile(value, probs = 0.25, na.rm = T), uCI_50 = quantile(value, probs = 0.75, na.rm = T))

#---------------#
#  COLLATE DATA #
#---------------#

## Extract time-dependent harvest mortality hazard rates
Tmax <- 23
time <- c(1:Tmax)
  
mH.data <- data.frame(Year = time+1996, median = rep(NA, length(time)), lCI = rep(NA, length(time)), uCI = rep(NA, length(time)))
   
for(t in 1:(Tmax-1)){
	
  	mHa <- out.mat[, paste("mH[1, ", t, "]", sep = "")]
    mH.data[t,2] <- quantile(mHa, 0.5)
    mH.data[t,3] <- quantile(mHa, 0.05)
    mH.data[t,4] <- quantile(mHa, 0.95)
   
}


## Add Age-at-harvest data  
mH.data$AaH_Data <- c(colSums(IPM.data$C), NA)

## Add different measures of harvest effort
mH.data$NoTrappers <- effData$NoTrappersTotal
mH.data$NoTrappersSuccess <- effData$NoTrappersSuccess
mH.data$PropSuccess <- effData$NoTrappersSuccess/effData$NoTrappersTotal
mH.data$NoAreas <- effData$NoAreas
mH.data$NoTraps <- effData$NoTraps
mH.data$NoTrapDays <- effData$NoTrapDays


#---------------------------#
#  INVESTIGATE CORRELATIONS #
#---------------------------#

# Basic correlation coefficients
cor(mH.data[,c(2,5:11)], use = 'complete.obs')

# Pairplot
pdf('Harvest_Correlations_Pairs.pdf', width = 8.5, height = 8)
ggpairs(mH.data[,c(2,5:11)], lower = list(continuous=wrap("points", color = 'blue', alpha = 0.5)),
  diag = list(continuous=wrap("barDiag", color = NA, fill = 'blue'))) + theme_bw() + theme(panel.grid = element_blank())
dev.off()

# Correlation plot (with significance)
pdf('Harvest_Correlations_Corrs.pdf', width = 8.5, height = 8)
corrplot(cor(mH.data[,c(2,5:11)], use = 'complete.obs'), type = 'upper', tl.col="black")
corrplot(cor(mH.data[,c(2,5:11)], use = 'complete.obs'), sig.level = 0.05, insig = 'blank', type = 'upper', p.mat = cor.mtest(mH.data[,c(2,5:11)])$p, tl.col="black")
dev.off()

# RESULTS:
# -> Estimates of harvest mortality are positively correlated with all measures of harvest effort, and significantly so for the AaH data, Number of trappers, and Number of areas trapped
# -> The correlations of harvest mortality with effort data are stronger than the correlations of Aah data with effort data, indicating that the estimates may be closer to the "truth" than the raw data
# -> There are also positive correlation among many of the effort measures
# -> Notably, the proportion of successful trappers is NEGATIVELY correlated with the number of trappers, areas, traps, and trap days, indicating that harvest may be density-dependent!


#------------------#
# PLOT TIME SERIES #
#------------------#

## Re-arrange data
mH.data2 <- melt(mH.data[,-c(3,4)], id.vars = c('Year'))
mH.data2$lCI <- c(mH.data$lCI, rep(NA, nrow(mH.data2)-length(mH.data$lCI)))
mH.data2$uCI <- c(mH.data$uCI, rep(NA, nrow(mH.data2)-length(mH.data$lCI)))

## Make a second data frame with standardized effort values
mH.data.sc <- mH.data
mH.data.sc$AaH_Data <- scale(mH.data$AaH_Data)
mH.data.sc$NoTrappers <- scale(mH.data$NoTrappers)
mH.data.sc$NoTrappersSuccess <- scale(mH.data$NoTrappersSuccess)
mH.data.sc$PropSuccess <- scale(mH.data$PropSuccess)
mH.data.sc$NoAreas <- scale(mH.data$NoAreas)
mH.data.sc$NoTraps <- scale(mH.data$NoTraps)
mH.data.sc$NoTrapDays <- scale(mH.data$NoTrapDays)

mH.data2.sc <- melt(mH.data.sc[,-c(3,4)], id.vars = c('Year'))
mH.data2.sc$lCI <- c(mH.data.sc$lCI, rep(NA, nrow(mH.data2.sc)-length(mH.data.sc$lCI)))
mH.data2.sc$uCI <- c(mH.data.sc$uCI, rep(NA, nrow(mH.data2.sc)-length(mH.data.sc$lCI)))


## Plot harvest mortality hazard rate
plot.mH <- ggplot(mH.data, aes(x = Year, y = median)) + geom_line(color = viridis(7)[1]) + geom_ribbon(aes(ymin = lCI, ymax = uCI), color = NA, alpha = 0.5, fill = viridis(7)[1]) + ylab('Estimated harvest mortality') + scale_x_continuous(breaks = c(1997:2019)) + theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), axis.text.x = element_blank(), legend.position = 'none', axis.title.x = element_blank())
plot.mH

## Plot absolute effort - Multiple panels
plot.E1 <- ggplot(subset(mH.data2, variable%in%c('AaH_Data', 'NoTrappers', 'NoAreas', 'NoTrapDays')), aes(x = Year, y = value)) + geom_line(color = magma(5)[3], alpha = 0.5) + geom_point(color = magma(5)[3]) + ylab('Standardized number') + scale_x_continuous(breaks = c(1997:2019)) + facet_wrap(~variable, scales = 'free_y', ncol = 1) + theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))
plot.E1

## Plot relative effort - One panel
plot.E2 <- ggplot(subset(mH.data2.sc, variable%in%c('AaH_Data', 'NoTrappers', 'NoAreas', 'NoTrapDays')), aes(x = Year, y = value)) + geom_line(aes(color = variable), alpha = 0.5) + geom_point(aes(color = variable, shape = variable)) + ylab('Standardized number') + scale_x_continuous(breaks = c(1997:2019)) + scale_color_manual(values = c(viridis(3)[2], magma(10)[2], magma(10)[6], magma(10)[8]), labels = c('Foxes harvested', 'Trappers', 'Trapping areas', 'Trap days')) + scale_shape_manual(values = c(19, 15, 17, 18), labels = c('Foxes harvested', 'Trappers', 'Trapping areas', 'Trap days')) + theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5), legend.title = element_blank(), legend.position = c(0.85, 0.7))
plot.E2


# Combined plots - Several panels
gA <- ggplotGrob(plot.mH)
gB <- ggplotGrob(plot.E1)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth)
gB$widths[2:5] <- as.list(maxWidth)
 
pdf("HarvestMortality_VS_Data_1.pdf", width = 6.5, height = 8)
grid.arrange(gA, gB, ncol=1, heights = c(1, 2))
dev.off()


# Combined plots - 2 panels
gA <- ggplotGrob(plot.mH)
gB <- ggplotGrob(plot.E2)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth)
gB$widths[2:5] <- as.list(maxWidth)
 
pdf("HarvestMortality_VS_Data_2.pdf", width = 6.5, height = 5)
grid.arrange(gA, gB, ncol=1, heights = c(1, 1.2))
dev.off()

# --> Functionality confirmed. 