#############################################################
#### INTEGRATED POPULATION MODEL FOR SVALBARD ARCTIC FOX ####
#############################################################

library(coda)
library(reshape2)
library(ggplot2)
library(viridis)
library(plyr)
library(gridExtra)

## Load workspace containing all data
#load('200228_AF_IPM_Data.RData')
IPM.data <- readRDS("AF_IPM_Data.rds") 

## Load posterior samples
# load('200429_AF_IPM_VersionB4.RData')
AF.IPM <- readRDS("AF_IPM.rds")
ls()

## Re-arrange data
out.mat <- as.matrix(AF.IPM)
data <- melt(out.mat)
colnames(data) <- c('index', 'parameter', 'value')

## Summarise posterior into median and 50% and 90% CI
data.sum <- ddply(data, .(parameter), summarise, median = median(value, na.rm = T), lCI_90 = quantile(value, probs = 0.05, na.rm = T), uCI_90 = quantile(value, probs = 0.95, na.rm = T), lCI_50 = quantile(value, probs = 0.25, na.rm = T), uCI_50 = quantile(value, probs = 0.75, na.rm = T))

#-----------------------------#
# PLOT POPULATION SIZE ~ TIME #
#-----------------------------#

## Subset data
popN <- paste('N.tot[', c(1:23), ']', sep = '')
data.Ntot <- subset(data.sum, parameter%in%popN)

## Add time
data.Ntot$indexT <- c(1,10:19, 2, 20:23, 3:9)
data.Ntot$Year <- data.Ntot$indexT+1996
data.Ntot <- data.Ntot[order(data.Ntot$Year),]

## Add harvest data
data.Ntot$HarvestData <- c(colSums(IPM.data$C), NA)

## Plot - Estimate with 95% CI
axis <- scale_x_continuous(breaks = c(1997:2019), limits = range(data.Ntot$Year))

plot.Ntot <- ggplot(data.Ntot, aes(x = Year, y = median)) + 
              geom_line(color = '#5C566B') + 
              geom_ribbon(aes(ymin = lCI_90, ymax = uCI_90), alpha = 0.5, fill = '#5C566B') + 
              ylab('Population size (females)') + 
              scale_x_continuous(breaks = c(1997:2019), limits = range(data.Ntot$Year)) + 
              theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))


## Plot - Estimate with 95% CI + harvest data
plot.Htot <- ggplot(data.Ntot, aes(x = Year, y = HarvestData)) + 
              #geom_point(color = viridis(3)[2]) + 
              #geom_line(color = viridis(3)[2], alpha = 0.5) + 
              geom_bar(stat = 'identity', fill = viridis(3)[2]) +
              ylab('Harvest (females)') + 
              scale_x_continuous(breaks = c(1997:2019), limits = range(data.Ntot$Year)) + 
              theme_bw() + theme(panel.grid.minor = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), panel.grid.major.y = element_blank())
#plot.Htot

aligned <- cowplot::align_plots(plot.Htot, plot.Ntot, align = "v", axis = 'b')
grid.arrange(aligned[[1]], aligned[[2]], nrow = 2, heights = c(0.4, 1))


pdf("Ntot_time_Cdata_Rev.pdf", width = 6.5, height = 5)
#grid.arrange(plot.Htot, plot.Ntot, nrow = 2, heights = c(0.4, 1))
grid.arrange(aligned[[1]], aligned[[2]], nrow = 2, heights = c(0.4, 1))
dev.off()

#------------------------------------------#
# LITTER SIZE AND LOCAL RECRUITMENT ~ TIME #
#------------------------------------------#

## Subset data
popR <- paste('R.tot[', c(1:23), ']', sep = '')
data.Rtot <- subset(data.sum, parameter%in%c(popR))

## Add time
data.Rtot$indexT <- c(1:23)
data.Rtot$Year <- data.Rtot$indexT+1996
data.Rtot <- data.Rtot[order(data.Rtot$Year),]

## Add den survey data
dataDS <- data.frame(indexT = IPM.data$DS_year, NoPups = IPM.data$NoPups)

data.Rtot$DenOcc <- c(IPM.data$NoOcc/IPM.data$NoMon)
data.Rtot$SumPups <- ddply(dataDS, .(indexT), summarise, SumPups = sum(NoPups))$SumPups

#--------------------------------------#
# LOCAL RECRUITS VS. IMMIGRANTS ~ TIME #
#--------------------------------------#

## Subset data
popImm <- paste('Imm[', c(2:23), ']', sep = '') 
data.Imm <- subset(data.sum, parameter%in%popImm)

## Add time
data.Imm$indexT <- c(2:23)
data.Imm$Year <- data.Imm$indexT+1996
data.Imm <- data.Imm[order(data.Imm$Year),]

## Merge with recruitment data
data.Imm <- rbind(data.Imm, data.Rtot[,1:8])

## Add factorial variable
data.Imm$Origin <- c(rep('Immigrant', 22), rep('Local recruit', 23))

## Plot - Estimate with 95% CI
pdf("RtotvsImm_time_Rev.pdf", width = 8, height = 4)
ggplot(data.Imm, aes(x = Year, y = median, group = Origin, color = Origin, fill = Origin)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = lCI_90, ymax = uCI_90), alpha = 0.3, color = NA) + 
  ylab('Number of females') + scale_x_continuous(breaks = c(1997:2019)) + 
  scale_color_manual(values = viridis(5)[c(4,2)]) + scale_fill_manual(values = viridis(5)[c(4,2)]) + 
  theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))
dev.off()

# --> Functionality confirmed. 