#############################################################
#### INTEGRATED POPULATION MODEL FOR SVALBARD ARCTIC FOX ####
#############################################################

library(coda)
library(reshape)
library(ggplot2)
library(viridis)
library(plyr)
library(gridExtra)

## Load posterior samples
load('200429_AF_IPM_VersionB4_Bmonitor.RData')
ls()

## Re-arrange data
out.mat <- as.matrix(AF.IPM.varB4)
data <- melt(out.mat)
colnames(data) <- c('index', 'parameter', 'value')

## Summarise posterior into median and 50% and 90% CI
data.sum <- ddply(data, .(parameter), summarise, median = median(value, na.rm = T), lCI_90 = quantile(value, probs = 0.05, na.rm = T), uCI_90 = quantile(value, probs = 0.95, na.rm = T), lCI_50 = quantile(value, probs = 0.25, na.rm = T), uCI_50 = quantile(value, probs = 0.75, na.rm = T))


#--------------------------------------#
# BREEDING POPULATION STRUCTURE ~ TIME #
#--------------------------------------#

## Subset data
data.Bat <- data.sum[49:140,]

## Add time
data.Bat$indexT <- rep(c(1, 10:19, 2, 20:23, 3:9), 4)
data.Bat$Year <- data.Bat$indexT+1996
data.Bat$AgeClass <- rep(c('1', '2', '3', '4+'), each = 23)
data.Bat <- data.Bat[order(data.Bat$Year),]

p1 <- ggplot(data.Bat, aes(x = Year, y = median, group = AgeClass, color = AgeClass, fill = AgeClass)) + geom_area(position = 'stack') + ylab('Breeding population size') + scale_x_continuous(breaks = c(1997:2019)) + scale_color_manual(values = cividis(5)[2:5]) + scale_fill_manual(values = cividis(5)[2:5]) + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))

pdf("Bat_time_stack.pdf", width = 6.5, height = 4)
p1
dev.off()

p2 <- ggplot(data.Bat, aes(x = Year, y = median, group = AgeClass, color = AgeClass, fill = AgeClass)) + geom_area(position = 'fill') + ylab('Age proportion of breeding population') + scale_x_continuous(breaks = c(1997:2019)) + scale_color_manual(values = cividis(5)[2:5]) + scale_fill_manual(values = cividis(5)[2:5]) + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))
pdf("Bat_time_stack_relative.pdf", width = 6.5, height = 4)
p2
dev.off()


gA <- ggplotGrob(p1 + theme(axis.title.x = element_blank(), axis.text.x = element_blank()))
gB <- ggplotGrob(p2)

pdf("Bat_time_stack_combined.pdf", width = 6.5, height = 6)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))
dev.off()



