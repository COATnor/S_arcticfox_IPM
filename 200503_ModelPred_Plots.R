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

write.csv(data.sum, 'PosteriorSummary.csv')

## Summary stats

# Harvest in second period
quantile(exp(log(out.mat[,'Mu.mH[2]']) + out.mat[,'betaHP.mH']), probs = c(0.5, 0.05, 0.95))
quantile(exp(log(out.mat[,'Mu.mH[1]']) + out.mat[,'betaHP.mH']), probs = c(0.5, 0.05, 0.95))

# % reduction from first to second period
quantile(1-(exp(log(out.mat[,'Mu.mH[2]']) + out.mat[,'betaHP.mH'])/out.mat[,'Mu.mH[2]']), probs = c(0.5, 0.05, 0.95))
quantile(1-(exp(log(out.mat[,'Mu.mH[1]']) + out.mat[,'betaHP.mH'])/out.mat[,'Mu.mH[1]']), probs = c(0.5, 0.05, 0.95))

# Survival in both periods
quantile(exp(-(out.mat[,'Mu.mH[2]'] + out.mat[,'Mu.mO[2]'])), probs = c(0.5, 0.05, 0.95))
quantile(exp(-(out.mat[,'Mu.mH[1]'] + out.mat[,'Mu.mO[1]'])), probs = c(0.5, 0.05, 0.95))

quantile(exp(-(exp(log(out.mat[,'Mu.mH[2]']) + out.mat[,'betaHP.mH']) + out.mat[,'Mu.mO[2]'])), probs = c(0.5, 0.05, 0.95))
quantile(exp(-(exp(log(out.mat[,'Mu.mH[1]']) + out.mat[,'betaHP.mH']) + out.mat[,'Mu.mO[1]'])), probs = c(0.5, 0.05, 0.95))

# Age-dependent fetus numbers

# Age class 1 (1-2 yr old)
quantile(exp(log(out.mat[,'mean.rho']) + out.mat[,'a.eff1']*2), probs = c(0.5, 0.05, 0.95))

# Age class 4+(older than 4)
quantile(exp(log(out.mat[,'mean.rho']) + out.mat[,'a.eff1']*5), probs = c(0.5, 0.05, 0.95))


# Average population size
mean.Ntot <- rep(NA, nrow(out.mat))
for(i in 1:nrow(out.mat)){
	Ntot <- out.mat[i, paste("N.tot[", 1:23, "]", sep = "")]
	mean.Ntot[i] <- mean(Ntot)
}
quantile(Ntot, probs = c(0.05, 0.5, 0.95))

# Population growth rate
lambda <- matrix(NA, nrow = 22, ncol = 3)
for(t in 1:22){
	Ntot <- out.mat[, paste("N.tot[", t, "]", sep = "")]
	Ntot_next <- out.mat[, paste("N.tot[", t+1, "]", sep = "")]
	lambda[t,] <- quantile(Ntot_next/Ntot, probs = c(0.5, 0.05, 0.95), na.rm = TRUE)
}

# Average population growth rate and SD
mean.lambda <- SD.lambda <- rep(NA, nrow(out.mat))

for(i in 1:nrow(out.mat)){
	
	lambda <- rep(NA, 22)
	for(t in 1:22){
		lambda[t] <- out.mat[i, paste("N.tot[", t+1, "]", sep = "")]/out.mat[i, paste("N.tot[", t, "]", sep = "")]
	}
	
	mean.lambda[i] <- mean(lambda)
	SD.lambda[i] <- sd(lambda)
}
quantile(mean.lambda, probs = c(0.05, 0.5, 0.95), na.rm = TRUE)
quantile(SD.lambda, probs = c(0.05, 0.5, 0.95), na.rm = TRUE)

# Population proportions
mean.p1 <- mean.p2 <- mean.p3 <- mean.p4 <- mean.p5 <- rep(NA, nrow(out.mat))
for(i in 1:nrow(out.mat)){
	p1 <- out.mat[i, paste("N[1, ", 1:23, "]", sep = "")]/out.mat[i, paste("N.tot[", 1:23, "]", sep = "")]
	p2 <- out.mat[i, paste("N[2, ", 1:23, "]", sep = "")]/out.mat[i, paste("N.tot[", 1:23, "]", sep = "")]
	p3 <- out.mat[i, paste("N[3, ", 1:23, "]", sep = "")]/out.mat[i, paste("N.tot[", 1:23, "]", sep = "")]
	p4 <- out.mat[i, paste("N[4, ", 1:23, "]", sep = "")]/out.mat[i, paste("N.tot[", 1:23, "]", sep = "")]
	p5 <- out.mat[i, paste("N[5, ", 1:23, "]", sep = "")]/out.mat[i, paste("N.tot[", 1:23, "]", sep = "")]
	
	mean.p1[i] <- mean(p1)
	mean.p2[i] <- mean(p2)
	mean.p3[i] <- mean(p3)
	mean.p4[i] <- mean(p4)
	mean.p5[i] <- mean(p5)
}
quantile(mean.p1, probs = c(0.05, 0.5, 0.95), na.rm = TRUE)
quantile(mean.p2, probs = c(0.05, 0.5, 0.95), na.rm = TRUE)
quantile(mean.p3, probs = c(0.05, 0.5, 0.95), na.rm = TRUE)
quantile(mean.p4, probs = c(0.05, 0.5, 0.95), na.rm = TRUE)
quantile(mean.p5, probs = c(0.05, 0.5, 0.95), na.rm = TRUE)

#------------------------------------------------#
# ASSESS CI OVERLAP WITH 0 FOR ALL FIXED EFFECTS #
#------------------------------------------------#

## Effects with 90% CI not overlapping 0
# - HPeriod -> mH
# - SeaIce -> mO
# - Trend -> Psi

## Effects with 50% CI not overlapping 0
# - GooseRep -> juv. mO
# - RdCarcass -> m0, mO, Psi
# - SeaIce -> Psi, rho
# - Trend -> m0, rho

## Effects with 50% CI overlapping 0
# - RdCarass -> rho
# - SeaIce -> m0
# - Trend -> mO


#-----------------------------#
# PLOT POPULATION SIZE ~ TIME #
#-----------------------------#

# ## Subset data
# popN <- paste('N.tot[', c(1:23), ']', sep = '')
# data.Ntot <- subset(data.sum, parameter%in%popN)
# 
# ## Add time
# data.Ntot$indexT <- c(1,10:19, 2, 20:23, 3:9)
# data.Ntot$Year <- data.Ntot$indexT+1996
# data.Ntot <- data.Ntot[order(data.Ntot$Year),]
# 
# ## Add harvest data
# data.Ntot$HarvestData <- c(colSums(IPM.data$C), NA)
# 
# ## Plot - Estimate with 95% CI
# plot.Ntot <- ggplot(data.Ntot, aes(x = Year, y = median)) + geom_line(color = '#5C566B') + geom_ribbon(aes(ymin = lCI_90, ymax = uCI_90), alpha = 0.5, fill = '#5C566B') + ylab('Total population size') + scale_x_continuous(breaks = c(1997:2019)) + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))
# 
# pdf("Ntot_time.pdf", width = 6.5, height = 4)
# plot.Ntot
# dev.off()
# 
# ## Plot - Estimate with 95% CI + harvest data
# plot.Htot <- ggplot(data.Ntot, aes(x = Year, y = HarvestData)) + geom_point(color = viridis(3)[2]) + geom_line(color = viridis(3)[2], alpha = 0.5) + ylab('# harvested') + scale_x_continuous(breaks = c(1997:2019)) + theme_bw() + theme(panel.grid.minor = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), panel.grid.major.y = element_blank())
# 
# pdf("Ntot_time_Cdata.pdf", width = 6.5, height = 5)
# grid.arrange(plot.Htot, plot.Ntot, nrow = 2, heights = c(0.4, 1))
# dev.off()

# NOTE: Use revised version from 210117_ModelPred_Plots_RevUpdate.R


#-----------------------------#
# POPULATION STRUCTURE ~ TIME #
#-----------------------------#

## Subset data
N.params <- paste0("N[", rep(c(1:5), each = 23), ", ", rep(c(1:23), 5), "]")
data.Nat <- subset(data.sum, parameter %in% N.params)

## Add time
data.Nat$indexT <- rep(c(1:23), each = 5)
data.Nat$Year <- data.Nat$indexT+1996
data.Nat$AgeClass <- rep(c('0', '1', '2', '3', '4+'), 23)
data.Nat <- data.Nat[order(data.Nat$Year),]

## Add harvest data
data.Nat$HarvestData <- c(as.vector(IPM.data$C), rep(NA, 5))

## Plot - Lines with 90% CIs
pdf("Nat_time.pdf", width = 6.5, height = 4)
ggplot(data.Nat, aes(x = Year, y = median, group = AgeClass, color = AgeClass, fill = AgeClass)) + geom_line() + geom_ribbon(aes(ymin = lCI_90, ymax = uCI_90), alpha = 0.3, color = NA) + ylab('Population size by age class') + scale_x_continuous(breaks = c(1997:2019)) + scale_color_viridis(discrete = T, option = 'cividis') + scale_fill_viridis(discrete = T, option = 'cividis') + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))
dev.off()

## Plot - Stacked lines (without uncertainty)
pdf("Nat_time_stack.pdf", width = 6.5, height = 4)
ggplot(data.Nat, aes(x = Year, y = median, group = AgeClass, color = AgeClass, fill = AgeClass)) + geom_area(position = 'stack') + ylab('Population size by age class') + scale_x_continuous(breaks = c(1997:2019)) + scale_color_viridis(discrete = T, option = 'cividis') + scale_fill_viridis(discrete = T, option = 'cividis') + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))
dev.off()

pdf("Nat_time_stack_relative.pdf", width = 6.5, height = 4)
ggplot(data.Nat, aes(x = Year, y = median, group = AgeClass, color = AgeClass, fill = AgeClass)) + geom_area(position = 'fill') + ylab('Population size by age class') + scale_x_continuous(breaks = c(1997:2019)) + scale_color_viridis(discrete = T, option = 'cividis') + scale_fill_viridis(discrete = T, option = 'cividis') + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))
dev.off()


## Plot - Lines with 90% CIs and harvest data (separated by age class)

plot.H1 <- ggplot(subset(data.Nat, AgeClass == '0'), aes(x = Year, y = HarvestData)) + geom_point(color = cividis(5)[1]) + geom_line(color = cividis(5)[1], alpha = 0.5) + ylab('# harvested') + scale_x_continuous(breaks = c(1997:2019)) + theme_bw() + theme(panel.grid.minor = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), panel.grid.major.y = element_blank())

plot.H2 <- ggplot(subset(data.Nat, AgeClass == '1'), aes(x = Year, y = HarvestData)) + geom_point(color = cividis(5)[2]) + geom_line(color = cividis(5)[2], alpha = 0.5) + ylab('# harvested') + scale_x_continuous(breaks = c(1997:2019)) + theme_bw() + theme(panel.grid.minor = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), panel.grid.major.y = element_blank())

plot.H3 <- ggplot(subset(data.Nat, AgeClass == '2'), aes(x = Year, y = HarvestData)) + geom_point(color = cividis(5)[3]) + geom_line(color = cividis(5)[3], alpha = 0.5) + ylab('# harvested') + scale_x_continuous(breaks = c(1997:2019)) + theme_bw() + theme(panel.grid.minor = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), panel.grid.major.y = element_blank())

plot.H4 <- ggplot(subset(data.Nat, AgeClass == '3'), aes(x = Year, y = HarvestData)) + geom_point(color = cividis(5)[4]) + geom_line(color = cividis(5)[4], alpha = 0.5) + ylab('# harvested') + scale_x_continuous(breaks = c(1997:2019)) + theme_bw() + theme(panel.grid.minor = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), panel.grid.major.y = element_blank())

plot.H5 <- ggplot(subset(data.Nat, AgeClass == '4+'), aes(x = Year, y = HarvestData)) + geom_point(color = cividis(5)[5]) + geom_line(color = cividis(5)[5], alpha = 0.5) + ylab('# harvested') + scale_x_continuous(breaks = c(1997:2019)) + theme_bw() + theme(panel.grid.minor = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), panel.grid.major.y = element_blank())

plot.N1 <- ggplot(subset(data.Nat, AgeClass =='0'), aes(x = Year, y = median)) + geom_line(color = cividis(5)[1]) + geom_ribbon(aes(ymin = lCI_90, ymax = uCI_90), alpha = 0.3, color = NA, fill = cividis(5)[1]) + ylab('Number in age class 0') + scale_x_continuous(breaks = c(1997:2019)) + scale_color_viridis(discrete = T, option = 'cividis') + scale_fill_viridis(discrete = T, option = 'cividis') + theme_bw() + theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))

plot.N2 <- ggplot(subset(data.Nat, AgeClass =='1'), aes(x = Year, y = median)) + geom_line(color = cividis(5)[2]) + geom_ribbon(aes(ymin = lCI_90, ymax = uCI_90), alpha = 0.3, color = NA, fill = cividis(5)[2]) + ylab('Number in age class 1') + scale_x_continuous(breaks = c(1997:2019)) + scale_color_viridis(discrete = T, option = 'cividis') + scale_fill_viridis(discrete = T, option = 'cividis') + theme_bw() + theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))

plot.N3 <- ggplot(subset(data.Nat, AgeClass =='2'), aes(x = Year, y = median)) + geom_line(color = cividis(5)[3]) + geom_ribbon(aes(ymin = lCI_90, ymax = uCI_90), alpha = 0.3, color = NA, fill = cividis(5)[3]) + ylab('Number in age class 2') + scale_x_continuous(breaks = c(1997:2019)) + scale_color_viridis(discrete = T, option = 'cividis') + scale_fill_viridis(discrete = T, option = 'cividis') + theme_bw() + theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))

plot.N4 <- ggplot(subset(data.Nat, AgeClass =='3'), aes(x = Year, y = median)) + geom_line(color = cividis(5)[4]) + geom_ribbon(aes(ymin = lCI_90, ymax = uCI_90), alpha = 0.3, color = NA, fill = cividis(5)[4]) + ylab('Number in age class 3') + scale_x_continuous(breaks = c(1997:2019)) + scale_color_viridis(discrete = T, option = 'cividis') + scale_fill_viridis(discrete = T, option = 'cividis') + theme_bw() + theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))

plot.N5 <- ggplot(subset(data.Nat, AgeClass =='4+'), aes(x = Year, y = median)) + geom_line(color = cividis(5)[5]) + geom_ribbon(aes(ymin = lCI_90, ymax = uCI_90), alpha = 0.3, color = NA, fill = cividis(5)[5]) + ylab('Number in age class 4+') + scale_x_continuous(breaks = c(1997:2019)) + scale_color_viridis(discrete = T, option = 'cividis') + scale_fill_viridis(discrete = T, option = 'cividis') + theme_bw() + theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))


pdf("Nat_time_Cdata.pdf", width = 6.5, height = 5)
grid.arrange(plot.H1, plot.N1, nrow = 2, heights = c(0.4, 1))
grid.arrange(plot.H2, plot.N2, nrow = 2, heights = c(0.4, 1))
grid.arrange(plot.H3, plot.N3, nrow = 2, heights = c(0.4, 1))
grid.arrange(plot.H4, plot.N4, nrow = 2, heights = c(0.4, 1))
grid.arrange(plot.H5, plot.N5, nrow = 2, heights = c(0.4, 1))
dev.off()


#--------------------------------------#
# BREEDING POPULATION STRUCTURE ~ TIME #
#--------------------------------------#

# ## Subset data
# data.Nat <- data.sum[338:452,]
# 
# ## Add time
# data.Nat$indexT <- rep(c(1, 10:19, 2, 20:23, 3:9), 5)
# data.Nat$Year <- data.Ntot$indexT+1996
# data.Nat$AgeClass <- rep(c('1', '2', '3', '4', '5+'), each = 23)
# data.Nat <- data.Nat[order(data.Nat$Year),]
# 
# ## Add harvest data
# data.Nat$HarvestData <- c(as.vector(IPM.data$C), rep(NA, 5))
# 
# 
# pdf("Bat_time_stack_relative.pdf", width = 6.5, height = 4)
# ggplot(data.Nat, aes(x = Year, y = median, group = AgeClass, color = AgeClass, fill = AgeClass)) + geom_area(position = 'fill') + ylab('Age proportion of breeding population') + scale_x_continuous(breaks = c(1997:2019)) + scale_color_viridis(discrete = T, option = 'cividis') + scale_fill_viridis(discrete = T, option = 'cividis') + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))
# dev.off()

# NOTE: Use revised version in 200503_BreedingPop_Composition.R

#------------------------------------------#
# LITTER SIZE AND LOCAL RECRUITMENT ~ TIME #
#------------------------------------------#

## Subset data
popR <- paste('R.tot[', c(1:23), ']', sep = '')
data.Rtot <- subset(data.sum, parameter%in%c(popR))

popLS <- paste('meanLS[', c(1:23), ']', sep = '')
data.LS <- subset(data.sum, parameter%in%c(popLS))

## Add time
data.Rtot$indexT <- c(1:23)
data.Rtot$Year <- data.Rtot$indexT+1996
data.Rtot <- data.Rtot[order(data.Rtot$Year),]

data.LS$indexT <- c(1:23)
data.LS$Year <- data.LS$indexT+1996
data.LS <- data.LS[order(data.LS$Year),]


## Add den survey data
dataDS <- data.frame(indexT = IPM.data$DS_year, NoPups = IPM.data$NoPups)

data.Rtot$DenOcc <- c(IPM.data$NoOcc/IPM.data$NoMon)
data.Rtot$SumPups <- ddply(dataDS, .(indexT), summarise, SumPups = sum(NoPups))$SumPups

dataDS.sum <- ddply(dataDS, .(indexT), summarise, AvgPups = mean(NoPups), SdPups = sd(NoPups))
data.LS$AvgPups <- dataDS.sum$AvgPups
data.LS$SdPups <- dataDS.sum$SdPups


## Plot - Recruitment estimate with 95% CI and den occupancy and productivity data
plot.Rtot <- ggplot(data.Rtot, aes(x = Year, y = median)) + geom_line(color = '#5C566B') + geom_ribbon(aes(ymin = lCI_90, ymax = uCI_90), alpha = 0.5, fill = '#5C566B') + ylab('Number of female pups') + scale_x_continuous(breaks = c(1997:2019)) + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5)) + geom_point(aes(x = Year, y = SumPups/2), color = magma(5)[3]) + geom_line(aes(x = Year, y = SumPups/2), color = magma(5)[3], alpha = 0.5) 

plot.DenOcc <- ggplot(data.Rtot, aes(x = Year, y = DenOcc)) + geom_point(color = magma(5)[4]) + geom_line(color = magma(5)[4], alpha = 0.5) + ylab('Den occupancy') + scale_x_continuous(breaks = c(1997:2019)) + theme_bw() + theme(panel.grid.minor = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), panel.grid.major.y = element_blank())


pdf("Rtot_time.pdf", width = 6.5, height = 4)
plot.Rtot
dev.off()

pdf("Rtot_time_DSdata.pdf", width = 6.5, height = 4)
grid.arrange(plot.DenOcc, plot.Rtot, nrow = 2, heights = c(0.4, 1))
dev.off()

## Plot - Average litter size with 95% CI and den survey data (incl. SD)
pdf("LS_time_DSdata.pdf", width = 6.5, height = 4)
ggplot(data.LS, aes(x = Year, y = median)) + geom_line(color = '#5C566B') + geom_ribbon(aes(ymin = lCI_90, ymax = uCI_90), alpha = 0.5, fill = '#5C566B') + ylab('Average litter size') + scale_x_continuous(breaks = c(1997:2019)) + theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5)) + geom_point(aes(x = Year, y = AvgPups), color = magma(5)[3]) + geom_errorbar(aes(ymin = AvgPups-SdPups, ymax = AvgPups+SdPups), color = magma(5)[3], alpha = 0.3)
dev.off()


#--------------------------------------#
# LOCAL RECRUITS VS. IMMIGRANTS ~ TIME #
#--------------------------------------#

# ## Subset data
# popImm <- paste('Imm[', c(2:23), ']', sep = '') 
# data.Imm <- subset(data.sum, parameter%in%popImm)
# 
# ## Add time
# data.Imm$indexT <- c(10:19, 2, 20:23, 3:9)
# data.Imm$Year <- data.Imm$indexT+1996
# data.Imm <- data.Imm[order(data.Imm$Year),]
# 
# ## Merge with recruitment data
# data.Imm <- rbind(data.Imm, data.Rtot[,1:8])
# 
# ## Add factorial variable
# data.Imm$Origin <- c(rep('Immigrant', 22), rep('Local recruit', 23))
# 
# ## Plot - Estimate with 95% CI
# pdf("RtotvsImm_time.pdf", width = 8, height = 4)
# ggplot(data.Imm, aes(x = Year, y = median, group = Origin, color = Origin, fill = Origin)) + geom_line() + geom_ribbon(aes(ymin = lCI_90, ymax = uCI_90), alpha = 0.3, color = NA) + ylab('Age class 1 individuals') + scale_x_continuous(breaks = c(1997:2019)) + scale_color_manual(values = viridis(5)[c(4,2)]) + scale_fill_manual(values = viridis(5)[c(4,2)]) + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))
# dev.off()

# NOTE: Use revised version from 210117_ModelPred_Plots_RevUpdate.R


#-----------------------------#
# JUVENILE VS. ADULT SURVIVAL #
#-----------------------------#

## Organise data for juveniles and adults
mS.data <- data.frame(param = rep(rep(c('S','mH','mO'), each = nrow(out.mat)), 2),
                     stage = rep(c('Adult (1+)','Juvenile'), each = nrow(out.mat)*3), 
                     value = c(exp(-(out.mat[,'Mu.mH[1]']+out.mat[,'Mu.mO[1]'])), out.mat[,'Mu.mH[1]'], out.mat[,'Mu.mO[1]'], exp(-(out.mat[,'Mu.mH[2]']+out.mat[,'Mu.mO[2]'])), out.mat[,'Mu.mH[2]'], out.mat[,'Mu.mO[2]']))


# Plot posterior distributions for both life stages
p0 <- ggplot(subset(mS.data, param == 'S'), aes(x = value, group = stage)) + geom_density(aes(fill = stage, color = stage), alpha = 0.3) + 
  ggtitle('Survival probability') + ylab('Density') + xlab('Estimate') + xlim(0, 0.7) +
  scale_color_manual(values = viridis(7)[c(1,5)]) + scale_fill_manual(values = viridis(7)[c(1,5)]) +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), axis.title.x = element_blank(), legend.title = element_blank(), axis.text.x = element_text(size = 12), plot.title = element_text(face = 'bold'))
p1 <- p0 + theme(legend.position = 'none')
p1

p2 <- ggplot(subset(mS.data, param == 'mH'), aes(x = value, group = stage)) + geom_density(aes(fill = stage, color = stage), alpha = 0.3) + 
  ggtitle('Harvest mortality hazard rate') + ylab('Density') + xlab('Estimate') + xlim(0, 1.5) + 
  scale_color_manual(values = viridis(7)[c(1,5)]) + scale_fill_manual(values = viridis(7)[c(1,5)]) +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), legend.position = 'none', axis.title.x = element_blank(), axis.text.x = element_text(size = 12), plot.title = element_text(face = 'bold'))
p2

p3 <- ggplot(subset(mS.data, param == 'mO'), aes(x = value, group = stage)) + geom_density(aes(fill = stage, color = stage), alpha = 0.3) + 
  ggtitle('Natural mortality hazard rate') + ylab('Density') + xlab('Estimate') + xlim(0, 1.5) + 
  scale_color_manual(values = viridis(7)[c(1,5)]) + scale_fill_manual(values = viridis(7)[c(1,5)]) +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), legend.position = 'none', axis.title.x = element_blank(), axis.text.x = element_text(size = 12), plot.title = element_text(face = 'bold'))
p3

pdf('Mortality&Survival_age.pdf', width = 8.27, height = 5.83) 
grid.arrange(p1 + theme(legend.position = c(0.91, 0.6), legend.title = element_blank(), legend.text = element_text(size = 12)), p2, p3, nrow = 3)
dev.off()

pdf('Survival_age.pdf', width = 8.27/1.5, height = 5.83/1.5) 
p0 + theme(legend.position = c(0.18, 0.85), legend.title = element_blank(), legend.text = element_text(size = 12))
dev.off()


#----------------------------------------------#
# PREGNANCY RATE & PLACENTAL SCARS ~ AGE CLASS #
#----------------------------------------------#

## Making age-dependent predictions for reproductive parameters
ageC <- c(seq(1,2, length.out = 20), seq(2,3, length.out = 20), seq(3,4, length.out = 20), seq(4,5, length.out = 20))
ageC <- ageC[-c(20,40,60)]

rep.predict <- function(MCMC.mat){
  
  Psi.pred <- rho.pred <- data.frame(age = ageC, median = rep(NA, length(ageC)), lCI = rep(NA, length(ageC)), uCI = rep(NA, length(ageC)))
  
  for(i in 1:length(ageC)){
    eta <- plogis(MCMC.mat[,"par.b"]*(MCMC.mat[,"par.c"] - ageC[i]))
    Psi <- MCMC.mat[,"par.a"]*eta
    
    Psi.pred[i,2] <- quantile(Psi, 0.5)
    Psi.pred[i,3] <- quantile(Psi, 0.05)
    Psi.pred[i,4] <- quantile(Psi, 0.95)
    
    rho <- exp(log(MCMC.mat[,'mean.rho']) + MCMC.mat[,'a.eff1']*ageC[i])
    
    rho.pred[i,2] <- quantile(rho, 0.5)
    rho.pred[i,3] <- quantile(rho, 0.05)
    rho.pred[i,4] <- quantile(rho, 0.95)
  }
  return(list(Psi = Psi.pred, rho = rho.pred))
}

RepPred <- rep.predict(out.mat)

## Plotting: Pregnancy rate
pdf('PregnancyRate_age.pdf', width = 5, height = 3.5) 
ggplot(RepPred$Psi, aes(x = age, y = median)) + geom_line(color = viridis(7)[3]) + geom_ribbon(aes(ymin = lCI, ymax = uCI), color = NA, fill = viridis(7)[3], alpha = 0.5) + ylab('Pregnancy rate') + xlab('Age class') + scale_x_continuous(breaks = c(1:5), labels = c('0','1','2','3','4+'), limits = c(2,5)) + theme_bw() + theme(panel.grid = element_blank())
dev.off()

## Plotting: Pregnancy rate with raw data
P2.data <- data.frame(Preg = IPM.data$P2, Age = IPM.data$P2_age)
P2.data.sum <- ddply(P2.data, .(Age), summarise, Mean.Preg = mean(Preg))

pdf('PregnancyRate_age_P2Data.pdf', width = 5, height = 3.5) 
ggplot(RepPred$Psi, aes(x = age, y = median)) + geom_line(color = viridis(7)[3]) + geom_ribbon(aes(ymin = lCI, ymax = uCI), color = NA, fill = viridis(7)[3], alpha = 0.5) + ylab('Pregnancy rate') + xlab('Age class') + geom_jitter(data = P2.data, aes(x = Age, y = Preg), color = magma(5)[3], alpha = 0.2, height = 0.1, width = 0.2) + geom_point(data = P2.data.sum, aes(x = Age, y = Mean.Preg), color = magma(5)[3], shape = 17, size = 2) + 
scale_x_continuous(breaks = c(1:5), labels = c('0','1','2','3','4+'), limits = c(1.8,5.2)) + 
theme_bw() + theme(panel.grid = element_blank())
dev.off()



## Plotting: Number of placental scars
P1.data <- data.frame(rho = IPM.data$P1, Age = IPM.data$P1_age)
P1.data.sum <- ddply(P1.data, .(Age), summarise, Mean.rho = mean(rho))

pdf('PlacentalScars_age.pdf', width = 5, height = 3.5) 
ggplot(RepPred$rho, aes(x = age, y = median)) + geom_line(color = viridis(7)[2]) + geom_ribbon(aes(ymin = lCI, ymax = uCI), color = NA, fill = viridis(7)[2], alpha = 0.5) + ylab('Fetus number') + xlab('Age class') + scale_x_continuous(breaks = c(1:5), labels = c('0','1','2','3','4+'), limits = c(2,5)) + theme_bw() + theme(panel.grid = element_blank())
dev.off()

pdf('PlacentalScars_age_P1data.pdf', width = 5, height = 3.5) 
ggplot(subset(RepPred$rho, age > 1), aes(x = age, y = median)) + geom_line(color = viridis(7)[2]) + geom_ribbon(aes(ymin = lCI, ymax = uCI), color = NA, fill = viridis(7)[2], alpha = 0.5) + ylab('Fetus number') + xlab('Age class') + geom_point(data = P1.data, aes(x = Age, y = rho), color = magma(5)[3], alpha = 0.3, position = position_jitter(w = 0.2, h = 0)) + geom_point(data = P1.data.sum, aes(x = Age, y = Mean.rho), color = magma(5)[3], shape = 17, size = 2) +
scale_x_continuous(breaks = c(2:5), labels = c('1','2','3','4+'), limits = c(1.8,5.2)) + theme_bw() + theme(panel.grid = element_blank())
dev.off()


#--------------------#
# VITAL RATES ~ TIME #
#--------------------#

## Predicting mortality, survival, pregnancy rate, and placental scars over time
VRtime.predict <- function(MCMC.mat, Tmax, RepAgeClass){
  
  time <- c(1:Tmax)
  
  Psi.pred <- rho.pred <- S0.pred <- data.frame(Year = time+1996, median = rep(NA, length(time)), lCI = rep(NA, length(time)), uCI = rep(NA, length(time)))
  
  mH.pred.j <- mO.pred.j <- S.pred.j <- data.frame(Year = time+1996, median = rep(NA, length(time)), lCI = rep(NA, length(time)), uCI = rep(NA, length(time)), LifeStage = 'Juvenile')
  
  mH.pred.a <- mO.pred.a <- S.pred.a <- data.frame(Year = time+1996, median = rep(NA, length(time)), lCI = rep(NA, length(time)), uCI = rep(NA, length(time)), LifeStage = 'Adult (1+)')
  
  
  # Mortality and survival
  for(t in 1:(Tmax-1)){
  
  	mHj <- exp(log(MCMC.mat[,"Mu.mH[2]"]) + MCMC.mat[,"betaHP.mH"]*c(rep(0, 12), rep(1, 11))[t] + MCMC.mat[,paste("epsilon.mH[",t,"]",sep ='')])
  	mHa <- exp(log(MCMC.mat[,"Mu.mH[1]"]) + MCMC.mat[,"betaHP.mH"]*c(rep(0, 12), rep(1, 11))[t] + MCMC.mat[,paste("epsilon.mH[",t,"]",sep ='')])
  	
    mOj <- exp(log(MCMC.mat[,"Mu.mO[2]"]) + MCMC.mat[,"betaT.mO"]*IPM.data$MeanWinterTemp[t+1] + MCMC.mat[,"betaRC.mO"]*IPM.data$RdCarcass[t+1] + MCMC.mat[,"betaGR.mO"]*c(IPM.data$GooseJuvProp[-23],0)[t] + MCMC.mat[,"betaSI.mO"]*IPM.data$JanJunSeaIceIsfj[t+1])
    mOa <- exp(log(MCMC.mat[,"Mu.mO[1]"]) + MCMC.mat[,"betaT.mO"]*IPM.data$MeanWinterTemp[t+1] + MCMC.mat[,"betaRC.mO"]*IPM.data$RdCarcass[t+1] + MCMC.mat[,"betaSI.mO"]*IPM.data$JanJunSeaIceIsfj[t+1])
    
    Sj <- exp(-(mHj + mOj))
    Sa <- exp(-(mHa + mOa))
    
    mH.pred.j[t,2] <- quantile(mHj, 0.5)
    mH.pred.j[t,3] <- quantile(mHj, 0.05)
    mH.pred.j[t,4] <- quantile(mHj, 0.95)
    
    mH.pred.a[t,2] <- quantile(mHa, 0.5)
    mH.pred.a[t,3] <- quantile(mHa, 0.05)
    mH.pred.a[t,4] <- quantile(mHa, 0.95)

    mO.pred.j[t,2] <- quantile(mOj, 0.5)
    mO.pred.j[t,3] <- quantile(mOj, 0.05)
    mO.pred.j[t,4] <- quantile(mOj, 0.95)
    
    mO.pred.a[t,2] <- quantile(mOa, 0.5)
    mO.pred.a[t,3] <- quantile(mOa, 0.05)
    mO.pred.a[t,4] <- quantile(mOa, 0.95) 
    
    S.pred.j[t,2] <- quantile(Sj, 0.5)
    S.pred.j[t,3] <- quantile(Sj, 0.05)
    S.pred.j[t,4] <- quantile(Sj, 0.95)
    
    S.pred.a[t,2] <- quantile(Sa, 0.5)
    S.pred.a[t,3] <- quantile(Sa, 0.05)
    S.pred.a[t,4] <- quantile(Sa, 0.95)      
    
  }
  
  
  mH.pred <- rbind(mH.pred.j, mH.pred.a)
  mO.pred <- rbind(mO.pred.j, mO.pred.a)
  S.pred <- rbind(S.pred.j, S.pred.a)
  
  
  # Pregnancy rate, placental scars, & denning survival
  for(t in 1:Tmax){
  	
    logit.eta <- MCMC.mat[,"par.b"]*(MCMC.mat[,"par.c"] - RepAgeClass) + MCMC.mat[,"betaY.Psi"]*t + MCMC.mat[,"betaRC.Psi"]*IPM.data$RdCarcass[t] +  MCMC.mat[,"betaT.Psi"]*IPM.data$MeanWinterTemp[t] + MCMC.mat[,"betaSI.Psi"]*IPM.data$JanJunSeaIceIsfj[t] + MCMC.mat[,paste("epsilon.Psi[",t,"]",sep ='')]
    
    eta <- plogis(logit.eta)
    
    Psi <- MCMC.mat[,"par.a"]*eta
    
    Psi.pred[t,2] <- quantile(Psi, 0.5)
    Psi.pred[t,3] <- quantile(Psi, 0.05)
    Psi.pred[t,4] <- quantile(Psi, 0.95)
    
    
    log.rho <- log(MCMC.mat[,'mean.rho']) + MCMC.mat[,'a.eff1']*RepAgeClass + MCMC.mat[,"betaSI.rho"]*IPM.data$JanJunSeaIceIsfj[t] + MCMC.mat[,"betaRC.rho"]*IPM.data$RdCarcass[t] + MCMC.mat[,"betaT.rho"]*IPM.data$MeanWinterTemp[t] + MCMC.mat[,paste("epsilon.rho[",t,"]",sep ='')]
    
    rho <- exp(log.rho)
    
    rho.pred[t,2] <- quantile(rho, 0.5)
    rho.pred[t,3] <- quantile(rho, 0.05)
    rho.pred[t,4] <- quantile(rho, 0.95)
    
    log.m0 <- -log(MCMC.mat[,'S0']) + MCMC.mat[,"betaSI.m0"]*IPM.data$JanJunSeaIceIsfj[t] + MCMC.mat[,"betaRC.m0"]*IPM.data$RdCarcass[t] + MCMC.mat[,paste("epsilon.m0[",t,"]",sep ='')]
    
    S0 <- exp(-exp(log.m0))
    
    S0.pred[t,2] <- quantile(S0, 0.5)
    S0.pred[t,3] <- quantile(S0, 0.05)
    S0.pred[t,4] <- quantile(S0, 0.95)
  }
  
  return(list(Psi = Psi.pred, rho = rho.pred, mH = mH.pred, mO = mO.pred, S = S.pred, S0 = S0.pred))
}
#VRPred <- VRtime.predict(out.mat, 23, 3)

## Arranging time-dependent vital rates
VRtime.arrange <- function(MCMC.mat, Tmax, RepAgeClass){
  
  time <- c(1:Tmax)
  
  Psi.pred <- rho.pred <- S0.pred <- data.frame(Year = time+1996, median = rep(NA, length(time)), lCI = rep(NA, length(time)), uCI = rep(NA, length(time)))
  
  mH.pred.j <- mO.pred.j <- S.pred.j <- data.frame(Year = time+1996, median = rep(NA, length(time)), lCI = rep(NA, length(time)), uCI = rep(NA, length(time)), LifeStage = 'Juvenile')
  
  mH.pred.a <- mO.pred.a <- S.pred.a <- data.frame(Year = time+1996, median = rep(NA, length(time)), lCI = rep(NA, length(time)), uCI = rep(NA, length(time)), LifeStage = 'Adult (1+)')
  
  
  # Mortality and survival
  for(t in 1:(Tmax-1)){
  
  	mHj <- MCMC.mat[, paste("mH[2, ", t, "]", sep = "")]
  	mHa <- MCMC.mat[, paste("mH[1, ", t, "]", sep = "")]
  	
    mOj <- MCMC.mat[, paste("mO[2, ", t, "]", sep = "")]
    mOa <- MCMC.mat[, paste("mO[1, ", t, "]", sep = "")]
    
    Sj <- exp(-(mHj + mOj))
    Sa <- exp(-(mHa + mOa))
    
    mH.pred.j[t,2] <- quantile(mHj, 0.5)
    mH.pred.j[t,3] <- quantile(mHj, 0.05)
    mH.pred.j[t,4] <- quantile(mHj, 0.95)
    
    mH.pred.a[t,2] <- quantile(mHa, 0.5)
    mH.pred.a[t,3] <- quantile(mHa, 0.05)
    mH.pred.a[t,4] <- quantile(mHa, 0.95)

    mO.pred.j[t,2] <- quantile(mOj, 0.5)
    mO.pred.j[t,3] <- quantile(mOj, 0.05)
    mO.pred.j[t,4] <- quantile(mOj, 0.95)
    
    mO.pred.a[t,2] <- quantile(mOa, 0.5)
    mO.pred.a[t,3] <- quantile(mOa, 0.05)
    mO.pred.a[t,4] <- quantile(mOa, 0.95) 
    
    S.pred.j[t,2] <- quantile(Sj, 0.5)
    S.pred.j[t,3] <- quantile(Sj, 0.05)
    S.pred.j[t,4] <- quantile(Sj, 0.95)
    
    S.pred.a[t,2] <- quantile(Sa, 0.5)
    S.pred.a[t,3] <- quantile(Sa, 0.05)
    S.pred.a[t,4] <- quantile(Sa, 0.95)      
    
  }
  
  
  mH.pred <- rbind(mH.pred.j, mH.pred.a)
  mO.pred <- rbind(mO.pred.j, mO.pred.a)
  S.pred <- rbind(S.pred.j, S.pred.a)
  
  
  # Pregnancy rate, placental scars, & denning survival
  for(t in 1:Tmax){
    
    Psi <- MCMC.mat[, paste("Psi[", RepAgeClass, ", ", t, "]", sep = "")]
    
    Psi.pred[t,2] <- quantile(Psi, 0.5)
    Psi.pred[t,3] <- quantile(Psi, 0.05)
    Psi.pred[t,4] <- quantile(Psi, 0.95)
    
    rho <- MCMC.mat[, paste("rho[", RepAgeClass, ", ", t, "]", sep = "")]
    
    rho.pred[t,2] <- quantile(rho, 0.5)
    rho.pred[t,3] <- quantile(rho, 0.05)
    rho.pred[t,4] <- quantile(rho, 0.95)
    
    S0 <- exp(-MCMC.mat[, paste("m0t[", t, "]", sep = "")])
    
    S0.pred[t,2] <- quantile(S0, 0.5)
    S0.pred[t,3] <- quantile(S0, 0.05)
    S0.pred[t,4] <- quantile(S0, 0.95)
  }
  
  return(list(Psi = Psi.pred, rho = rho.pred, mH = mH.pred, mO = mO.pred, S = S.pred, S0 = S0.pred))
}

VRPred <- VRtime.arrange(out.mat, 23, 3)

## Plotting: Pregnancy rate
pdf('PregnancyRate_time.pdf', width = 6, height = 3.5) 
ggplot(VRPred$Psi, aes(x = Year, y = median)) + geom_line(color = viridis(7)[3]) + geom_ribbon(aes(ymin = lCI, ymax = uCI), color = NA, fill = viridis(7)[3], alpha = 0.5) + ylab('Pregnancy rate (age class 2)') + scale_x_continuous(breaks = c(1997:2019)) + theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))
dev.off()

## Plotting: Number of placental scars
pdf('PlacentalScars_time.pdf', width = 6, height = 3.5) 
ggplot(VRPred$rho, aes(x = Year, y = median)) + geom_line(color = viridis(7)[2]) + geom_ribbon(aes(ymin = lCI, ymax = uCI), color = NA, fill = viridis(7)[2], alpha = 0.5) + ylab('Fetus number (age class 2)') + scale_x_continuous(breaks = c(1997:2019)) + theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))
dev.off()

## Plotting: Denning survival
pdf('DenningSurvival_time.pdf', width = 6, height = 3.5) 
ggplot(VRPred$S0, aes(x = Year, y = median)) + geom_line(color = viridis(7)[4]) + geom_ribbon(aes(ymin = lCI, ymax = uCI), color = NA, fill = viridis(7)[4], alpha = 0.5) + ylab('Denning survival') + scale_x_continuous(breaks = c(1997:2019)) + theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))
dev.off()

## Plotting: Survival and mortality
p1 <- ggplot(VRPred$S, aes(x = Year, y = median)) + geom_line(aes(color = LifeStage)) + geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = LifeStage), color = NA, alpha = 0.5) + ylab('Estimate') + ggtitle('Survival probability') + scale_x_continuous(breaks = c(1997:2019)) + scale_color_manual(values = viridis(7)[c(5,1)]) + scale_fill_manual(values = viridis(7)[c(5,1)]) + theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), axis.text.x = element_blank(), legend.position = 'none', axis.title.x = element_blank(), plot.title = element_text(face = 'bold'))

p2 <- ggplot(VRPred$mH, aes(x = Year, y = median)) + geom_line(aes(color = LifeStage)) + geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = LifeStage), color = NA, alpha = 0.5) + ylab('Estimate') + ggtitle('Harvest mortality hazard rate') + scale_x_continuous(breaks = c(1997:2019)) + scale_color_manual(values = viridis(7)[c(5,1)]) + scale_fill_manual(values = viridis(7)[c(5,1)]) + theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), axis.text.x = element_blank(), legend.position = 'none', axis.title.x = element_blank(), plot.title = element_text(face = 'bold'))

p3 <- ggplot(VRPred$mO, aes(x = Year, y = median)) + geom_line(aes(color = LifeStage)) + geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = LifeStage), color = NA, alpha = 0.5) + ylab('Estimate') + ggtitle('Natural mortality hazard rate') + scale_x_continuous(breaks = c(1997:2019)) + scale_color_manual(values = viridis(7)[c(5,1)]) + scale_fill_manual(values = viridis(7)[c(5,1)]) + theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5, size = 12), legend.position = 'none', axis.title.x = element_blank(), plot.title = element_text(face = 'bold'))


pdf('Mortality&Survival_time.pdf', width = 8, height = 7) 
grid.arrange(p1, p2 + theme(legend.position = c(0.85, 0.75)), p3, nrow = 3)
dev.off()


#----------------------------------------#
# VITAL RATES ~ ENVIRONMENTAL COVARIATES #
#----------------------------------------#

## Making vital rate predictions depening on (standardized) env. covariates
SDcov <- seq(-2.5, 2.5, length.out = 100)

cov.predict <- function(MCMC.mat, effect, RepAgeClass){
  
  # Prepare data frames
  Psi.pred <- rho.pred <- mO.pred <- S0.pred <- data.frame(SDcov = SDcov, median = rep(NA, length(SDcov)), lCI = rep(NA, length(SDcov)), uCI = rep(NA, length(SDcov)), Covariate = effect)
  
  Psi.pred$VitalRate <- 'Pregnancy rate'
  rho.pred$VitalRate <- 'Fetus number'
  mO.pred$VitalRate <- 'Natural mortality'
  S0.pred$VitalRate <- 'Denning survival'
  
  # Set effect sizes  
  if(effect == 'RdCarcass'){
  	beta.Psi <- MCMC.mat[,"betaRC.Psi"]
  	beta.rho <- MCMC.mat[,"betaRC.rho"]
  	beta.mO <- MCMC.mat[,"betaRC.mO"]
  	beta.m0 <- MCMC.mat[,"betaRC.m0"]
  }

  if(effect == 'SeaIce'){
  	beta.Psi <- MCMC.mat[,"betaSI.Psi"]
  	beta.rho <- MCMC.mat[,"betaSI.rho"]
  	beta.mO <- MCMC.mat[,"betaSI.mO"]
  	beta.m0 <- MCMC.mat[,"betaSI.m0"]
  }
  
  if(effect == 'Goose'){
  	beta.Psi <- 0
  	beta.rho <- 0
  	beta.mO <- MCMC.mat[,"betaG.mO"]
  	beta.m0 <- 0
  }
  
  
  # Make predictions 
  for(i in 1:length(SDcov)){
  	
  	# Pregnacy rate
    eta <- plogis(MCMC.mat[,"par.b"]*(MCMC.mat[,"par.c"] - RepAgeClass) + beta.Psi*SDcov[i])
    Psi <- MCMC.mat[,"par.a"]*eta
    
    Psi.pred[i,2] <- quantile(Psi, 0.5)
    Psi.pred[i,3] <- quantile(Psi, 0.05)
    Psi.pred[i,4] <- quantile(Psi, 0.95)
    
    # Placental scars
    rho <- exp(log(MCMC.mat[,'mean.rho']) + MCMC.mat[,'a.eff1']*RepAgeClass + beta.rho*SDcov[i])
    
    rho.pred[i,2] <- quantile(rho, 0.5)
    rho.pred[i,3] <- quantile(rho, 0.05)
    rho.pred[i,4] <- quantile(rho, 0.95)
    
    # Natural mortality
    if(effect == 'Goose'){
    	mO <- exp(log(MCMC.mat[,'Mu.mO[2]']) + beta.mO*SDcov[i])
    }else{
    	mO <- exp(log(MCMC.mat[,'Mu.mO[1]']) + beta.mO*SDcov[i])
    }
    
    mO.pred[i,2] <- quantile(mO, 0.5)
    mO.pred[i,3] <- quantile(mO, 0.05)
    mO.pred[i,4] <- quantile(mO, 0.95)
    
    # Denning survival
    m0 <- exp(log(-log(MCMC.mat[,'S0'])) + beta.m0*SDcov[i])
    S0 <- exp(-m0)
    
    S0.pred[i,2] <- quantile(S0, 0.5)
    S0.pred[i,3] <- quantile(S0, 0.05)
    S0.pred[i,4] <- quantile(S0, 0.95)
    
  }
  
  results <- rbind(Psi.pred, rho.pred, mO.pred, S0.pred)
  
  return(results)
}


CarcassEff <- cov.predict(out.mat, 'RdCarcass', 3)
SeaIceEff <- cov.predict(out.mat, 'SeaIce', 3)
GooseEff <- cov.predict(out.mat, 'Goose', 3)


## Plotting: Carcass Effects
pdf('CarcassEffects.pdf', width = 6, height = 5) 
ggplot(CarcassEff, aes(x = SDcov, y = median)) + geom_line(aes(color = VitalRate)) + geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = VitalRate), color = NA, alpha = 0.5) + ylab('Estimate') + xlab('Covariate value') + theme_bw() + theme(panel.grid = element_blank(), legend.position = 'none', plot.title = element_text(hjust = 0.5)) + scale_color_manual(values = c(viridis(7)[6], viridis(7)[2], viridis(7)[1], viridis(7)[3])) + scale_fill_manual(values = c(viridis(7)[6], viridis(7)[2], viridis(7)[1], viridis(7)[3])) + facet_wrap(~VitalRate, scales = 'free_y') + ggtitle('Effects of reindeer carcasses on vital rates')
dev.off()

## Plotting: Sea Ice Effects
pdf('SeaIceEffects.pdf', width = 6, height = 5) 
ggplot(SeaIceEff, aes(x = SDcov, y = median)) + geom_line(aes(color = VitalRate)) + geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = VitalRate), color = NA, alpha = 0.5) + ylab('Estimate') + xlab('Covariate value') + theme_bw() + theme(panel.grid = element_blank(), legend.position = 'none', plot.title = element_text(hjust = 0.5)) + scale_color_manual(values = c(viridis(7)[6], viridis(7)[2], viridis(7)[1], viridis(7)[3])) + scale_fill_manual(values = c(viridis(7)[6], viridis(7)[2], viridis(7)[1], viridis(7)[3])) + facet_wrap(~VitalRate, scales = 'free_y') + ggtitle('Effects of sea ice on vital rates')
dev.off()

## Plotting: Goose Effects
pdf('GooseEffects.pdf', width = 6, height = 5) 
ggplot(GooseEff, aes(x = SDcov, y = median)) + geom_line(aes(color = VitalRate)) + geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = VitalRate), color = NA, alpha = 0.5) + ylab('Estimate') + xlab('Covariate value') + theme_bw() + theme(panel.grid = element_blank(), legend.position = 'none', plot.title = element_text(hjust = 0.5)) + scale_color_manual(values = c('grey95', 'grey95', viridis(7)[5], 'grey95')) + scale_fill_manual(values = c('grey95', 'grey95', viridis(7)[5], 'grey95')) + facet_wrap(~VitalRate, scales = 'free_y') + ggtitle('Effects of goose reproduction on vital rates')
dev.off()


## Plotting: Carcass and sea ice effects
JointEff <- rbind(SeaIceEff, CarcassEff)

pdf('SeaIce&CarcassEffects.pdf', width = 7, height = 5) 
ggplot(JointEff, aes(x = SDcov, y = median, group = Covariate)) + geom_line(aes(color = Covariate, linetype = Covariate)) + geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = Covariate), color = NA, alpha = 0.3) + ylab('Estimate') + xlab('Standardized covariate value') + theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5)) + scale_color_manual(values = c(cividis(4)[1], cividis(4)[3])) + scale_fill_manual(values = c(cividis(4)[1], cividis(4)[3])) + facet_wrap(~VitalRate, scales = 'free_y')
dev.off()

# --> Functionality confirmed. 