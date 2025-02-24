#############################################################
#### INTEGRATED POPULATION MODEL FOR SVALBARD ARCTIC FOX ####
#############################################################

library(coda)
library(reshape)
library(ggplot2)
library(viridis)
library(plyr)
library(gridExtra)

## Load workspace containing all data
load('200228_AF_IPM_Data.RData')


## Load posterior samples
load('200429_AF_IPM_VersionB4.RData')
ls()

## Re-arrange data
out.mat <- as.matrix(AF.IPM.varB4)
data <- melt(out.mat)
colnames(data) <- c('index', 'parameter', 'value')

#---------------------------------------------------------------#
# FUNCTIONS FOR CALCULATION OF PROPORTION OF VARIANCE EXPLAINED #
#---------------------------------------------------------------#

## Juvenile background mortality
varexpl.mOj = function(MCMC.mat, Tmax, i){
	
  var.SeaIce = var(MCMC.mat[i, "betaSI.mO"]*IPM.data$JanJunSeaIceIsfj[2:Tmax])
  var.RdCarcass = var(MCMC.mat[i, "betaRC.mO"]*IPM.data$RdCarcass[2:Tmax])
  var.Goose = var(MCMC.mat[i, "betaG.mO"]*IPM.data$GooseJuvProp[1:(Tmax-1)])
  var.Trend = var(MCMC.mat[i, "betaY.mO"]*(c(1:(Tmax-1))-12))
  var.Random = MCMC.mat[i, "sigma.mO"]^2
  
  var.Total = var.SeaIce + var.RdCarcass + var.Goose + var.Trend + var.Random
  
  output <- data.frame(sample = i, 
  SeaIce = var.SeaIce/var.Total,
  RdCarcass = var.RdCarcass/var.Total,
  Goose = var.Goose/var.Total,
  Trend = var.Trend/var.Total,
  Random = var.Random/var.Total
  )
  
  rownames(output) <- NULL
  
  return(output)
}


## Adult background mortality
varexpl.mOa = function(MCMC.mat, Tmax, i){
	
  var.SeaIce = var(MCMC.mat[i, "betaSI.mO"]*IPM.data$JanJunSeaIceIsfj[2:Tmax])
  var.RdCarcass = var(MCMC.mat[i, "betaRC.mO"]*IPM.data$RdCarcass[2:Tmax])
  var.Trend = var(MCMC.mat[i, "betaY.mO"]*(c(1:(Tmax-1))-12))
  var.Random = MCMC.mat[i, "sigma.mO"]^2
  
  var.Total = var.SeaIce + var.RdCarcass + var.Trend + var.Random
  
  output <- data.frame(sample = i, 
  SeaIce = var.SeaIce/var.Total,
  RdCarcass = var.RdCarcass/var.Total,
  Trend = var.Trend/var.Total,
  Random = var.Random/var.Total
  )
  
  rownames(output) <- NULL
  
  return(output)
}

## Pregnancy rate
varexpl.Psi = function(MCMC.mat, Tmax, i){
	
  var.SeaIce = var(MCMC.mat[i, "betaSI.Psi"]*IPM.data$JanJunSeaIceIsfj[1:Tmax])
  var.RdCarcass = var(MCMC.mat[i, "betaRC.Psi"]*IPM.data$RdCarcass[1:Tmax])
  var.Trend = var(MCMC.mat[i, "betaY.Psi"]*(c(1:Tmax)-12))
  var.Random = MCMC.mat[i, "sigma.Psi"]^2
  
  var.Total = var.SeaIce + var.RdCarcass + var.Trend + var.Random
  
  output <- data.frame(sample = i, 
  SeaIce = var.SeaIce/var.Total,
  RdCarcass = var.RdCarcass/var.Total,
  Trend = var.Trend/var.Total,
  Random = var.Random/var.Total
  )
  
  rownames(output) <- NULL
  
  return(output)
}

## Fetus number rate
varexpl.rho = function(MCMC.mat, Tmax, i){
	
  var.SeaIce = var(MCMC.mat[i, "betaSI.rho"]*IPM.data$JanJunSeaIceIsfj[1:Tmax])
  var.RdCarcass = var(MCMC.mat[i, "betaRC.rho"]*IPM.data$RdCarcass[1:Tmax])
  var.Trend = var(MCMC.mat[i, "betaY.rho"]*(c(1:Tmax)-12))
  var.Random = MCMC.mat[i, "sigma.rho"]^2
  
  var.Total = var.SeaIce + var.RdCarcass + var.Trend + var.Random
  
  output <- data.frame(sample = i, 
  SeaIce = var.SeaIce/var.Total,
  RdCarcass = var.RdCarcass/var.Total,
  Trend = var.Trend/var.Total,
  Random = var.Random/var.Total
  )
  
  rownames(output) <- NULL
  
  return(output)
}

## Denning mortality
varexpl.m0 = function(MCMC.mat, Tmax, i){
	
  var.SeaIce = var(MCMC.mat[i, "betaSI.m0"]*IPM.data$JanJunSeaIceIsfj[1:Tmax])
  var.RdCarcass = var(MCMC.mat[i, "betaRC.m0"]*IPM.data$RdCarcass[1:Tmax])
  var.Trend = var(MCMC.mat[i, "betaY.m0"]*(c(1:Tmax)-12))
  var.Random = MCMC.mat[i, "sigma.m0"]^2
  
  var.Total = var.SeaIce + var.RdCarcass + var.Trend + var.Random
  
  output <- data.frame(sample = i,
  SeaIce = var.SeaIce/var.Total,
  RdCarcass = var.RdCarcass/var.Total,
  Trend = var.Trend/var.Total,
  Random = var.Random/var.Total
  )
  
  rownames(output) <- NULL
  
  return(output)
}


#---------------------------------------------#
# CALCULATE PROPORTIONS OF VARIANCE EXPLAINED #
#---------------------------------------------#

## Running the functions
mOj.props <- do.call("rbind", sapply(1:nrow(out.mat), FUN = function(i) varexpl.mOj(out.mat, 23, i), simplify = FALSE))

mOa.props <- do.call("rbind", sapply(1:nrow(out.mat), FUN = function(i) varexpl.mOa(out.mat, 23, i), simplify = FALSE))

Psi.props <- do.call("rbind", sapply(1:nrow(out.mat), FUN = function(i) varexpl.Psi(out.mat, 23, i), simplify = FALSE))

rho.props <- do.call("rbind", sapply(1:nrow(out.mat), FUN = function(i) varexpl.rho(out.mat, 23, i), simplify = FALSE))

m0.props <- do.call("rbind", sapply(1:nrow(out.mat), FUN = function(i) varexpl.m0(out.mat, 23, i), simplify = FALSE))


## Re-arrange and merge data
mOj.data <- melt(mOj.props, id.vars = 'sample')
mOj.data$Parameter <- 'Juvenile natural mortality'

mOa.data <- melt(mOa.props, id.vars = 'sample')
mOa.data$Parameter <- 'Adult natural mortality'

Psi.data <- melt(Psi.props, id.vars = 'sample')
Psi.data$Parameter <- 'Pregnancy rate'

rho.data <- melt(rho.props, id.vars = 'sample')
rho.data$Parameter <- 'Fetus number'

m0.data <- melt(m0.props, id.vars = 'sample')
m0.data$Parameter <- 'Denning mortality'

results <- rbind(mOj.data, mOa.data, Psi.data, rho.data, m0.data)

## Save results
save(results, file = '200503_VarianceDecomposition.RData')


#--------------#
# PLOT RESULTS #
#--------------#

## Multipanel plots of all vital rates

# Overlapping densities
ggplot(results, aes(x = value, group = variable)) + geom_density(aes(fill = variable, color = variable), alpha = 0.3) + scale_fill_viridis(discrete = T) + scale_color_viridis(discrete = T) + facet_wrap(~Parameter, ncol = 1, scales = 'free') + theme_bw() + theme(panel.grid = element_blank())

# Violins
ggplot(results, aes(x = variable, y = value, group = variable)) + geom_violin(aes(fill = variable, color = variable), alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + xlab('') + scale_fill_viridis(discrete = T) + scale_color_viridis(discrete = T) + facet_wrap(~Parameter, ncol = 1, scales = 'free_y') + theme_bw() + theme(legend.position = 'none', panel.grid = element_blank()) 

# Half-Violins 
ggplot(results, aes(x = variable, y = value, group = variable)) + geom_violinhalf(aes(fill = variable, color = variable), alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + xlab('') + scale_fill_viridis(discrete = T) + scale_color_viridis(discrete = T) + facet_wrap(~Parameter, ncol = 1, scales = 'free_y') + theme_bw() + theme(legend.position = 'none', panel.grid = element_blank()) 

# Half-Violins - flipped
ggplot(results, aes(x = variable, y = value, group = variable)) + geom_violinhalf(aes(fill = variable, color = variable), alpha = 0.5, scale = 'width', trim = T) + xlab('') + scale_fill_viridis(discrete = T) + scale_color_viridis(discrete = T) + facet_wrap(~Parameter, nrow = 2, scales = 'free_y') + theme_bw() + theme(panel.grid = element_blank(), axis.text.y = element_blank()) + coord_flip()


# --> The results here are really not very conclusive, as there is very heavy skew towards 0 for everything (except RV in pregnancy rate)
# --> The best way to visualize this may be using horizontal half-violins, but this plot should not be more than a supplementary figure

## Final plot
results$Parameter <- factor(results$Parameter, levels = c('Juvenile natural mortality', 'Adult natural mortality', 'Denning mortality', 'Pregnancy rate', 'Fetus number'))

pdf('200503_VarianceDecomposition_Violins.pdf', width = 8*0.8, height = 5*0.8)
ggplot(results, aes(x = variable, y = value, group = variable)) + geom_violinhalf(aes(fill = variable), color = NA, scale = 'width', trim = T) + xlab('') + scale_fill_viridis(discrete = T, labels = c('Sea ice', 'Reindeer carcasses', 'Goose reproduction', 'Time trend', 'Random year effects')) + facet_wrap(~Parameter, nrow = 2, scales = 'free_y') + ylab('Proportion of variance explained') + theme_bw() + theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = c(0.8, 0.25), legend.title = element_blank()) + coord_flip() + guides(fill = guide_legend(reverse=T))
dev.off()


#-------------------#
# EXTRACT QUANTILES #
#-------------------#

results.sum <- ddply(results, .(Parameter, variable), summarise, 
median = median(value), 
q_0.05 = quantile(value, probs = 0.05), 
q_0.95 = quantile(value, probs = 0.95),
q_0.25 = quantile(value, probs = 0.25),
q_0.75 = quantile(value, probs = 0.75))

# Based on medians:
# mO: RV > SeaIce > RdCarcass > (Goose >) Trend
# Psi: RV > Trend > RdCarcass > SeaIce > WinterTemp
# rho: RV > Trend > SeaIce > WinterTemp > RdCarcass
# m0: RV > RdCarcass > Trend = SeaIce

