library(reshape2)
library(ggplot2)
library(ggridges)
library(viridis)
library(patchwork)
library(plyr)

## Load LTRE results
load('AF_LTRE_Fixed.RData')

str(results) # LTRE with contributions of survival probabilities and mortality hazard rates

###########################
#### REFORMATTING DATA ####
###########################

## Adding summarised contributions to the data frame
results2 <- results
results2$Annual_mortality <- results$cont_mHj+results$cont_mOj+results$cont_mHa+results$cont_mOa
results2$Denning_mortality <- results$cont_m0
results2$Annual_survival <- results$cont_Sj+results$cont_Sa
results2$Denning_survival <- results$cont_S0
results2$Pregnancy_rate <- results$cont_Psi2+results$cont_Psi3+results$cont_Psi4+results$cont_Psi5
results2$Fetus_number <- results$cont_rho2+results$cont_rho3+results$cont_rho4+results$cont_rho5
results2$Population_structure <- results$cont_N_1+results$cont_N_2+results$cont_N_3+results$cont_N_4+results$cont_N_5
results2$Immigration <- results$cont_Imm

results2 <- results2[,23:32]

## Transforming data frames into vertical format
LTRE.data <- melt(results, id.vars = c('t1', 't2'))
LTRE.data.sum <- melt(results2, id.vars = c('t1', 't2'))

## Add year
LTRE.data$Year <- LTRE.data$t1 + 1996
LTRE.data.sum$Year <- LTRE.data.sum$t1 + 1996

## Make summaries of posteriors (medians + 95% CIs)
LTRE.quant <- ddply(LTRE.data, .(Year, variable), summarise, median = median(value), lCI = quantile(value, probs = 0.025), uCI = quantile(value, probs = 0.975))

LTRE.quant.sum <- ddply(LTRE.data.sum, .(Year, variable), summarise, median = median(value), lCI = quantile(value, probs = 0.025), uCI = quantile(value, probs = 0.975))


## Make subsets without survival
LTRE.quant.H <- subset(LTRE.quant,!(variable%in%c('cont_Sj', 'cont_Sa', 'cont_S0'))) 
LTRE.quant.sumH <- subset(LTRE.quant.sum,!(variable%in%c('Annual_survival', 'Denning_survival'))) 

#################################################################
#### STACKED BAR (AND LINE) PLOTS - HAZARD RATES, SUMMARIES #####
#################################################################

## Re-order and re-name factor levels
LTRE.quant.sumH$variable <- factor(LTRE.quant.sumH$variable, levels = c('Annual_mortality', 'Denning_mortality', 'Pregnancy_rate', 'Fetus_number', 'Population_structure', 'Immigration'), labels = c('Annual mortality', 'Denning mortality', 'Pregnancy rate', 'Fetus number', 'Population structure', 'Immigration'))

## Plot for summarized contributions (unscaled)
p.bar <- ggplot(LTRE.quant.sumH, aes(x = Year, y = median, group = variable)) + geom_bar(aes(fill = variable, color = variable), stat = 'identity', position = 'stack') + ylab('Contribution') + scale_fill_viridis(discrete = T) + scale_color_viridis(discrete = T)  + geom_hline(yintercept = 0, linetype = 'dotted') + theme_bw() + theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
p.bar

pdf('FixedLTRE_Bars.pdf', height = 5, width = 8.3)
p.bar
dev.off()

## Plot for summarized contributions (scaled)
p.bar.sc <- ggplot(LTRE.quant.sumH, aes(x = Year, y = abs(median), group = variable)) + geom_bar(aes(fill = variable, color = variable), stat = 'identity', position = 'fill') + ylab('Absolute contribution') + scale_fill_viridis(discrete = T) + scale_color_viridis(discrete = T)  + theme_bw() + theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
p.bar.sc

pdf('FixedLTRE_ScaledBars.pdf', height = 5, width = 8.3)
p.bar.sc
dev.off()

## A variation: stacked lines for summarized contributions (scaled)
p.stacklines.sc <- ggplot(LTRE.quant.sumH, aes(x = Year, y = abs(median), group = variable)) + geom_area(aes(fill = variable, color = variable), stat = 'identity', position = 'fill') + ylab('Absolute contribution') + scale_fill_viridis(discrete = T) + scale_color_viridis(discrete = T)  + theme_bw() + theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
p.stacklines.sc

pdf('FixedLTRE_StackedLines.pdf', height = 5, width = 8.3)
p.stacklines.sc
dev.off()


## Plot for mortality components
m.params <- c('cont_mHj', 'cont_mOj', 'cont_mHa', 'cont_mOa', 'cont_m0')
m.labels <- expression(m[j]^H, m[j]^O, m[a]^H, m[a]^O, m[0])
m.colors <- c(viridis(6)[1], alpha(viridis(6)[1], 0.75), alpha(viridis(6)[1], 0.5), alpha(viridis(6)[1], 0.25), viridis(6)[2])

p.stacklines.m <- ggplot(subset(LTRE.quant.H, variable%in%m.params), aes(x = Year, y = abs(median), group = variable)) + geom_area(aes(fill = variable), color = NA, stat = 'identity', position = 'fill') + ylab('Absolute contribution') + 
scale_fill_manual(values = m.colors, labels = m.labels) + scale_color_manual(values = m.colors, labels = m.labels) + theme_bw() + theme(legend.title = element_blank(), panel.grid = element_blank()) + labs(tag = 'a)') 
p.stacklines.m


## Plot for regnancy rates
Psi.params <- c('cont_Psi2', 'cont_Psi3', 'cont_Psi4', 'cont_Psi5')
Psi.labels <- expression(Psi[1], Psi[2], Psi[3], Psi["4+"])
Psi.colors <- c(viridis(6)[3], alpha(viridis(6)[3], 0.75), alpha(viridis(6)[3], 0.5), alpha(viridis(6)[3], 0.25))

p.stacklines.Psi <- ggplot(subset(LTRE.quant.H, variable%in%Psi.params), aes(x = Year, y = abs(median), group = variable)) + geom_area(aes(fill = variable), color = NA, stat = 'identity', position = 'fill') + ylab('Absolute contribution') + 
scale_fill_manual(values = Psi.colors, labels = Psi.labels) + scale_color_manual(values = Psi.colors, labels = Psi.labels) + theme_bw() + theme(legend.title = element_blank(), panel.grid = element_blank()) + labs(tag = 'b)') 
p.stacklines.Psi


## Plot for fetus numbers
rho.params <- c('cont_rho2', 'cont_rho3', 'cont_rho4', 'cont_rho5')
rho.labels <- expression(rho[1], rho[2], rho[3], rho["4+"])
rho.colors <- c(viridis(6)[4], alpha(viridis(6)[4], 0.75), alpha(viridis(6)[4], 0.5), alpha(viridis(6)[4], 0.25))

p.stacklines.rho <- ggplot(subset(LTRE.quant.H, variable%in%rho.params), aes(x = Year, y = abs(median), group = variable)) + geom_area(aes(fill = variable), color = NA, stat = 'identity', position = 'fill') + ylab('Absolute contribution') + 
  scale_fill_manual(values = rho.colors, labels = rho.labels) + scale_color_manual(values = rho.colors, labels = rho.labels) + theme_bw() + theme(legend.title = element_blank(), panel.grid = element_blank()) + labs(tag = 'c)') 
p.stacklines.rho


## Plot for population structure and immigration
N.params <- c('cont_N_1', 'cont_N_2', 'cont_N_3', 'cont_N_4', 'cont_N_5', 'cont_Imm')
N.labels <- expression(N[0], N[1], N[2], N[3], N["4+"], I)
N.colors <- c(viridis(6)[5], alpha(viridis(6)[5], 0.75), alpha(viridis(6)[5], 0.5), alpha(viridis(6)[5], 0.25), alpha(viridis(6)[5], 0.05), viridis(6)[6])

p.stacklines.N <- ggplot(subset(LTRE.quant.H, variable%in%N.params), aes(x = Year, y = abs(median), group = variable)) + geom_area(aes(fill = variable), color = NA, stat = 'identity', position = 'fill') + ylab('Absolute contribution') + 
scale_fill_manual(values = N.colors, labels = N.labels) + scale_color_manual(values = N.colors, labels = N.labels) + theme_bw() + theme(legend.title = element_blank(), panel.grid = element_blank()) + labs(tag = 'd)') 
p.stacklines.N


## Plot separate panels together
pdf('FixedLTRE_StackedLines_Groups.pdf', width = 8.27, height = 11.69*0.9)
p.stacklines.m / p.stacklines.Psi / p.stacklines.rho / p.stacklines.N
dev.off()


#####################################################################
#### MULTIPANEL LINE AND RIBBON PLOTS - HAZARD RATES, SUMMARIES #####
#####################################################################

## Plot for summarized contributions
p.ribbons <- ggplot(LTRE.quant.sumH, aes(x = Year, y = median, group = variable)) + geom_line(aes(color = variable)) + geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = variable), alpha = 0.3) + ylab('Contribution') + scale_fill_viridis(discrete = T) + scale_color_viridis(discrete = T) + geom_hline(yintercept = 0, linetype = 'dotted', size = 0.25) + theme_bw() + theme(legend.position = 'none', panel.grid = element_blank()) + facet_wrap(~variable, ncol = 1)
p.ribbons

pdf('FixedLTRE_Ribbons.pdf', width = 8.27*0.6, height = 11.69*0.8)
p.ribbons
dev.off()


## Plot for mortality components
m.params <- c('cont_mHj', 'cont_mOj', 'cont_mHa', 'cont_mOa', 'cont_m0')

p.ribbons.m <- ggplot(subset(LTRE.quant.H, variable%in%m.params), aes(x = Year, y = median, group = variable)) + geom_line(aes(color = variable)) + geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = variable), alpha = 0.3) + ylab('Contribution') + xlab('') + scale_fill_manual(values = c(rep(viridis(6)[1], 4), viridis(6)[2])) + scale_color_manual(values = c(rep(viridis(6)[1], 4), viridis(6)[2])) + geom_hline(yintercept = 0, linetype = 'dotted') + theme_bw() + theme(legend.position = 'none', panel.grid = element_blank()) + facet_wrap(~variable, ncol = 1)
p.ribbons.m

## Plot for pregnancy rates
Psi.params <- c('cont_Psi2', 'cont_Psi3', 'cont_Psi4', 'cont_Psi5')

p.ribbons.Psi <- ggplot(subset(LTRE.quant.H, variable%in%Psi.params), aes(x = Year, y = median, group = variable)) + geom_line(color = viridis(6)[3]) + geom_ribbon(aes(ymin = lCI, ymax = uCI), fill = viridis(6)[3], alpha = 0.3) + ylab('Contribution') + xlab('') + geom_hline(yintercept = 0, linetype = 'dotted') + theme_bw() + theme(legend.position = 'none', panel.grid = element_blank()) + facet_wrap(~variable, ncol = 1)
p.ribbons.Psi


## Plot for fetus numbers
rho.params <- c('cont_rho2', 'cont_rho3', 'cont_rho4', 'cont_rho5')

p.ribbons.rho <- ggplot(subset(LTRE.quant.H, variable%in%rho.params), aes(x = Year, y = median, group = variable)) + geom_line(color = viridis(6)[4]) + geom_ribbon(aes(ymin = lCI, ymax = uCI), fill = viridis(6)[4], alpha = 0.3) + ylab('Contribution') + xlab('') + geom_hline(yintercept = 0, linetype = 'dotted') + theme_bw() + theme(legend.position = 'none', panel.grid = element_blank()) + facet_wrap(~variable, ncol = 1)
p.ribbons.rho

## Plot for population structure and immigration
N.params <- c('cont_N_1', 'cont_N_2', 'cont_N_3', 'cont_N_4', 'cont_N_5', 'cont_Imm')

p.ribbons.N <- ggplot(subset(LTRE.quant.H, variable%in%N.params), aes(x = Year, y = median, group = variable)) + geom_line(aes(color = variable)) + geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = variable), alpha = 0.3) + ylab('Contribution') + xlab('') + scale_fill_manual(values = c(rep(viridis(6)[5], 5), viridis(6)[6])) + scale_color_manual(values = c(rep(viridis(6)[5], 5), viridis(6)[6])) + geom_hline(yintercept = 0, linetype = 'dotted') + theme_bw() + theme(legend.position = 'none', panel.grid = element_blank()) + facet_wrap(~variable, ncol = 1)
p.ribbons.N

# --> Functionality confirmed. 

