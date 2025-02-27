library(reshape2)
library(ggplot2)
library(ggridges)
library(viridis)
library(patchwork)

## Load LTRE results
load('AF_LTRE_Random.RData')

str(LTRE.S.data) # LTRE with contributions of survival probabilities
str(LTRE.H.data) # LTRE with contributions of mortality hazard rates


###################################
#### SUMMARISING CONTRIBUTIONS ####
###################################

## Make summarised contributions

# With survival
LTRE.S.sum <- data.frame(
type = rep(c('Annual survival', 'Denning survival', 'Pregnancy rate', 'Fetus number', 'Population structure', 'Immigration'), each = dim(contS.all)[1]), 
contribution = c(contS.all[,'cont_Sj']+contS.all[,'cont_Sa'],
				 contS.all[,'cont_S0'],
				 contS.all[,'cont_Psi2']+contS.all[,'cont_Psi3']+contS.all[,'cont_Psi4']+contS.all[,'cont_Psi5'],
				 contS.all[,'cont_rho2']+contS.all[,'cont_rho3']+contS.all[,'cont_rho4']+contS.all[,'cont_rho5'],
				 contS.all[,'cont_N_1']+contS.all[,'cont_N_2']+contS.all[,'cont_N_3']+contS.all[,'cont_N_4']+contS.all[,'cont_N_5'],
				 contS.all[,'cont_Imm']))


# With mortality
LTRE.H.sum <- data.frame(
  type = rep(c('Annual mortality', 'Denning mortality', 'Pregnancy rate', 'Fetus number', 'Population structure', 'Immigration'), each = dim(contH.all)[1]), 
  contribution = c(contH.all[,'cont_mHj']+contH.all[,'cont_mOj']+contH.all[,'cont_mHa']+contH.all[,'cont_mOa'],
                   contH.all[,'cont_m0'],
                   contH.all[,'cont_Psi2']+contH.all[,'cont_Psi3']+contH.all[,'cont_Psi4']+contH.all[,'cont_Psi5'],
                   contH.all[,'cont_rho2']+contH.all[,'cont_rho3']+contH.all[,'cont_rho4']+contH.all[,'cont_rho5'],
                   contH.all[,'cont_N_1']+contH.all[,'cont_N_2']+contH.all[,'cont_N_3']+contH.all[,'cont_N_4']+contH.all[,'cont_N_5'],
                   contH.all[,'cont_Imm']))


##########################################
#### MULTI-PANEL PLOT - HAZARD RATES #####
##########################################

## Re-order factor levels
LTRE.H.sum$type <- factor(LTRE.H.sum$type, levels = c('Annual mortality', 'Denning mortality', 'Pregnancy rate', 'Fetus number', 'Population structure', 'Immigration'))

## Plot for summarized contributions
addline_format <- function(x,...){
    gsub('\\s','\n',x)
}

p.summary <- ggplot(LTRE.H.sum, aes(x = type, y = contribution, group = type)) + geom_violin(aes(fill = type, color = type), alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + ylim(-0.025, 0.125) + ylab('Contribution') + xlab('') + scale_fill_viridis(discrete = T) + scale_color_viridis(discrete = T) + scale_x_discrete(labels = addline_format(c('Annual mortality', 'Denning mortality', 'Pregnancy rate', 'Fetus number', 'Population structure', 'Immigration'))) + geom_hline(yintercept = 0, linetype = 'dotted') + theme_bw() + theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12)) + labs(tag = 'a)') 
p.summary

## Plot for mortality components
m.params <- c('cont_mHj', 'cont_mOj', 'cont_mHa', 'cont_mOa', 'cont_m0')

p.m <- ggplot(subset(LTRE.H.data, variable%in%m.params), aes(x = variable, y = value, group = variable)) + geom_violin(aes(fill = variable, color = variable), alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + ylab('Contribution') + xlab('') + scale_fill_manual(values = c(rep(viridis(6)[1], 4), viridis(6)[2])) + scale_color_manual(values = c(rep(viridis(6)[1], 4), viridis(6)[2]))  + scale_x_discrete(labels = expression(m[j]^H, m[j]^O, m[a]^H, m[a]^O, m[0])) + geom_hline(yintercept = 0, linetype = 'dotted') + theme_bw() + theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12)) + labs(tag = 'b)') 
p.m

## Plot for pregnancy rates
Psi.params <- c('cont_Psi2', 'cont_Psi3', 'cont_Psi4', 'cont_Psi5')

p.Psi <- ggplot(subset(LTRE.H.data, variable%in%Psi.params), aes(x = variable, y = value, group = variable)) + geom_violin(fill = viridis(6)[3], color = viridis(6)[3], alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + ylab('Contribution') + xlab('') + scale_x_discrete(labels = expression(Psi[1], Psi[2], Psi[3], Psi["4+"])) + geom_hline(yintercept = 0, linetype = 'dotted') + theme_bw() + theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12)) + labs(tag = 'c)') 
p.Psi

## Plot for fetus numbers
rho.params <- c('cont_rho2', 'cont_rho3', 'cont_rho4', 'cont_rho5')

p.rho <- ggplot(subset(LTRE.H.data, variable%in%rho.params), aes(x = variable, y = value, group = variable)) + geom_violin(fill = viridis(6)[4], color = viridis(6)[4], alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + ylab('Contribution') + xlab('') + scale_x_discrete(labels = expression(rho[1], rho[2], rho[3], rho["4+"])) + geom_hline(yintercept = 0, linetype = 'dotted') + theme_bw() + theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12)) + labs(tag = 'd)') 
p.rho


## Plot for population structure and immigration
N.params <- c('cont_N_1', 'cont_N_2', 'cont_N_3', 'cont_N_4', 'cont_N_5', 'cont_Imm')

p.N <- ggplot(subset(LTRE.H.data, variable%in%N.params), aes(x = variable, y = value, group = variable)) + geom_violin(aes(fill = variable, color = variable), alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + ylab('Contribution') + xlab('') + scale_fill_manual(values = c(rep(viridis(6)[5], 5), viridis(6)[6])) + scale_color_manual(values = c(rep(viridis(6)[5], 5), viridis(6)[6]))  + scale_x_discrete(labels = expression(N[0], N[1], N[2], N[3], N["4+"], I)) + geom_hline(yintercept = 0, linetype = 'dotted') + theme_bw() + theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12)) + labs(tag = 'e)') 
p.N

## Join plots together (using patchwork package)
pdf('RandomLTRE_Violins.pdf', width = 8.27, height = 11.69)
p.summary / (p.m + ylim(-0.02, 0.08) | p.Psi + ylim(-0.02, 0.08)) / (p.rho | p.N + ylim(-0.02, 0.08))
dev.off()

# --> Functionality confirmed. 