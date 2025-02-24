library(reshape)
library(ggplot2)
library(ggridges)
library(viridis)
library(patchwork)

## Load LTRE results
load('200503_AF_LTRE_Random.RData')

str(elas.data) 

##################################
#### SUMMARISING ELASTICITIES ####
##################################

## Make summarised contributions

# With survival
elas.S.sum <- data.frame(
  type = rep(c('Annual survival', 'Denning survival', 'Pregnancy rate', 'Fetus number', 'Population structure', 'Immigration'), each = dim(elas.all)[1]), 
  elasticity = c(elas.all[,'elas_Sj']+elas.all[,'elas_Sa'],
                   elas.all[,'elas_S0'],
                   elas.all[,'elas_Psi2']+elas.all[,'elas_Psi3']+elas.all[,'elas_Psi4']+elas.all[,'elas_Psi5'],
                   elas.all[,'elas_rho2']+elas.all[,'elas_rho3']+elas.all[,'elas_rho4']+elas.all[,'elas_rho5'],
                   elas.all[,'elas_N_1']+elas.all[,'elas_N_2']+elas.all[,'elas_N_3']+elas.all[,'elas_N_4']+elas.all[,'elas_N_5'],
                   elas.all[,'elas_Imm']))


# With mortality
elas.H.sum <- data.frame(
  type = rep(c('Annual mortality', 'Denning mortality', 'Pregnancy rate', 'Fetus number', 'Population structure', 'Immigration'), each = dim(elas.all)[1]), 
  elasticity = c(elas.all[,'elas_mHj']+elas.all[,'elas_mOj']+elas.all[,'elas_mHa']+elas.all[,'elas_mOa'],
                   elas.all[,'elas_m0'],
                   elas.all[,'elas_Psi2']+elas.all[,'elas_Psi3']+elas.all[,'elas_Psi4']+elas.all[,'elas_Psi5'],
                   elas.all[,'elas_rho2']+elas.all[,'elas_rho3']+elas.all[,'elas_rho4']+elas.all[,'elas_rho5'],
                   elas.all[,'elas_N_1']+elas.all[,'elas_N_2']+elas.all[,'elas_N_3']+elas.all[,'elas_N_4']+elas.all[,'elas_N_5'],
                   elas.all[,'elas_Imm']))


##########################################
#### MULTI-PANEL PLOT - HAZARD RATES #####
##########################################

## Re-order factor levels
elas.H.sum$type <- factor(elas.H.sum$type, levels = c('Annual mortality', 'Denning mortality', 'Pregnancy rate', 'Fetus number', 'Population structure', 'Immigration'))

## Plot for summarized contributions
addline_format <- function(x,...){
    gsub('\\s','\n',x)
}

p.summary <- ggplot(elas.H.sum, aes(x = type, y = abs(elasticity), group = type)) + geom_violin(aes(fill = type, color = type), alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + ylab('Absolute elasticity') + xlab('') + scale_fill_viridis(discrete = T) + scale_color_viridis(discrete = T) + scale_x_discrete(labels = addline_format(c('Annual mortality', 'Denning mortality', 'Pregnancy rate', 'Fetus number', 'Population structure', 'Immigration'))) + theme_bw() + theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12)) + labs(tag = 'a)') 
p.summary

## Plot for mortality components
m.params <- c('elas_mHj', 'elas_mOj', 'elas_mHa', 'elas_mOa', 'elas_m0')

p.m <- ggplot(subset(elas.data, variable%in%m.params), aes(x = variable, y = value, group = variable)) + geom_violin(aes(fill = variable, color = variable), alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + ylab('Elasticity') + xlab('') + scale_fill_manual(values = c(rep(viridis(6)[1], 4), viridis(6)[2])) + scale_color_manual(values = c(rep(viridis(6)[1], 4), viridis(6)[2]))  + scale_x_discrete(labels = expression(m[j]^H, m[j]^O, m[a]^H, m[a]^O, m[0])) + theme_bw() + theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12)) + labs(tag = 'b)') 
p.m

## Plot for pregnancy rates
Psi.params <- c('elas_Psi2', 'elas_Psi3', 'elas_Psi4', 'elas_Psi5')

p.Psi <- ggplot(subset(elas.data, variable%in%Psi.params), aes(x = variable, y = value, group = variable)) + geom_violin(fill = viridis(6)[3], color = viridis(6)[3], alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + ylab('Elasticity') + xlab('') + scale_x_discrete(labels = expression(Psi[1], Psi[2], Psi[3], Psi["4+"])) + theme_bw() + theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12)) + labs(tag = 'c)') 
p.Psi

## Plot for fetus numbers
rho.params <- c('elas_rho2', 'elas_rho3', 'elas_rho4', 'elas_rho5')

p.rho <- ggplot(subset(elas.data, variable%in%rho.params), aes(x = variable, y = value, group = variable)) + geom_violin(fill = viridis(6)[4], color = viridis(6)[4], alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + ylab('Elasticity') + xlab('') + scale_x_discrete(labels = expression(rho[1], rho[2], rho[3], rho["4+"])) + theme_bw() + theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12)) + labs(tag = 'd)') 
p.rho


## Plot for population structure and immigration
N.params <- c('elas_N_1', 'elas_N_2', 'elas_N_3', 'elas_N_4', 'elas_N_5', 'elas_Imm')

p.N <- ggplot(subset(elas.data, variable%in%N.params), aes(x = variable, y = value, group = variable)) + geom_violin(aes(fill = variable, color = variable), alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + ylab('Elasticity') + xlab('') + scale_fill_manual(values = c(rep(viridis(6)[5], 5), viridis(6)[6])) + scale_color_manual(values = c(rep(viridis(6)[5], 5), viridis(6)[6]))  + scale_x_discrete(labels = expression(N[0], N[1], N[2], N[3], N["4+"], I)) + theme_bw() + theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12)) + labs(tag = 'e)') 
p.N

## Join plots together (using patchwork package)
pdf('Elasticity_Violins.pdf', width = 8.27, height = 11.69)
p.summary / (p.m | p.Psi) / (p.rho | p.N)
dev.off()

