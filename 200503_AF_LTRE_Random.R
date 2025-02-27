library(coda)
library(matrixStats)
library(reshape2)
library(ggplot2)

## Load MCMC samples
load("Data/ArcticFox_IPM_PosteriorSamples.RData")
str(MCMCsamples)

## Transform samples into matrix (all chains combined)
out.mat <- as.matrix(MCMCsamples)
str(out.mat)

## Take a subset of samples for testing
#out.mat <- out.mat[1:10,]

#sample.index <- seq(1, dim(out.mat)[1], by = 32)
#out.mat <- out.mat[sample.index,]

## Set sample and year number
nosamples <- dim(out.mat)[1]
noyears <- 23
  
  
###############
#### SETUP ####
###############

## Prepare matrices to rearrange samples - Vital rates & population sizes

# Time-varying vital rates
Sj <- mHj <- mOj <- matrix(NA, nrow = nosamples, ncol = noyears)
Sa <- mHa <- mOa <- matrix(NA, nrow = nosamples, ncol = noyears)
S0 <- m0 <- matrix(NA, nrow = nosamples, ncol = noyears)
Psi2 <- Psi3 <- Psi4 <- Psi5 <- matrix(NA, nrow = nosamples, ncol = noyears)
rho2 <- rho3 <- rho4 <- rho5 <- matrix(NA, nrow = nosamples, ncol = noyears)

# Time-varying population sizes and growth rates
N_1 <- matrix(NA, nrow = nosamples, ncol = noyears)
N_2 <- matrix(NA, nrow = nosamples, ncol = noyears)
N_3 <- matrix(NA, nrow = nosamples, ncol = noyears)
N_4 <- matrix(NA, nrow = nosamples, ncol = noyears)
N_5 <- matrix(NA, nrow = nosamples, ncol = noyears)
N_tot <- matrix(NA, nrow = nosamples, ncol = noyears)
lambda <- matrix(NA, nrow = nosamples, ncol = noyears)

# Time-varying immigrant numbers
Imm <- matrix(NA, nrow = nosamples, ncol = noyears)


## Fill samples into vectors and matrices
for(i in 1:nosamples){
  
  for(t in 1:noyears){
  	
  	# Time-varying vital rates   
    Psi2[i,t] <- out.mat[i, paste("Psi[2, ", t, "]", sep = "")]
    Psi3[i,t] <- out.mat[i, paste("Psi[3, ", t, "]", sep = "")]
    Psi4[i,t] <- out.mat[i, paste("Psi[4, ", t, "]", sep = "")]
    Psi5[i,t] <- out.mat[i, paste("Psi[5, ", t, "]", sep = "")]
      
    rho2[i,t] <- out.mat[i, paste("rho[2, ", t, "]", sep = "")]
    rho3[i,t] <- out.mat[i, paste("rho[3, ", t, "]", sep = "")]
    rho4[i,t] <- out.mat[i, paste("rho[4, ", t, "]", sep = "")]
    rho5[i,t] <- out.mat[i, paste("rho[5, ", t, "]", sep = "")]	
  	
  	m0[i,t] <- out.mat[i, paste("m0t[", t, "]", sep = "")]
    S0[i,t] <- exp(-m0[i,t])
    
    # Time-varying population sizes
    N_1[i,t] <- out.mat[i, paste("N[1, ", t, "]", sep="")]
    N_2[i,t] <- out.mat[i, paste("N[2, ", t, "]", sep="")]
    N_3[i,t] <- out.mat[i, paste("N[3, ", t, "]", sep="")]
    N_4[i,t] <- out.mat[i, paste("N[4, ", t, "]", sep="")]
    N_5[i,t] <- out.mat[i, paste("N[5, ", t, "]", sep="")]
    N_tot[i,t] <- out.mat[i, paste("N.tot[", t, "]", sep="")]
    
    # Time-varying immigrant numbers
	Imm[i,t] <- out.mat[i, paste("Imm[", t, "]", sep="")]
  }
  
  for(t in 1:(noyears-1)){
   
    # Time-varying vital rates
  	mHj[i,t] <- out.mat[i, paste("mH[2, ", t, "]", sep = "")]
  	mOj[i,t] <- out.mat[i, paste("mO[2, ", t, "]", sep = "")] 
  	Sj[i,t] <- exp(-(mHj[i,t] + mOj[i,t]))
  	
  	mHa[i,t] <- out.mat[i, paste("mH[1, ", t, "]", sep = "")]
  	mOa[i,t] <- out.mat[i, paste("mO[1, ", t, "]", sep = "")] 
  	Sa[i,t] <- exp(-(mHa[i,t] + mOa[i,t]))

    
    # Population growth rate
    lambda[i,t] <- N_tot[i,t+1]/N_tot[i,t]
  }
}


## Make time-average population sizes
N_1_mean <- rowMeans(N_1[,1:22])
N_2_mean <- rowMeans(N_2[,1:22])
N_3_mean <- rowMeans(N_3[,1:22])
N_4_mean <- rowMeans(N_4[,1:22])
N_5_mean <- rowMeans(N_5[,1:22])
N_tot_mean <- rowMeans(N_tot[,1:22])

n_1_mean <- rowMeans(N_1)/rowMeans(N_tot)
n_2_mean <- rowMeans(N_2)/rowMeans(N_tot)
n_3_mean <- rowMeans(N_3)/rowMeans(N_tot)
n_4_mean <- rowMeans(N_4)/rowMeans(N_tot)
n_5_mean <- rowMeans(N_5)/rowMeans(N_tot)

## Make average immigrant numbers
Imm_mean <- rowMeans(Imm[,2:23])

## Make time-average vital rates
Sj_mean <- rowMeans(Sj, na.rm = T)
mHj_mean <- rowMeans(mHj, na.rm = T)
mOj_mean <- rowMeans(mOj, na.rm = T)

Sa_mean <- rowMeans(Sa, na.rm = T)
mHa_mean <- rowMeans(mHa, na.rm = T)
mOa_mean <- rowMeans(mOa, na.rm = T)

S0_mean <- rowMeans(S0[,2:23], na.rm = T)
m0_mean <- rowMeans(m0[,2:23], na.rm = T)

Psi2_mean <- rowMeans(Psi2[,2:23], na.rm = T)
Psi3_mean <- rowMeans(Psi3[,2:23], na.rm = T)
Psi4_mean <- rowMeans(Psi4[,2:23], na.rm = T)
Psi5_mean <- rowMeans(Psi5[,2:23], na.rm = T)

rho2_mean <- rowMeans(rho2[,2:23], na.rm = T)
rho3_mean <- rowMeans(rho3[,2:23], na.rm = T)
rho4_mean <- rowMeans(rho4[,2:23], na.rm = T)
rho5_mean <- rowMeans(rho5[,2:23], na.rm = T)

# NOTE: I am averaging across indeces 1:22 for juvenile & adult survival/mortality (and population sizes), and 2:23 for all reproductive parameters (and immigration) because in the matrix, the former appear with index t and the latter with index t+1

## Make time-average population growth rate
lambda_mean <- rowMeans(lambda, na.rm = T)


################################################
#### CALCULATION OF TRANSIENT SENSITIVITIES ####
################################################

## Calculate transient sensitivities for vital rates and population size/structure (evaluated at the temporal mean)

sens_Sj <- rep(NA, nosamples)
sens_mHj <- rep(NA, nosamples)
sens_mOj <- rep(NA, nosamples)

sens_Sa <- rep(NA, nosamples)
sens_mHa <- rep(NA, nosamples)
sens_mOa <- rep(NA, nosamples)

sens_S0 <- rep(NA, nosamples)
sens_m0 <- rep(NA, nosamples)

sens_Psi2 <- rep(NA, nosamples)
sens_Psi3 <- rep(NA, nosamples)
sens_Psi4 <- rep(NA, nosamples)
sens_Psi5 <- rep(NA, nosamples)

sens_rho2 <- rep(NA, nosamples)
sens_rho3 <- rep(NA, nosamples)
sens_rho4 <- rep(NA, nosamples)
sens_rho5 <- rep(NA, nosamples)

sens_N_1 <- rep(NA, nosamples)
sens_N_2 <- rep(NA, nosamples)
sens_N_3 <- rep(NA, nosamples)
sens_N_4 <- rep(NA, nosamples)
sens_N_5 <- rep(NA, nosamples)

sens_Imm <- rep(NA, nosamples)


for(i in 1:nosamples){
  
	sens_Sj[i] <- (0.5*Psi2_mean[i]*rho2_mean[i]*S0_mean[i] + 1)*n_1_mean[i]
	
	sens_mHj[i] <- -exp(-(mHj_mean[i] + mOj_mean[i]))*sens_Sj[i]
	sens_mOj[i] <- sens_mHj[i]


	sens_Sa[i] <- (0.5*Psi3_mean[i]*rho3_mean[i]*S0_mean[i] + 1)*n_2_mean[i] + (0.5*Psi4_mean[i]*rho4_mean[i]*S0_mean[i] + 1)*n_3_mean[i] + (0.5*Psi5_mean[i]*rho5_mean[i]*S0_mean[i] + 1)*(n_4_mean[i] + n_5_mean[i])
	
	sens_mHa[i] <- -exp(-(mHa_mean[i] + mOa_mean[i]))*sens_Sj[i]
	sens_mOa[i] <- sens_mHa[i] 


	sens_S0[i] <- 0.5*(Psi2_mean[i]*rho2_mean[i]*n_1_mean[i] + Psi3_mean[i]*rho3_mean[i]*n_2_mean[i] + Psi4_mean[i]*rho4_mean[i]*n_3_mean[i] + Psi5_mean[i]*rho5_mean[i]*(n_4_mean[i] + n_5_mean[i]))
	
	sens_m0[i] <- -exp(-m0_mean[i])*sens_S0[i]


	sens_Psi2[i] <- 0.5*S0_mean[i]*rho2_mean[i]*n_1_mean[i]
	
	sens_Psi3[i] <- 0.5*S0_mean[i]*rho3_mean[i]*n_2_mean[i]
	
	sens_Psi4[i] <- 0.5*S0_mean[i]*rho4_mean[i]*n_3_mean[i]
	
	sens_Psi5[i] <- 0.5*S0_mean[i]*rho5_mean[i]*(n_4_mean[i] + n_5_mean[i])


	sens_rho2[i] <- 0.5*S0_mean[i]*Psi2_mean[i]*n_1_mean[i]
	
	sens_rho3[i] <- 0.5*S0_mean[i]*Psi3_mean[i]*n_2_mean[i]
	
	sens_rho4[i] <- 0.5*S0_mean[i]*Psi4_mean[i]*n_3_mean[i]
	
	sens_rho5[i] <- 0.5*S0_mean[i]*Psi5_mean[i]*(n_4_mean[i] + n_5_mean[i])


	sens_N_1[i] <- (Sj_mean[i]*(0.5*Psi2_mean[i]*rho2_mean[i]*S0_mean[i] + 1) - lambda_mean[i]) / N_tot_mean[i]
	
	sens_N_2[i] <- (Sa_mean[i]*(0.5*Psi3_mean[i]*rho3_mean[i]*S0_mean[i] + 1) - lambda_mean[i]) / N_tot_mean[i]
	
	sens_N_3[i] <- (Sa_mean[i]*(0.5*Psi4_mean[i]*rho4_mean[i]*S0_mean[i] + 1) - lambda_mean[i]) / N_tot_mean[i]
	
	sens_N_4[i] <- (Sa_mean[i]*(0.5*Psi5_mean[i]*rho5_mean[i]*S0_mean[i] + 1) - lambda_mean[i]) / N_tot_mean[i]
	
	sens_N_5[i] <- sens_N_4[i]
	
	
	sens_Imm[i] <- 1 / N_tot_mean[i]
 
}


## Get posterior summaries for transient sensitivities
quantile(sens_Sj, probs = c(0.05, 0.5, 0.95))
quantile(sens_mHj, probs = c(0.05, 0.5, 0.95))
quantile(sens_mOj, probs = c(0.05, 0.5, 0.95))

quantile(sens_Sa, probs = c(0.05, 0.5, 0.95))
quantile(sens_mHa, probs = c(0.05, 0.5, 0.95))
quantile(sens_mOa, probs = c(0.05, 0.5, 0.95))

quantile(sens_S0, probs = c(0.05, 0.5, 0.95))
quantile(sens_m0, probs = c(0.05, 0.5, 0.95))

quantile(sens_Psi2, probs = c(0.05, 0.5, 0.95))
quantile(sens_Psi3, probs = c(0.05, 0.5, 0.95))
quantile(sens_Psi4, probs = c(0.05, 0.5, 0.95))
quantile(sens_Psi5, probs = c(0.05, 0.5, 0.95))
quantile(sens_Psi2+sens_Psi3+sens_Psi4+sens_Psi5, probs = c(0.05, 0.5, 0.95))

quantile(sens_rho2, probs = c(0.05, 0.5, 0.95))
quantile(sens_rho3, probs = c(0.05, 0.5, 0.95))
quantile(sens_rho4, probs = c(0.05, 0.5, 0.95))
quantile(sens_rho5, probs = c(0.05, 0.5, 0.95))
quantile(sens_rho2+sens_rho3+sens_rho4+sens_rho5, probs = c(0.05, 0.5, 0.95))

quantile(sens_N_1, probs = c(0.05, 0.5, 0.95))
quantile(sens_N_2, probs = c(0.05, 0.5, 0.95))
quantile(sens_N_3, probs = c(0.05, 0.5, 0.95))
quantile(sens_N_4, probs = c(0.05, 0.5, 0.95))
quantile(sens_N_5, probs = c(0.05, 0.5, 0.95))

quantile(sens_Imm, probs = c(0.05, 0.5, 0.95))

## Plot sensitivities
sens.all <- cbind(sens_Sj, sens_mHj, sens_mOj, 
				  sens_Sa, sens_mHa, sens_mOa, 
				  sens_S0, sens_m0, 
				  sens_Psi2, sens_Psi3, sens_Psi4, sens_Psi5,
				  sens_rho2, sens_rho3, sens_rho4, sens_rho4,
				  sens_N_1, sens_N_2, sens_N_3, sens_N_4, sens_N_5,
				  sens_Imm)

sens.data <- melt(as.data.frame(sens.all))

ggplot(subset(sens.data, value != 0), aes(x = variable, y = value, group = variable)) + geom_boxplot()
ggplot(subset(sens.data, value != 0), aes(x = variable, y = abs(value), group = variable)) + geom_boxplot()


###############################################
#### CALCULATION OF TRANSIENT ELASTICITIES ####
###############################################

## Calculate transient elasiticities for vital rates and population size/structure (evaluated at the temporal mean)

elas_Sj <- sens_Sj*(Sj_mean/lambda_mean)
elas_mHj <- sens_mHj*(mHj_mean/lambda_mean)
elas_mOj <- sens_mOj*(mOj_mean/lambda_mean)

elas_Sa <- sens_Sa*(Sa_mean/lambda_mean)
elas_mHa <- sens_mHa*(mHa_mean/lambda_mean)
elas_mOa <- sens_mOa*(mOa_mean/lambda_mean)

elas_S0 <- sens_S0*(S0_mean/lambda_mean)
elas_m0 <- sens_m0*(m0_mean/lambda_mean)

elas_Psi2 <- sens_Psi2*(Psi2_mean/lambda_mean)
elas_Psi3 <- sens_Psi3*(Psi3_mean/lambda_mean)
elas_Psi4 <- sens_Psi4*(Psi4_mean/lambda_mean)
elas_Psi5 <- sens_Psi5*(Psi5_mean/lambda_mean)

elas_rho2 <- sens_rho2*(rho2_mean/lambda_mean)
elas_rho3 <- sens_rho3*(rho3_mean/lambda_mean)
elas_rho4 <- sens_rho4*(rho4_mean/lambda_mean)
elas_rho5 <- sens_rho5*(rho5_mean/lambda_mean)

elas_N_1 <- sens_N_1*(N_1_mean/lambda_mean)
elas_N_2 <- sens_N_2*(N_2_mean/lambda_mean)
elas_N_3 <- sens_N_3*(N_3_mean/lambda_mean)
elas_N_4 <- sens_N_4*(N_4_mean/lambda_mean)
elas_N_5 <- sens_N_5*(N_5_mean/lambda_mean)

elas_Imm <- sens_Imm*(Imm_mean/lambda_mean)


## Get posterior summaries for transient elasticities
quantile(elas_Sj, probs = c(0.05, 0.5, 0.95))
quantile(elas_mHj, probs = c(0.05, 0.5, 0.95))
quantile(elas_mOj, probs = c(0.05, 0.5, 0.95))

quantile(elas_Sa, probs = c(0.05, 0.5, 0.95))
quantile(elas_mHa, probs = c(0.05, 0.5, 0.95))
quantile(elas_mOa, probs = c(0.05, 0.5, 0.95))

quantile(elas_S0, probs = c(0.05, 0.5, 0.95))
quantile(elas_m0, probs = c(0.05, 0.5, 0.95))

quantile(elas_Psi2, probs = c(0.05, 0.5, 0.95))
quantile(elas_Psi3, probs = c(0.05, 0.5, 0.95))
quantile(elas_Psi4, probs = c(0.05, 0.5, 0.95))
quantile(elas_Psi5, probs = c(0.05, 0.5, 0.95))
quantile(sum(elas_Psi2, elas_Psi3, elas_Psi4, elas_Psi5), probs = c(0.05, 0.5, 0.95))

quantile(elas_rho2, probs = c(0.05, 0.5, 0.95))
quantile(elas_rho3, probs = c(0.05, 0.5, 0.95))
quantile(elas_rho4, probs = c(0.05, 0.5, 0.95))
quantile(elas_rho5, probs = c(0.05, 0.5, 0.95))
quantile(sum(elas_rho2, elas_rho3, elas_rho4, elas_rho5), probs = c(0.05, 0.5, 0.95))

quantile(elas_N_1, probs = c(0.05, 0.5, 0.95))
quantile(elas_N_2, probs = c(0.05, 0.5, 0.95))
quantile(elas_N_3, probs = c(0.05, 0.5, 0.95))
quantile(elas_N_4, probs = c(0.05, 0.5, 0.95))
quantile(elas_N_5, probs = c(0.05, 0.5, 0.95))

quantile(elas_Imm, probs = c(0.05, 0.5, 0.95))

## Plot elasticities
elas.all <- cbind(elas_Sj, elas_mHj, elas_mOj, 
                  elas_Sa, elas_mHa, elas_mOa, 
                  elas_S0, elas_m0, 
                  elas_Psi2, elas_Psi3, elas_Psi4, elas_Psi5,
                  elas_rho2, elas_rho3, elas_rho4, elas_rho5,
                  elas_N_1, elas_N_2, elas_N_3, elas_N_4, elas_N_5,
                  elas_Imm)

elas.data <- melt(as.data.frame(elas.all))

ggplot(subset(elas.data, value != 0), aes(x = variable, y = value, group = variable)) + geom_boxplot()
ggplot(subset(elas.data, value != 0), aes(x = variable, y = abs(value), group = variable)) + geom_boxplot()

###############################################################
#### CALCULATE LTRE CONTRIBUTIONS - SURVIVAL PROBABILITIES ####
###############################################################

## Prepare vectors to store results
cont_Sj <- rep(NA, nosamples)

cont_Sa <- rep(NA, nosamples)

cont_S0 <- rep(NA, nosamples)

cont_Psi2 <- rep(NA, nosamples)
cont_Psi3 <- rep(NA, nosamples)
cont_Psi4 <- rep(NA, nosamples)
cont_Psi5 <- rep(NA, nosamples)

cont_rho2 <- rep(NA, nosamples)
cont_rho3 <- rep(NA, nosamples)
cont_rho4 <- rep(NA, nosamples)
cont_rho5 <- rep(NA, nosamples)

cont_N_1 <- rep(NA, nosamples)
cont_N_2 <- rep(NA, nosamples)
cont_N_3 <- rep(NA, nosamples)
cont_N_4 <- rep(NA, nosamples)
cont_N_5 <- rep(NA, nosamples)

cont_Imm <- rep(NA, nosamples)

est_var <- matrix(NA, nrow = 17, ncol = nosamples)
est_covar <- matrix(NA, nrow = 17, ncol = nosamples)


## Calculate LTRE contributions (random design LTRE)

for(i in 1:nosamples){
  
  ## Make matrix of vital rates and population sizes
  dp_stoch <- cbind(Sj[i,1:(noyears-1)],
					Sa[i,1:(noyears-1)],
  					S0[i,2:noyears],
  					Psi2[i,2:noyears],
  					Psi3[i,2:noyears],
  					Psi4[i,2:noyears],
  					Psi5[i,2:noyears],
  					rho2[i,2:noyears],
  					rho3[i,2:noyears],
  					rho4[i,2:noyears],
  					Psi5[i,2:noyears],
                    N_1[i,1:(noyears-1)],
                    N_2[i,1:(noyears-1)],
                    N_3[i,1:(noyears-1)],
                    N_4[i,1:(noyears-1)],
                    N_5[i,1:(noyears-1)],
                    Imm[i,2:noyears]
                    )
  
  ## Derive process variances and covariances
  dp_varcov <- var(dp_stoch)
  
  ## Save total estimated (co)variance per parameter
  est_var[,i] <- diag(dp_varcov)
  est_covar[,i] <- rowSums(dp_varcov)
  
  ## Make a vector of sensitivities
  sensvec <- c(	sens_Sj[i],
				sens_Sa[i],
				sens_S0[i],
				sens_Psi2[i],
				sens_Psi3[i],
				sens_Psi4[i],
				sens_Psi5[i],
				sens_rho2[i],
				sens_rho3[i],
				sens_rho4[i],
				sens_rho5[i],
				sens_N_1[i],
				sens_N_2[i],
				sens_N_3[i],
				sens_N_4[i],
				sens_N_5[i],
				sens_Imm[i]
			)
  
  ## Calculate demographic contributions
  # NOTE: Here we multiply sensitivities and (co)variances
  
  cont.mat <- matrix(NA, nrow = length(sensvec), ncol = length(sensvec))
  for(k in 1:length(sensvec)){
    for(m in 1:length(sensvec)){
      cont.mat[k,m] <- dp_varcov[k,m]*sensvec[k]*sensvec[m]
    }
  }
  
  ## Summarise contributions (sum of variances and covariances)
  cont <- rowSums(cont.mat)
  
  cont_Sj[i] <- cont[1]
  cont_Sa[i] <- cont[2]
  cont_S0[i] <- cont[3]
  cont_Psi2[i] <- cont[4]
  cont_Psi3[i] <- cont[5]
  cont_Psi4[i] <- cont[6]
  cont_Psi5[i] <- cont[7]
  cont_rho2[i] <- cont[8]
  cont_rho3[i] <- cont[9]
  cont_rho4[i] <- cont[10]
  cont_rho5[i] <- cont[11]
  cont_N_1[i] <- cont[12]
  cont_N_2[i] <- cont[13]
  cont_N_3[i] <- cont[14]
  cont_N_4[i] <- cont[15]
  cont_N_5[i] <- cont[16]
  cont_Imm[i] <- cont[17]
}

## Check sum of contributions agains variance in lambda
contS.all <- cbind(cont_Sj, cont_Sa, cont_S0,
				  cont_Psi2, cont_Psi3, cont_Psi4, cont_Psi5,
				  cont_rho2, cont_rho3, cont_rho4, cont_rho5,
				  cont_N_1, cont_N_2, cont_N_3, cont_N_4, cont_N_5,
				  cont_Imm)

contS.data <- as.data.frame(contS.all)

total.contS <- rowSums(contS.all)
quantile(total.contS, probs = c(0.05, 0.5, 0.995))

tempvar_lambda <- rowVars(lambda[,1:(noyears-1)])
quantile(tempvar_lambda, probs = c(0.05, 0.5, 0.995))
# --> The variation explained by LTRE contributions and the variation in lambda are quite close
# --> Strangely, the former is larger than the latter, but that may have to do with the fact that the LTRE analysis ignores demographic stochasticity

## Compare contributions
colMeans(contS.all)

## Plot posteriors for LTRE contributions
LTRE.S.data <- melt(contS.data)

ggplot(subset(LTRE.S.data, value != 0), aes(x = variable, y = value, group = variable)) + geom_boxplot()
ggplot(subset(LTRE.S.data, value != 0), aes(x = variable, y = abs(value), group = variable)) + geom_boxplot()
# --> This indicates the largest contributions come from juvenile and adult survival, as well as immigration


###############################################################
#### CALCULATE LTRE CONTRIBUTIONS - MORTALITY HAZARD RATES ####
###############################################################

## Prepare vectors to store results
cont_mHj <- rep(NA, nosamples)
cont_mOj <- rep(NA, nosamples)

cont_mHa <- rep(NA, nosamples)
cont_mOa <- rep(NA, nosamples)

cont_m0 <- rep(NA, nosamples)

cont_Psi2 <- rep(NA, nosamples)
cont_Psi3 <- rep(NA, nosamples)
cont_Psi4 <- rep(NA, nosamples)
cont_Psi5 <- rep(NA, nosamples)

cont_rho2 <- rep(NA, nosamples)
cont_rho3 <- rep(NA, nosamples)
cont_rho4 <- rep(NA, nosamples)
cont_rho5 <- rep(NA, nosamples)

cont_N_1 <- rep(NA, nosamples)
cont_N_2 <- rep(NA, nosamples)
cont_N_3 <- rep(NA, nosamples)
cont_N_4 <- rep(NA, nosamples)
cont_N_5 <- rep(NA, nosamples)

cont_Imm <- rep(NA, nosamples)

est_var <- matrix(NA, nrow = 19, ncol = nosamples)
est_covar <- matrix(NA, nrow = 19, ncol = nosamples)


## Calculate LTRE contributions (random design LTRE)

for(i in 1:nosamples){
  
  ## Make matrix of vital rates and population sizes
  dp_stoch <- cbind(mHj[i,1:(noyears-1)],
  					mOj[i,1:(noyears-1)],
  					mHa[i,1:(noyears-1)],
  					mOa[i,1:(noyears-1)],
  					m0[i,2:noyears],
  					Psi2[i,2:noyears],
  					Psi3[i,2:noyears],
  					Psi4[i,2:noyears],
  					Psi5[i,2:noyears],
  					rho2[i,2:noyears],
  					rho3[i,2:noyears],
  					rho4[i,2:noyears],
  					Psi5[i,2:noyears],
                    N_1[i,1:(noyears-1)],
                    N_2[i,1:(noyears-1)],
                    N_3[i,1:(noyears-1)],
                    N_4[i,1:(noyears-1)],
                    N_5[i,1:(noyears-1)],
                    Imm[i,2:noyears]
                    )
  
  ## Derive process variances and covariances
  dp_varcov <- var(dp_stoch)
  
  ## Save total estimated (co)variance per parameter
  est_var[,i] <- diag(dp_varcov)
  est_covar[,i] <- rowSums(dp_varcov)
  
  ## Make a vector of sensitivities
  sensvec <- c(	sens_mHj[i],
				sens_mOj[i],
				sens_mHa[i],
				sens_mOa[i],
				sens_m0[i],
				sens_Psi2[i],
				sens_Psi3[i],
				sens_Psi4[i],
				sens_Psi5[i],
				sens_rho2[i],
				sens_rho3[i],
				sens_rho4[i],
				sens_rho5[i],
				sens_N_1[i],
				sens_N_2[i],
				sens_N_3[i],
				sens_N_4[i],
				sens_N_5[i],
				sens_Imm[i]
			)
  
  ## Calculate demographic contributions
  # NOTE: Here we multiply sensitivities and (co)variances
  
  cont.mat <- matrix(NA, nrow = length(sensvec), ncol = length(sensvec))
  for(k in 1:length(sensvec)){
    for(m in 1:length(sensvec)){
      cont.mat[k,m] <- dp_varcov[k,m]*sensvec[k]*sensvec[m]
    }
  }
  
  ## Summarise contributions (sum of variances and covariances)
  cont <- rowSums(cont.mat)
  
  cont_mHj[i] <- cont[1]
  cont_mOj[i] <- cont[2]
  cont_mHa[i] <- cont[3]
  cont_mOa[i] <- cont[4]
  cont_m0[i] <- cont[5]
  cont_Psi2[i] <- cont[6]
  cont_Psi3[i] <- cont[7]
  cont_Psi4[i] <- cont[8]
  cont_Psi5[i] <- cont[9]
  cont_rho2[i] <- cont[10]
  cont_rho3[i] <- cont[11]
  cont_rho4[i] <- cont[12]
  cont_rho5[i] <- cont[13]
  cont_N_1[i] <- cont[14]
  cont_N_2[i] <- cont[15]
  cont_N_3[i] <- cont[16]
  cont_N_4[i] <- cont[17]
  cont_N_5[i] <- cont[18]
  cont_Imm[i] <- cont[19]
}

## Check sum of contributions agains variance in lambda
contH.all <- cbind(cont_mHj, cont_mOj, cont_mHa, cont_mOa, cont_m0, 
				  cont_Psi2, cont_Psi3, cont_Psi4, cont_Psi5,
				  cont_rho2, cont_rho3, cont_rho4, cont_rho5,
				  cont_N_1, cont_N_2, cont_N_3, cont_N_4, cont_N_5,
				  cont_Imm)

contH.data <- as.data.frame(contH.all)

total.contH <- rowSums(contH.all)
quantile(total.contH, probs = c(0.025, 0.5, 0.975))

tempvar_lambda <- rowVars(lambda[,1:(noyears-1)])
quantile(tempvar_lambda, probs = c(0.025, 0.5, 0.975))
# --> With this one, there is less "overestimation" of variation explained by LTRE contributions (closer to observed variation in lambda)

## Compare contributions
colMeans(contH.all)

## Plot posteriors for LTRE contributions
LTRE.H.data <- melt(contH.data)

ggplot(subset(LTRE.H.data, value != 0), aes(x = variable, y = value, group = variable)) + geom_boxplot()
ggplot(subset(LTRE.H.data, value != 0), aes(x = variable, y = abs(value), group = variable)) + geom_boxplot()
# --> This now shows that changes in natural survival contributed more than changes in harvest
# --> It also paints a more even picture where pregnancy rate is similarly important to cause-specific mortality 

save.image("AF_LTRE_Random.RData")

# --> Functionality confirmed. 