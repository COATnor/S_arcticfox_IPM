library(coda)
library(matrixStats)
library(reshape)

## Load MCMC samples
load('ArcticFox_IPM_PosteriorSamples.RData')
str(MCMCsamples)
# = posterior samples in coda format

## Transform samples into matrix (all chains combined)
out.mat <- as.matrix(MCMCsamples)
str(out.mat)

## Set sample and year number
nosamples <- dim(out.mat)[1]
noyears <- 23
  
  
########################################
#### FUNCTION FOR FIXED DESIGN LTRE ####
########################################


fixedD.LTRE = function(t1, t2){
	
	#-------------------------------------------------
	# STEP 1) Organise samples for relevant quantities
	#-------------------------------------------------
	
	## Extract vital rates for relevant time-steps (all samples)
	mHj <- cbind(out.mat[, paste("mH[2, ", t1, "]", sep = "")], out.mat[, paste("mH[2, ", t2, "]", sep = "")])
	mOj <- cbind(out.mat[, paste("mO[2, ", t1, "]", sep = "")], out.mat[, paste("mO[2, ", t2, "]", sep = "")])
	Sj <- exp(-(mHj + mOj))
 
	mHa <- cbind(out.mat[, paste("mH[1, ", t1, "]", sep = "")], out.mat[, paste("mH[1, ", t2, "]", sep = "")])
	mOa <- cbind(out.mat[, paste("mO[1, ", t1, "]", sep = "")], out.mat[, paste("mO[1, ", t2, "]", sep = "")])
	Sa <- exp(-(mHa + mOa))
	 
	m0 <- cbind(out.mat[, paste("m0t[", t1+1, "]", sep = "")], out.mat[, paste("m0t[", t2+1, "]", sep = "")])
	S0 <- exp(-m0)
	
	
	Psi2 <- cbind(out.mat[, paste("Psi[2, ", t1+1, "]", sep = "")], out.mat[, paste("Psi[2, ", t2+1, "]", sep = "")])
	Psi3 <- cbind(out.mat[, paste("Psi[3, ", t1+1, "]", sep = "")], out.mat[, paste("Psi[3, ", t2+1, "]", sep = "")])
	Psi4 <- cbind(out.mat[, paste("Psi[4, ", t1+1, "]", sep = "")], out.mat[, paste("Psi[4, ", t2+1, "]", sep = "")])
	Psi5 <- cbind(out.mat[, paste("Psi[5, ", t1+1, "]", sep = "")], out.mat[, paste("Psi[5, ", t2+1, "]", sep = "")])
	
	rho2 <- cbind(out.mat[, paste("rho[2, ", t1+1, "]", sep = "")], out.mat[, paste("rho[2, ", t2+1, "]", sep = "")])
	rho3 <- cbind(out.mat[, paste("rho[3, ", t1+1, "]", sep = "")], out.mat[, paste("rho[3, ", t2+1, "]", sep = "")])
	rho4 <- cbind(out.mat[, paste("rho[4, ", t1+1, "]", sep = "")], out.mat[, paste("rho[4, ", t2+1, "]", sep = "")])
	rho5 <- cbind(out.mat[, paste("rho[5, ", t1+1, "]", sep = "")], out.mat[, paste("rho[5, ", t2+1, "]", sep = "")])


	## Extract population sizes for relevant time-steps (all samples)
	# NOTE: Here, I am also taking the additional time-step t2+1 (for calculating lambda at t2)
	N_1 <- cbind(out.mat[, paste("N[1, ", t1, "]", sep="")], out.mat[, paste("N[1, ", t2, "]", sep="")], out.mat[, paste("N[1, ", t2+1, "]", sep="")])
	N_2 <- cbind(out.mat[, paste("N[2, ", t1, "]", sep="")], out.mat[, paste("N[2, ", t2, "]", sep="")], out.mat[, paste("N[2, ", t2+1, "]", sep="")])
	N_3 <- cbind(out.mat[, paste("N[3, ", t1, "]", sep="")], out.mat[, paste("N[3, ", t2, "]", sep="")], out.mat[, paste("N[3, ", t2+1, "]", sep="")])
	N_4 <- cbind(out.mat[, paste("N[4, ", t1, "]", sep="")], out.mat[, paste("N[4, ", t2, "]", sep="")], out.mat[, paste("N[4, ", t2+1, "]", sep="")])
	N_5 <- cbind(out.mat[, paste("N[5, ", t1, "]", sep="")], out.mat[, paste("N[5, ", t2, "]", sep="")], out.mat[, paste("N[5, ", t2+1, "]", sep="")])
	
	N_tot <- N_1 + N_2 + N_3 + N_4 + N_5
	

	## Extract immigrant numbers for relevant time-steps (all samples)
	Imm <- cbind(out.mat[, paste("Imm[", t1, "]", sep="")], out.mat[, paste("Imm[", t2, "]", sep="")])

	
	## Calculate population proportions and growth rates
	n_1 <- N_1[,1:2] / N_tot[,1:2]
	n_2 <- N_2[,1:2] / N_tot[,1:2]
	n_3 <- N_3[,1:2] / N_tot[,1:2]
	n_4 <- N_4[,1:2] / N_tot[,1:2]
	n_5 <- N_5[,1:2] / N_tot[,1:2]
	
	lambda <- cbind(N_tot[,2]/N_tot[,1], N_tot[,3]/N_tot[,2])



	#------------------------------------------------
	# STEP 2) Calculate derivatives (for mean values)
	#------------------------------------------------

	sens_Sj <- (0.5*rowMeans(Psi2)*rowMeans(rho2)*rowMeans(S0) + 1)*rowMeans(n_1)
	sens_mHj <- -exp(-(rowMeans(mHj) + rowMeans(mOj)))*sens_Sj
	sens_mOj <- sens_mHj

	sens_Sa <- (0.5*rowMeans(Psi3)*rowMeans(rho3)*rowMeans(S0) + 1)*rowMeans(n_2) + (0.5*rowMeans(Psi4)*rowMeans(rho4)*rowMeans(S0) + 1)*rowMeans(n_3) + (0.5*rowMeans(Psi5)*rowMeans(rho5)*rowMeans(S0) + 1)*rowMeans(n_4 + n_5)
	sens_mHa <- -exp(-(rowMeans(mHa) + rowMeans(mOa)))*sens_Sj
	sens_mOa <- sens_mHa 

	sens_S0 <- 0.5*(rowMeans(Psi2)*rowMeans(rho2)*rowMeans(n_1) + rowMeans(Psi3)*rowMeans(rho3)*rowMeans(n_2) + rowMeans(Psi4)*rowMeans(rho4)*rowMeans(n_3) + rowMeans(Psi5)*rowMeans(rho5)*rowMeans(n_4 + n_5))
	sens_m0 <- -exp(-rowMeans(m0))*sens_S0

	sens_Psi2 <- 0.5*rowMeans(S0)*rowMeans(rho2)*rowMeans(n_1)
	sens_Psi3 <- 0.5*rowMeans(S0)*rowMeans(rho3)*rowMeans(n_2)
	sens_Psi4 <- 0.5*rowMeans(S0)*rowMeans(rho4)*rowMeans(n_3)
	sens_Psi5 <- 0.5*rowMeans(S0)*rowMeans(rho5)*rowMeans(n_4 + n_5)

	sens_rho2 <- 0.5*rowMeans(S0)*rowMeans(Psi2)*rowMeans(n_1)
	sens_rho3 <- 0.5*rowMeans(S0)*rowMeans(Psi3)*rowMeans(n_2)
	sens_rho4 <- 0.5*rowMeans(S0)*rowMeans(Psi4)*rowMeans(n_3)
	sens_rho5 <- 0.5*rowMeans(S0)*rowMeans(Psi5)*rowMeans(n_4 + n_5)

	sens_N_1 <- (rowMeans(Sj)*(0.5*rowMeans(Psi2)*rowMeans(rho2)*rowMeans(S0) + 1) - rowMeans(lambda)) / rowMeans(N_tot[,1:2])
	sens_N_2 <- (rowMeans(Sa)*(0.5*rowMeans(Psi3)*rowMeans(rho3)*rowMeans(S0) + 1) - rowMeans(lambda)) / rowMeans(N_tot[,1:2])
	sens_N_3 <- (rowMeans(Sa)*(0.5*rowMeans(Psi4)*rowMeans(rho4)*rowMeans(S0) + 1) - rowMeans(lambda)) / rowMeans(N_tot[,1:2])
	sens_N_4 <- (rowMeans(Sa)*(0.5*rowMeans(Psi5)*rowMeans(rho5)*rowMeans(S0) + 1) - rowMeans(lambda)) / rowMeans(N_tot[,1:2])
	sens_N_5 <- sens_N_4


	sens_Imm <- 1 / rowMeans(N_tot[,1:2])


	#----------------------------------------------------
	# STEP 3) Calculate contributions (fixed design LTRE)
	#----------------------------------------------------

	cont_Sj <- (Sj[,2] - Sj[,1])*sens_Sj
	cont_mHj <- (mHj[,2] - mHj[,1])*sens_mHj
	cont_mOj <- (mOj[,2] - mOj[,1])*sens_mOj
	
	cont_Sa <- (Sa[,2]-Sa[,1])*sens_Sa
	cont_mHa <- (mHa[,2] - mHa[,1])*sens_mHa
	cont_mOa <- (mOa[,2] - mOa[,1])*sens_mOa

	cont_S0 <- (S0[,2] - S0[,1])*sens_S0
	cont_m0 <- (m0[,2] - m0[,1])*sens_m0

	cont_Psi2 <- (Psi2[,2] - Psi2[,1])*sens_Psi2
	cont_Psi3 <- (Psi3[,2] - Psi3[,1])*sens_Psi3
	cont_Psi4 <- (Psi4[,2] - Psi4[,1])*sens_Psi4
	cont_Psi5 <- (Psi5[,2] - Psi5[,1])*sens_Psi5

	cont_rho2 <- (rho2[,2] - rho2[,1])*sens_rho2
	cont_rho3 <- (rho3[,2] - rho3[,1])*sens_rho3
	cont_rho4 <- (rho4[,2] - rho4[,1])*sens_rho4
	cont_rho5 <- (rho5[,2] - rho5[,1])*sens_rho5

	cont_N_1 <- (N_1[,2] - N_1[,1])*sens_N_1
	cont_N_2 <- (N_2[,2] - N_2[,1])*sens_N_2
	cont_N_3 <- (N_3[,2] - N_3[,1])*sens_N_3
	cont_N_4 <- (N_4[,2] - N_4[,1])*sens_N_4
	cont_N_5 <- (N_5[,2] - N_5[,1])*sens_N_5

	cont_Imm <- (Imm[,2] - Imm[,1])*sens_Imm


	#-----------------------------------
	# STEP 4) Arrange and output results
	#-----------------------------------
	
	## Arrange contributions in a data frame
	cont.all <- cbind(cont_Sj, cont_mHj, cont_mOj,
					  cont_Sa, cont_mHa, cont_mOa,
					  cont_S0, cont_m0,
				  	  cont_Psi2, cont_Psi3, cont_Psi4, cont_Psi5,
				  	  cont_rho2, cont_rho3, cont_rho4, cont_rho5,
				  	  cont_N_1, cont_N_2, cont_N_3, cont_N_4, cont_N_5,
				  	  cont_Imm)

	cont.data <- as.data.frame(cont.all)
	
	## Add information on which time-steps were compared (t1 & t2)
	cont.data$t1 <- t1
	cont.data$t2 <- t2
	
	## Return data
	return(cont.data)
}


#################################################################
#### RUNNING FIXED DESIGN LTRE FOR ALL SUBSEQUENT TIME-STEPS ####
#################################################################

## Make a list of pairs of matrices / time-steps to compare
t1 <- c(1:21)
t2 <- t1 + 1
years <- cbind(t1, t2)

## Run analysis for all pairs of matrices / time-steps
results.fixed <- do.call("rbind", sapply(1:nrow(years), FUN = function(x) fixedD.LTRE(t1 = years[x,1], t2 = years[x,2]), simplify = FALSE))


