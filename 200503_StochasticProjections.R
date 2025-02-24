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
load('200429_AF_IPM_VersionB4.RData')

out.mat <- as.matrix(AF.IPM.varB4)


#-------------------------------------------#
# GENERATE AND STORE YEAR-SPECIFIC MATRICES #
#-------------------------------------------#

## Set maximum time-step in data
Tmax <- 23

## Prepare an array to store projection matrices
mat.Y <- array(0, dim = c(5, 5, Tmax-1))

## Prepare a matrix to store initial population sizes
N0.Y <- matrix(NA, nrow = 5, ncol = Tmax-1)

## Prepare a vector to store immigrant numbers
Imm.Y <- rep(NA, Tmax-1)

## Extract year-specific matrices and population sizes based on posterior means
for(t in 1:(Tmax-1)){
	
	# Extract vital rates
	Sj <- median(exp(-(out.mat[, paste("mH[2, ", t, "]", sep = "")] + out.mat[, paste("mO[2, ", t, "]", sep = "")])))
  	Sa <- median(exp(-(out.mat[, paste("mH[1, ", t, "]", sep = "")] + out.mat[, paste("mO[1, ", t, "]", sep = "")])))
  
  	Psi2 <- median(out.mat[, paste("Psi[2, ", t+1, "]", sep = "")])
  	Psi3 <- median(out.mat[, paste("Psi[3, ", t+1, "]", sep = "")])
  	Psi4 <- median(out.mat[, paste("Psi[4, ", t+1, "]", sep = "")])
  	Psi5 <- median(out.mat[, paste("Psi[5, ", t+1, "]", sep = "")])
  	
  	rho2 <- median(out.mat[, paste("rho[2, ", t+1, "]", sep = "")])
  	rho3 <- median(out.mat[, paste("rho[3, ", t+1, "]", sep = "")])
  	rho4 <- median(out.mat[, paste("rho[4, ", t+1, "]", sep = "")])
  	rho5 <- median(out.mat[, paste("rho[5, ", t+1, "]", sep = "")])
  	
  	S0 <- median(exp(-out.mat[, paste("m0t[", t+1, "]", sep = "")]))
	
	
	# Write projection matrix
	mat.Y[1, 1, t] <- Sj*0.5*Psi2*rho2*S0
	mat.Y[1, 2, t] <- Sa*0.5*Psi3*rho3*S0
	mat.Y[1, 3, t] <- Sa*0.5*Psi4*rho4*S0
	mat.Y[1, 4, t] <- Sa*0.5*Psi5*rho5*S0
	mat.Y[1, 5, t] <- Sa*0.5*Psi5*rho5*S0
	
	mat.Y[2, 1, t] <- Sj
	
	mat.Y[3, 2, t] <- Sa
	
	mat.Y[4, 3, t] <- Sa
	
	mat.Y[5, 4, t] <- Sa
	mat.Y[5, 5, t] <- Sa

	
	# Extract starting population size
	N0.Y[1, t] <- median(out.mat[, paste("N[1, ", t, "]", sep = "")])
	N0.Y[2, t] <- median(out.mat[, paste("N[2, ", t, "]", sep = "")])
	N0.Y[3, t] <- median(out.mat[, paste("N[3, ", t, "]", sep = "")])
	N0.Y[4, t] <- median(out.mat[, paste("N[4, ", t, "]", sep = "")])
	N0.Y[5, t] <- median(out.mat[, paste("N[5, ", t, "]", sep = "")])
	
	# Extract immigrant numbers
	Imm.Y[t] <- median(out.mat[, paste("Imm[", t+1, "]", sep = "")])
}


#--------------------------------------------------------------------------#
# FUNCTION FOR STOCHASTIC PROJECTION USING A PRE-DEFINED SEQUENCE OF YEARS #
#--------------------------------------------------------------------------#

stoch.Proj = function(YearSeq){
	
	# Determine projection length
	Tmax <- length(YearSeq)
	
	# Make population vector
	N <- matrix(NA, nrow = 5, ncol = Tmax+1)
	
	# Set starting population size
	N[,1] <- N0.Y[,YearSeq[1]]
	
	# Project population
	for(t in 1:Tmax){
		
		N[,t+1] <- (mat.Y[,,YearSeq[t]]%*%N[,t]) + c(Imm.Y[YearSeq[t]], 0, 0, 0, 0)
	}
	
	return(N)	
}

#--------------------------------------------------------------------------#
# FUNCTION TO RUN A SERIES OF PROJECTIONS AND RETURN RESULTS AS DATA FRAME #
#--------------------------------------------------------------------------#

stoch.Sims = function(i, data.Tmax, sim.Tmax){ # i = simulation number
	
	# Make a random sequence of years
	YearSeq <- sample(c(1:(data.Tmax-1)), sim.Tmax, replace = T)
	
	# Population projection
	N <- stoch.Proj(YearSeq)
	
	# Arrange projection results in a data frame
	output <- data.frame(
	SimNo = i, 
	SimYear = c(1:(length(YearSeq)+1)), 
	Ntot = colSums(N),
	N1 = N[1,],
	N2 = N[2,],
	N3 = N[3,],
	N4 = N[4,],
	N5 = N[5,],
	p1 = N[1,]/colSums(N),
	p2 = N[2,]/colSums(N),
	p3 = N[3,]/colSums(N),
	p4 = N[4,]/colSums(N),
	p5 = N[5,]/colSums(N)
	)
	
	# Return results
	return(output)
}


#------------------------------------------#
# RUNNING SIMULATIONS AND PLOTTING RESULTS #
#------------------------------------------#

## Set number of replicates to simulate
SimNoMax <- 200

## Set number of years to simulate for
SimYearMax <- 100

## Running simulations
sim.results <- do.call("rbind", sapply(1:SimNoMax, FUN = function(i) stoch.Sims(i, data.Tmax = 23, sim.Tmax = SimYearMax), simplify = FALSE))

## Plot simulations over time

# Total population size
pdf("Simulation_PopSize.pdf", width = 6.5, height = 4)
ggplot(sim.results, aes(x = SimYear, y = Ntot, group = SimNo)) + geom_line(color = magma(8)[4], alpha = 0.3, size = 0.3) + geom_hline(aes(yintercept = 0), linetype = 'dotted', size = 0.3) + ylab('Population Size') + xlab('Simulation Year') + theme_bw() + theme(panel.grid = element_blank())
dev.off()

pdf("Simulation_LogPopSize.pdf", width = 6.5, height = 4)
ggplot(sim.results, aes(x = SimYear, y = log(Ntot), group = SimNo)) + geom_line(color = magma(8)[4], alpha = 0.3, size = 0.3) + ylab('Log Population Size') + xlab('Simulation Year') + theme_bw() + theme(panel.grid = element_blank())
dev.off()


# Population age structure
p1 <- ggplot(sim.results, aes(x = SimYear, y = p1, group = SimNo)) + geom_line(color = cividis(5)[1], alpha = 0.3, size = 0.3) + ylim(0, max(sim.results$p1)) + ylab('Proportion') + xlab('Simulation Year') + theme_bw() + theme(panel.grid = element_blank()) + ggtitle('Age class 1')

p2 <- ggplot(sim.results, aes(x = SimYear, y = p2, group = SimNo)) + geom_line(color = cividis(5)[2], alpha = 0.3, size = 0.3) + ylim(0, max(sim.results$p1)) + ylab('Proportion') + xlab('Simulation Year') + theme_bw() + theme(panel.grid = element_blank()) + ggtitle('Age class 2')

p3 <- ggplot(sim.results, aes(x = SimYear, y = p3, group = SimNo)) + geom_line(color = cividis(5)[3], alpha = 0.3, size = 0.3) + ylim(0, max(sim.results$p1)) + ylab('Proportion') + xlab('Simulation Year') + theme_bw() + theme(panel.grid = element_blank()) + ggtitle('Age class 3')

p4 <- ggplot(sim.results, aes(x = SimYear, y = p4, group = SimNo)) + geom_line(color = cividis(5)[4], alpha = 0.3, size = 0.3) + ylim(0, max(sim.results$p1)) + ylab('Proportion') + xlab('Simulation Year') + theme_bw() + theme(panel.grid = element_blank()) + ggtitle('Age class 4')

p5 <- ggplot(sim.results, aes(x = SimYear, y = p5, group = SimNo)) + geom_line(color = cividis(5)[5], alpha = 0.3, size = 0.3) + ylim(0, max(sim.results$p1)) + ylab('Proportion') + xlab('Simulation Year') + theme_bw() + theme(panel.grid = element_blank()) + ggtitle('Age class 5')


pdf("Simulation_PopStructure.pdf", width = 6.5, height = 8)
grid.arrange(
p1 + theme(axis.title.x = element_blank(), axis.text.x = element_blank()), 
p2 + theme(axis.title.x = element_blank(), axis.text.x = element_blank()), 
p3 + theme(axis.title.x = element_blank(), axis.text.x = element_blank()), 
p4 + theme(axis.title.x = element_blank(), axis.text.x = element_blank()), 
p5, ncol = 1, heights = c(1, 1, 1, 1, 1.25))
dev.off()

