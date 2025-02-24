library(reshape)
library(plyr)

## Load LTRE results
load('200503_AF_LTRE_Random.RData')

str(LTRE.H.data) # LTRE with contributions of mortality hazard rates
str(est_var)
str(est_covar)

## Summarize estimated variances / covariances
var.sum <- t(apply(est_var, 1, function(x) quantile(x, probs = c(0.5, 0.05, 0.95))))
covar.sum <- t(apply(est_covar, 1, function(x) quantile(x, probs = c(0.5, 0.05, 0.95))))

## Summarize absolute LTRE contributions
cont.sum <- ddply(LTRE.H.data, .(variable), summarise, 
Q50 = quantile(abs(value), probs = 0.5), Q05 = quantile(abs(value), probs = 0.05), Q95 = quantile(abs(value), probs = 0.95))


## Organise results into a data frame
results <- data.frame(
	Parameter = c("mHj", "mOj", "mHa", "mOa", "m0", paste0("Psi", 1:4), paste0("rho", 1:4), paste0("N", 0:4), "I"),
	VarianceEst = paste0(round(var.sum[,1], digits = 2), " [", round(var.sum[,2], digits = 2), ", ", round(var.sum[,3], digits = 2), "]"),
	CovarianceEst = paste0(round(covar.sum[,1], digits = 2), " [", round(covar.sum[,2], digits = 3), ", ", round(covar.sum[,3], digits = 2), "]"),
	LTRECont = paste0(round(cont.sum[,2], digits = 3), " [", round(cont.sum[,3], digits = 3), ", ", round(cont.sum[,4], digits = 3), "]")
)

## Save results
write.csv(results, '200504_VarEst_VS_LTRECont.csv', row.names = F)
