## ----message=FALSE------------------------------------------------------------------------------------------------------------------
#######################################################################################
# Data preparation
#######################################################################################

rm(list=ls())

# Required packages
library(deSolve)
library(survival)
library(ggplot2)
#library(devtools)
#install_github("FJRubio67/HazReg")
library(HazReg)
library(Rtwalk)
library(knitr)
library(spBayes)
library(bshazard)


## ----include=FALSE------------------------------------------------------------------------------------------------------------------
source("/Users/FJRubio/Dropbox/ODESurv/ODESurv/Codes/routines/routines.R")
#source("C:/Users/Javier/Dropbox/ODESurv/Codes/routines/routines.R")


## ----eval=FALSE---------------------------------------------------------------------------------------------------------------------
## source("routines.R")


## -----------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Data preparation
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


head(rotterdam)

dim(rotterdam)

# Kaplan-Meier estimator for the survival times
km <- survfit(Surv(rotterdam$dtime/365.24, rotterdam$death) ~ 1)

plot(km$time, km$surv, type = "l", col = "black", lwd = 2, lty = 1, 
     ylim = c(0,1), xlab = "Time", ylab = "Survival")



# New data frame: logical status, time in years, survival times sorted
df <- data.frame(time = rotterdam$dtime, status = rotterdam$death)
df$status <- as.logical(rotterdam$death)
df$time <- df$time/365.24

df <- df[order(df$time),]

# Required quantities
status <- as.logical(df$status)
t_obs <- df$time[status]
survtimes <- df$time


## -----------------------------------------------------------------------------------------------------------------------------------
#==================================================================================================
# Maximum Likelihood Analysis
#==================================================================================================


#--------------------------------------------------------------------------------------------------
# Fitting a Weibull distribution
#--------------------------------------------------------------------------------------------------

# Initial value
initW <- c(0,0)

# Optimisation step
OPTW <- GHMLE(initW, survtimes, status, hstr = "baseline", dist = "Weibull", method = "nlminb", maxit = 10000)

MLEW <- exp(OPTW$OPT$par)

# Fitted Weibull hazard
fithw <- Vectorize(function(t) hweibull( t, exp(OPTW$OPT$par[1]), exp(OPTW$OPT$par[2]) ) )

# Fitted Weibull cumulative hazard
fitchw <- Vectorize(function(t) chweibull( t, exp(OPTW$OPT$par[1]), exp(OPTW$OPT$par[2]) ) )

# Fitted Weibull survival
fitsw <- Vectorize(function(t) exp(-chweibull( t, exp(OPTW$OPT$par[1]), exp(OPTW$OPT$par[2]) ) ))

# AIC
AICW <- 2*OPTW$OPT$objective + 2*length(OPTW$OPT$par)

# BIC
BICW <- 2*OPTW$OPT$objective + length(OPTW$OPT$par)*log(length(survtimes))


## -----------------------------------------------------------------------------------------------------------------------------------
#==================================================================================================
# Maximum Likelihood Analysis
#==================================================================================================

#--------------------------------------------------------------------------------------------------
# Fitting a PGW distribution
#--------------------------------------------------------------------------------------------------

# Initial value
initPGW <- c(0,0,0)

# Optimisation step
OPTPGW <- GHMLE(initPGW, survtimes, status, hstr = "baseline", dist = "PGW", 
                method = "nlminb", maxit = 10000)

MLEPGW <- exp(OPTPGW$OPT$par)

# Fitted Weibull hazard
fithpgw <- Vectorize(function(t) hpgw( t, MLEPGW[1], MLEPGW[2], MLEPGW[3] ) )

# Fitted Weibull cumulative hazard
fitchpgw <- Vectorize(function(t) chpgw( t, MLEPGW[1], MLEPGW[2], MLEPGW[3] ) )

# Fitted Weibull survival
fitspgw <- Vectorize(function(t) exp(-chpgw( t, MLEPGW[1], MLEPGW[2], MLEPGW[3] ) ))

# AIC
AICPGW <- 2*OPTPGW$OPT$objective + 2*length(OPTPGW$OPT$par)

# BIC
BICPGW <- 2*OPTPGW$OPT$objective + length(OPTPGW$OPT$par)*log(length(survtimes))



## -----------------------------------------------------------------------------------------------------------------------------------
#==================================================================================================
# MLE Analysis
#==================================================================================================

#--------------------------------------------------------------------------------------------------
# Hazard-Response ODE model for the hazard function: Solver solution
#--------------------------------------------------------------------------------------------------

# Initial point
initHR <- log(c(lambda = 2.5, kappa = 0.1, alpha = 2, beta = 2.))
h0 = 1e-2
q0 = 1e-6

# Optimisation step
OPTHR <- nlminb(initHR , log_likHR, control = list(iter.max = 10000)) 

OPTHR

MLEHR <- c(exp(OPTHR$par[1:4]))

# AIC
AICHR <- 2*OPTHR$objective + 2*length(OPTHR$par)

# BIC
BICHR <- 2*OPTHR$objective + length(OPTHR$par)*log(length(survtimes))


## -----------------------------------------------------------------------------------------------------------------------------------
# AIC comparison
AICs <- c(AICW, AICPGW, AICHR)

AICs

# Best model (Hazard-Response)
which.min(AICs)

# BIC comparison
BICs <- c(BICW, BICPGW, BICHR)

BICs

# Best model (Hazard-Response)
which.min(BICs)

# Fitted hazard and cumulative hazard: Solver
paramsHR  <- c(lambda = exp(OPTHR$par[1]), kappa = exp(OPTHR$par[2]), alpha = exp(OPTHR$par[3]), beta = exp(OPTHR$par[4]))
initHR0      <- c(h = h0, q = q0, H = 0 )
times = seq(0,20,by = 0.1)
out <- ode(initHR0, times, hazmodHR, paramsHR, method = "lsode", jacfunc = jacODE, jactype = "fullusr")


# Comparison: hazard functions
plot(as.matrix(out[,c(1,2)]), xlim = c(0,max(times)), ylim = c(0,0.125), type = "l", lwd = 2,
     xlab = "Time", ylab = "Hazard", main = "Hazard-Response ODE", cex.axis = 1.5, cex.lab = 1.5)
curve(fithw, 0, max(times), lwd= 2, lty = 2, col = "gray", add = TRUE)
curve(fithpgw, 0, max(times), lwd= 2, lty = 3, col = "gray", add = TRUE)
legend("topleft", legend = c("HR Solver","Weibull", "PGW"), lty = c(1,2,3), 
       lwd = c(2,2,2), col = c("black","gray","gray"))

# Comparison: hazard vs. treament
plot(as.matrix(out[,c(1,2)]), xlim = c(0,max(times)), ylim = c(0,0.125), type = "l", lwd = 2,
     xlab = "Time", ylab = "Hazard", main = "Hazard-Response ODE", cex.axis = 1.5, cex.lab = 1.5)
points(as.matrix(out[,c(1,3)]), xlim = c(0,max(times)), ylim = c(0,0.125), type = "l", lwd = 2, lty = 2,
     xlab = "Time", ylab = "Hazard", main = "Hazard-Response ODE", cex.axis = 1.5, cex.lab = 1.5)
abline(h = MLEHR[2], lwd = 2, col = "red")

# Comparison: cumulative hazard functions
plot(as.matrix(out[,c(1,4)]), xlim = c(0,max(times)), ylim = c(0,2.25), type = "l", lwd = 2,
     xlab = "Time", ylab = "Cumulative Hazard", main = "Hazard-Response ODE", cex.axis = 1.5, cex.lab = 1.5)
curve(fitchw, 0, max(times), lwd= 2, lty = 2, col = "gray", add = TRUE)
curve(fitchpgw, 0, max(times), lwd= 2, lty = 3, col = "gray", add = TRUE)
legend("topleft", legend = c("HR Solver","Weibull","PGW"), lty = c(1,2,3), 
       lwd = c(2,2,2), col = c("black","gray","gray"))

# Comparison: survival functions
plot(as.vector(out[,1]),exp(-as.vector(out[,4])), xlim = c(0,max(times)), ylim = c(0,1), type = "l", lwd = 2,
     xlab = "Time", ylab = "Survival", main = "Hazard-Response ODE", cex.axis = 1.5, cex.lab = 1.5)
curve(fitsw, 0, max(times), lwd= 2, lty = 2, col = "gray", add = TRUE)
curve(fitspgw, 0, max(times), lwd= 2, lty = 3, col = "gray", add = TRUE)
points(km$time, km$surv, type = "l", col = "green", lwd = 2, lty = 1)
legend("topright", legend = c("HR Solver","Weibull","PGW","KM"), lty = c(1,2,3,1), 
       lwd = c(2,2,2,2), col = c("black","gray","gray","green"))



## -----------------------------------------------------------------------------------------------------------------------------------
#==================================================================================================
# Bayesian Analysis
#==================================================================================================

#--------------------------------------------------------------------------------------------------
# Priors
#--------------------------------------------------------------------------------------------------

par(mfrow = c(1,2))
p_sigmaW <- Vectorize(function(t) dgamma(t, shape = 2, scale = 2))
curve(p_sigmaW,0,15, n = 1000, xlab = ~sigma, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2)

p_nuW <- Vectorize(function(t) dgamma(t, shape = 2, scale = 2))
curve(p_nuW,0,15, n = 1000, xlab = ~nu, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2)


par(mfrow = c(1,1))

#--------------------------------------------------------------------------------------------------
# Weibull model for the hazard function
#--------------------------------------------------------------------------------------------------

# Support
SupportW <- function(x) {   TRUE }

# Random initial points
X0W <- function(x) { OPTW$OPT$par + runif(2,-0.01,0.01) }

# twalk for analytic solution
set.seed(1234)
infoW <- Runtwalk( dim=2,  Tr=110000,  Obj=log_postW, Supp=SupportW, 
                   x0=X0W(), xp0=X0W(), PlotLogPost = FALSE) 


# Posterior sample after burn-in and thinning
ind=seq(10000,110000,100) 

# Summaries
summW <- apply(exp(infoW$output[ind,]),2,summary)
colnames(summW) <- c("sigma","nu")
kable(summW, digits = 3)

# KDEs
sigmapW <- exp(infoW$output[,1][ind])
nupW <- exp(infoW$output[,2][ind])

plot(density(sigmapW), main = "", xlab = expression(sigma), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
curve(p_sigmaW,13,17, n = 1000, xlab = ~lambda, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE)


plot(density(nupW), main = "", xlab = expression(nu), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
curve(p_nuW,1.1,1.5, n = 1000, xlab = ~kappa, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE)


## -----------------------------------------------------------------------------------------------------------------------------------
#==================================================================================================
# Bayesian Analysis
#==================================================================================================

#--------------------------------------------------------------------------------------------------
# Priors
#--------------------------------------------------------------------------------------------------

par(mfrow = c(1,3))
p_sigmaPGW <- Vectorize(function(t) dgamma(t, shape = 2, scale = 2))
curve(p_sigmaPGW,0,15, n = 1000, xlab = ~sigma, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2)

p_nuPGW <- Vectorize(function(t) dgamma(t, shape = 2, scale = 2))
curve(p_nuPGW,0,15, n = 1000, xlab = ~nu, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2)

p_gammaPGW <- Vectorize(function(t) dgamma(t, shape = 2, scale = 2))
curve(p_gammaPGW,0,15, n = 1000, xlab = ~gamma, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2)

par(mfrow = c(1,1))

#--------------------------------------------------------------------------------------------------
# Weibull model for the hazard function
#--------------------------------------------------------------------------------------------------

# Support
SupportPGW <- function(x) {   TRUE }

# Random initial points
X0PGW <- function(x) { OPTPGW$OPT$par + runif(3,-0.01,0.01) }


# twalk for analytic solution
set.seed(1234)
infoPGW <- Runtwalk( dim=3,  Tr=110000,  Obj=log_postPGW, Supp=SupportPGW, 
                     x0=X0PGW(), xp0=X0PGW(), PlotLogPost = FALSE) 


# Posterior sample after burn-in and thinning
ind=seq(10000,110000,100) 

# Summaries
summPGW <- apply(exp(infoPGW$output[ind,]),2,summary)
colnames(summPGW) <- c("sigma","nu","gamma")
kable(summPGW, digits = 3)

# KDEs
sigmapPGW <- exp(infoPGW$output[,1][ind])
nupPGW <- exp(infoPGW$output[,2][ind])
gammapPGW <- exp(infoPGW$output[,3][ind])

plot(density(sigmapPGW), main = "", xlab = expression(sigma), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
curve(p_sigmaPGW,2,5, n = 1000, xlab = ~lambda, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE)


plot(density(nupPGW), main = "", xlab = expression(nu), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
curve(p_nuPGW,1.5,2.7, n = 1000, xlab = ~kappa, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE)

plot(density(gammapPGW), main = "", xlab = expression(gamma), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
curve(p_gammaPGW,0,15, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE)


## -----------------------------------------------------------------------------------------------------------------------------------
#==================================================================================================
# Bayesian Analysis
#==================================================================================================

#--------------------------------------------------------------------------------------------------
# Priors
#--------------------------------------------------------------------------------------------------

par(mfrow = c(2,2))
p_lambdaHR <- Vectorize(function(t) dgamma(t, shape = 2, scale = 2))
curve(p_lambdaHR,0,15, n = 1000, xlab = ~lambda, ylab = "Prior Density",
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2)

p_kappaHR <- Vectorize(function(t) dgamma(t, shape = 2, scale = 2))
curve(p_kappaHR,0,15, n = 1000, xlab = ~kappa, ylab = "Prior Density",
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2)

p_alphaHR <- Vectorize(function(t) dgamma(t, shape = 2, scale = 2))
curve(p_alphaHR,0,15, n = 1000, xlab = ~alpha, ylab = "Prior Density",
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2)

p_betaHR <- Vectorize(function(t) dgamma(t, shape = 2, scale = 2))
curve(p_betaHR,0,15, n = 1000, xlab = ~beta, ylab = "Prior Density",
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2)

par(mfrow = c(1,1))


#--------------------------------------------------------------------------------------------------
# Hazard-Response ODE model for the hazard function: Solver solution
#--------------------------------------------------------------------------------------------------


n.batch <- 1100
batch.length <- 100

lp <- function(par) -log_postHR(par)

inits <- OPTHR$par


## ----include=FALSE------------------------------------------------------------------------------------------------------------------
set.seed(1234)
infoHR <- adaptMetropGibbs(ltd=lp, starting=inits, accept.rate=0.44, batch=n.batch, 
                           batch.length=batch.length, report=100, verbose=FALSE)


## ----eval=FALSE---------------------------------------------------------------------------------------------------------------------
## set.seed(1234)
## infoHR <- adaptMetropGibbs(ltd=lp, starting=inits, accept.rate=0.44, batch=n.batch,
##                            batch.length=batch.length, report=100, verbose=FALSE)


## -----------------------------------------------------------------------------------------------------------------------------------
chainHR <- infoHR$p.theta.samples[,1:4]

# Burning and thinning the chain
burn <- 1e4
thin <- 100
NS <- n.batch*batch.length
ind <- seq(burn,NS,thin)



# Summaries
summHR <- apply(exp(chainHR[ind,1:4]),2,summary)
colnames(summHR) <- c("lambda","kappa","alpha","beta")
kable(summHR, digits = 3)

# KDEs
lambdapHR <- exp(chainHR[,1][ind])
kappapHR <- exp(chainHR[,2][ind])
alphapHR <- exp(chainHR[,3][ind])
betapHR <- exp(chainHR[,4][ind])


plot(density(lambdapHR), main = "", xlab = expression(lambda), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
curve(p_lambdaHR,0,3, n = 1000, xlab = ~lambda, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE)


plot(density(kappapHR), main = "", xlab = expression(kappa), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
curve(p_kappaHR,0.05,0.3, n = 1000, xlab = ~kappa, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE)

plot(density(alphapHR), main = "", xlab = expression(alpha), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
curve(p_alphaHR,0,10, n = 1000, xlab = ~alpha, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE)

plot(density(betapHR), main = "", xlab = expression(beta), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
curve(p_betaHR,0,15, n = 1000, xlab = ~beta, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE)



## -----------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# Predictive Hazard functions
#---------------------------------------------------------------------------------------

# Predictive Weibull hazard
predhW <- Vectorize(function(t){
  num <- den <- temp <- vector()
  for(i in 1:length(ind)){
    num[i] <- exp(-chweibull( t, sigmapW[i], nupW[i]))*hweibull( t, sigmapW[i], nupW[i])
    den[i] <- exp(-chweibull( t, sigmapW[i], nupW[i]))
  }
  return(mean(num)/mean(den))
})

# Predictive PGW
predhPGW <- Vectorize(function(t){
  num <- den <- temp <- vector()
  for(i in 1:length(ind)) num[i] <- exp(-chpgw( t, sigmapPGW[i], nupPGW[i], gammapPGW[i]))*
      hpgw( t, sigmapPGW[i], nupPGW[i], gammapPGW[i])
  for(i in 1:length(ind)) den[i] <- exp(-chpgw( t, sigmapPGW[i], nupPGW[i], gammapPGW[i]))
  return(mean(num)/mean(den))
})



# Creating the credible envelopes
tvec <- seq(0,20,by = 0.01)
ntvec <- length(tvec)

# Weibull
hCIW <- matrix(0, ncol = ntvec, nrow = length(ind))

for(j in 1:length(ind)){
  for(k in 1:ntvec){
    hCIW[j,k ] <- hweibull( tvec[k],sigmapW[j], nupW[j]) 
  }
} 

hW <-  predhW(tvec)


hCIWL <- apply(hCIW, 2, ql)
hCIWU <- apply(hCIW, 2, qu)

# PGW
hCIPGW <- matrix(0, ncol = ntvec, nrow = length(ind))

for(j in 1:length(ind)){
  for(k in 1:ntvec){
    hCIPGW[j,k ] <- hpgw( tvec[k],sigmapPGW[j], nupPGW[j], gammapPGW[j]) 
  }
} 

hPGWs <-  predhPGW(tvec)

hCIPGWL <- apply(hCIPGW, 2, ql)
hCIPGWU <- apply(hCIPGW, 2, qu)

# Hazard-Response

hCIHR <- matrix(0, ncol = ntvec, nrow = length(ind))
chCIHR <- matrix(0, ncol = ntvec, nrow = length(ind))
SCIHR <- matrix(0, ncol = ntvec, nrow = length(ind))
qCIHR <- matrix(0, ncol = ntvec, nrow = length(ind))

for(j in 1:length(ind)){
  paramsHRj  <- c(lambda = lambdapHR[j], kappa = kappapHR[j], alpha = alphapHR[j],
                  beta = betapHR[j])
  initHRj      <- c(h = h0, q = q0, H = 0 )
  outj <- ode(initHRj, tvec, hazmodHR, paramsHRj, method = "lsode", jacfunc = jacODE, jactype = "fullusr")
  
  hCIHR[j, ] <- outj[,2]
  qCIHR[j, ] <- outj[,3]
  chCIHR[j, ] <- outj[,4]
  SCIHR[j, ] <- exp(-outj[,4])
  
} 


numHR <- hCIHR*exp(-chCIHR)
denHR <- exp(-chCIHR)

hpredHR <- colMeans(numHR)/colMeans(denHR)


hCIHRL <- apply(hCIHR, 2, ql)
hCIHRU <- apply(hCIHR, 2, qu)

qCIHRL <- apply(qCIHR, 2, ql)
qCIHRU <- apply(qCIHR, 2, qu)


# Plots

# Weibull
plot(tvec,  hW, type = "l", ylim = c(0,0.15), xlab = "Time", ylab = "Predictive Hazard", 
     cex.axis = 1.5, cex.lab = 1.5, lwd =1)
points(tvec,  hCIWL, col = "gray", type = "l")
points(tvec,  hCIWU, col = "gray", type = "l")
polygon(c(tvec, rev(tvec)), c(hCIWL[order(tvec)], rev(hCIWU[order(tvec)])),
        col = "gray", border = NA)
points(tvec,  hW,type = "l", col = "black", lwd = 1)

# PGW
points(tvec,  hPGWs, type = "l", xlab = "Time", ylab = "Predictive Hazard", 
       cex.axis = 1.5, cex.lab = 1.5, lwd =1, lty = 2)
points(tvec,  hCIPGWL, col = "gray", type = "l")
points(tvec,  hCIPGWU, col = "gray", type = "l")
polygon(c(tvec, rev(tvec)), c(hCIPGWL[order(tvec)], rev(hCIPGWU[order(tvec)])),
        col = "gray", border = NA)
points(tvec,  hPGWs,type = "l", col = "black", lwd = 1, lty =2)

# Hazard-Response
points(tvec,  hpredHR, type = "l", xlab = "Time", ylab = "Predictive Hazard", 
       cex.axis = 1.5, cex.lab = 1.5, lwd =1, lty = 3)
points(tvec,  hCIHRL, col = "gray", type = "l")
points(tvec,  hCIHRU, col = "gray", type = "l")
polygon(c(tvec, rev(tvec)), c(hCIHRL[order(tvec)], rev(hCIHRU[order(tvec)])),
        col = "gray", border = NA)
points(tvec,  hpredHR,type = "l", col = "black", lwd = 1, lty =3)
points(tvec,  hPGWs,type = "l", col = "black", lwd = 1, lty =2)
legend("bottomright", legend = c("Weibull", "PGW", "HR"), lty = c(1,2,3), 
       lwd = c(2,2,2), col = c("black","black","black"))





## -----------------------------------------------------------------------------------------------------------------------------------
# Hazard vs. Response
plot(tvec,  hpredHR, type = "l", xlab = "Time", ylab = "Predictive Hazard", 
       cex.axis = 1.5, cex.lab = 1.5, lwd =2, lty = 1, ylim = c(0,0.1))

qpredHR <- colMeans(qCIHR)

points(tvec,  qpredHR, type = "l", xlab = "Time", ylab = "Predictive Hazard", 
       cex.axis = 1.5, cex.lab = 1.5, lwd =2, lty = 2)

legend("bottomright", legend = c("Hazard","Response"), lty = c(1,2), 
       lwd = c(2,2), col = c("black","black"))


## -----------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------- 
# Predictive Survival functions 
#--------------------------------------------------------------------------------------- 

# Predictive Weibull survival 
predsW <- Vectorize(function(t){ 
  temp <- vector() 
  for(i in 1:length(ind)) temp[i] <- exp(-chweibull(t,sigmapW[i], nupW[i]) )  
  return(mean(temp)) 
}) 

# Predictive PGW survival 
predsPGW <- Vectorize(function(t){ 
  temp <- vector() 
  for(i in 1:length(ind)) temp[i] <- exp(-chpgw( t,sigmapPGW[i], nupPGW[i], gammapPGW[i]) )  
  return(mean(temp)) 
}) 


# Predictive Hazard-Response survival: Solver solution 
predsHR <- colMeans(denHR)


# Weibull 
SCIW <- matrix(0, ncol = ntvec, nrow = length(ind)) 

for(j in 1:length(ind)){ 
  for(k in 1:ntvec){ 
    SCIW[j,k ] <- exp(-chweibull( tvec[k],sigmapW[j], nupW[j])) 
  } 
}

SW <-  predsW(tvec) 

SCIWL <- apply(SCIW, 2, ql) 
SCIWU <- apply(SCIW, 2, qu) 

# PGW

SCIPGW <- matrix(0, ncol = ntvec, nrow = length(ind)) 

for(j in 1:length(ind)){ 
  for(k in 1:ntvec){ 
    SCIPGW[j,k ] <- exp(-chpgw( tvec[k],sigmapPGW[j], nupPGW[j], gammapPGW[j])) 
  } 
}

SPGW <-  predsPGW(tvec) 

SCIPGWL <- apply(SCIPGW, 2, ql) 
SCIPGWU <- apply(SCIPGW, 2, qu) 


# Hazard-Response 
SpredCIHR <- colMeans(denHR) 

SCIHRL <- apply(SCIHR, 2, ql) 
SCIHRU <- apply(SCIHR, 2, qu) 

# Plots (no credible envelopes) 

# Weibull 
plot(tvec,  SW, type = "l", ylim = c(0,1), xlab = "Time", ylab = "Predictive Survival",  
     cex.axis = 1.5, cex.lab = 1.5, lwd =1) 

# PGW 
points(tvec,  SPGW, type = "l", ylim = c(0,1), xlab = "Time", ylab = "Predictive Survival",  
       cex.axis = 1.5, cex.lab = 1.5, lwd =1, lty = 2) 

# Hazard-Response 
points(tvec,  SpredCIHR, type = "l", ylim = c(0,1), xlab = "Time", ylab = "Predictive Survival",  
       cex.axis = 1.5, cex.lab = 1.5, lwd =1, lty = 3) 

legend("topright", legend = c("Weibull", "PGW", "HR"), lty = c(1,2,3),  
       lwd = c(2,2,2), col = c("black","black","black")) 


# Plots (with credible envelopes) 

# Weibull 
plot(tvec,  SW, type = "l", ylim = c(0,1), xlab = "Time", ylab = "Predictive Survival",  
     cex.axis = 1.5, cex.lab = 1.5, lwd =1) 
points(tvec,  SCIWL, col = "gray", type = "l") 
points(tvec,  SCIWU, col = "gray", type = "l") 
polygon(c(tvec, rev(tvec)), c(SCIWL[order(tvec)], rev(SCIWU[order(tvec)])), 
        col = "gray", border = NA) 
points(tvec,  SW,type = "l", col = "black", lwd = 1) 

# PGW 
points(tvec,  SPGW, type = "l", ylim = c(0,1), xlab = "Time", ylab = "Predictive Survival",  
       cex.axis = 1.5, cex.lab = 1.5, lwd =1, lty = 2) 
points(tvec,  SCIPGWL, col = "gray", type = "l") 
points(tvec,  SCIPGWU, col = "gray", type = "l") 
polygon(c(tvec, rev(tvec)), c(SCIPGWL[order(tvec)], rev(SCIPGWU[order(tvec)])), 
        col = "gray", border = NA) 
points(tvec,  SPGW,type = "l", col = "black", lwd = 1, lty =2) 

# Hazard-Response 
points(tvec,  SpredCIHR, type = "l", ylim = c(0,1), xlab = "Time", ylab = "Predictive Survival",  
       cex.axis = 1.5, cex.lab = 1.5, lwd =1, lty = 3) 
points(tvec,  SCIHRL, col = "gray", type = "l") 
points(tvec,  SCIHRU, col = "gray", type = "l") 
polygon(c(tvec, rev(tvec)), c(SCIHRL[order(tvec)], rev(SCIHRU[order(tvec)])), 
        col = "gray", border = NA) 
points(tvec,  SpredCIHR,type = "l", col = "black", lwd = 1, lty =3) 

points(tvec,  SPGW,type = "l", col = "black", lwd = 1, lty =2) 

legend("topright", legend = c("Weibull", "PGW", "HR"), lty = c(1,2,3),  
       lwd = c(2,2,2), col = c("black","black","black")) 



## ----message=FALSE------------------------------------------------------------------------------------------------------------------
# fitting the estimator
fit <- bshazard(Surv(df$time, df$status) ~ 1, data = df, nbin = 100, degree = 3, verbose = FALSE)


## -----------------------------------------------------------------------------------------------------------------------------------
# Hazard vs. Response
plot(tvec,  hpredHR, type = "l", xlab = "Time", ylab = "Hazard", 
       cex.axis = 1.5, cex.lab = 1.5, lwd =2, lty = 1, ylim = c(0,0.125))

# bshazard with confidence intervals
points(fit$time, fit$hazard, type='l', lwd = 2, ylim = c(0,0.15), col = "darkgray")
lines(fit$time, fit$lower.ci, lty = 2, lwd = 1, col = "gray")
lines(fit$time, fit$upper.ci, lty = 2, lwd = 1, col = "gray")

# Hazard vs. Response
points(tvec,  hpredHR, type = "l", xlab = "Time", ylab = "Hazard", 
     cex.axis = 1.5, cex.lab = 1.5, lwd =2, lty = 1, ylim = c(0,0.15))



## -----------------------------------------------------------------------------------------------------------------------------------

# Log-posterior (unknown initial conditions)
lpX <- function(par) -log_postHRX(par)

initsX <- c(OPTHR$par,h0,q0)


## ----include=FALSE----------------------------------------------------------------------------------------
set.seed(1234)
infoHRX <- adaptMetropGibbs(ltd=lpX, starting=initsX, accept.rate=0.44, batch=n.batch, 
                           batch.length=batch.length, report=100, verbose=FALSE)


chainHRX <- infoHRX$p.theta.samples[,1:6]


# Summaries
summHRX <- apply(exp(chainHRX[ind,1:6]),2,summary)
colnames(summHRX) <- c("lambda","kappa","alpha","beta","h0","q0")
kable(summHRX, digits = 3)

# KDEs
lambdapHRX <- exp(chainHRX[,1][ind])
kappapHRX <- exp(chainHRX[,2][ind])
alphapHRX <- exp(chainHRX[,3][ind])
betapHRX <- exp(chainHRX[,4][ind])
h0pHRX <- exp(chainHRX[,5][ind])
q0pHRX <- exp(chainHRX[,6][ind])


plot(density(lambdapHRX), main = "", xlab = expression(lambda), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
curve(p_lambdaHR,0,3, n = 1000, xlab = ~lambda, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE)


plot(density(kappapHRX), main = "", xlab = expression(kappa), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
curve(p_kappaHR,0.05,0.3, n = 1000, xlab = ~kappa, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE)

plot(density(alphapHRX), main = "", xlab = expression(alpha), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
curve(p_alphaHR,0,10, n = 1000, xlab = ~alpha, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE)

plot(density(betapHRX), main = "", xlab = expression(beta), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
curve(p_betaHR,0,15, n = 1000, xlab = ~beta, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE)

p_h0HR <- Vectorize(function(t) dgamma(t, shape = 2, scale = 5e-3))
plot(density(h0pHRX), main = "", xlab = expression(h[0]), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
curve(p_h0HR,0,0.05, n = 1000, xlab = ~beta, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE)

p_q0HR <- Vectorize(function(t) dgamma(t, shape = 2, scale = 5e-7))
plot(density(q0pHRX), main = "", xlab = expression(q[0]), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
curve(p_q0HR,0,7e-6, n = 1000, xlab = ~beta, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE)


# Comparison with the posteriors obtained by fixing the initial conditions

plot(density(lambdapHR), main = "", xlab = expression(lambda), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
points(density(lambdapHRX), lwd = 2, type = "l", lty = 2)
legend("topright", legend = c("fixed","prior"), lty = c(1,2), lwd = c(2,2))


plot(density(kappapHR), main = "", xlab = expression(kappa), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2, ylim = c(0,30))
points(density(kappapHRX), lwd = 2, type = "l", lty = 2)
legend("topright", legend = c("fixed","prior"), lty = c(1,2), lwd = c(2,2))


plot(density(alphapHR), main = "", xlab = expression(alpha), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
points(density(alphapHRX), lwd = 2, type = "l", lty = 2)
legend("topright", legend = c("fixed","prior"), lty = c(1,2), lwd = c(2,2))

plot(density(betapHR), main = "", xlab = expression(beta), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
points(density(betapHRX), lwd = 2, type = "l", lty = 2)
legend("topright", legend = c("fixed","prior"), lty = c(1,2), lwd = c(2,2))


