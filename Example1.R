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
library(spBayesSurv)


## ----include=FALSE------------------------------------------------------------------------------------------------------------------
#source("/Users/FJRubio/Dropbox/ODESurv/ODESurv/Codes/routines/routines.R")
#source("/Users/javierrubio/Dropbox/ODESurv/ODESurv/Codes/routines/routines.R")
source("C:/Users/Javier/Dropbox/ODESurv/ODESurv/Codes/routines/routines.R")


## ----eval=FALSE---------------------------------------------------------------------------------------------------------------------
## source("routines.R")


## -----------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Data preparation
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

data(LeukSurv)
head(LeukSurv)

# New data frame: logical status, time in years, survival times sorted
df <- data.frame(time = LeukSurv$time, status = LeukSurv$cens)
df$status <- as.logical(LeukSurv$cens)
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
# Logistic ODE model for the hazard function: Analytic solution
#--------------------------------------------------------------------------------------------------

# Initial point
initL <- c(log(0.1),log(2),log(1))

# Optimisation step
OPTL <- nlminb(initL , log_likL, control = list(iter.max = 10000)) 
OPTL

# Fitted Logistic hazard: Analytic solution
fithazL <- Vectorize(function(t) hlogisode( t,exp(OPTL$par[1]), exp(OPTL$par[2]), exp(OPTL$par[3]) ) )

# Fitted Logistic cumulative hazard: Analytic solution
fitchazL <- Vectorize(function(t) chlogisode( t,exp(OPTL$par[1]), exp(OPTL$par[2]), exp(OPTL$par[3]) ) )

# Fitted Logistic survival: Analytic solution
fitsL <- Vectorize(function(t) exp(-chlogisode( t,exp(OPTL$par[1]), exp(OPTL$par[2]), exp(OPTL$par[3]) ) ))


# AIC
AICL <- 2*OPTL$objective + 2*length(OPTL$par)

# AIC
BICL <- 2*OPTL$objective + length(OPTL$par)*log(length(survtimes))


## -----------------------------------------------------------------------------------------------------------------------------------

# AIC comparison
AICs <- c(AICW, AICL)

AICs

# Best model (Logistic)
which.min(AICs)

# BIC comparison
BICs <- c(BICW, BICL)

BICs

# Best model (Logistic)
which.min(BICs)

# Comparison: hazard functions
curve(fithw, 0, max(survtimes), lwd= 2, lty = 1, col = "black", ylim = c(0,3.5), n = 1000,
      xlab = "Time", ylab = "Hazard", main = "", cex.axis = 1.5, cex.lab = 1.5)
curve(fithazL, 0, max(survtimes), lwd= 2, lty = 2, add = TRUE, n = 1000)

legend("topright", legend = c("Weibull", "Logistic"), lty = c(1,2), 
       lwd = c(2,2), col = c("black","black"))

# Comparison: cumulative hazard functions
curve(fitchw, 0, max(survtimes), lwd= 2, lty = 1, col = "black", ylim = c(0,3.5), n = 1000,
      xlab = "Time", ylab = "Cumulative Hazard", main = "", cex.axis = 1.5, cex.lab = 1.5)
curve(fitchazL, 0, max(survtimes), lwd= 2, lty = 2, add = TRUE, n = 1000)

legend("topleft", legend = c("Weibull", "Logistic"), lty = c(1,2), 
       lwd = c(2,2), col = c("black","black","black"))

# Comparison: survival functions

# Kaplan-Meier estimator 
km <- survfit(Surv(survtimes, status) ~ 1)


curve(fitsw, 0, max(survtimes), lwd= 2, lty = 1, col = "black", ylim = c(0,1), n = 1000,
      xlab = "Time", ylab = "Survival", main = "", cex.axis = 1.5, cex.lab = 1.5)
curve(fitsL, 0, max(survtimes), lwd= 2, lty = 2, add = TRUE, n = 1000)
points(km$time, km$surv, type = "l", col = "gray", lwd = 2, lty = 1)


legend("topright", legend = c("Weibull", "Logistic", "Kaplan-Meier"), lty = c(1,2,1), 
       lwd = c(2,2,2), col = c("black","black","gray"))



## -----------------------------------------------------------------------------------------------------------------------------------
#==================================================================================================
# Bayesian Analysis
#==================================================================================================

#--------------------------------------------------------------------------------------------------
# Priors
#--------------------------------------------------------------------------------------------------

par(mfrow = c(1,2))
p_sigma <- Vectorize(function(t) dgamma(t, shape = 2, scale = 2))
curve(p_sigma,0,15, n = 1000, xlab = ~sigma, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2)

p_nu <- Vectorize(function(t) dgamma(t, shape = 2, scale = 2))
curve(p_nu,0,15, n = 1000, xlab = ~nu, ylab = "Prior Density", 
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
sigmap <- exp(infoW$output[,1][ind])
nup <- exp(infoW$output[,2][ind])

plot(density(sigmap), main = "", xlab = expression(sigma), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
curve(p_sigma,0.75,1.5, n = 1000, xlab = ~lambda, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE)


plot(density(nup), main = "", xlab = expression(nu), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
curve(p_nu,0.40,0.6, n = 1000, xlab = ~kappa, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE)


## -----------------------------------------------------------------------------------------------------------------------------------
#==================================================================================================
# Bayesian Analysis
#==================================================================================================

#--------------------------------------------------------------------------------------------------
# Priors
#--------------------------------------------------------------------------------------------------

par(mfrow = c(1,3))
p_lambda <- Vectorize(function(t) dgamma(t, shape = 2, scale = 2))
curve(p_lambda,0,15, n = 1000, xlab = ~lambda, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2)

p_kappa <- Vectorize(function(t) dgamma(t, shape = 2, scale = 2))
curve(p_kappa,0,15, n = 1000, xlab = ~kappa, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2)

p_h0 <- Vectorize(function(t) dgamma(t, shape = 2, scale = 2))
curve(p_h0,0,15, n = 1000, xlab = expression(h[0]), ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2)

par(mfrow = c(1,1))


#--------------------------------------------------------------------------------------------------
# Logistic ODE model for the hazard function: Analytic solution
#--------------------------------------------------------------------------------------------------

# Support
SupportL <- function(x) {   TRUE }

# Random initial points
X0L <- function(x) { OPTL$par + runif(3,-0.01,0.01) }

# twalk for analytic solution
set.seed(1234)
infoL <- Runtwalk( dim=3,  Tr=110000,  Obj=log_postL, Supp=SupportL, 
                   x0=X0L(), xp0=X0L(),PlotLogPost = FALSE) 


# Posterior sample after burn-in and thinning
ind=seq(10000,110000,100) 

# Summaries
summL <- apply(exp(infoL$output[ind,]),2,summary)
colnames(summL) <- c("lambda","kappa","h_0")
kable(summL, digits = 3)

# KDEs
lambdap <- exp(infoL$output[,1][ind])
kappap <- exp(infoL$output[,2][ind])
h0p <- exp(infoL$output[,3][ind])

plot(density(lambdap), main = "", xlab = expression(lambda), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
curve(p_lambda,0,0.5, n = 1000, xlab = ~lambda, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE)


plot(density(kappap), main = "", xlab = expression(kappa), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
curve(p_kappa,0,0.15, n = 1000, xlab = ~kappa, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE)


plot(density(h0p), main = "", xlab = "h_0", ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
curve(p_h0,2,7, n = 1000, xlab = expression(h[0]), ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE)



## -----------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# Predictive Hazard functions
#---------------------------------------------------------------------------------------

# Predictive Weibull hazard: Analytic solution
predhW <- Vectorize(function(t){
  num <- den <- temp <- vector()
  for(i in 1:length(ind)){
    num[i] <- exp(-chweibull( t, sigmap[i], nup[i]))*hweibull( t, sigmap[i], nup[i])
    den[i] <- exp(-chweibull( t, sigmap[i], nup[i]))
  }
  return(mean(num)/mean(den))
})

# Predictive Logistic hazard: Analytic solution
predhL <- Vectorize(function(t){
  num <- den <- temp <- vector()
  for(i in 1:length(ind)) num[i] <- exp(-chlogisode( t,lambdap[i], kappap[i], h0p[i]))*hlogisode( t,lambdap[i], kappap[i], h0p[i])
  for(i in 1:length(ind)) den[i] <- exp(-chlogisode( t,lambdap[i], kappap[i], h0p[i]))
  return(mean(num)/mean(den))
})

# Creating the credible envelopes
tvec <- seq(0,14,by = 0.01)
ntvec <- length(tvec)

# Weibull
hCIW <- matrix(0, ncol = ntvec, nrow = length(ind))

for(j in 1:length(ind)){
  for(k in 1:ntvec){
    hCIW[j,k ] <- hweibull( tvec[k],sigmap[j], nup[j]) 
  }
} 

hW <-  predhW(tvec)


hCIWL <- apply(hCIW, 2, ql)
hCIWU <- apply(hCIW, 2, qu)

# Logistic
hCIL <- matrix(0, ncol = ntvec, nrow = length(ind))

for(j in 1:length(ind)){
  for(k in 1:ntvec){
    hCIL[j,k ] <- hlogisode( tvec[k],lambdap[j], kappap[j], h0p[j]) 
  }
} 

hL <-  predhL(tvec)

hCILL <- apply(hCIL, 2, ql)
hCILU <- apply(hCIL, 2, qu)

# Plots

# Weibull
plot(tvec,  hW, type = "l", ylim = c(0,4), xlab = "Time", ylab = "Predictive Hazard", 
     cex.axis = 1.5, cex.lab = 1.5, lwd =1)
points(tvec,  hCIWL, col = "gray", type = "l")
points(tvec,  hCIWU, col = "gray", type = "l")
polygon(c(tvec, rev(tvec)), c(hCIWL[order(tvec)], rev(hCIWU[order(tvec)])),
        col = "gray", border = NA)
points(tvec,  hW,type = "l", col = "black", lwd = 1)

# Logistic
points(tvec,  hL, type = "l", ylim = c(0,4), xlab = "Time", ylab = "Predictive Hazard", 
     cex.axis = 1.5, cex.lab = 1.5, lwd =1, lty = 2)
points(tvec,  hCILL, col = "gray", type = "l")
points(tvec,  hCILU, col = "gray", type = "l")
polygon(c(tvec, rev(tvec)), c(hCILL[order(tvec)], rev(hCILU[order(tvec)])),
        col = "gray", border = NA)
points(tvec,  hL,type = "l", col = "black", lwd = 1, lty =2)


legend("topright", legend = c("Weibull", "Logistic"), lty = c(1,2), 
       lwd = c(2,2), col = c("black","black"))


## -----------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# Predictive Survival functions
#---------------------------------------------------------------------------------------

# Predictive Weibull survival
predsW <- Vectorize(function(t){
  temp <- vector()
  for(i in 1:length(ind)) temp[i] <- exp(-chweibull(t,sigmap[i], nup[i]) ) 
  return(mean(temp))
})

# Predictive Logistic survival: Analytic solution
predsL <- Vectorize(function(t){
  temp <- vector()
  for(i in 1:length(ind)) temp[i] <- exp(-chlogisode( t,lambdap[i], kappap[i], h0p[i]) ) 
  return(mean(temp))
})



# Weibull
SCIW <- matrix(0, ncol = ntvec, nrow = length(ind))

for(j in 1:length(ind)){
  for(k in 1:ntvec){
    SCIW[j,k ] <- exp(-chweibull( tvec[k],sigmap[j], nup[j]))
  }
} 

SW <-  predsW(tvec)


SCIWL <- apply(SCIW, 2, ql)
SCIWU <- apply(SCIW, 2, qu)

# Logistic
SCIL <- matrix(0, ncol = ntvec, nrow = length(ind))

for(j in 1:length(ind)){
  for(k in 1:ntvec){
    SCIL[j,k ] <- exp(-chlogisode( tvec[k],lambdap[j], kappap[j], h0p[j])) 
  }
} 

SL <-  predsL(tvec)

SCILL <- apply(SCIL, 2, ql)
SCILU <- apply(SCIL, 2, qu)


# Plots (no credible envelopes)

# Weibull
plot(tvec,  SW, type = "l", ylim = c(0,1), xlab = "Time", ylab = "Predictive Survival", 
     cex.axis = 1.5, cex.lab = 1.5, lwd =1)

# Logistic
points(tvec,  SL, type = "l", ylim = c(0,1), xlab = "Time", ylab = "Predictive Survival", 
       cex.axis = 1.5, cex.lab = 1.5, lwd =1, lty = 2)


legend("topright", legend = c("Weibull", "Logistic"), lty = c(1,2), 
       lwd = c(2,2), col = c("black","black"))


# Plots (with credible envelopes)

# Weibull
plot(tvec,  SW, type = "l", ylim = c(0,1), xlab = "Time", ylab = "Predictive Survival", 
     cex.axis = 1.5, cex.lab = 1.5, lwd =1)
points(tvec,  SCIWL, col = "gray", type = "l")
points(tvec,  SCIWU, col = "gray", type = "l")
polygon(c(tvec, rev(tvec)), c(SCIWL[order(tvec)], rev(SCIWU[order(tvec)])),
        col = "gray", border = NA)
points(tvec,  SW,type = "l", col = "black", lwd = 1)

# Logistic
points(tvec,  SL, type = "l", ylim = c(0,1), xlab = "Time", ylab = "Predictive Survival", 
       cex.axis = 1.5, cex.lab = 1.5, lwd =1, lty = 2)
points(tvec,  SCILL, col = "gray", type = "l")
points(tvec,  SCILU, col = "gray", type = "l")
polygon(c(tvec, rev(tvec)), c(SCILL[order(tvec)], rev(SCILU[order(tvec)])),
        col = "gray", border = NA)
points(tvec,  SL,type = "l", col = "black", lwd = 1, lty =2)

legend("topright", legend = c("Weibull", "Logistic"), lty = c(1,2), 
       lwd = c(2,2), col = c("black","black"))


