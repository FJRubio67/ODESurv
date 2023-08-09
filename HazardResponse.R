## ---- include=FALSE---------------------------------------------------------------------------------------
# Required packages
library(deSolve)
library(demography)
library(survival)


## ----include=FALSE----------------------------------------------------------------------------------------
source("/Users/FJRubio/Dropbox/ODESurv/ODESurv/Codes/routines/routines.R")


## ----eval=FALSE-------------------------------------------------------------------------------------------
## source("routines.R")


## ---------------------------------------------------------------------------------------------------------

# Parameter values
lambda <- 1.8
kappa <- 0.1
alpha <- 6
beta <- 4.8
h0 <- 1e-2
q0 <- 1e-6

# True parameter values for plots
trueparf <- c(lambda, kappa, alpha, beta, kappa, alpha, h0, q0)
# True parameter values for simulations
truepar <- c(lambda, kappa, alpha, beta)

# Hazard function: Solver
paramsHR  <-  c(lambda, kappa, alpha, beta)
initHT0      <- c(h = h0, q = q0, H = 0 )
times = seq(0,36,by = 0.1)
out <- ode(initHT0, times, hazmodHR, paramsHR, method = "lsode", jacfunc = jacODE, jactype = "fullusr")


# Hazard function
plot(as.matrix(out[,c(1,2)]), xlim = c(0,max(times)), ylim = c(0,0.125), 
     type = "l", lwd = 2,
     xlab = "Time", ylab = "Hazard", main = "Hazard-Response ODE", 
     cex.axis = 1.5, cex.lab = 1.5)
points(as.matrix(out[,c(1,3)]), lty=2, type="l", lwd = 2)
legend("bottomright", legend=c("Hazard","Response"), lwd = c(2,2), lty = c(1,2))

# Survival function
plot(as.vector(out[,1]),exp(-as.vector(out[,4])), xlim = c(0,max(times)), ylim = c(0,1), type = "l", lwd = 2,
     xlab = "Time", ylab = "Survival", main = "Hazard-Response ODE", cex.axis = 1.5, cex.lab = 1.5)



## ---------------------------------------------------------------------------------------------------------

# Approximate simulation
set.seed(1234)
sim <- simHR(n = 5000, par = truepar, tmin = 0, tmax = 150, step = 1e-3)$sim

# Censoring time
cens <- 20

# Vital status
status <- ifelse(sim < cens, 1, 0)

mean(status)

# Observed times
time <- ifelse(sim < cens, sim, cens)

# Kaplan-Meier estimator for the simulated times
kmsim <- survfit(Surv(time, status) ~ 1)

# Comparison
plot(as.vector(out[,1]),exp(-as.vector(out[,4])), xlim = c(0,max(times)), ylim = c(0,1), type = "l", lwd = 2,
     xlab = "Time", ylab = "Survival", main = "Hazard-Response ODE", cex.axis = 1.5, cex.lab = 1.5)
points(kmsim$time, kmsim$surv, type = "l", col = "gray", lwd = 2, lty = 1, ylim = c(0,1),
     xlab = "Time", ylab = "Survival")
points(kmsim$time, kmsim$surv, type = "l", col = "gray", lwd = 2, lty = 2)
legend("topright", legend = c("True model","Kaplan-Meier"), col = c("black","gray"), lwd = c(2,2), lty = c(1,2))


